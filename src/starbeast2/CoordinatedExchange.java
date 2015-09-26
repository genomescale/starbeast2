// Based on Exchange.java in BEAST2, which in turn was based on ExchangeOperator.java in BEAST.

package starbeast2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;

import com.google.common.collect.SetMultimap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
* @author Remco Bouckaert
* @author Alexei Drummond
* @author Huw Ogilvie
* @author Andrew Rambaut
 */

@Description("Implements the co-ordinated species and gene tree operator described in Yang & Rannala 2015. "
        + "This performs a narrow exchange operation, then prunes-and-regrafts gene tree nodes made invalid "
        + "by the operation onto a valid contemporary branch (without changing node heights). "
        + "See http://doi.org/10.1093/molbev/msu279 for full details.")
public class CoordinatedExchange extends Operator {
    public Input<MultispeciesCoalescent> mscInput = new Input<MultispeciesCoalescent>("multispeciesCoalescent", "Multispecies coalescent (gene trees contained within a species tree).");

    private Node brother;
    private Node parent;
    private Node uncle;

    private List<SortedMap<Node, Node>> forwardMovedNodes;
    private SetMultimap<Integer, Node> forwardGraftNodes;

    @Override
    public void initAndValidate() {
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final MultispeciesCoalescent msc = mscInput.get();

        double fLogHastingsRatio = narrow(msc);

        // only rearrange gene trees if the species tree has changed
        if (fLogHastingsRatio != Double.NEGATIVE_INFINITY) {
            fLogHastingsRatio += rearrangeGeneTrees(msc);
        }

        return fLogHastingsRatio;
    }

    private int isg(final Node n) {
      return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    private int sisg(final Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }

    public double narrow(final MultispeciesCoalescent msc) {
        final Tree speciesTree = msc.getSpeciesTree();
        final int nInternalNodes = speciesTree.getInternalNodeCount();
        if (nInternalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        Node grandparent = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        while (grandparent.getLeft().isLeaf() && grandparent.getRight().isLeaf()) {
            grandparent = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        }

        parent = grandparent.getLeft();
        uncle = grandparent.getRight();
        if (parent.getHeight() < uncle.getHeight()) {
            parent = grandparent.getRight();
            uncle = grandparent.getLeft();
        }

        if (parent.isLeaf()) {
            return Double.NEGATIVE_INFINITY;
        }

        brother = (Randomizer.nextBoolean() ? parent.getLeft() : parent.getRight());

        // must be done before any changes are made to the gene or species trees
        forwardMovedNodes = msc.getMovedChildren(brother);
        forwardGraftNodes = msc.getGraftBranches(uncle);

        int validGP = 0;
        {
            for(int i = nInternalNodes + 1; i < 1 + 2*nInternalNodes; ++i) {
                validGP += isg(speciesTree.getNode(i));
            }
        }

        final int c2 = sisg(parent) + sisg(uncle);

        exchangeNodes(parent, grandparent, brother, uncle);

        final int validGPafter = validGP - c2 + sisg(parent) + sisg(uncle);

        return Math.log((float)validGP/validGPafter);
    }

    // for testing purposes
    public void manipulateSpeciesTree(Node argBrother) {
        brother = argBrother;
        parent = brother.getParent();

        final MultispeciesCoalescent msc = mscInput.get();
        final Node grandparent = parent.getParent();
        if (grandparent.getLeft().getNr() == parent.getNr()) {
            uncle = grandparent.getRight();
        } else {
            uncle = grandparent.getLeft();
        }

        forwardMovedNodes = msc.getMovedChildren(brother);
        forwardGraftNodes = msc.getGraftBranches(uncle);

        exchangeNodes(parent, grandparent, brother, uncle);
    }

    public double rearrangeGeneTrees(final MultispeciesCoalescent msc) {
        final List<GeneTreeWithinSpeciesTree> geneTrees = msc.getGeneTrees();
        final int nGeneTrees = geneTrees.size();
        final List<Map<Node, Integer>> forwardGraftCounts = new ArrayList<>();

        for (int j = 0; j < nGeneTrees; j++) {
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            final Map<Node, Integer> jForwardGraftCounts = new HashMap<>();

            final Set<Node> jForwardGraftNodes = forwardGraftNodes.get(j);
            final SortedMap<Node, Node> jForwardMovedNodes = forwardMovedNodes.get(j);

            for (Entry<Node, Node> nodeEntry: jForwardMovedNodes.entrySet()) {
                final Node movedNode = nodeEntry.getKey();
                final Node disownedChild = nodeEntry.getValue();
                final double movedNodeHeight = movedNode.getHeight();

                final List<Node> validGraftBranches = new ArrayList<>();
                int forwardGraftCount = 0;
                for (Node potentialGraft: jForwardGraftNodes) {
                    final double potentialGraftBottom = potentialGraft.getHeight();
                    double potentialGraftTop;

                    if (potentialGraft.isRoot()) {
                        potentialGraftTop = Double.POSITIVE_INFINITY;
                    } else {
                        potentialGraftTop = potentialGraft.getParent().getHeight();
                    }

                    if (movedNodeHeight >= potentialGraftBottom && movedNodeHeight < potentialGraftTop) {
                        forwardGraftCount++;
                        validGraftBranches.add(potentialGraft);
                    }
                }

                // no compatible branches to graft this node on to
                // this only occurs when there is missing data and the gene tree root is in the "parent" branch
                if (forwardGraftCount == 0) {
                    return Double.NEGATIVE_INFINITY;
                } else {
                    final Node chosenGraft = validGraftBranches.get(Randomizer.nextInt(forwardGraftCount));
                    pruneAndRegraft(movedNode, chosenGraft, disownedChild);
                    movedNode.makeDirty(Tree.IS_FILTHY);
                    jForwardGraftCounts.put(movedNode, forwardGraftCount);
                }
            }

            assert checkTreeSanity(geneTree.getRoot());

            forwardGraftCounts.add(jForwardGraftCounts);
        }

        assert msc.computeCoalescentTimes(); // this move should always preserve gene-tree-within-species-tree compatibility

        SetMultimap<Integer, Node> reverseGraftNodes = msc.getGraftBranches(brother);

        double logHastingsRatio = 0.0;
        for (int j = 0; j < nGeneTrees; j++) {
            final Map<Node, Integer> jForwardGraftCounts = forwardGraftCounts.get(j);
            final Set<Node> jReverseGraftNodes = reverseGraftNodes.get(j);
            for (Node movedNode: jForwardGraftCounts.keySet()) {
                final double movedNodeHeight = movedNode.getHeight();
                final int forwardGraftCount = jForwardGraftCounts.get(movedNode);
                int reverseGraftCount = 0;

                for (Node potentialGraft: jReverseGraftNodes) {
                    final double potentialGraftBottom = potentialGraft.getHeight();
                    double potentialGraftTop;
                    if (potentialGraft.isRoot()) {
                        potentialGraftTop = Double.POSITIVE_INFINITY;
                    } else {
                        potentialGraftTop = potentialGraft.getParent().getHeight();
                    }

                    if (movedNodeHeight >= potentialGraftBottom && movedNodeHeight < potentialGraftTop) {
                        reverseGraftCount++;
                    }
                }

                logHastingsRatio += Math.log(forwardGraftCount) - Math.log(reverseGraftCount);
            }
        }

        return logHastingsRatio;
    }

    private boolean checkTreeSanity(Node node) {
        List<Node> children = node.getChildren();
        for (Node childNode: children) {
            assert childNode.getParent() == node;
            assert childNode.getHeight() <= node.getHeight();
            if (!node.isLeaf()) {
                checkTreeSanity(childNode);
            }
        }

        if (node.isLeaf()) {
            assert children.size() == 0;
        } else {
            assert children.size() == 2;
        }

        return true;
    }

    // removes nodeToMove from the segmented line between its disownedChild and oldParent
    // reattaches it between the newChild node and the parent of the newChild node
    // does not change any node heights
    protected static void pruneAndRegraft(final Node nodeToMove, final Node newChild, final Node disownedChild) {
        final Node oldParent = nodeToMove.getParent();
        final Node newParent = newChild.getParent();
        
        assert oldParent != null;
        assert newParent != null;

        oldParent.addChild(disownedChild);
        nodeToMove.removeChild(disownedChild);
        nodeToMove.addChild(newChild);
        newParent.removeChild(newChild);

        if (oldParent != newParent) {
            oldParent.removeChild(nodeToMove);
            newParent.addChild(nodeToMove);
        }

        disownedChild.makeDirty(Tree.IS_FILTHY);
        nodeToMove.makeDirty(Tree.IS_FILTHY);
        newChild.makeDirty(Tree.IS_FILTHY);
        oldParent.makeDirty(Tree.IS_FILTHY);
        newParent.makeDirty(Tree.IS_FILTHY);
    }

    protected static void exchangeNodes(final Node parent, final Node grandparent, final Node brother, final Node uncle) {
        grandparent.removeChild(uncle);
        parent.addChild(uncle);

        parent.removeChild(brother);
        grandparent.addChild(brother);

        parent.makeDirty(Tree.IS_FILTHY);
        grandparent.makeDirty(Tree.IS_FILTHY);        
        brother.makeDirty(Tree.IS_FILTHY);
        uncle.makeDirty(Tree.IS_FILTHY);
    }
}
