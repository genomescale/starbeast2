// Based on Exchange.java in BEAST2, which in turn was based on ExchangeOperator.java in BEAST.

package starbeast2;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;
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

    private int parentBranchNumber;
    private int uncleBranchNumber;
    private int brotherBranchNumber;

    private ListMultimap<Integer, Node> forwardParentalNodes;
    private SetMultimap<Integer, Node> forwardFraternalNodes;
    private SetMultimap<Integer, Node> forwardAvuncularNodes;

    final static NodeHeightComparator nhc = new NodeHeightComparator();

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

        Node parent = grandparent.getLeft();
        Node uncle = grandparent.getRight();
        if (parent.getHeight() < uncle.getHeight()) {
            parent = grandparent.getRight();
            uncle = grandparent.getLeft();
        }

        if (parent.isLeaf()) {
            return Double.NEGATIVE_INFINITY;
        }

        Node brother = (Randomizer.nextBoolean() ? parent.getLeft() : parent.getRight());

        parentBranchNumber = parent.getNr();
        uncleBranchNumber = uncle.getNr();
        brotherBranchNumber = brother.getNr();

        // must be done before any changes are made to the gene or species trees
        forwardParentalNodes = msc.getBranchNodes(parentBranchNumber);
        forwardFraternalNodes = msc.getAssociatedNodes(brotherBranchNumber, parentBranchNumber);
        forwardAvuncularNodes = msc.getAssociatedNodes(uncleBranchNumber, uncleBranchNumber);

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
        final MultispeciesCoalescent msc = mscInput.get();
        final Node brother = argBrother;
        final Node parent = brother.getParent();
        final Node grandparent = parent.getParent();
        Node uncle;
        if (grandparent.getLeft().getNr() == parent.getNr()) {
            uncle = grandparent.getRight();
        } else {
            uncle = grandparent.getLeft();
        }

        brotherBranchNumber = brother.getNr();
        parentBranchNumber = parent.getNr();
        uncleBranchNumber = uncle.getNr();

        forwardParentalNodes = msc.getBranchNodes(parentBranchNumber);
        forwardFraternalNodes = msc.getAssociatedNodes(brotherBranchNumber, parentBranchNumber);
        forwardAvuncularNodes = msc.getAssociatedNodes(uncleBranchNumber, uncleBranchNumber);

        exchangeNodes(parent, grandparent, brother, uncle);
    }

    public double rearrangeGeneTrees(final MultispeciesCoalescent msc) {
        final List<GeneTreeWithinSpeciesTree> geneTrees = msc.getGeneTrees();
        final int nGeneTrees = geneTrees.size();
        final List<Multiset<Node>> forwardGraftCounts = new ArrayList<Multiset<Node>>();

        for (int j = 0; j < nGeneTrees; j++) {
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            final Multiset<Node> jForwardGraftCounts = HashMultiset.create();

            final Set<Node> fraternalNodes = forwardFraternalNodes.get(j);
            final Set<Node> avuncularNodes = forwardAvuncularNodes.get(j);

            // must evaluate gene tree nodes oldest to youngest (in node age)
            // otherwise destination graft branches will be severed below subsequent gene nodes
            final List<Node> parentBranchNodes = forwardParentalNodes.get(j);
            parentBranchNodes.sort(nhc);

            for (final Node geneTreeNode: parentBranchNodes) {
                final Node leftChildNode = geneTreeNode.getLeft();
                final Node rightChildNode = geneTreeNode.getRight();
                // final boolean leftContainsBrother = geneTree.checkOverlap(leftChildNode, brotherBranchNumber);
                // final boolean rightContainsBrother = geneTree.checkOverlap(rightChildNode, brotherBranchNumber);
                final boolean leftDescendsViaBrother = fraternalNodes.contains(leftChildNode);
                final boolean rightDescendsViaBrother = fraternalNodes.contains(rightChildNode);

                // if exactly one gene tree node child branch exclusively descends via the "sister"
                // then this gene tree node needs to be pruned and reattached
                if (leftDescendsViaBrother ^ rightDescendsViaBrother) {
                    final double nodeHeight = geneTreeNode.getHeight();

                    final List<Node> validGraftBranches = new ArrayList<>();
                    for (Node potentialGraft: avuncularNodes) {
                        final double potentialGraftBottom = potentialGraft.getHeight();
                        double potentialGraftTop;
                        if (potentialGraft.isRoot()) {
                            potentialGraftTop = Double.POSITIVE_INFINITY;
                        } else {
                            potentialGraftTop = potentialGraft.getParent().getHeight();
                        }
                        if (nodeHeight >= potentialGraftBottom && nodeHeight <= potentialGraftTop) {
                            jForwardGraftCounts.add(geneTreeNode);
                            validGraftBranches.add(potentialGraft);
                        }
                    }
                    // there should always be a compatible branch...
                    // ASSUMING NO MISSING DATA (so this assert is not valid for production starbeast2)
                    final int forwardGraftCount = jForwardGraftCounts.count(geneTreeNode);
                    assert forwardGraftCount != 0;
                    if (forwardGraftCount == 0) { // no compatible branches to graft this node on to
                        return Double.NEGATIVE_INFINITY;
                    } else {
                        final Node chosenNode = validGraftBranches.get(Randomizer.nextInt(forwardGraftCount));

                        // this gene tree node will be grafted to the branch defined by "chosenNode"
                        if (leftDescendsViaBrother) { // the left child will be disowned and reattached directly to its grandparent
                            pruneAndRegraft(geneTreeNode, chosenNode, leftChildNode);
                            rightChildNode.makeDirty(Tree.IS_FILTHY);
                        } else { // the right child will be disowned and reattached directly to its grandparent
                            pruneAndRegraft(geneTreeNode, chosenNode, rightChildNode);
                            leftChildNode.makeDirty(Tree.IS_FILTHY);
                        }
                    }
                }
            }

            assert checkTreeSanity(geneTree.getRoot());

            forwardGraftCounts.add(jForwardGraftCounts);
        }

        final boolean allCompatible = msc.computeCoalescentTimes(); // rebuild the gene-trees-within-species-tree to account for the tree changes
        assert allCompatible; // this move should always preserve gene-tree-within-species-tree compatibility

        final SetMultimap<Integer, Node> reverseAvuncularNodes = msc.getAssociatedNodes(brotherBranchNumber, brotherBranchNumber);
        double logHastingsRatio = 0.0;
        for (int j = 0; j < nGeneTrees; j++) {
            final Multiset<Node> jForwardGraftCounts = forwardGraftCounts.get(j);
            System.out.println(jForwardGraftCounts.size());
            final Set<Node> avuncularNodes = reverseAvuncularNodes.get(j);
            for (Node geneTreeNode: jForwardGraftCounts.elementSet()) {
                final double geneTreeNodeHeight = geneTreeNode.getHeight();
                final int forwardGraftCount = jForwardGraftCounts.count(geneTreeNode);
                int reverseGraftCount = 0;

                for (Node potentialGraftBranch: avuncularNodes) {
                    final double lowerHeight = potentialGraftBranch.getHeight();
                    double upperHeight;
                    if (potentialGraftBranch.isRoot()) {
                        upperHeight = Double.POSITIVE_INFINITY;
                    } else {
                        upperHeight = potentialGraftBranch.getParent().getHeight();
                    }

                    if (geneTreeNodeHeight >= lowerHeight && geneTreeNodeHeight <= upperHeight) {
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
    // TODO get this working when the source or destination nodes are the root node
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
