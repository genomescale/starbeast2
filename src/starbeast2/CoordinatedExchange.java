// Based on Exchange.java in BEAST2, which in turn was based on ExchangeOperator.java in BEAST.

package starbeast2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
    public Input<MultispeciesCoalescent> mscInput = new Input<MultispeciesCoalescent>("speciescoalescent", "Multispecies coalescent (gene trees contained within a species tree).");

    private Node parent;
    private Node grandparent;
    private Node brother;
    private Node uncle;

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

        grandparent = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        while (grandparent.getLeft().isLeaf() && grandparent.getRight().isLeaf()) {
            grandparent = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        }

        parent = grandparent.getLeft();
        uncle = grandparent.getRight();
        if (parent.getHeight() < uncle.getHeight()) {
            parent = grandparent.getRight();
            uncle = grandparent.getLeft();
        }

        if( parent.isLeaf() ) {
            // tree with dated tips
            return Double.NEGATIVE_INFINITY;
        }

        int validGP = 0;
        {
            for(int i = nInternalNodes + 1; i < 1 + 2*nInternalNodes; ++i) {
                validGP += isg(speciesTree.getNode(i));
            }
        }

        final int c2 = sisg(parent) + sisg(uncle);

        brother = (Randomizer.nextBoolean() ? parent.getLeft() : parent.getRight());
        exchangeNodes(parent, grandparent, brother, uncle);

        final int validGPafter = validGP - c2 + sisg(parent) + sisg(uncle);

        return Math.log((float)validGP/validGPafter);
    }

    // for testing purposes
    public void manipulateSpeciesTree(Node argBrother) {
        brother = argBrother;
        parent = brother.getParent();
        grandparent = parent.getParent();
        if (grandparent.getLeft().getNr() == parent.getNr()) {
            uncle = grandparent.getRight();
        } else {
            uncle = grandparent.getLeft();
        }

        exchangeNodes(parent, grandparent, brother, uncle);
    }

    public double rearrangeGeneTrees(final MultispeciesCoalescent msc) {
        final List<GeneTreeWithinSpeciesTree> geneTrees = msc.getGeneTrees();
        final int nGeneTrees = geneTrees.size();
        final List<Map<Integer, Integer>> forwardGraftCounts = new ArrayList<Map<Integer, Integer>>();

        final int parentBranchNumber = parent.getNr();
        final int uncleBranchNumber = uncle.getNr();
        final int brotherBranchNumber = brother.getNr();

        for (int j = 0; j < nGeneTrees; j++) {
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            final Map<Integer, Integer> jForwardGraftCounts = new HashMap<Integer, Integer>();
            final Set<Node> parentBranchNodes = geneTree.branchNodeMap.get(parentBranchNumber);
            for (Node geneTreeNode: parentBranchNodes) {
                final Node leftChildNode = geneTreeNode.getLeft();
                final Node rightChildNode = geneTreeNode.getRight();
                final boolean leftContainsBrother = geneTree.lineageOverlap.get(leftChildNode.getNr()).contains(brotherBranchNumber);
                final boolean rightContainsBrother = geneTree.lineageOverlap.get(rightChildNode.getNr()).contains(brotherBranchNumber);

                // if exactly one gene tree node child branch exclusively descends via the "sister"
                // then this gene tree node needs to be pruned and reattached
                if (leftContainsBrother ^ rightContainsBrother) {
                    final double nodeHeight = geneTreeNode.getHeight();
                    final List<Node> validGraftBranches = findGraftBranches(geneTree, uncleBranchNumber, nodeHeight);
                    final int forwardGraftCount = validGraftBranches.size();
                    if (forwardGraftCount == 0) { // no compatible branches to graft this node on to
                        return Double.NEGATIVE_INFINITY;
                    } else {
                        jForwardGraftCounts.put(geneTreeNode.getNr(), forwardGraftCount);
                        final Node chosenNode = validGraftBranches.get(Randomizer.nextInt(forwardGraftCount));

                        // this gene tree node will be grafted to the branch defined by "chosenNode"
                        if (leftContainsBrother) { // the left child will be disowned and reattached directly to its grandparent
                            pruneAndRegraft(geneTreeNode, chosenNode, leftChildNode);
                            rightChildNode.makeDirty(Tree.IS_FILTHY);
                        } else { // the right child will be disowned and reattached directly to its grandparent
                            pruneAndRegraft(geneTreeNode, chosenNode, rightChildNode);
                            leftChildNode.makeDirty(Tree.IS_FILTHY);
                        }
                    }
                }
            }

            forwardGraftCounts.add(jForwardGraftCounts);
        }

        msc.computeCoalescentTimes(); // rebuild the gene-trees-within-species-tree to account for the tree changes

        double logHastingsRatio = 0.0;

        for (int j = 0; j < nGeneTrees; j++) {
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            final Map<Integer, Integer> geneTreeForwardGraftCounts = forwardGraftCounts.get(j);
            for (int geneTreeNodeNumber: geneTreeForwardGraftCounts.keySet()) {
                final double nodeHeight = geneTree.getNodeHeight(geneTreeNodeNumber);
                final List<Node> reverseGraftBranches = findGraftBranches(geneTree, brotherBranchNumber, nodeHeight);

                final int forwardGraftCount = geneTreeForwardGraftCounts.get(geneTreeNodeNumber);
                final int reverseGraftCount = reverseGraftBranches.size();

                logHastingsRatio += Math.log((float) forwardGraftCount / reverseGraftCount);
            }
        }

        return logHastingsRatio;
    }

    private static List<Node> findGraftBranches(final GeneTreeWithinSpeciesTree geneTree, final int uncleNumber, final double reattachHeight) {
        // identify branches in the former "uncle" (now "brother") which overlap with the height of the node to be moved
        final Set<Node> potentialGraftBranches = geneTree.branchOverlap.get(uncleNumber);
        final List<Node> validGraftBranches = new ArrayList<Node>();
        for (Node potentialGraftBranch: potentialGraftBranches) {
            final Node graftParent = potentialGraftBranch.getParent();

            final double tipwardHeight = potentialGraftBranch.getHeight();
            final double rootwardHeight = graftParent.getHeight();

            if ((reattachHeight >= tipwardHeight) && (reattachHeight <= rootwardHeight)) {
                validGraftBranches.add(potentialGraftBranch);
            }
        }

        return validGraftBranches;
    }

    // removes nodeSPR from the segmented line between its disownedChild and oldParent
    // reattaches it between the newChild node and the parent of the newChild node
    // does not change any node heights
    protected static void pruneAndRegraft(final Node nodeSPR, final Node newChild, final Node disownedChild) {
        final Node oldParent = nodeSPR.getParent();
        final Node newParent = newChild.getParent();

        oldParent.removeChild(nodeSPR);
        oldParent.addChild(disownedChild);
        oldParent.makeDirty(Tree.IS_FILTHY);

        newParent.removeChild(newChild);
        newParent.addChild(nodeSPR);
        newParent.makeDirty(Tree.IS_FILTHY);
        
        nodeSPR.removeChild(disownedChild);
        nodeSPR.addChild(newChild);
        nodeSPR.makeDirty(Tree.IS_FILTHY);

        newChild.makeDirty(Tree.IS_FILTHY);
        disownedChild.makeDirty(Tree.IS_FILTHY);
    }

    protected static void exchangeNodes(final Node parent, final Node grandparent, final Node brother, final Node uncle) {
        parent.removeChild(brother);
        parent.addChild(uncle);
        parent.makeDirty(Tree.IS_FILTHY);

        grandparent.removeChild(uncle);
        grandparent.addChild(brother);
        grandparent.makeDirty(Tree.IS_FILTHY);
        
        brother.makeDirty(Tree.IS_FILTHY);
        uncle.makeDirty(Tree.IS_FILTHY);
    }
}
