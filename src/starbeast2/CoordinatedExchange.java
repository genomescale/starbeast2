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
* @author Huw Ogilvie
 */

@Description("Implements the co-ordinated species and gene tree operator described in Yang & Rannala 2015. "
        + "This performs a narrow exchange operation, then prunes-and-regrafts gene tree nodes made invalid "
        + "by the operation onto a valid contemporary branch (without changing node heights). "
        + "See http://doi.org/10.1093/molbev/msu279 for full details.")
public class CoordinatedExchange extends Operator {
    public Input<MultispeciesCoalescent> mscInput = new Input<MultispeciesCoalescent>("speciescoalescent", "Multispecies coalescent (gene trees contained within a species tree).");
    MultispeciesCoalescent multispeciesCoalescent;

    @Override
    public void initAndValidate() {
    }

    private int isg(final Node n) {
      return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    private int sisg(final Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }

    // Based on the "narrow" method of the BEASTv2.3 Exchange operator.
    @Override
    public double proposal() {
        final MultispeciesCoalescent msc = mscInput.get();
        final Tree speciesTree = msc.getSpeciesTree();
        final List<GeneTreeWithinSpeciesTree> geneTrees = msc.getGeneTrees();
        final int nGeneTrees = geneTrees.size();
        final List<Map<Integer, Integer>> forwardGraftCounts = new ArrayList<Map<Integer, Integer>>();
        msc.computeCoalescentTimes(); // needs doing before any gene or species tree changes

        final int nInternalNodes = speciesTree.getInternalNodeCount();
        if (nInternalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        Node grandparent = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        while (grandparent.getLeft().isLeaf() && grandparent.getRight().isLeaf()) {
            grandparent = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        }

        Node parent = grandparent.getLeft();
        Node uncle  = grandparent.getRight();
        if (parent.getHeight() < uncle.getHeight()) {
            parent = grandparent.getRight();
            uncle  = grandparent.getLeft();
        }

        int validGP = 0;
        for(int i = nInternalNodes + 1; i < 1 + 2*nInternalNodes; ++i) {
            validGP += isg(speciesTree.getNode(i));
        }

        final int c2 = sisg(parent) + sisg(uncle);

        Node brother;
        if (Randomizer.nextBoolean()) {
            brother = parent.getLeft();
        } else {
            brother = parent.getRight();
        }

        replace(parent, brother, uncle);
        replace(grandparent, uncle, brother);

        final int validGPafter = validGP - c2 + sisg(parent) + sisg(uncle);
        double fLogHastingsRatio = Math.log((float) validGP / validGPafter);

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
                final boolean leftContainsBrother = geneTree.lineageOverlap.get(leftChildNode.getNr()).contains(brother);
                final boolean rightContainsBrother = geneTree.lineageOverlap.get(rightChildNode.getNr()).contains(brother);

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
                        } else { // the right child will be disowned and reattached directly to its grandparent
                            pruneAndRegraft(geneTreeNode, chosenNode, rightChildNode);
                        }
                    }
                }
            }
            forwardGraftCounts.add(jForwardGraftCounts);
        }

        msc.computeCoalescentTimes(); // rebuild the gene-trees-within-species-tree to account for the tree changes

        for (int j = 0; j < nGeneTrees; j++) {
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            final Map<Integer, Integer> jForwardGraftCounts = forwardGraftCounts.get(j);
            for (int geneTreeNodeNumber: jForwardGraftCounts.keySet()) {
                final double nodeHeight = geneTree.getNodeHeight(geneTreeNodeNumber);
                final List<Node> reverseGraftBranches = findGraftBranches(geneTree, brotherBranchNumber, nodeHeight);

                final int forwardGraftCount = jForwardGraftCounts.get(geneTreeNodeNumber);
                final int reverseGraftCount = reverseGraftBranches.size();

                fLogHastingsRatio += Math.log((float) forwardGraftCount / reverseGraftCount);
            }
        }

        return fLogHastingsRatio;
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
    private static void pruneAndRegraft(final Node nodeSPR, final Node newChild, final Node disownedChild) {
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
    }

    // Copied with minimal modifications from the BEASTv2.3 TreeOperator class.
    private static void replace(final Node node, final Node child, final Node replacement) {
        node.removeChild(child);
        node.addChild(replacement);
        node.makeDirty(Tree.IS_FILTHY);
        replacement.makeDirty(Tree.IS_FILTHY);
    }
}
