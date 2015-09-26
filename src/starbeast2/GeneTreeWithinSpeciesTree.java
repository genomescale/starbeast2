package starbeast2;

import java.util.Arrays;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
 */

public class GeneTreeWithinSpeciesTree extends TreeDistribution {
    public Input<Double> ploidyInput = new Input<Double>("ploidy", "Ploidy (copy number) for this gene, typically a whole number or half (default is 2).", 2.0);
    protected double ploidy;

    private int geneTreeLeafNodeCount;
    private int geneTreeNodeCount;

    final protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create(); // the coalescent event times for this gene tree for all species tree branches
    final protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create(); // the number of lineages at the tipward end of each branch

    protected int[] geneNodeSpeciesAssignment;

    public void initAndValidate() throws Exception {
        if (treeInput.get() == null) {
            final String thisID = getID();
            if (thisID == null) {
                throw new Exception(String.format("A gene tree input named \"tree\" is required but was not supplied."));
            } else {
                throw new Exception(String.format("A gene tree input named \"tree\" is required but was not supplied for \"%s\".", thisID));
            }
        }

        ploidy = ploidyInput.get();

        geneTreeNodeCount = treeInput.get().getNodeCount();
        geneNodeSpeciesAssignment = new int[geneTreeNodeCount];

        geneTreeLeafNodeCount = treeInput.get().getLeafNodeCount();
    }

    protected boolean computeCoalescentTimes(TreeInterface speciesTree, Map<String, Integer> tipNumberMap) {
        final TreeInterface geneTree = treeInput.get();

        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
        Arrays.fill(geneNodeSpeciesAssignment, -1); // -1 means no species assignment for that gene tree node has been made yet

        coalescentLineageCounts.clear();
        coalescentTimes.clear();

        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
            final int speciesTreeLeafNumber = tipNumberMap.get(geneTreeLeafNode.getID());
            final Node speciesTreeLeafNode = speciesTree.getNode(speciesTreeLeafNumber);
            coalescentLineageCounts.add(speciesTreeLeafNumber);

            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
            final int firstCoalescenceNumber = firstCoalescenceNode.getNr();

            if (!recurseCoalescenceEvents(firstCoalescenceNode, firstCoalescenceNumber, speciesTreeLeafNode, speciesTreeLeafNumber)) {
                // this gene tree IS NOT compatible with the species tree
                return false;
            }
        }

        // the number of coalescent times should equal the number of internal gene tree nodes (each a coalescent event)
        assert coalescentTimes.size() == geneTreeNodeCount - geneTreeLeafNodeCount;

        // this gene tree IS compatible with the species tree
        return true;
    }

    private boolean recurseCoalescenceEvents(final Node geneTreeNode, final int geneTreeNodeNumber, final Node speciesTreeNode, final int speciesTreeNodeNumber) {
        final double geneTreeNodeHeight = geneTreeNode.getHeight();

        // check if the next coalescence event occurs in an ancestral branch
        if (!speciesTreeNode.isRoot()) {
            final Node speciesTreeParentNode = speciesTreeNode.getParent();
            if (geneTreeNodeHeight >= speciesTreeParentNode.getHeight()) {
                final int speciesTreeParentNodeNumber = speciesTreeParentNode.getNr();
                coalescentLineageCounts.add(speciesTreeParentNodeNumber);

                return recurseCoalescenceEvents(geneTreeNode, geneTreeNodeNumber, speciesTreeParentNode, speciesTreeParentNodeNumber);
            }
        }

        // this code executes if the next coalescence event occurs within the current branch
        final int existingSpeciesAssignment = geneNodeSpeciesAssignment[geneTreeNodeNumber];
        if (existingSpeciesAssignment == -1) {
            geneNodeSpeciesAssignment[geneTreeNodeNumber] = speciesTreeNodeNumber;
            coalescentTimes.put(speciesTreeNodeNumber, geneTreeNodeHeight);
            final Node geneTreeParentNode = geneTreeNode.getParent();
            if (geneTreeParentNode == null) {
                // this is the root of the gene tree and no incompatibilities were detected
                return true;
            } else {
                // if this is not the root of the gene tree, check the next coalescence
                final int geneTreeParentNumber = geneTreeParentNode.getNr();
                return recurseCoalescenceEvents(geneTreeParentNode, geneTreeParentNumber, speciesTreeNode, speciesTreeNodeNumber);
            }
        } else if (existingSpeciesAssignment == speciesTreeNodeNumber) {
            return true; // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
        } else {
            return false; // this gene tree IS NOT compatible with the species tree
        }
    }

    // identify nodes to be moved as part of a coordinated exchange move
    protected boolean findMovedChildren(Node geneTreeNode, SortedMap<Node, Node> movedNodes, Set<String> brotherDescendants, double lowerHeight, double upperHeight) {
        if (geneTreeNode.isLeaf()) {
            final String descendantName = geneTreeNode.getID();
            return brotherDescendants.contains(descendantName);
        }

        final Node leftChild = geneTreeNode.getLeft();
        final Node rightChild = geneTreeNode.getRight();

        final boolean leftOverlapsBrother = findMovedChildren(leftChild, movedNodes, brotherDescendants, lowerHeight, upperHeight);
        final boolean rightOverlapsBrother = findMovedChildren(rightChild, movedNodes, brotherDescendants, lowerHeight, upperHeight);

        final double nodeHeight = geneTreeNode.getHeight();
        if (nodeHeight >= lowerHeight && nodeHeight < upperHeight) {
            if (leftOverlapsBrother && !rightOverlapsBrother) {
                movedNodes.put(geneTreeNode, leftChild);
            } else if (!leftOverlapsBrother && rightOverlapsBrother) {
                movedNodes.put(geneTreeNode, rightChild);
            }
        }

        return leftOverlapsBrother || rightOverlapsBrother;
    }

    // identify nodes that can serve as graft branches as part of a coordinated exchange move
    protected boolean findGraftBranches(Node geneTreeNode, Set<Node> graftNodes, Set<String> branchDescendants) {
        if (geneTreeNode.isLeaf()) {
            final String descendantName = geneTreeNode.getID();
            return branchDescendants.contains(descendantName);
        }

        final Node leftChild = geneTreeNode.getLeft();
        final Node rightChild = geneTreeNode.getRight();
        final boolean leftOverlaps = findGraftBranches(leftChild, graftNodes, branchDescendants);
        final boolean rightOverlaps = findGraftBranches(rightChild, graftNodes, branchDescendants);

        // subtree defined by a child node overlaps species subtree defined by branch
        if (leftOverlaps || rightOverlaps) {
            if (leftOverlaps) {
                graftNodes.add(leftChild);
            }

            if (rightOverlaps) {
                graftNodes.add(rightChild);
            }

            return true;
        }

        return false;
    }

    protected Node getRoot() {
        return treeInput.get().getRoot();
    }

    protected double getTreeHeight() {
        return treeInput.get().getRoot().getHeight();
    }
}
