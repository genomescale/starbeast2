package starbeast2;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
    private int speciesTreeNodeCount;
    private int[] geneNodeSpeciesAssignment;

    protected int[] coalescentLineageCounts; // the number of lineages at the tipward end of each branch
    protected Node[] geneTreeNodesArray;
    final protected List<List<Double>> coalescentTimes = new ArrayList<List<Double>>(); // the coalescent event times for this gene tree for all species tree branches
    final protected List<Set<Integer>> lineageOverlap = new ArrayList<Set<Integer>>(); // species tree branches that overlap with each gene tree branch or descendant gene tree branches 
    final protected List<Set<Integer>> branchOverlap = new ArrayList<Set<Integer>>(); // gene tree branches that overlap with each species tree branch
    final protected List<Set<Integer>> branchNodeMap = new ArrayList<Set<Integer>>(); // gene tree nodes within each species tree branch

    public double getTreeHeight() {
        final Node geneTreeRootNode = treeInput.get().getRoot();
        final double geneTreeHeight = geneTreeRootNode.getHeight();
        
        return geneTreeHeight;
    }

    public double getNodeHeight(int nodeNumber) {
        final Node arbitraryNode = treeInput.get().getNode(nodeNumber);
        final double nodeHeight = arbitraryNode.getHeight();
        
        return nodeHeight;
    }

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
        for (int geneNodeNumber = 0; geneNodeNumber < geneTreeNodeCount; geneNodeNumber++) {
            lineageOverlap.add(new HashSet<Integer>());
        }

        geneTreeNodesArray = treeInput.get().getNodesAsArray();
        geneTreeLeafNodeCount = treeInput.get().getLeafNodeCount();
    }

    public void initCoalescentArrays(TreeInterface speciesTree) {
        speciesTreeNodeCount = speciesTree.getNodeCount();
        coalescentLineageCounts = new int[speciesTreeNodeCount];

        for (int speciesNodeNumber = 0; speciesNodeNumber < speciesTreeNodeCount; speciesNodeNumber++) {
            coalescentTimes.add(new ArrayList<Double>());
            branchOverlap.add(new HashSet<Integer>());
            branchNodeMap.add(new HashSet<Integer>());
        }
    }

    public boolean computeCoalescentTimes(TreeInterface speciesTree, HashMap<String, Integer> tipNumberMap) {
        final TreeInterface geneTree = treeInput.get();

        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
        Arrays.fill(coalescentLineageCounts, 0);
        Arrays.fill(geneNodeSpeciesAssignment, -1); // -1 means no species assignment for that gene tree node has been made yet

        // rebuild each list of per-branch coalescent times for the same reason
        for (int i = 0; i < speciesTreeNodeCount; i++) {
            coalescentTimes.get(i).clear();
            branchNodeMap.get(i).clear();
            branchOverlap.get(i).clear();
        }

        // rebuild each set of gene tree node inheriting species for the same reason
        for (int i = 0; i < geneTreeNodeCount; i++) {
            lineageOverlap.get(i).clear();
        }

        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
            final int speciesTreeLeafNumber = tipNumberMap.get(geneTreeLeafNode.getID());
            final Node speciesTreeLeafNode = speciesTree.getNode(speciesTreeLeafNumber);
            branchNodeMap.get(speciesTreeLeafNumber).add(geneTreeLeafNumber);

            coalescentLineageCounts[speciesTreeLeafNumber]++;

            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
            final int firstCoalescenceNumber = firstCoalescenceNode.getNr();

            if (!recurseCoalescenceEvents(firstCoalescenceNode, firstCoalescenceNumber, geneTreeLeafNode, geneTreeLeafNumber, speciesTreeLeafNode, speciesTreeLeafNumber)) {
                // this gene tree IS NOT compatible with the species tree
                return false;
            }
        }

        // begin sanity check
        int coalescenceCount = 0;
        for (List<Double> times: coalescentTimes) {
            coalescenceCount += times.size();
        }

        assert coalescenceCount == geneTreeNodeCount - geneTreeLeafNodeCount;
        // end sanity check

        // this gene tree IS compatible with the species tree
        return true;
    }

    private boolean recurseCoalescenceEvents(final Node geneTreeNode, final int geneTreeNodeNumber, final Node geneTreeChildNode, final int geneTreeChildNumber, final Node speciesTreeNode, final int speciesTreeNodeNumber) {
        branchOverlap.get(speciesTreeNodeNumber).add(geneTreeChildNumber);
        lineageOverlap.get(geneTreeChildNumber).add(speciesTreeNodeNumber);

        final double geneTreeNodeHeight = geneTreeNode.getHeight();
        final Node speciesTreeParentNode = speciesTreeNode.getParent();

        double speciesTreeParentHeight;
        if (speciesTreeParentNode == null) { // current node is the root node
            speciesTreeParentHeight = Double.POSITIVE_INFINITY; // there are no ancestoral branches to consider
        } else {
            speciesTreeParentHeight = speciesTreeParentNode.getHeight();
        }

        if (geneTreeNodeHeight < speciesTreeParentHeight) { // this gene coalescence occurs on the species tree current branch
            lineageOverlap.get(geneTreeNodeNumber).addAll(lineageOverlap.get(geneTreeChildNumber));
            final int existingSpeciesAssignment = geneNodeSpeciesAssignment[geneTreeNodeNumber];
            if (existingSpeciesAssignment == -1) {
                branchNodeMap.get(speciesTreeNodeNumber).add(geneTreeNodeNumber);
                geneNodeSpeciesAssignment[geneTreeNodeNumber] = speciesTreeNodeNumber;
                coalescentTimes.get(speciesTreeNodeNumber).add(geneTreeNodeHeight);
                final Node geneTreeParentNode = geneTreeNode.getParent();
                if (geneTreeParentNode == null) {
                    // this is the root of the gene tree and no incompatibilities were detected
                    return true;
                } else {
                    // if this is not the root of the gene tree, check the next coalescence
                    final int geneTreeParentNumber = geneTreeParentNode.getNr();
                    return recurseCoalescenceEvents(geneTreeParentNode, geneTreeParentNumber, geneTreeNode, geneTreeNodeNumber, speciesTreeNode, speciesTreeNodeNumber);
                }
            } else if (existingSpeciesAssignment == speciesTreeNodeNumber) {
                return true; // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
            } else {
                return false; // this gene tree IS NOT compatible with the species tree
            }
        } else { // this gene coalescence occurs on an ancestral branch
            final int speciesTreeParentNodeNumber = speciesTreeParentNode.getNr();
            coalescentLineageCounts[speciesTreeParentNodeNumber]++;

            return recurseCoalescenceEvents(geneTreeNode, geneTreeNodeNumber, geneTreeChildNode, geneTreeChildNumber, speciesTreeParentNode, speciesTreeParentNodeNumber);
        }
    }

    public Node[] getNodes() {
        return treeInput.get().getNodesAsArray();
    }

    public Node getRoot() {
        return treeInput.get().getRoot();
    }
}
