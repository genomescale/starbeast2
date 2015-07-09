package msc;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
 */

public class GeneTreeWithinSpeciesTree extends TreeDistribution {
    public Input<Double> ploidyInput = new Input<Double>("ploidy", "Ploidy (copy number) for this gene, typically a whole number or half (default is 2).", 2.0);
    public double ploidy;

    private List<Node> geneTreeLeafNodes;
    private int geneTreeNodeCount;
    private int speciesTreeNodeCount;
    private int[] geneTreeSpeciesAssignment;
    
    public int[] coalescentEventCounts; // this is typically given the pronumeral "k"
    public int[] coalescentLineageCounts; // this is typically given the pronumeral "n"
    public List<List<Double>> coalescentTimes;

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

        geneTreeLeafNodes = treeInput.get().getRoot().getAllLeafNodes();
        geneTreeNodeCount = treeInput.get().getNodeCount();
        geneTreeSpeciesAssignment = new int[geneTreeNodeCount];
        
        coalescentTimes = new ArrayList<List<Double>>();
    }

    public void initCoalescentArrays(TreeInterface speciesTree) {
        speciesTreeNodeCount = speciesTree.getNodeCount();
        coalescentEventCounts = new int[speciesTreeNodeCount];
        coalescentLineageCounts = new int[speciesTreeNodeCount];

        coalescentTimes = new ArrayList<List<Double>>();
        for (int i = 0; i < speciesTreeNodeCount; i++) {
            coalescentTimes.add(new ArrayList<Double>());
        }
    }

    public boolean computeCoalescentTimes(TreeInterface speciesTree, HashMap<String, Integer> tipNumberMap) {
        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
        Arrays.fill(coalescentEventCounts, 0);
        Arrays.fill(coalescentLineageCounts, 0);
        Arrays.fill(geneTreeSpeciesAssignment, -1); // -1 means no species assignment for that gene tree node has been made yet

        // rebuild each list of per-branch coalescent times for the same reason
        for (int i = 0; i < speciesTreeNodeCount; i++) {
            coalescentTimes.get(i).clear();
        }

        for (final Node geneTreeLeafNode: geneTreeLeafNodes) {
            final int speciesTreeLeafNumber = tipNumberMap.get(geneTreeLeafNode.getID());
            final Node speciesTreeLeafNode = speciesTree.getNode(speciesTreeLeafNumber);

            coalescentLineageCounts[speciesTreeLeafNumber]++;

            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
            final int firstCoalescenceNumber = firstCoalescenceNode.getNr();

            if (!recurseCoalescenceEvents(firstCoalescenceNode, speciesTreeLeafNode, firstCoalescenceNumber, speciesTreeLeafNumber)) {
                // this gene tree IS NOT compatible with the species tree
                return false;
            }
        }

        // this gene tree IS compatible with the species tree
        return true;
    }

    private boolean recurseCoalescenceEvents(Node geneTreeNode, Node speciesTreeNode, final int geneTreeNodeNumber, final int speciesTreeNodeNumber) {
        final double geneTreeNodeHeight = geneTreeNode.getHeight();

        final Node speciesTreeParentNode = speciesTreeNode.getParent();

        double speciesTreeParentHeight;
        if (speciesTreeParentNode == null) { // current node is the root node
            speciesTreeParentHeight = Double.POSITIVE_INFINITY; // there are no ancestoral branches to consider
        } else {
            speciesTreeParentHeight = speciesTreeParentNode.getHeight();
        }

        if (geneTreeNodeHeight < speciesTreeParentHeight) { // this gene coalescence occurs on the species tree current branch
            final int existingSpeciesAssignment = geneTreeSpeciesAssignment[geneTreeNodeNumber];
            if (existingSpeciesAssignment == -1) {
                geneTreeSpeciesAssignment[geneTreeNodeNumber] = speciesTreeNodeNumber;
                coalescentEventCounts[speciesTreeNodeNumber]++;
                coalescentTimes.get(speciesTreeNodeNumber).add(geneTreeNodeHeight);

                final Node geneTreeParentNode = geneTreeNode.getParent();
                if (geneTreeParentNode != null) {
                    // if this is not the root of the gene tree, check the next coalescence
                    final int geneTreeParentNodeNumber = geneTreeParentNode.getNr();
                    return recurseCoalescenceEvents(geneTreeParentNode, speciesTreeNode, geneTreeParentNodeNumber, speciesTreeNodeNumber);
                } else {
                    return true; // this is the root of the gene tree and no incompatibilities were detected
                }
            } else if (existingSpeciesAssignment == speciesTreeNodeNumber) {
                return true; // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
            } else {
                return false; // this gene tree IS NOT compatible with the species tree
            }
        } else { // this gene coalescence occurs on an ancestral branch
            final int speciesTreeParentNodeNumber = speciesTreeParentNode.getNr();
            coalescentLineageCounts[speciesTreeParentNodeNumber]++;

            return recurseCoalescenceEvents(geneTreeNode, speciesTreeParentNode, geneTreeNodeNumber, speciesTreeParentNodeNumber);
        }
    }
}
