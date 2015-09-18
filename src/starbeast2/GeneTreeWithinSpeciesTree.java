package starbeast2;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.SetMultimap;

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
    final protected ListMultimap<Integer, Node> branchNodeMap = ArrayListMultimap.create(); // gene tree nodes within each species tree branch
    final protected SetMultimap<Integer, Node> speciesGeneDescent = HashMultimap.create(); // gene tree nodes within each species tree branch and descendants

    protected int[] coalescentLineageCounts; // the number of lineages at the tipward end of each branch
    protected int[] geneNodeSpeciesAssignment;

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

        geneTreeLeafNodeCount = treeInput.get().getLeafNodeCount();
    }

    public boolean computeCoalescentTimes(TreeInterface speciesTree, HashMap<String, Integer> tipNumberMap) {
        final TreeInterface geneTree = treeInput.get();

        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
        coalescentLineageCounts = new int[speciesTree.getNodeCount()];
        Arrays.fill(geneNodeSpeciesAssignment, -1); // -1 means no species assignment for that gene tree node has been made yet

        coalescentTimes.clear();
        branchNodeMap.clear();
        speciesGeneDescent.clear();

        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
            final int speciesTreeLeafNumber = tipNumberMap.get(geneTreeLeafNode.getID());
            final Node speciesTreeLeafNode = speciesTree.getNode(speciesTreeLeafNumber);
            branchNodeMap.put(speciesTreeLeafNumber, geneTreeLeafNode);

            coalescentLineageCounts[speciesTreeLeafNumber]++;
            speciesGeneDescent.put(speciesTreeLeafNumber, geneTreeLeafNode);

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
        final Node speciesTreeParentNode = speciesTreeNode.getParent();

        double speciesTreeParentHeight;
        if (speciesTreeParentNode == null) { // current node is the root node
            speciesTreeParentHeight = Double.POSITIVE_INFINITY; // there are no ancestral branches to consider
        } else {
            speciesTreeParentHeight = speciesTreeParentNode.getHeight();
        }

        if (geneTreeNodeHeight < speciesTreeParentHeight) { // this gene coalescence occurs on the species tree current branch
            speciesGeneDescent.put(speciesTreeNodeNumber, geneTreeNode);

            final int existingSpeciesAssignment = geneNodeSpeciesAssignment[geneTreeNodeNumber];
            if (existingSpeciesAssignment == -1) {
                branchNodeMap.put(speciesTreeNodeNumber, geneTreeNode);
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
        } else { // this gene coalescence occurs on an ancestral branch
            final int speciesTreeParentNodeNumber = speciesTreeParentNode.getNr();
            final Set<Node> descendantGeneNodes = speciesGeneDescent.get(speciesTreeNodeNumber);
            speciesGeneDescent.putAll(speciesTreeParentNodeNumber, descendantGeneNodes);
            coalescentLineageCounts[speciesTreeParentNodeNumber]++;

            return recurseCoalescenceEvents(geneTreeNode, geneTreeNodeNumber, speciesTreeParentNode, speciesTreeParentNodeNumber);
        }
    }

    protected void subtreeOverlap(final Node geneTreeNode, final int speciesTreeNodeNumber) {
        
    }

    public Node[] getNodes() {
        final TreeInterface geneTree = treeInput.get();
        return geneTree.getNodesAsArray();
    }

    public Node getRoot() {
        return treeInput.get().getRoot();
    }
}
