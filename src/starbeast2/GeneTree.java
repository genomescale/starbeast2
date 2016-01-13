package starbeast2;

import java.util.Arrays;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.HashMultiset;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
 */

public class GeneTree extends CalculationNode {
    public Input<SpeciesNetwork> speciesNetworkInput = new Input<>("speciesNetwork", "Species network for embedding the gene tree.", Validate.REQUIRED);
    public Input<Tree> geneTreeInput = new Input<>("geneTree", "Gene tree embedded in the Species network.", Validate.REQUIRED);
    public Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for this gene (default is 2).", 2.0);
    protected double ploidy;

    private int geneTreeLeafNodeCount;
    private int geneTreeNodeCount;
    private double[][] speciesOccupancy;
    private boolean needsUpdate;

    // the coalescent event times for this gene tree for all species network branches
    final protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create();
    final protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create(); // the number of lineages at the tipward end of each branch

    protected int[] geneNodeSpeciesAssignment;
    protected boolean geneTreeCompatible;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = true;
        return true;
    }

    public void restore() {
        needsUpdate = true;
        super.restore();
    }

    public void initAndValidate() throws Exception {
        ploidy = ploidyInput.get();

        geneTreeLeafNodeCount = geneTreeInput.get().getLeafNodeCount();
        geneTreeNodeCount = geneTreeInput.get().getNodeCount();
        geneNodeSpeciesAssignment = new int[geneTreeNodeCount];

        needsUpdate = true;
    }

    protected boolean computeCoalescentTimes() {
        if (needsUpdate) {
            update();
        }

        // the number of coalescent times should equal the number of internal gene tree nodes (each a coalescent event)
        if (geneTreeCompatible) {
            assert coalescentTimes.size() == geneTreeNodeCount - geneTreeLeafNodeCount;
            // this gene tree IS compatible with the species tree
            return true;
        } else {
            return false;
        }
    }

    void update() {
        final SpeciesNetwork mapping = speciesNetworkInput.get();
        final Network speciesNetwork = mapping.getNetwork();
        final Map<String, Integer> tipNumberMap = mapping.getTipNumberMap();

        final int speciesNetworkNodeCount = speciesNetwork.getNodeCount();
        speciesOccupancy = new double[geneTreeNodeCount][2*speciesNetworkNodeCount];

        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
        Arrays.fill(geneNodeSpeciesAssignment, -1); // -1 means no species assignment for that gene tree node has been made yet

        coalescentLineageCounts.clear();
        coalescentTimes.clear();

        final TreeInterface geneTree = geneTreeInput.get();
        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
            final int speciesNetworkLeafNumber = tipNumberMap.get(geneTreeLeafNode.getID());
            final NetworkNode speciesNetworkLeafNode = speciesNetwork.getNode(speciesNetworkLeafNumber);
            coalescentLineageCounts.add(2 * speciesNetworkLeafNumber);

            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
            final int firstCoalescenceNumber = firstCoalescenceNode.getNr();
            final double lastHeight = 0.0;

            if (!recurseCoalescenceEvents(geneTreeLeafNumber, lastHeight, firstCoalescenceNode, firstCoalescenceNumber, speciesNetworkLeafNode, 2*speciesNetworkLeafNumber)) {
                // this gene tree IS NOT compatible with the species tree
                geneTreeCompatible = false;
                needsUpdate = false;
                return;
            }
        }

        geneTreeCompatible = true;
        needsUpdate = false;
    }

    private boolean recurseCoalescenceEvents(final int lastGeneTreeNodeNumber, final double lastHeight, final Node geneTreeNode, final int geneTreeNodeNumber, final NetworkNode speciesNetworkNode, final int speciesNetworkNodeNumber) {
        final double geneTreeNodeHeight = geneTreeNode.getHeight();

        // check if the next coalescence event occurs in an ancestral branch
        if (!speciesNetworkNode.isRoot()) {
            final NetworkNode speciesNetworkParentNode = speciesNetworkNode.getParent(); // (geneTreeNodeNumber);
            final double speciesNetworkParentHeight = speciesNetworkParentNode.getHeight();
            if (geneTreeNodeHeight >= speciesNetworkParentHeight) {
                speciesOccupancy[lastGeneTreeNodeNumber][speciesNetworkNodeNumber] = speciesNetworkParentHeight - lastHeight;
                final int speciesNetworkParentNodeNumber = 2 * speciesNetworkParentNode.getNr(); // + speciesNetworkParentNode.getOffset();
                coalescentLineageCounts.add(speciesNetworkParentNodeNumber);

                return recurseCoalescenceEvents(lastGeneTreeNodeNumber, speciesNetworkParentHeight, geneTreeNode, geneTreeNodeNumber, speciesNetworkParentNode, speciesNetworkParentNodeNumber);
            }
        }

        // this code executes if the next coalescence event occurs within the current branch
        speciesOccupancy[lastGeneTreeNodeNumber][speciesNetworkNodeNumber] = geneTreeNodeHeight - lastHeight;
        final int existingSpeciesAssignment = geneNodeSpeciesAssignment[geneTreeNodeNumber];
        if (existingSpeciesAssignment == -1) {
            geneNodeSpeciesAssignment[geneTreeNodeNumber] = speciesNetworkNodeNumber;
            coalescentTimes.put(speciesNetworkNodeNumber, geneTreeNodeHeight);
            final Node nextGeneTreeNode = geneTreeNode.getParent();
            if (nextGeneTreeNode == null) {
                // this is the root of the gene tree and no incompatibilities were detected
                return true;
            } else {
                // if this is not the root of the gene tree, check the subsequent (back in time) coalescence event
                final int nextGeneTreeNodeNumber = nextGeneTreeNode.getNr();
                return recurseCoalescenceEvents(geneTreeNodeNumber, geneTreeNodeHeight, nextGeneTreeNode, nextGeneTreeNodeNumber, speciesNetworkNode, speciesNetworkNodeNumber);
            }
        } else {
            // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
            // otherwise, this gene tree IS NOT compatible with the species tree
            return (existingSpeciesAssignment == speciesNetworkNodeNumber);
        }
    }

    public double[] getOccupancy(Node node) {
        if (needsUpdate) {
            update();
        }

        final int geneTreeNodeNumber = node.getNr();
        return speciesOccupancy[geneTreeNodeNumber];
    }
}
