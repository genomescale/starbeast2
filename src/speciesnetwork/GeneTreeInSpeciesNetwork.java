package speciesnetwork;

import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.List;

import beast.core.parameter.IntegerParameterList;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.HashMultiset;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameterList;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import speciesnetwork.SpeciesNetwork.traversal;

/**
* @author Huw Ogilvie
 */

public class GeneTreeInSpeciesNetwork extends CalculationNode {
    public Input<SpeciesNetwork> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network for embedding the gene tree.", Validate.REQUIRED);
    public Input<Tree> geneTreeInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", Validate.REQUIRED);
    public Input<IntegerParameterList> embeddingInput =
            new Input<>("embedding", "Map of gene tree traversal within the species network.", Validate.REQUIRED);
    public Input<Double> ploidyInput =
            new Input<>("ploidy", "Ploidy (copy number) for this gene (default is 2).", 2.0);
    protected double ploidy;

    private traversal[][] traversalMatrix;

    private int geneTreeLeafNodeCount;
    private int geneTreeNodeCount;
    private int speciesLeafNodeCount;
    private int reticulationNodeOffset;
    private int speciesBranchCount;
    private boolean needsUpdate;

    protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create(); // the coalescent event times for this gene tree for all species tree branches
    protected ListMultimap<Integer, Double> storedCoalescentTimes = ArrayListMultimap.create(); // the coalescent event times for this gene tree for all species tree branches
    protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create(); // the number of lineages at the tipward end of each branch
    protected Multiset<Integer> storedCoalescentLineageCounts = HashMultiset.create(); // the number of lineages at the tipward end of each branch

    protected int[] geneNodeBranchAssignment;
    protected int[] storedGeneNodeBranchAssignment;
    protected double[][] speciesOccupancy;
    protected double[][] storedSpeciesOccupancy;
    protected boolean geneTreeCompatible;
    protected boolean storedGeneTreeCompatible;

    /**
     * 2d array list of mapping gene tree into the species network
     */
    protected IntegerParameterList treeMappingToNetworkList;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = embeddingInput.isDirty() || geneTreeInput.isDirty() || speciesNetworkInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
        storedCoalescentTimes.clear();
        storedCoalescentLineageCounts.clear();

        storedCoalescentTimes.putAll(coalescentTimes);
        storedCoalescentLineageCounts.addAll(coalescentLineageCounts);

        storedSpeciesOccupancy = new double[speciesOccupancy.length][speciesOccupancy[0].length];
        System.arraycopy(geneNodeBranchAssignment, 0, storedGeneNodeBranchAssignment, 0, geneNodeBranchAssignment.length);
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);

        storedGeneTreeCompatible = geneTreeCompatible;

        super.store();
    }

    @Override
    public void restore() {
        ListMultimap<Integer, Double> tmpCoalescentTimes = coalescentTimes;
        Multiset<Integer> tmpCoalescentLineageCounts = coalescentLineageCounts;
        int[] tmpGeneNodeSpeciesAssignment = geneNodeBranchAssignment;
        double[][] tmpSpeciesOccupancy = speciesOccupancy;
        boolean tmpGeneTreeCompatible = geneTreeCompatible;

        coalescentTimes = storedCoalescentTimes;
        coalescentLineageCounts = storedCoalescentLineageCounts;
        speciesOccupancy = storedSpeciesOccupancy;
        geneNodeBranchAssignment = storedGeneNodeBranchAssignment;
        geneTreeCompatible = storedGeneTreeCompatible;

        storedCoalescentTimes = tmpCoalescentTimes;
        storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
        storedSpeciesOccupancy = tmpSpeciesOccupancy;
        storedGeneNodeBranchAssignment = tmpGeneNodeSpeciesAssignment;
        storedGeneTreeCompatible = tmpGeneTreeCompatible;

        super.restore();
    }

    public void initAndValidate() throws Exception {
        ploidy = ploidyInput.get();

        geneTreeNodeCount = geneTreeInput.get().getNodeCount();
        geneNodeBranchAssignment = new int[geneTreeNodeCount];
        storedGeneNodeBranchAssignment = new int[geneTreeNodeCount];

        geneTreeLeafNodeCount = geneTreeInput.get().getLeafNodeCount();

        // generate map of species tree tip node names to node numbers
        final SpeciesNetwork speciesNetwork = speciesNetworkInput.get();
        final HashMap<String, Integer> speciesNumberMap = new HashMap<>();

        for (NetworkNode leafNode: speciesNetwork.getNetwork().getLeafNodes()) {
            final String speciesName = leafNode.getID();
            final int speciesNumber = leafNode.getNr();

            speciesNumberMap.put(speciesName, speciesNumber);
        }

        geneTreeCompatible = false;
        storedGeneTreeCompatible = false;
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
        final Network speciesNetwork = speciesNetworkInput.get().getNetwork();
        final TreeInterface geneTree = geneTreeInput.get();
        final Map<String, Integer> tipNumberMap = speciesNetworkInput.get().getTipNumberMap();

        int speciesNodeCount = speciesNetwork.getNodeCount();
        int reticulationNodeCount = speciesNetwork.getReticulationNodeCount();
        speciesLeafNodeCount = speciesNetwork.getLeafNodeCount();
        // equivalent to the number of leaf nodes plus number of internal speciation nodes
        reticulationNodeOffset = speciesNodeCount - (reticulationNodeCount + 1);
        // each reticulation node has two branches
        speciesBranchCount = speciesNodeCount + reticulationNodeCount;
        // traversalNodeCount excludes the leaves and the root
        final int traversalNodeCount = speciesNodeCount - (speciesLeafNodeCount + 1);
        final int geneTreeNodeCount = speciesNetwork.getNodeCount();

        traversalMatrix = new traversal[geneTreeNodeCount - 1][traversalNodeCount];
        speciesOccupancy = new double[geneTreeNodeCount][speciesBranchCount];

        final IntegerParameterList embedding = embeddingInput.get();
        for (int i = 0; i < traversalNodeCount - 1; i++) {
            for (int j = 0; j < geneTreeNodeCount - 1; j++) {
                switch (embedding.get(i).getValue(j)) {
                case -1:
                    traversalMatrix[j][i] = traversal.NEITHER;
                    break;
                case 0:
                    traversalMatrix[j][i] = traversal.LEFT;
                    break;
                case 1:
                    traversalMatrix[j][i] = traversal.RIGHT;
                    break;
                }
            }
        }

        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
        Arrays.fill(geneNodeBranchAssignment, -1);
        // -1 means no species assignment for that gene tree node has been made yet

        coalescentLineageCounts.clear();
        coalescentTimes.clear();

        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
            final int speciesLeafNodeNumber = tipNumberMap.get(geneTreeLeafNode.getID());
            final NetworkNode speciesLeafNode = speciesNetwork.getNode(speciesLeafNodeNumber);
            final traversal speciesOrientation = speciesLeafNode.getOrientation();
            coalescentLineageCounts.add(speciesLeafNodeNumber);

            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
            final double lastHeight = 0.0;
            if (!recurseCoalescenceEvents(geneTreeLeafNumber, lastHeight, firstCoalescenceNode,
                                          speciesLeafNode, speciesLeafNodeNumber, speciesOrientation)) {
                // this gene tree IS NOT compatible with the species tree
                geneTreeCompatible = false;
                needsUpdate = false;
                return;
            }
        }

        geneTreeCompatible = true;
        needsUpdate = false;
    }

    private boolean recurseCoalescenceEvents(final int lastGeneTreeNodeNumber, final double lastHeight, final Node geneTreeNode,
                                             final NetworkNode speciesNode, final int speciesBranchNumber, final traversal orientation) {
        // check if the next coalescence event occurs in an ancestral branch
        if (!speciesNode.isRoot()) {
            final NetworkNode speciesParentNode = (orientation == traversal.LEFT) ? speciesNode.getLeftParent() : speciesNode.getRightParent();
            final int speciesParentNodeNumber = speciesParentNode.getNr();
            final traversal nextOrientation = traversalMatrix[lastGeneTreeNodeNumber][speciesParentNodeNumber - speciesLeafNodeCount];
            if (nextOrientation != traversal.NEITHER) { // this gene lineage traverses through the species parent node (left or right)
                final double speciesParentNodeHeight = speciesParentNode.getHeight();
                int speciesParentBranchNumber;
                if (speciesParentNode.isRoot()) {
                    speciesParentBranchNumber = speciesBranchCount - 1;
                } else if (speciesParentNode.isReticulation()) {
                    final int reticulationNumber = speciesParentNodeNumber - reticulationNodeOffset;
                    speciesParentBranchNumber = reticulationNodeOffset + (reticulationNumber * 2);
                    if (orientation == traversal.RIGHT) speciesParentBranchNumber++;
                } else { // leaf or internal speciation node
                    speciesParentBranchNumber = speciesParentNodeNumber;
                }
                speciesOccupancy[lastGeneTreeNodeNumber][speciesBranchNumber] = speciesParentNodeHeight - lastHeight;
                coalescentLineageCounts.add(speciesParentBranchNumber);
                return recurseCoalescenceEvents(lastGeneTreeNodeNumber, speciesParentNodeHeight, geneTreeNode,
                        speciesParentNode, speciesParentBranchNumber, nextOrientation);
            }
        }

        // this code executes if the next coalescence event occurs within the current branch
        final double geneTreeNodeHeight = geneTreeNode.getHeight();
        final int geneTreeNodeNumber = geneTreeNode.getNr();
        final int existingBranchAssignment = geneNodeBranchAssignment[geneTreeNodeNumber];
        speciesOccupancy[lastGeneTreeNodeNumber][speciesBranchNumber] = geneTreeNodeHeight - lastHeight;
        if (existingBranchAssignment == -1) {
            geneNodeBranchAssignment[geneTreeNodeNumber] = speciesBranchNumber;
            coalescentTimes.put(speciesBranchNumber, geneTreeNodeHeight);
            final Node nextGeneTreeNode = geneTreeNode.getParent();
            if (nextGeneTreeNode == null) {
                // this is the root of the gene tree and no incompatibilities were detected
                return true;
            } else {
                // if this is not the root of the gene tree, check the subsequent (back in time) coalescence event
                return recurseCoalescenceEvents(geneTreeNodeNumber, geneTreeNodeHeight, nextGeneTreeNode,
                                                speciesNode, speciesBranchNumber, orientation);
            }
        } else {
            // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
            // return false if this gene tree IS NOT compatible with the species network
            return existingBranchAssignment == speciesBranchNumber;
        }
    }

    public double[][] getSpeciesOccupancy() {
        if (needsUpdate) update();

        return speciesOccupancy;
    }

    protected IntegerParameterList getTreeMappingToNetwork() {
        return treeMappingToNetworkList;
    }

    protected Tree getTree() {
        return geneTreeInput.get();
    }

    protected Node getRoot() {
        return geneTreeInput.get().getRoot();
    }

    protected double getTreeHeight() {
        return geneTreeInput.get().getRoot().getHeight();
    }

    /**
     * @return the first tip node which is descendant of
     * @param gTreeNode
     * this can be in Tree.java as gTreeNode.getGeneTreeTipDescendant()
     */
    public Node getGeneNodeDescendantTip(Node gTreeNode) {
        final TreeInterface geneTree = geneTreeInput.get();
        final List<Node> gTreeTips = geneTree.getExternalNodes();  // tips
        for (Node tip : gTreeTips) {
            Node node = tip;
            while(node != null && !node.equals(gTreeNode)) {
                node = node.getParent();
            }
            if (node != null)
                return tip;  // find you!
        }
        return null;  // looped all the tips but nothing found
    }
}
