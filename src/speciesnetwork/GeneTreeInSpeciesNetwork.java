package speciesnetwork;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.HashMultiset;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
* @author Chi Zhang
 */

public class GeneTreeInSpeciesNetwork extends CalculationNode {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network for embedding the gene tree.", Validate.REQUIRED);
    public Input<Tree> geneTreeInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", Validate.REQUIRED);
    public Input<IntegerParameter> embeddingInput =
            new Input<>("embedding", "Map of gene tree traversal within the species network.", Validate.REQUIRED);
    public Input<Double> ploidyInput =
            new Input<>("ploidy", "Ploidy (copy number) for this gene (default is 2).", 2.0);
    protected double ploidy;

    private int geneTreeNodeCount;
    private int speciesLeafNodeCount;
    private int speciesBranchCount;
    private boolean needsUpdate;
    private IntegerParameter embedding;

    protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create(); // the coalescent event times for this gene tree for all species tree branches
    protected ListMultimap<Integer, Double> storedCoalescentTimes = ArrayListMultimap.create(); // the coalescent event times for this gene tree for all species tree branches
    protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create(); // the number of lineages at the tipward end of each branch
    protected Multiset<Integer> storedCoalescentLineageCounts = HashMultiset.create(); // the number of lineages at the tipward end of each branch

    protected double[][] speciesOccupancy;
    protected double[][] storedSpeciesOccupancy;

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
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);

        super.store();
    }

    @Override
    public void restore() {
        ListMultimap<Integer, Double> tmpCoalescentTimes = coalescentTimes;
        Multiset<Integer> tmpCoalescentLineageCounts = coalescentLineageCounts;
        double[][] tmpSpeciesOccupancy = speciesOccupancy;

        coalescentTimes = storedCoalescentTimes;
        coalescentLineageCounts = storedCoalescentLineageCounts;
        speciesOccupancy = storedSpeciesOccupancy;

        storedCoalescentTimes = tmpCoalescentTimes;
        storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
        storedSpeciesOccupancy = tmpSpeciesOccupancy;

        super.restore();
    }

    public void initAndValidate() throws Exception {
        ploidy = ploidyInput.get();
        geneTreeNodeCount = geneTreeInput.get().getNodeCount();
        needsUpdate = true;
    }

    protected void computeCoalescentTimes() {
        if (needsUpdate) {
            update();
        }
    }

    void update() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final TreeInterface geneTree = geneTreeInput.get();
        embedding = embeddingInput.get();

        speciesLeafNodeCount = speciesNetwork.getLeafNodeCount();
        // each reticulation node has two branches
        speciesBranchCount = speciesNetwork.getBranchCount();
        speciesOccupancy = new double[geneTreeNodeCount][speciesBranchCount];

        // reset coalescent arrays as these values need to be recomputed after any changes to the species or gene tree
        coalescentLineageCounts.clear();
        coalescentTimes.clear();

        final Node geneTreeRoot = geneTree.getRoot();
        final NetworkNode speciesNetworkRoot = speciesNetwork.getRoot();
        final int speciesRootBranchNumber = speciesBranchCount - 1;
        recurseCoalescentEvents(geneTreeRoot, speciesNetworkRoot, speciesRootBranchNumber, Double.POSITIVE_INFINITY);
        needsUpdate = false;
    }

    // forward in time recursion, unlike StarBEAST 2
    private void recurseCoalescentEvents(final Node geneTreeNode, final NetworkNode speciesNetworkNode, final int speciesBranchNumber, final double lastHeight) {
        final double geneNodeHeight = geneTreeNode.getHeight();
        final double speciesNodeHeight = speciesNetworkNode.getHeight();
        final int geneTreeNodeNumber = geneTreeNode.getNr();

        // check if coalescence node occurs in a descendant species network branch
        if (geneNodeHeight < speciesNodeHeight) {
            speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight - speciesNodeHeight;
            coalescentLineageCounts.add(speciesBranchNumber);
            final int traversalNodeNumber = speciesNetworkNode.getNr() - speciesLeafNodeCount;
            final int traversalDirection = embedding.getMatrixValue(traversalNodeNumber, geneTreeNodeNumber);
            final NetworkNode nextSpeciesNode = (traversalDirection == 0) ? speciesNetworkNode.getLeftChild() : speciesNetworkNode.getRightChild();
            final int nextSpeciesBranchNumber = (traversalDirection == 0) ? nextSpeciesNode.getLeftBranchNumber() : nextSpeciesNode.getRightBranchNumber();
            recurseCoalescentEvents(geneTreeNode, nextSpeciesNode, nextSpeciesBranchNumber, speciesNodeHeight);
        } else if (geneTreeNode.isLeaf()) { // assumes tip node heights are always zero
            speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight;
            coalescentLineageCounts.add(speciesBranchNumber);
        } else {
            speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight - geneNodeHeight;
            coalescentTimes.put(speciesBranchNumber, geneNodeHeight);
            final Node leftChild = geneTreeNode.getLeft();
            final Node rightChild = geneTreeNode.getRight();
            recurseCoalescentEvents(leftChild, speciesNetworkNode, speciesBranchNumber, geneNodeHeight);
            recurseCoalescentEvents(rightChild, speciesNetworkNode, speciesBranchNumber, geneNodeHeight);
        }
    }

    public double[][] getSpeciesOccupancy() throws Exception {
        if (needsUpdate) update();

        return speciesOccupancy;
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
