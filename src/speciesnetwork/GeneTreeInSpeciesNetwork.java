package speciesnetwork;

import java.util.Comparator;
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

    private boolean needsUpdate;
    private IntegerParameter embedding;
    private int speciesLeafNodeCount;

    // the coalescent times of this gene tree for all species branches
    protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create();
    protected ListMultimap<Integer, Double> storedCoalescentTimes = ArrayListMultimap.create();
    // the number of lineages at the tipward end of each species branch
    protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create();
    protected Multiset<Integer> storedCoalescentLineageCounts = HashMultiset.create();

    protected double[][] speciesOccupancy;
    protected double[][] storedSpeciesOccupancy;
    protected double logGammaSum;
    protected double storedLogGammaSum;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = embeddingInput.isDirty() ||
                      geneTreeInput.isDirty() || speciesNetworkInput.isDirty();
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
        
        storedLogGammaSum = logGammaSum;

        super.store();
    }

    @Override
    public void restore() {
        ListMultimap<Integer, Double> tmpCoalescentTimes = coalescentTimes;
        Multiset<Integer> tmpCoalescentLineageCounts = coalescentLineageCounts;
        double[][] tmpSpeciesOccupancy = speciesOccupancy;
        double tmpLogGammaSum = logGammaSum;

        coalescentTimes = storedCoalescentTimes;
        coalescentLineageCounts = storedCoalescentLineageCounts;
        speciesOccupancy = storedSpeciesOccupancy;
        logGammaSum = storedLogGammaSum;

        storedCoalescentTimes = tmpCoalescentTimes;
        storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
        storedSpeciesOccupancy = tmpSpeciesOccupancy;
        storedLogGammaSum = tmpLogGammaSum;

        super.restore();
    }

    public void initAndValidate() {
        ploidy = ploidyInput.get();
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
        logGammaSum = 0.0;

        final int geneTreeNodeCount = geneTree.getNodeCount();
        final int speciesBranchCount = (speciesNetwork.getNodeCount() * 2) - 1;
        speciesOccupancy = new double[geneTreeNodeCount][speciesBranchCount];
        speciesLeafNodeCount = speciesNetwork.getLeafNodeCount();

        // reset coalescent arrays as these values need to be recomputed after any changes to the species or gene tree
        coalescentLineageCounts.clear();
        coalescentTimes.clear();

        final Node geneTreeRoot = geneTree.getRoot();
        final NetworkNode speciesNetworkRoot = speciesNetwork.getRoot();
        final int speciesRootBranchNumber = speciesBranchCount - 1;
        try {
            recurseCoalescentEvents(geneTreeRoot, speciesNetworkRoot, speciesRootBranchNumber, Double.POSITIVE_INFINITY);
        } catch (Exception e) {
            e.printStackTrace();
        }
        needsUpdate = false;
    }

    // forward in time recursion, unlike StarBEAST 2
    private void recurseCoalescentEvents(final Node geneTreeNode, final NetworkNode speciesNetworkNode,
                                         final int speciesBranchNumber, final double lastHeight) throws Exception {
        final double geneNodeHeight = geneTreeNode.getHeight();
        final double speciesNodeHeight = speciesNetworkNode.getHeight();
        final int geneTreeNodeNumber = geneTreeNode.getNr();

        // check if coalescent node occurs in a descendant species network branch
        if (geneNodeHeight < speciesNodeHeight) {
            speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight - speciesNodeHeight;
            coalescentLineageCounts.add(speciesBranchNumber);
            final int traversalNodeNumber = speciesNetworkNode.getNr() - speciesLeafNodeCount;
            if (speciesNetworkNode.isReticulation()) {
                final double gammaP = speciesNetworkNode.inheritProb;
                if (speciesBranchNumber % 2 == 0) { // determine traversal direction backward in time
                    logGammaSum += Math.log(gammaP);
                } else {
                    logGammaSum += Math.log(1.0 - gammaP);
                }
            }
            // traversal direction forward in time
            final int nextSpeciesBranchNumber = embedding.getMatrixValue(traversalNodeNumber, geneTreeNodeNumber);
            final NetworkNode nextSpeciesNode = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNumber);
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

    /* unambiguous, so only identical nodes are considered equal
     * if height is equal, consider distance from root (depth)
     * if depth is equal, consider assigned node number
     */
    final class NodeHeightComparator implements Comparator<Node> {
        final int lessThan = -1;
        final int greaterThan = 1;

        @Override
        public int compare(Node nodeA, Node nodeB) {
            final double heightA = nodeA.getHeight();
            final double heightB = nodeB.getHeight();
            if (heightA == heightB) {
                final int depthA = calculateNodeDepth(nodeA);
                final int depthB = calculateNodeDepth(nodeB);
                if (depthA == depthB) {
                    final int nodeNumberA = nodeA.getNr();
                    final int nodeNumberB = nodeB.getNr();
                    if (nodeNumberA == nodeNumberB) return 0;
                    return nodeNumberA > nodeNumberB ? greaterThan : lessThan;
                }
                return depthA > depthB ? greaterThan : lessThan;
            }
            return heightA > heightB ? greaterThan : lessThan;
        }

        private int calculateNodeDepth(Node node) {
            if (node.isRoot()) return 0;
            return calculateNodeDepth(node.getParent()) - 1;
        }
    }
}
