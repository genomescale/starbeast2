package starbeast2;

import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

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
    public Input<SpeciesNetwork> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network for embedding the gene tree.", Validate.REQUIRED);
    public Input<Tree> geneTreeInput =
            new Input<>("geneTree", "Gene tree embedded in the Species network.", Validate.REQUIRED);
    public Input<Double> ploidyInput =
            new Input<>("ploidy", "Ploidy (copy number) for this gene (default is 2).", 2.0);
    protected double ploidy;

    private int geneTreeLeafNodeCount;
    private int geneTreeNodeCount;
    private boolean needsUpdate;

    protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create(); // the coalescent event times for this gene tree for all species tree branches
    protected ListMultimap<Integer, Double> storedCoalescentTimes = ArrayListMultimap.create(); // the coalescent event times for this gene tree for all species tree branches
    protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create(); // the number of lineages at the tipward end of each branch
    protected Multiset<Integer> storedCoalescentLineageCounts = HashMultiset.create(); // the number of lineages at the tipward end of each branch

    protected int[] geneNodeSpeciesAssignment;
    protected int[] storedGeneNodeSpeciesAssignment;
    protected double[][] speciesOccupancy;
    protected double[][] storedSpeciesOccupancy;
    protected boolean geneTreeCompatible;
    protected boolean storedGeneTreeCompatible;

    /**
     * gene tree lineage inheritance direction
     * Tracing backward in time, true -> left parent and false -> right parent.
     * The length of the boolean array (1st dimension) equals to the number of gene tips.
     * The length of the list (2nd dimension) equals to the number of reticulation (hybridization) nodes.
     * There have methods to insert/delete a list element (boolean array) when adding/deleting a reticulation event,
     * and change the booleans if the gene tree lineages change the ancestral population at a reticulation event.
     */
    public List<boolean[]> lineageInheritance = new ArrayList<>();

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = geneTreeInput.isDirty() || speciesNetworkInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
        storedCoalescentTimes.clear();
        storedCoalescentLineageCounts.clear();

        storedCoalescentTimes.putAll(coalescentTimes);
        storedCoalescentLineageCounts.addAll(coalescentLineageCounts);

        storedSpeciesOccupancy = new double[speciesOccupancy.length][speciesOccupancy[0].length];
        System.arraycopy(geneNodeSpeciesAssignment, 0, storedGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment.length);
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);

        storedGeneTreeCompatible = geneTreeCompatible;

        super.store();
    }

    @Override
    public void restore() {
        ListMultimap<Integer, Double> tmpCoalescentTimes = coalescentTimes;
        Multiset<Integer> tmpCoalescentLineageCounts = coalescentLineageCounts;
        int[] tmpGeneNodeSpeciesAssignment = geneNodeSpeciesAssignment;
        double[][] tmpSpeciesOccupancy = speciesOccupancy;
        boolean tmpGeneTreeCompatible = geneTreeCompatible;

        coalescentTimes = storedCoalescentTimes;
        coalescentLineageCounts = storedCoalescentLineageCounts;
        speciesOccupancy = storedSpeciesOccupancy;
        geneNodeSpeciesAssignment = storedGeneNodeSpeciesAssignment;
        geneTreeCompatible = storedGeneTreeCompatible;

        storedCoalescentTimes = tmpCoalescentTimes;
        storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
        storedSpeciesOccupancy = tmpSpeciesOccupancy;
        storedGeneNodeSpeciesAssignment = tmpGeneNodeSpeciesAssignment;
        storedGeneTreeCompatible = tmpGeneTreeCompatible;

        super.restore();
    }

    public void initAndValidate() throws Exception {
        ploidy = ploidyInput.get();

        geneTreeNodeCount = geneTreeInput.get().getNodeCount();
        geneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        storedGeneNodeSpeciesAssignment = new int[geneTreeNodeCount];

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
        final Map<String, Integer> tipNumberMap = speciesNetworkInput.get().getTipNumberMap();

        final int speciesNetworkNodeCount = speciesNetwork.getNodeCount();
        speciesOccupancy = new double[geneTreeNodeCount][2*speciesNetworkNodeCount];

        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
        Arrays.fill(geneNodeSpeciesAssignment, -1);
        // -1 means no species assignment for that gene tree node has been made yet

        coalescentLineageCounts.clear();
        coalescentTimes.clear();

        final TreeInterface geneTree = geneTreeInput.get();
        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
            final int speciesNetworkLeafNumber = tipNumberMap.get(geneTreeLeafNode.getID());
            final NetworkNode speciesNetworkLeafNode = speciesNetwork.getNode(speciesNetworkLeafNumber);
            coalescentLineageCounts.add(2 * speciesNetworkLeafNumber);

            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
            // final int firstCoalescenceNumber = firstCoalescenceNode.getNr();
            final double lastHeight = 0.0;
            if (!recurseCoalescenceEvents(geneTreeLeafNumber, lastHeight, firstCoalescenceNode,
                                            speciesNetworkLeafNode, 2*speciesNetworkLeafNumber)) {
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
                                             final NetworkNode speciesNetworkNode, final int speciesNetworkNodeNumber) {
        final double geneTreeNodeHeight = geneTreeNode.getHeight();
        final int geneTreeNodeNumber = geneTreeNode.getNr();

        // check if the next coalescence event occurs in an ancestral branch
        if (!speciesNetworkNode.isRoot()) {
            final NetworkNode speciesNetworkParentNode = getNetworkParentNode(speciesNetworkNode, geneTreeNode);
            final double speciesNetworkParentHeight = speciesNetworkParentNode.getHeight();
            if (geneTreeNodeHeight >= speciesNetworkParentHeight) {
                speciesOccupancy[lastGeneTreeNodeNumber][speciesNetworkNodeNumber] = speciesNetworkParentHeight - lastHeight;
                final int speciesNetworkParentNodeNumber = 2 * speciesNetworkParentNode.getNr()
                                                             + getOffset(speciesNetworkParentNode, geneTreeNode);
                coalescentLineageCounts.add(speciesNetworkParentNodeNumber);
                return recurseCoalescenceEvents(lastGeneTreeNodeNumber, speciesNetworkParentHeight, geneTreeNode,
                                                speciesNetworkParentNode, speciesNetworkParentNodeNumber);
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
                return recurseCoalescenceEvents(geneTreeNodeNumber, geneTreeNodeHeight, nextGeneTreeNode,
                                                speciesNetworkNode, speciesNetworkNodeNumber);
            }
        } else {
            // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
            // return false if this gene tree IS NOT compatible with the species tree
            return (existingSpeciesAssignment == speciesNetworkNodeNumber);
        }
    }

    public double[][] getSpeciesOccupancy() {
        if (needsUpdate) update();

        return speciesOccupancy;
    }

    /* public double[] getOccupancy(Node node) {
        if (needsUpdate) {
            update();
        }

        final int geneTreeNodeNumber = node.getNr();
        return speciesOccupancy[geneTreeNodeNumber];
    } */

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
     * @param networkNode the current network node
     * @param gTreeNode the gene tree node (coalescent event)
     * @return the parent network node toward which the gene tree node goes
     * the method is here instead of in NetworkNode class because it needs information of the gene tree node
     */
    public NetworkNode getNetworkParentNode(NetworkNode networkNode, Node gTreeNode) {
        if (!networkNode.isReticulation()) {
            if (networkNode.getLeftParent() != null)
                return networkNode.getLeftParent();
            else if (networkNode.getRightParent() != null)
                return networkNode.getRightParent();
            else
                return null; // networkNode is root
        } else {  // networkNode is a reticulation node
            // find the tip node which is descendant of gTreeNode, then its number
            final int gTreeTipNr = getGeneTreeTipDescendant(gTreeNode).getNr();
            // find the boolean of inheritance (gTreeTipNr'th array element)
            if (lineageInheritance.get(networkNode.getReticulateNr())[gTreeTipNr])
                return networkNode.getLeftParent();
            else
                return networkNode.getRightParent();
        }
    }

    /**
     * @param networkNode the current network node
     * @param gTreeNode the gene tree node (coalescent event)
     * @return 1 if the gene tree node goes to right parent at reticulation node, 0 otherwise
     */
    protected int getOffset(NetworkNode networkNode, Node gTreeNode) {
        if (!networkNode.isReticulation()) {
            return 0;
        } else {  // networkNode is a reticulation node
            // find the tip node which is descendant of gTreeNode, then its number
            final int gTreeTipNr = getGeneTreeTipDescendant(gTreeNode).getNr();
            // find the boolean of inheritance (gTreeTipNr'th array element)
            if (lineageInheritance.get(networkNode.getReticulateNr())[gTreeTipNr])
                return 0;
            else
                return 1;
        }
    }

    /**
     * @return the first tip node which is descendant of
     * @param gTreeNode
     * this can be in Tree.java as gTreeNode.getGeneTreeTipDescendant()
     */
    private Node getGeneTreeTipDescendant(Node gTreeNode) {
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
