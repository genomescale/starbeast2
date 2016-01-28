package starbeast2;

import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.List;

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
     * The length of the boolean array (2nd dimension) equals to the number of gene tips.
     * The length of the list (1st dimension) equals to the number of reticulation (hybridization) nodes.
     * There have methods to insert/delete a list element (boolean array) when adding/deleting a reticulation event,
     * and change the booleans if the gene tree lineages change the ancestral population at a reticulation event.
     */
    // public List<boolean[]> lineageInheritance = new ArrayList<>();

    /**
     * This data structure is replaced by a list of ancestral species nodes for each gene tree node.
     * The information probably can be stored in metaData in Node.java, but is left here at the moment.
     * The 1st dimension is for the gene tree nodes/lineages.
     * The 2nd dimension is for the network nodes that each gene tree lineage has traversed.
     */
    public ListMultimap<Node, NetworkNode> traversedNetworkNodes = ArrayListMultimap.create();

    /**
     * The 1st dimension is for the gene tree nodes/lineages. For each gene tree node/lineage,
     * the 2nd dimension is for the directions at the network nodes that the lineage goes.
     */
    public ListMultimap<Node, Boolean> traversedInheritances = ArrayListMultimap.create();

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

        /* traversedNetworkNodes.clear();
        traversedInheritances.clear();
        for (int i = 0; i < geneTreeNodeCount; i++) {
            traversedNetworkNodes.add(new ArrayList<>());
            traversedInheritances.add(new ArrayList<>());
        } */

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
            if (!recurseCoalescenceEvents(geneTreeLeafNode, lastHeight, firstCoalescenceNode,
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

    private boolean recurseCoalescenceEvents(final Node lastGeneTreeNode, final double lastHeight, final Node geneTreeNode,
                                             final NetworkNode speciesNetworkNode, final int speciesNetworkPopNumber) {
        final double geneTreeNodeHeight = geneTreeNode.getHeight();
        final int geneTreeNodeNumber = geneTreeNode.getNr();
        final int lastGeneTreeNodeNumber = lastGeneTreeNode.getNr();

        // check if the next coalescence event occurs in an ancestral branch
        if (!speciesNetworkNode.isRoot()) {
            final NetworkNode speciesNetworkParentNode = speciesNetworkNode.getParent(lastGeneTreeNode);
            final double speciesNetworkParentHeight = speciesNetworkParentNode.getHeight();
            if (geneTreeNodeHeight >= speciesNetworkParentHeight) {
                speciesOccupancy[lastGeneTreeNodeNumber][speciesNetworkPopNumber] = speciesNetworkParentHeight - lastHeight;
                final int speciesNetworkParentNodeNumber = 2 * speciesNetworkParentNode.getNr()
                                                             + getOffset(speciesNetworkParentNode, lastGeneTreeNode);
                coalescentLineageCounts.add(speciesNetworkParentNodeNumber);
                return recurseCoalescenceEvents(lastGeneTreeNode, speciesNetworkParentHeight, geneTreeNode,
                                                speciesNetworkParentNode, speciesNetworkParentNodeNumber);
            }
        }

        // this code executes if the next coalescence event occurs within the current branch
        speciesOccupancy[lastGeneTreeNodeNumber][speciesNetworkPopNumber] = geneTreeNodeHeight - lastHeight;
        final int existingSpeciesAssignment = geneNodeSpeciesAssignment[geneTreeNodeNumber];
        if (existingSpeciesAssignment == -1) {
            geneNodeSpeciesAssignment[geneTreeNodeNumber] = speciesNetworkPopNumber;
            coalescentTimes.put(speciesNetworkPopNumber, geneTreeNodeHeight);
            final Node nextGeneTreeNode = geneTreeNode.getParent();
            if (nextGeneTreeNode == null) {
                // this is the root of the gene tree and no incompatibilities were detected
                return true;
            } else {
                // if this is not the root of the gene tree, check the subsequent (back in time) coalescence event
                return recurseCoalescenceEvents(geneTreeNode, geneTreeNodeHeight, nextGeneTreeNode,
                                                speciesNetworkNode, speciesNetworkPopNumber);
            }
        } else {
            // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
            // return false if this gene tree IS NOT compatible with the species tree
            return existingSpeciesAssignment == speciesNetworkPopNumber;
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
            // find the networkNode that the lineage of gTreeNode has traversed, get its index
            final int index = traversedNetworkNodes.get(gTreeNode).indexOf(networkNode);
            // need to check the index is not out of bounds ???
            // find the boolean of inheritance: true -> left, false -> right
            if (traversedInheritances.get(gTreeNode).get(index))
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
            // find the networkNode that the lineage of gTreeNode has traversed, get its index
            final int index = traversedNetworkNodes.get(gTreeNode).indexOf(networkNode);
            // find the boolean of inheritance: true -> left, false -> right
            return (traversedInheritances.get(gTreeNode).get(index)) ? 0 : 1;
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
