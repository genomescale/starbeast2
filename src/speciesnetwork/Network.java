package speciesnetwork;

import java.io.PrintStream;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;

/**
 * Network class to replace Tree class
 * It includes both bifurcation node (in-degree 1, out-degree 2 or 0(tip)) and reticulation node (in degree 2, out-degree 1).
 * @author Chi Zhang
 * @author Huw Ogilvie
 */

@Description("Network representing reticulate evolution of species")
public class Network extends StateNode {
    final public Input<Network> networkInput =
            new Input<>("initial", "Network to start with.");
    final public Input<String> nodeTypeInput =
            new Input<>("nodetype", "Type of the node in the network.", NetworkNode.class.getName());
    final public Input<TaxonSet> taxonSetInput =
            new Input<>("taxonset", "Set of taxa at the leafs of the network.");

    /**
     * state of dirtiness of a node in the network
     * DIRTY means a property on the node has changed, but not the topology. "property" includes the node height
     *       and that branch length to its parent.
     * FILTHY means the node's parent or child has changed.
     */
    public static final int IS_CLEAN = 0, IS_DIRTY = 1, IS_FILTHY = 2;

    // number of nodes
    protected int nodeCount = -1;
    protected int speciationNodeCount = -1;
    protected int leafNodeCount = -1;
    protected int reticulationNodeCount = -1;

    protected NetworkNode root;
    protected NetworkNode storedRoot;

    /**
     * array of all nodes in the network
     */
    protected NetworkNode[] networkNodes = null;
    protected NetworkNode[] storedNetworkNodes = null;

    /**
     * array of taxa names for the nodes in the network
     * such that taxaNames[node.getNr()] == node.getID()
     */
    String[] taxaNames = null;

    @Override
    public void initAndValidate() {
        if (nodeCount < 0) {
            if (taxonSetInput.get() != null) {
                makeCaterpillar(0, 1, false);
            } else {
                // make dummy network with a single root node
                root = newNode();
                root.labelNr = 0;
                root.height = 0;
                root.network = this;
                nodeCount = 1;
                speciationNodeCount = reticulationNodeCount = 0;
                leafNodeCount = 1;
            }
        }

        initArrays();

        // ensure all nodes have their taxon names set up
        // needs more things here???
    }

    public void makeCaterpillar(final double minInternalHeight, final double step, final boolean finalize) {
        // make a caterpillar
        final List<String> taxa = taxonSetInput.get().asStringList();
        NetworkNode left = newNode();
        left.labelNr = 0;
        left.height = 0;
        left.setID(taxa.get(0));
        for (int i = 1; i < taxa.size(); i++) {
            final NetworkNode right = newNode();
            right.labelNr = i;
            right.height = 0;
            right.setID(taxa.get(i));
            final NetworkNode parent = newNode();
            parent.labelNr = taxa.size() + i - 1;
            parent.height = minInternalHeight + i * step;
            // left.parent = parent;
            parent.setLeftChild(left);
            // right.parent = parent;
            parent.setRightChild(right);
            left = parent;
        }
        root = left;
        leafNodeCount = taxa.size();
        nodeCount = leafNodeCount * 2 - 1;
        speciationNodeCount = leafNodeCount - 1;

        if (finalize) {
            initArrays();
        }
    }

    /**
     * constructors
     */
    public Network() {
    }

    // constructor overload
    public Network(final NetworkNode rootNode) {
        setRoot(rootNode);
        initArrays();
    }

    // construct a network from newick string -- will not automatically adjust tips to zero
    public Network(final String sNewick) {
    }

    public NetworkNode newNode() {
        try {
            return (NetworkNode) Class.forName(nodeTypeInput.get()).newInstance();
        } catch (Exception e) {
            throw new RuntimeException("Cannot create node of type " + nodeTypeInput.get() + ": " + e.getMessage());
        }
    }

    protected void initArrays() {
        // initialise network-as-array representation + its stored variant
        networkNodes = new NetworkNode[nodeCount];
        listNodes(root, networkNodes);
        storedNetworkNodes = new NetworkNode[nodeCount];
        final NetworkNode copy = root.copy();
        listNodes(copy, storedNetworkNodes);
    }

    public NetworkNode getRoot() {
        return root;
    }

    public void setRoot(final NetworkNode root) {
        this.root = root;
        nodeCount = this.root.getNodeCount();
        // ensure root is the last node in networkNodes
        if (networkNodes != null && root.labelNr != networkNodes.length - 1) {
            final int rootPos = networkNodes.length - 1;
            NetworkNode tmp = networkNodes[rootPos];
            networkNodes[rootPos] = root;
            networkNodes[root.labelNr] = tmp;
            tmp.labelNr = root.labelNr;
            networkNodes[rootPos].labelNr = rootPos;
        }
    }

    /**
     * @return the number of nodes
     */
    public int getNodeCount() {
        nodeCount = root.getNodeCount();
        return nodeCount;
    }

    public int getSpeciationNodeCount() {
        speciationNodeCount = root.getSpeciationNodeCount();
        return speciationNodeCount;
    }

    public int getLeafNodeCount() {
        if (leafNodeCount < 0)  // assuming no change after set
            leafNodeCount = root.getLeafNodeCount();
        return leafNodeCount;
    }

    public int getReticulationNodeCount() {
        reticulationNodeCount = root.getReticulationNodeCount();
        return reticulationNodeCount;
    }

    /**
     * @return the number of branches in the network, including the root branch
     */
    public int getBranchCount() {
        return getNodeCount() + getReticulationNodeCount();
    }

    /**
     * @return the number of branches at the given time
     */
    public int getBranchCount(double time) {
        int nB = 1;
        for (NetworkNode node : networkNodes) {
            if (node.getHeight() > time) {
                if (node.isReticulation()) nB--;
                else nB++;
            }
        }
        return  nB;
    }

    /**
     * @return a list of leaf nodes contained in this network
     */
    public List<NetworkNode> getLeafNodes() {
        final List<NetworkNode> lNodes = new ArrayList<>();
        for (final NetworkNode node : networkNodes) {
            if (node.isLeaf()) lNodes.add(node);
        }
        return lNodes;
    }

    /**
     * @return a list of internal nodes contained in this network
     */
    public List<NetworkNode> getInternalNodes() {
        final List<NetworkNode> iNodes = new ArrayList<>();
        for (final NetworkNode node : networkNodes) {
            if (!node.isLeaf()) iNodes.add(node);
        }
        return iNodes;
    }

    /**
     * @return a list of reticulation nodes contained in this network
     */
    public List<NetworkNode> getReticulationNodes() {
        final List<NetworkNode> iNodes = new ArrayList<>();
        for (final NetworkNode node : networkNodes) {
            if (node.isReticulation()) iNodes.add(node);
        }
        return iNodes;
    }

    public NetworkNode getNode(final int nr) {
        return networkNodes[nr];
    }

    public NetworkNode getNode(final String label) {
        for (NetworkNode netNode: networkNodes) {
            if (netNode.getID().equals(label)) return netNode;
        }
        return null;
    }

    /**
     * @return an array of all the nodes in this network
     */
    public NetworkNode[] getAllNodesAsArray() {
        return networkNodes;
    }

    /**
     * convert network to array representation
     */
    void listNodes(final NetworkNode node, final NetworkNode[] nodes) {
        nodes[node.getNr()] = node;
        node.network = this;  //(JH) I don't understand this code
        for (final NetworkNode child : node.getChildren()) {
            listNodes(child, nodes);
        }
    }

    public double getNetworkLength () {
        double length = 0;
        for (final NetworkNode node : root.getAllChildNodes()) {
            if(node.getLeftParent() != null) length += node.getLeftLength();
            if(node.getRightParent()!= null) length += node.getRightLength();
        }
        return length;
    }

    /**
     * @return an array of taxon names in order of their node numbers
     */
    public String[] getTaxaNames() {
        if (taxaNames == null || taxaNames.length == 0) {
            final TaxonSet taxonSet = taxonSetInput.get();
            if (taxonSet != null) {
                String[] array = new String[taxonSet.asStringList().size()];
                taxaNames = taxonSet.asStringList().toArray(array);
            } else {
                taxaNames = new String[getNodeCount()];
                collectTaxaNames(getRoot());
            }
        }

        // sanity check
        if (taxaNames.length == 1 && taxaNames[0] == null) {
        }
        return taxaNames;
    }
    private void collectTaxaNames(final NetworkNode node) {
        if (node.getID() != null) {
            taxaNames[node.getNr()] = node.getID();
        }
        if (node.isLeaf()) {
            if (node.getID() == null) {
                node.setID("node" + node.getNr());
            }
        } else {
            for (NetworkNode child : node.getChildren()) {
                collectTaxaNames(child);
            }
        }
    }

    public String toString() {
        return root.toString();
    }

    @Override
    public void setEverythingDirty(final boolean isDirty) {
        setSomethingIsDirty(isDirty);
        if (!isDirty) {
            for(NetworkNode node : networkNodes) {
                node.isDirty = IS_CLEAN;
            }
        } else {
            for(NetworkNode node : networkNodes) {
                node.isDirty = IS_FILTHY;
            }
        }
    }

    /**
     * @return a deep copy of this network
     */
    @Override
    public Network copy() {
        Network network = new Network();
        network.setID(getID());
        network.index = index;
        network.root = root.copy();
        network.nodeCount = nodeCount;
        network.speciationNodeCount = speciationNodeCount;
        network.leafNodeCount = leafNodeCount;
        network.reticulationNodeCount = reticulationNodeCount;
        return network;
    }

    /**
     * copy of all values into existing network
     */
    @Override
    public void assignTo(final StateNode other) {
        final Network network = (Network) other;
        final NetworkNode[] nodes = new NetworkNode[nodeCount];
        listNodes(network.root, nodes);
        network.setID(getID());
        network.index = index;
        root.assignTo(nodes);
        network.root = nodes[root.getNr()];
        network.nodeCount = nodeCount;
        network.speciationNodeCount = speciationNodeCount;
        network.leafNodeCount = leafNodeCount;
        network.reticulationNodeCount = reticulationNodeCount;
    }

    /**
     * copy of all values from existing network
     */
    @Override
    public void assignFrom(final StateNode other) {
        final Network network = (Network) other;
        final NetworkNode[] nodes = new NetworkNode[network.getNodeCount()];
        for (int i = 0; i < network.getNodeCount(); i++) {
            nodes[i] = newNode();
        }
        setID(network.getID());
        index = network.index;
        root = nodes[network.root.getNr()];
        root.assignFrom(nodes, network.root);
        root.leftParent = root.rightParent = null;
        nodeCount = network.nodeCount;
        speciationNodeCount = network.speciationNodeCount;
        leafNodeCount = network.leafNodeCount;
        reticulationNodeCount = network.reticulationNodeCount;

        initArrays();
    }

    /**
     * as assignFrom, but only copy network structure
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        final Network network = (Network) other;
        if (networkNodes == null) {
            initArrays();
        }
        root = networkNodes[network.root.getNr()];
        final NetworkNode[] otherNodes = network.networkNodes;
        final int rootNr = root.getNr();
        assignFrom(0, rootNr, otherNodes);
        root.height = otherNodes[rootNr].height;
        root.leftParent = root.rightParent = null;
        if (otherNodes[rootNr].getLeftChild() != null) {
            root.setLeftChild(networkNodes[otherNodes[rootNr].getLeftChild().getNr()]);
        } else {
            root.setLeftChild(null);
        }
        if (otherNodes[rootNr].getRightChild() != null) {
            root.setRightChild(networkNodes[otherNodes[rootNr].getRightChild().getNr()]);
        } else {
            root.setRightChild(null);
        }
        assignFrom(rootNr + 1, nodeCount, otherNodes);
    }
    /**
     * helper to assignFromFragile *
     */
    private void assignFrom(final int start, final int end, final NetworkNode[] otherNodes) {
        for (int i = start; i < end; i++) {
            NetworkNode sink = networkNodes[i];
            NetworkNode src = otherNodes[i];
            sink.height = src.height;
            if (src.getLeftParent() != null) {
                sink.setLeftParent(networkNodes[src.getLeftParent().getNr()]);
                if (src.getRightParent() != null)
                    sink.setRightParent(networkNodes[src.getRightParent().getNr()]);
                else
                    sink.setRightParent(null);
            }
            if (src.getLeftChild() != null) {
                sink.setLeftChild(networkNodes[src.getLeftChild().getNr()]);
                if (src.getRightChild() != null) {
                    sink.setRightChild(networkNodes[src.getRightChild().getNr()]);
                } else {
                    sink.setRightChild(null);
                }
            }
        }
    }

    @Override
    public int scale(final double scale) {
        root.scale(scale);
        return getSpeciationNodeCount() + getReticulationNodeCount();
    }

    /**
     * reconstruct tree from XML fragment in the form of a DOM node
     */
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        // parser???

        initArrays();
    }

    @Override
    protected void store() {
        storeNodes(0, nodeCount);
        storedRoot = storedNetworkNodes[root.getNr()];
    }

    /**
     * Stores nodes with index i, for start <= i < end
     * (i.e. including start but not including end)
     *
     * @param start the first index to be stored
     * @param end   nodes are stored up to but not including this index
     */
    private void storeNodes(final int start, final int end) {
        for (int i = start; i < end; i++) {
            final NetworkNode sink = storedNetworkNodes[i];
            final NetworkNode src = networkNodes[i];
            sink.height = src.height;
            if (src.leftParent != null)
                sink.leftParent = storedNetworkNodes[src.leftParent.getNr()];
            else
                sink.leftParent = null;
            if (src.rightParent != null)
                sink.rightParent = storedNetworkNodes[src.rightParent.getNr()];
            else
                sink.rightParent = null;
            if (src.leftChild != null)
                sink.leftChild = storedNetworkNodes[src.leftChild.getNr()];
            else
                sink.leftChild = null;
            if (src.rightChild != null)
                sink.rightChild = storedNetworkNodes[src.rightChild.getNr()];
            else
                sink.rightChild = null;
        }
    }

    @Override
    public void restore() {
        nodeCount = storedNetworkNodes.length;

        final NetworkNode[] tmp = networkNodes;
        networkNodes = storedNetworkNodes;
        storedNetworkNodes = tmp;
        root = networkNodes[storedRoot.getNr()];

        hasStartedEditing = false;

        for(NetworkNode n : networkNodes) {
            n.isDirty = Network.IS_CLEAN;
        }
        // postCache = null;
    }

    public static void printTranslate(NetworkNode node, PrintStream out, int nodeCount) {
    }

    public static void printTaxa(final NetworkNode node, final PrintStream out, final int nodeCount) {
        final List<String> translateLines = new ArrayList<>();
        printTranslate(node, out, nodeCount);
        Collections.sort(translateLines);
        for (String line : translateLines) {
            line = line.split("\\s+")[2];
            out.println("\t\t\t" + line.replace(',', ' '));
        }
    }

    @Override
    public void init(PrintStream out) {
        NetworkNode node = getRoot();
        out.println("#NEXUS\n");
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + getLeafNodeCount() + ";");
        out.println("\t\tTaxlabels");
        printTaxa(node, out, getNodeCount() / 2);
        out.println("\t\t\t;");
        out.println("End;");

        out.println("Begin networks;");
        out.println("\tTranslate");
        printTranslate(node, out, getNodeCount() / 2);
        out.print(";");
    }

    @Override
    public void log(int sample, PrintStream out) {
        // Network network = (Network) getCurrent();
        out.print("network STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        // final int[] dummy = new int[1];
        //??? final String newick = network.getRoot().toSortedNewick(dummy);
        // out.print(newick);
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        out.print("End;");
    }

    @Override
    public int getDimension() {
        return getNodeCount();
    }

    @Override
    public double getArrayValue() {
        return root.height;
    }

    @Override
    public double getArrayValue(int nr) {
        return networkNodes[nr].height;
    }

    public void addLeafNode(final NetworkNode newNode) {
        final NetworkNode[] tmp = new NetworkNode[nodeCount + 1];
        System.arraycopy(networkNodes, 0, tmp, 1, nodeCount);
        networkNodes = tmp;
        networkNodes[0] = newNode;
        newNode.network = this;
        leafNodeCount++;
        resetNodeNumbers();
    }

    public void addSpeciationNode(final NetworkNode newNode) {
        final NetworkNode[] tmp = new NetworkNode[nodeCount + 1];
        System.arraycopy(networkNodes, 0, tmp, 0, leafNodeCount);
        System.arraycopy(networkNodes, leafNodeCount, tmp, leafNodeCount + 1, speciationNodeCount + reticulationNodeCount);
        networkNodes = tmp;
        networkNodes[leafNodeCount] = newNode;
        newNode.network = this;
        speciationNodeCount++;
        resetNodeNumbers();
    }

    public void addReticulationNode(final NetworkNode newNode) {
        final NetworkNode[] tmp = new NetworkNode[nodeCount + 1];
        System.arraycopy(networkNodes, 0, tmp, 0, nodeCount - 1);
        tmp[nodeCount] = networkNodes[nodeCount - 1];
        networkNodes = tmp;
        networkNodes[nodeCount - 1] = newNode;
        newNode.network = this;
        reticulationNodeCount++;
        resetNodeNumbers();
    }

    private void resetNodeNumbers() {
        nodeCount = networkNodes.length;
        for (int i = 0; i < nodeCount; i++) {
            networkNodes[i].setNr(i);
        }
    }

    public void deleteSpeciationNode(NetworkNode topNode) {
        // TODO Auto-generated method stub
        
    }

    public void deleteReticulationNode(NetworkNode netNode) {
        // TODO Auto-generated method stub
        
    }
}
