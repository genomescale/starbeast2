package starbeast2;

import java.io.PrintStream;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;

/**
 * Network class to replace Tree class
 * It includes both bifurcation node (in-degree 1, out-degree 2 or 0(tip)) and reticulation node (in degree 2, out-degree 1).
 */

@Description("Network representing reticulate evolution of species")
public class Network extends StateNode {  //implements TreeInterface
    final public Input<Network> networkInput = new Input<>("initial", "network to start with");
    final public Input<String> nodeTypeInput = new Input<>("nodetype", "type of the node in the network", NetworkNode.class.getName());

    /**
     * state of dirtiness of a node in the tree
     * DIRTY means a property on the node has changed, but not the topology. "property" includes the node height
     *       and that branch length to its parent.
     * FILTHY means the node's parent or child has changed.
     */
    public static final int IS_CLEAN = 0, IS_DIRTY = 1, IS_FILTHY = 2;

    // number of nodes
    protected int nodeCount = -1;
    protected int internalNodeCount = -1;
    protected int leafNodeCount = -1;

    protected NetworkNode root;
    protected NetworkNode storedRoot;

    /**
     * array of all nodes in the network
     */
    protected NetworkNode[] networkNodes = null;
    protected NetworkNode[] storedNetworkNodes = null;

    @Override
    public void initAndValidate() throws Exception {
        // bla bla
    }


    public Network() {
    }

    public Network(final NetworkNode rootNode) {
        setRoot(rootNode);
        // initArrays();
    }

    public Network(final String sNewick) {
    }

    protected NetworkNode newNode() {
        try {
            return (NetworkNode) Class.forName(nodeTypeInput.get()).newInstance();
        } catch (Exception e) {
            throw new RuntimeException("Cannot create node of type " + nodeTypeInput.get() + ": " + e.getMessage());
        }
    }

    public NetworkNode getRoot() {
        return root;
    }

    public void setRoot(final NetworkNode root) {
        this.root = root;
        // ensure root is the last node in networkNodes
    }
    /**
     * @return the number of nodes
     */
    public int getNodeCount() {
        if (nodeCount < 0)
            nodeCount = this.root.getNodeCount();
        return nodeCount;
    }

    public int getInternalNodeCount() {
        if (internalNodeCount < 0) {
            internalNodeCount = root.getInternalNodeCount();
        }
        return internalNodeCount;
    }

    public int getLeafNodeCount() {
        if (leafNodeCount < 0) {
            leafNodeCount = root.getLeafNodeCount();
        }
        return leafNodeCount;
    }

    /**
     * @return a list of leaf nodes contained in this network
     */
    public List<NetworkNode> getLeafNodes() {
        final ArrayList<NetworkNode> lNodes = new ArrayList<>();
        for (final NetworkNode node : networkNodes) {
            if (node.isLeaf()) lNodes.add(node);
        }
        return lNodes;
    }

    public NetworkNode getNode(final int nr) {
        return networkNodes[nr];
    }

    public NetworkNode[] getNodesAsArray() {
        return networkNodes;
    }

    public String toString() {
        return root.toString();
    }

    @Override
    public void setEverythingDirty(final boolean isDirty) {
        setSomethingIsDirty(isDirty);
        if (!isDirty) {
            for(NetworkNode n : networkNodes) {
                n.isDirty = IS_CLEAN;
            }
        } else {
            for(NetworkNode n : networkNodes) {
                n.isDirty = IS_FILTHY;
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
        network.internalNodeCount = internalNodeCount;
        network.leafNodeCount = leafNodeCount;
        return network;
    }

    /**
     * copy of all values into existing network
     */
    @Override
    public void assignTo(final StateNode other) {
        final Network network = (Network) other;
        final NetworkNode[] nodes = new NetworkNode[nodeCount];
        // ??? listNodes(network.root, nodes);
        network.setID(getID());
        //tree.index = index;
        root.assignTo(nodes);
        network.root = nodes[root.getNr()];
        network.nodeCount = nodeCount;
        network.internalNodeCount = internalNodeCount;
        network.leafNodeCount = leafNodeCount;
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
        //index = tree.index;
        root = nodes[network.root.getNr()];
        root.assignFrom(nodes, network.root);
        root.parents = null;
        nodeCount = network.nodeCount;
        internalNodeCount = network.internalNodeCount;
        leafNodeCount = network.leafNodeCount;

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
        root.parents = null;
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
    public int scale(final double scale) throws Exception {
        root.scale(scale);
        return getInternalNodeCount(); // is this number correct??? should return nBifurcationNode+2*nReticulationNode?
    }

    /**
     * reconstruct tree from XML fragment in the form of a DOM node
     */
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        // initArrays();
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
            // ???
        }
    }

    @Override
    public void restore() {
        nodeCount = storedNetworkNodes.length;

        final NetworkNode[] tmp = storedNetworkNodes;
        storedNetworkNodes = networkNodes;
        networkNodes = tmp;
        root = networkNodes[storedRoot.getNr()];

        // leafNodeCount = root.getLeafNodeCount();

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
    public void init(PrintStream out) throws Exception {
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
        Network network = (Network) getCurrent();
        out.print("network STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        final int[] dummy = new int[1];
        final String newick = network.getRoot().toSortedNewick(dummy);
        out.print(newick);
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
}
