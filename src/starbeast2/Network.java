package starbeast2;

import java.io.*;
import java.util.*;

import beast.core.*;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import beast.evolution.tree.*;

/**
 * Network class to replace Tree class
 * It includes both bifurcation node (in-degree 1, out-degree 2 or 0(tip)) and reticulation node (in degree 2, out-degree 1).
 */

@Description("Network representing reticulate evolution of species")
public class Network extends StateNode {  //implements TreeInterface
    public Input<Network> networkInput = new Input<>("initial", "species network to start with");

    /**
     * state of dirtiness of a node in the tree
     * DIRTY means a property on the node has changed, but not the topology. "property" includes the node height
     *       and that branch length to its parent.
     * FILTHY means the node's parent or child has changed.
     */
    public static final int IS_CLEAN = 0, IS_DIRTY = 1, IS_FILTHY = 2;

    // number of nodes
    // protected int nodeCount = -1;

    protected NetworkNode root;

    /**
     * array of all nodes in the network *
     */
    protected NetworkNode[] networkNodes = null;

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
        return this.root.getNodeCount();
    }

    /**
     * @return a list of leaf nodes contained in this network
     */
    public List<NetworkNode> getLeafNodes() {
        final ArrayList<NetworkNode> lNodes = new ArrayList<>();
        for (final NetworkNode node: networkNodes) {
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
            for( Node n : m_nodes ) {
                n.isDirty = IS_CLEAN;
            }
            //  root.makeAllDirty(IS_CLEAN);
        } else {
            for( Node n : m_nodes ) {
                n.isDirty = IS_FILTHY;
            }
            //    root.makeAllDirty(IS_FILTHY);
        }
    }

    /**
     * @return a deep copy of this network.
     */
    @Override
    public Network copy() {
        Network network = new Network();
        network.setID(getID());
        network.index = index;
        network.root = root.copy();
        return network;
    }

    /**
     * copy of all values into existing tree *
     */
    @Override
    public void assignTo(final StateNode other) {
        final Tree tree = (Tree) other;
        final Node[] nodes = new Node[nodeCount];
        listNodes(tree.root, nodes);
        tree.setID(getID());
        //tree.index = index;
        root.assignTo(nodes);
        tree.root = nodes[root.getNr()];
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
    }

    /**
     * copy of all values from existing tree *
     */
    @Override
    public void assignFrom(final StateNode other) {
        final Tree tree = (Tree) other;
        final Node[] nodes = new Node[tree.getNodeCount()];//tree.getNodesAsArray();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            nodes[i] = newNode();
        }
        setID(tree.getID());
        //index = tree.index;
        root = nodes[tree.root.getNr()];
        root.assignFrom(nodes, tree.root);
        root.parent = null;
        nodeCount = tree.nodeCount;
        internalNodeCount = tree.internalNodeCount;
        leafNodeCount = tree.leafNodeCount;
        initArrays();
    }

    /**
     * as assignFrom, but only copy tree structure *
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        final Tree tree = (Tree) other;
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[tree.root.getNr()];
        final Node[] otherNodes = tree.m_nodes;
        final int rootNr = root.getNr();
        assignFrom(0, rootNr, otherNodes);
        root.height = otherNodes[rootNr].height;
        root.parent = null;
        if (otherNodes[rootNr].getLeft() != null) {
            root.setLeft(m_nodes[otherNodes[rootNr].getLeft().getNr()]);
        } else {
            root.setLeft(null);
        }
        if (otherNodes[rootNr].getRight() != null) {
            root.setRight(m_nodes[otherNodes[rootNr].getRight().getNr()]);
        } else {
            root.setRight(null);
        }
        assignFrom(rootNr + 1, nodeCount, otherNodes);
    }

    /**
     * helper to assignFromFragile *
     */
    private void assignFrom(final int start, final int end, final Node[] otherNodes) {
        for (int i = start; i < end; i++) {
            Node sink = m_nodes[i];
            Node src = otherNodes[i];
            sink.height = src.height;
            sink.parent = m_nodes[src.parent.getNr()];
            if (src.getLeft() != null) {
                sink.setLeft(m_nodes[src.getLeft().getNr()]);
                if (src.getRight() != null) {
                    sink.setRight(m_nodes[src.getRight().getNr()]);
                } else {
                    sink.setRight(null);
                }
            }
        }
    }

    @Override
    public int scale(final double scale) throws Exception {
        root.scale(scale);
        return getInternalNodeCount()- getDirectAncestorNodeCount();
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

        // this condition can only be true for sampled ancestor trees
        if (m_storedNodes.length != nodeCount) {
            final Node[] tmp = new Node[nodeCount];
            System.arraycopy(m_storedNodes, 0, tmp, 0, m_storedNodes.length - 1);
            if (nodeCount > m_storedNodes.length) {
                tmp[m_storedNodes.length - 1] = m_storedNodes[m_storedNodes.length - 1];
                tmp[nodeCount - 1] = newNode();
                tmp[nodeCount - 1].setNr(nodeCount - 1);
            }
            m_storedNodes = tmp;
        }

        storeNodes(0, nodeCount);
        storedRoot = m_storedNodes[root.getNr()];
    }


    /**
     * Stores nodes with index i, for start <= i < end
     * (i.e. including start but not including end)
     *
     * @param start the first index to be stored
     * @param end   nodes are stored up to but not including this index
     */
    private void storeNodes(final int start, final int end) {
        // Use direct members for speed (we are talking 5-7% or more from total time for large trees :)
        for (int i = start; i < end; i++) {
            final Node sink = m_storedNodes[i];
            final Node src = m_nodes[i];
            sink.height = src.height;

            if ( src.parent != null ) {
                sink.parent = m_storedNodes[src.parent.getNr()];
            } else {
                // currently only called in the case of sampled ancestor trees
                // where root node is not always last in the list
                sink.parent = null;
            }

            final List<Node> children = sink.children;
            final List<Node> srcChildren = src.children;

            if( children.size() == srcChildren.size() ) {
                // shave some more time by avoiding list clear and add
                for (int k = 0; k < children.size(); ++k) {
                    final Node srcChild = srcChildren.get(k);
                    // don't call addChild, which calls  setParent(..., true);
                    final Node c = m_storedNodes[srcChild.getNr()];
                    c.parent = sink;
                    children.set(k, c);
                }
            } else {
                children.clear();
                //sink.removeAllChildren(false);
                for (final Node srcChild : srcChildren) {
                    // don't call addChild, which calls  setParent(..., true);
                    final Node c = m_storedNodes[srcChild.getNr()];
                    c.parent = sink;
                    children.add(c);
                    //sink.addChild(c);
                }
            }
        }
    }

    @Override
    public void restore() {

        // necessary for sampled ancestor trees
        nodeCount = m_storedNodes.length;

        final Node[] tmp = m_storedNodes;
        m_storedNodes = m_nodes;
        m_nodes = tmp;
        root = m_nodes[storedRoot.getNr()];

        // necessary for sampled ancestor trees,
        // we have the nodes, no need for expensive recursion
        leafNodeCount = 0;
        for( Node n : m_nodes ) {
            leafNodeCount += n.isLeaf() ? 1 : 0;
        }

        //leafNodeCount = root.getLeafNodeCount();

        hasStartedEditing = false;

        for( Node n : m_nodes ) {
            n.isDirty = Tree.IS_CLEAN;
        }

        postCache = null;
    }

    @Override
    public void init(PrintStream out) throws Exception {
        Node node = getRoot();
        out.println("#NEXUS\n");
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + getLeafNodeCount() + ";");
        out.println("\t\tTaxlabels");
        printTaxa(node, out, getNodeCount() / 2);
        out.println("\t\t\t;");
        out.println("End;");

        out.println("Begin trees;");
        out.println("\tTranslate");
        printTranslate(node, out, getNodeCount() / 2);
        out.print(";");
    }

    @Override
    public void log(int sample, PrintStream out) {
        Tree tree = (Tree) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final int[] dummy = new int[1];
        final String newick = tree.getRoot().toSortedNewick(dummy);
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
