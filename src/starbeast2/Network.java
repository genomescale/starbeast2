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

    /**
     * counters of number of nodes
     */
    protected int nodeCount = -1;
    protected int bifurcationNodeCount = -1;
    protected int reticulationNodeCount = -1;
    protected int leafNodeCount = -1;

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

    /**
     * getters and setters
     *
     * @return the number of nodes
     */
    public int getNodeCount() {
        if (nodeCount < 0) {
            nodeCount = this.root.getNodeCount();
        }

        //System.out.println("nodeCount=" + nodeCount);
        return nodeCount;
    }

    /**
     * @return a list of leaf nodes contained in this network
     */
    public List<NetworkNode> getLeafNodes() {
        final ArrayList<NetworkNode> lNodes = new ArrayList<>();
        // for (int i = 0; i < getNodeCount(); i++) {
        //    final NetworkNode node = getNode(i);
        //    if (node.isLeaf()) externalNodes.add(node);
        // }
        return lNodes;
    }

    public NetworkNode getRoot() {
        return root;
    }

    public void setRoot(final NetworkNode root) {
        this.root = root;
        // ensure root is the last node in networkNodes
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

    /**
     * reconstruct tree from XML fragment in the form of a DOM node
     */
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        // initArrays();
    }

    @Override
    public double getArrayValue() {
        return root.height;
    }

    @Override
    public double getArrayValue(int iValue) {
        return networkNodes[iValue].height;
    }

}
