package starbeast2;

import java.util.*;

import beast.core.BEASTObject;
import beast.core.Description;

/**
 * NetworkNode that is like Node but has 2 parents and 2 children.
 * The right parent is null if it is a tree node.
 * The right child is null if it is a reticulation node.
 */

@Description("Network node for binary rooted network")
public class NetworkNode extends BEASTObject {

    /**
     * label nr of node
     */
    protected int labelNr;

    // whether the node has been visited, say by a recursive method
    protected boolean visited = false;

    /**
     * height of this node
     */
    protected double height = Double.MAX_VALUE;

    /**
     * list of children of this node
     */
    List<NetworkNode> children = new ArrayList<>();

    /**
     * list of parents of this node
     */
    List<NetworkNode> parents = new ArrayList<>();

    /**
     * status of this node after an operation is performed on the state
     */
    int isDirty = Network.IS_CLEAN;

    /**
     * meta-data contained in square brackets in Newick
     */
    public String metaDataString;

    /**
     * arbitrarily labeled metadata on this node
     */
    protected Map<String, Object> metaData = new TreeMap<>();

    /**
     * the network that this node is a part of
     */
    protected Network network;

    public NetworkNode() {
    }

    public NetworkNode(final String id) throws Exception {
        setID(id);
        initAndValidate();
    }

    @Override
    public void initAndValidate() throws Exception {
        // do nothing
    }

    public Network getNetwork() {
        return network;
    }

    public int getNr() {
        return labelNr;
    }

    public void setNr(final int nr) {
        labelNr = nr;
    }

    public double getHeight() {
        return height;
    }

    /**
     * A Node IS_DIRTY if its value (like height) has changed.
     * A Node IS_FILTHY if its parent or child has changed.
     * Otherwise the node IS_CLEAN.
     */
    public int isDirty() {
        return isDirty;
    }

    public void makeDirty(final int nDirty) {
        isDirty |= nDirty;
    }

    public int getParentCount() {
        return parents.size();
    }

    public int getChildCount() {
        return children.size();
    }

    /* parents */
    public NetworkNode getLeftParent() {
        if (parents.size() == 0) {
            return null;
        }
        return parents.get(0);
    }

    public NetworkNode getRightParent() {
        if (parents.size() <= 1) {
            return null;
        }
        return parents.get(1);
    }

    void setLeftParent(final NetworkNode leftParent) {
        // startEditing();
        if (parents.size() == 0) {
            parents.add(leftParent);
        } else {
            parents.set(0, leftParent);
        }
        isDirty = Network.IS_FILTHY;
    }

    void setRightParent(final NetworkNode rightParent) {
        // startEditing();
        switch (parents.size()) {
            case 0:
                parents.add(null);
            case 1:
                parents.add(rightParent);
                break;
            default:
                parents.set(1, rightParent);
                break;
        }
        isDirty = Network.IS_FILTHY;
    }

    /* children */
    public NetworkNode getLeftChild() {
        if (children.size() == 0) {
            return null;
        }
        return children.get(0);
    }

    public NetworkNode getRightChild() {
        if (children.size() <= 1) {
            return null;
        }
        return children.get(1);
    }

    public void setLeftChild(final NetworkNode leftChild) {
        if (children.size() == 0) {
            children.add(leftChild);
        } else {
            children.set(0, leftChild);
        }
    }

    public void setRightChild(final NetworkNode rightChild) {
        switch (children.size()) {
            case 0:
                children.add(null);
            case 1:
                children.add(rightChild);
                break;
            default:
                children.set(1, rightChild);
                break;
        }
    }

    /**
     * @return true if current node is root node
     */
    public boolean isRoot() {
        return parents.size() == 0;
    }
    /**
     * @return true if current node is leaf node
     */
    public boolean isLeaf() {
        return children.size() == 0;
    }

    public boolean isVisited() {return visited;}

    public void setVisited () {visited = true;}

    public void resetVisited () {visited = false;}

    /* reset the visited indicators */
    public void resetAllVisited () {
        for (final NetworkNode child : children) {
            child.resetVisited();
        }
        resetVisited ();
    }

    public int getNodeCount() {
        resetAllVisited ();
        return recurseNodeCount();
    }

    public int recurseNodeCount() {
        int nodes = 1;
        setVisited();
        for (final NetworkNode child : children) {
            if (!child.isVisited())
                nodes += child.recurseNodeCount();
        }
        return nodes;
    }

    public int getLeafNodeCount() {
        resetAllVisited ();
        return recurseLeafNodeCount();
    }

    public int recurseLeafNodeCount() {
        if (isLeaf()) return 1;
        int nodes = 0;
        setVisited();
        for (final NetworkNode child : children) {
            if (!child.isVisited())
                nodes += child.recurseLeafNodeCount();
        }
        return nodes;
    }

    public int getInternalNodeCount() {
        resetAllVisited ();
        return recurseInternalNodeCount();
    }

    public int recurseInternalNodeCount() {
        if (isLeaf()) return 0;
        int nodes = 1;
        setVisited();
        for (final NetworkNode child : children) {
            if (!child.isVisited())
                nodes += child.recurseInternalNodeCount();
        }
        return nodes;
    }

    /**
     * many other methods below
     */

}
