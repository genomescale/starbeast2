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
     * use even number for the left parent (0, 2, ...), odd number (+1) for the right (1, 3, ...)
     */
    protected int labelNr;

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
    protected Map<String, Object> metaData = new TreeMap<String, Object>();

    /**
     * the species network that this node is a part of
     */
    protected SpeciesNetwork speciesNetwork;

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

    public SpeciesNetwork getSpeciesNetwork() {
        return speciesNetwork;
    }

    public int getNr() {
        return labelNr;
    }

    public void setNr(final int iLabel) {
        labelNr = iLabel;
    }

    public double getHeight() {
        return height;
    }

    /**
     * A Node is IS_DIRTY if its value (like height) has changed.
     * A Node is IS_FILTHY if its parent or child has changed.
     * Otherwise the node is IS_CLEAN.
     */
    public int isDirty() {
        return isDirty;
    }

    public void makeDirty(final int nDirty) {
        isDirty |= nDirty;
    }

    /* This recursion needs revision: duplicated visit */
    public void makeAllDirty(final int nDirty) {
        isDirty = nDirty;
        if (!isLeaf()) {
            getLeftChild().makeAllDirty(nDirty);
            if (getRightChild() != null) {
                getRightChild().makeAllDirty(nDirty);
            }
        }
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

    /**
     * get all leaf node under this node, if this node is leaf then list.size() = 0.
     *
     * @return
     */
    public List<NetworkNode> getAllLeafNodes() {
        final List<NetworkNode> leafNodes = new ArrayList<>();
        if (!this.isLeaf()) getAllLeafNodes(leafNodes);
        return leafNodes;
    }

    // recursive
    public void getAllLeafNodes(final List<NetworkNode> leafNodes) {
        if (this.isLeaf()) {
            leafNodes.add(this);
        }
        for (NetworkNode child : children)
            child.getAllLeafNodes(leafNodes);
    }

    public int getNodeCount() {
        return getLeafNodeCount() + getBifurcationNodeCount() + getReticulationNodeCount();
    }

    /**
     * many other methods below
     */
    /*


    public String toString() {
    }

     */

}
