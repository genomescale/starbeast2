package starbeast2;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

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

    /**
     * height of this node
     */
    protected double height = Double.MAX_VALUE;

    /**
     * children of this node
     */
    NetworkNode leftChild;
    NetworkNode rightChild;

    /**
     * parents of this node
     */
    NetworkNode leftParent;
    NetworkNode rightParent;

    protected int nParents;
    protected int nChildren;

    private NetworkNode clone;

    // whether the node has been visited, say by a recursive method
    protected boolean visited = false;
    /**
     * status of this node after an operation is performed on the state
     */
    int isDirty = Network.IS_CLEAN;

    private void updateSizes() {
        nParents = 0;
        nChildren = 0;
        if (leftParent  != null) nParents++;
        if (rightParent != null) nParents++;
        if (leftChild   != null) nChildren++;
        if (rightChild  != null) nChildren++;
    }

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

    protected void startEditing() {
        if (network != null && network.getState() != null) {
            network.startEditing(null);
        }
    }

    public int getParentCount() {
        return nParents;
    }

    public int getChildCount() {
        return nChildren;
    }

    /* parents */
    public NetworkNode getLeftParent() {
        return leftParent;
    }

    public NetworkNode getRightParent() {
        return rightParent;
    }

    void setLeftParent(final NetworkNode newLeftParent) {
        startEditing();
        leftParent = newLeftParent;
        leftParent.leftChild = this;
        isDirty = Network.IS_FILTHY;
        updateSizes();
    }

    void setRightParent(final NetworkNode newRightParent) {
        startEditing();
        rightParent = newRightParent;
        rightParent.rightChild = this;
        isDirty = Network.IS_FILTHY;
        updateSizes();
    }

    /* children */
    public NetworkNode getLeftChild() {
        return leftChild;
    }

    public NetworkNode getRightChild() {
        return rightChild;
    }

    public void setLeftChild(final NetworkNode newLeftChild) {
        startEditing();
        leftChild = newLeftChild;
        leftChild.leftParent = this;
        updateSizes();
    }

    public void setRightChild(final NetworkNode newRightChild) {
        startEditing();
        rightChild = newRightChild;
        rightChild.rightParent = this;
        updateSizes();
    }

    /**
     * @return unmodifiable list of children of this node
     */
    public List<NetworkNode> getChildren() {
        List<NetworkNode> children = new ArrayList<>();
        if (leftChild != null) children.add(leftChild);
        if (rightChild != null) children.add(rightChild);
        return children;
    }

    public List<NetworkNode> getParents() {
        List<NetworkNode> parents = new ArrayList<>();
        if (leftParent != null) parents.add(leftParent);
        if (rightParent != null) parents.add(rightParent);
        return parents;
    }

    /**
     * @return true if current node is root node
     */
    public boolean isRoot() {
        return nParents == 0;
    }
    /**
     * @return true if current node is leaf node
     */
    public boolean isLeaf() {
        return nChildren == 0;
    }

    /**
     * @return true if current node is reticulation node
     */
    boolean isReticulation() {
        return nParents == 2;
    }

    /* get and (re)set the visited indicator */
    public boolean isVisited() {
        return visited;
    }
    public void setVisited() {
        visited = true;
    }
    public void resetVisited() {
        visited = false;
    }

    /* reset all the visited indicators */
    public void recursiveResetVisited() {
        resetVisited ();
        if (leftChild != null) leftChild.recursiveResetVisited();
        if (rightChild != null) rightChild.recursiveResetVisited();
    }

    public int getNodeCount() {
        recursiveResetVisited();
        return recurseNodeCount();
    }
    private int recurseNodeCount() {
        if (visited) return 0;

        int nodeCount = 1;
        if (leftChild != null) nodeCount += leftChild.recurseNodeCount();
        if (rightChild != null) nodeCount += rightChild.recurseNodeCount();

        setVisited();
        return nodeCount;
    }

    public int getLeafNodeCount() {
        recursiveResetVisited();
        return recurseLeafNodeCount();
    }
    private int recurseLeafNodeCount() {
        if (visited) {
            return 0;
        } else if (nChildren == 0) {
            return 1;
        }

        int nodeCount = 0;
        if (leftChild != null) nodeCount += leftChild.recurseLeafNodeCount();
        if (rightChild != null) nodeCount += rightChild.recurseLeafNodeCount();

        setVisited();
        return nodeCount;
    }

    public int getInternalNodeCount() {
        recursiveResetVisited();
        return recurseInternalNodeCount();
    }
    private int recurseInternalNodeCount() {
        if (visited || nChildren == 0) return 0;

        int nodeCount = 1;
        if (leftChild != null) nodeCount += leftChild.recurseInternalNodeCount();
        if (rightChild != null) nodeCount += rightChild.recurseInternalNodeCount();

        setVisited();
        return nodeCount;
    }

    /**
     * @return (deep) copy of node
     */
    public NetworkNode copy() {
        recursiveResetClones();
        return recursiveCopy();
    }

    private void recursiveResetClones() {
        clone = null;
        if (leftChild != null) leftChild.recursiveResetClones();
        if (rightChild != null) rightChild.recursiveResetClones();
    }

    private NetworkNode recursiveCopy() {
        if (clone == null) {
            final NetworkNode clone = new NetworkNode();
            clone.height = height;
            clone.labelNr = labelNr;
            clone.metaDataString = metaDataString;
            clone.metaData = new TreeMap<>(metaData);
            clone.setID(getID());
            if (leftChild != null) clone.setLeftChild(getLeftChild().recursiveCopy());
            if (rightChild != null) clone.setRightChild(getRightChild().recursiveCopy());
        }
        return clone;
    }

    /**
     * assign values to a network in array representation
     */
    public void assignTo (final NetworkNode[] nodes) {
        recursiveResetVisited();
        recursiveAssignTo(nodes);
    }
    private void recursiveAssignTo(final NetworkNode[] nodes) {
        if (!visited) {
            final NetworkNode node = nodes[getNr()];
            node.height = height;
            node.labelNr = labelNr;
            node.metaDataString = metaDataString;
            node.metaData = new TreeMap<>(metaData);
            node.setID(getID());
            setVisited();
            if (leftChild != null) {
                node.setLeftChild(nodes[getLeftChild().getNr()]);
                getLeftChild().recursiveAssignTo(nodes);
            }
            if (rightChild != null) {
                node.setRightChild(nodes[getRightChild().getNr()]);
                getRightChild().recursiveAssignTo(nodes);
            }
        }
    }

    /**
     * assign values from a network in array representation
     */
    public void assignFrom (final NetworkNode[] nodes, final NetworkNode node) {
        recursiveResetVisited();
        recursiveAssignFrom(nodes, node);
    }
    private void recursiveAssignFrom(final NetworkNode[] nodes, final NetworkNode node) {
        if (!visited) {
            height = node.height;
            labelNr = node.labelNr;
            metaDataString = node.metaDataString;
            metaData = new TreeMap<>(node.metaData);
            setID(node.getID());
            setVisited();
            if (node.getLeftChild() != null) {
                setLeftChild(nodes[node.getLeftChild().getNr()]);
                getLeftChild().recursiveAssignFrom(nodes, node.getLeftChild());
            }
            if (node.getRightChild() != null) {
                setRightChild(nodes[node.getRightChild().getNr()]);
                getRightChild().recursiveAssignFrom(nodes, node.getRightChild());
            }
        }
    }

    /**
     * scale height of this node and all its descendants
     * @param scale scale factor
     */
    public void scale (final double scale) throws Exception {
        recursiveResetVisited();
        recursiveScale(scale);
    }
    private void recursiveScale(final double scale) {
        startEditing();
        isDirty |= Network.IS_DIRTY;
        if (!isLeaf() && !visited) {
            height *= scale;
            setVisited();
            if (getLeftChild() != null)
                getLeftChild().recursiveScale(scale);
            if (getRightChild() != null)
                getRightChild().recursiveScale(scale);
            if (height < getLeftChild().height || height < getRightChild().height)
                throw new IllegalArgumentException("Scale gives negative branch length");
        }
    }

    /**
     * many other methods below
     */

}
