package speciesnetwork;

import java.util.*;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.evolution.tree.Node;

/**
 * NetworkNode that is like Node but has 2 parents and 2 children.
 * The right parent is null if it is a tree node.
 * The right child is null if it is a reticulation node.
 * @author Chi Zhang
 * @author Huw Ogilvie
 */

@Description("Network node for binary rooted network")
public class NetworkNode extends BEASTObject {

    /**
     * label nr of node
     */
    protected int labelNr;

    /**
     * inheritance probability associated with the left parent branch
     */
    protected double inheritProb;

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

    // a set of labels
    Set<String> labels = new TreeSet<>();

    private NetworkNode clone;

    // whether the node has been visited, say by a recursive method
    protected boolean visited = false;
    /**
     * status of this node after an operation is performed on the state
     */
    int isDirty = Network.IS_CLEAN;

    protected void updateSizes() {
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

    public NetworkNode(final String id) {
        setID(id);
        initAndValidate();
    }

    // instantiate a new network node with the same height, labels and metadata as a tree node
    // this does not copy the parents or children
    public NetworkNode(Node treeNode) {
        setID(treeNode.getID());
        height = treeNode.getHeight();
        metaDataString = treeNode.metaDataString;
        for (String metaDataKey: treeNode.getMetaDataNames()) {
            Object metaDataValue = treeNode.getMetaData(metaDataKey);
            metaData.put(metaDataKey, metaDataValue);
        }

        initAndValidate();
    }

    @Override
    public void initAndValidate() {
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
    
    public int getReticulationNumber() {
        if (leftParent == null || rightParent == null) throw new RuntimeException();
        final int reticulationNodeOffset = network.getNodeCount() - network.getReticulationNodeCount() - 1;
        return labelNr - reticulationNodeOffset;
    }

    public double getHeight() {
        return height;
    }

    public void setHeight(final double height) {
        startEditing();
        this.height = height;
        isDirty |= Network.IS_DIRTY;
        if (getLeftChild() != null)
            getLeftChild().isDirty |= Network.IS_DIRTY;
        if (getRightChild() != null)
            getRightChild().isDirty |= Network.IS_DIRTY;
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
        isDirty = Network.IS_FILTHY;
        updateSizes();
    }

    public void setRightChild(final NetworkNode newRightChild) {
        startEditing();
        rightChild = newRightChild;
        rightChild.rightParent = this;
        isDirty = Network.IS_FILTHY;
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
    public boolean isReticulation() {
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

    // returns total node count (leaf, internal including root) of subtree defined by this node
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

    public int getSpeciationNodeCount() {
        recursiveResetVisited();
        return recurseSpeciationNodeCount();
    }
    private int recurseSpeciationNodeCount() {
        if (visited) return 0;

        // don't count reticulation nodes
        int nodeCount = (leftChild != null && rightChild != null) ? 1 : 0;
        if (leftChild != null) nodeCount += leftChild.recurseSpeciationNodeCount();
        if (rightChild != null) nodeCount += rightChild.recurseSpeciationNodeCount();

        setVisited();
        return nodeCount;
    }

    public int getReticulationNodeCount() {
        recursiveResetVisited();
        return recurseReticulationNodeCount();
    }
    private int recurseReticulationNodeCount() {
        if (visited) return 0;

        // only count reticulation nodes
        int nodeCount = (leftParent != null && rightParent != null) ? 1 : 0;
        if (leftChild != null) nodeCount += leftChild.recurseReticulationNodeCount();
        if (rightChild != null) nodeCount += rightChild.recurseReticulationNodeCount();

        setVisited();
        return nodeCount;
    }

    /**
     * get all child node under this node, if this node is leaf then list.size() = 0.
     */
    public List<NetworkNode> getAllChildNodes() {
        final List<NetworkNode> childNodes = new ArrayList<>();
        recursiveResetVisited();
        if (!this.isLeaf()) getAllChildNodes(childNodes);
        return childNodes;
    }
    // recursive
    public void getAllChildNodes(final List<NetworkNode> childNodes) {
        if (visited) return;

        childNodes.add(this);
        setVisited();
        if (leftChild != null) leftChild.getAllChildNodes(childNodes);
        if (rightChild != null) rightChild.getAllChildNodes(childNodes);
    }

    /**
     * @return length of branch between this node and its parent
     */
    public final double getLeftLength() {
        if (isRoot() || getLeftParent() == null) return 0;
        else return getLeftParent().height - height;
    }

    public final double getRightLength() {
        if (isRoot() || getRightParent() == null) return 0;
        else return getRightParent().height - height;
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
            clone = new NetworkNode();
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
    public void scale (final double scale) {
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
            //if (height < getLeftChild().height || height < getRightChild().height)
            //    throw new IllegalArgumentException("Scale gives negative branch length");
        }
    }

    // swap the orientation of parents (so the left parent becomes the right parent and vice versa)
    public void reorientParents() {
        // swap parents
        NetworkNode tmpParent = leftParent;
        leftParent = rightParent;
        rightParent = tmpParent;
        
        // update parent nodes to reflect changes
        // cascade changes to parent of half-siblings if they exist
        if (leftParent != null) {
            NetworkNode halfSibling = leftParent.leftChild;
            leftParent.rightChild = null;
            if (halfSibling != null) halfSibling.reorientParents();
            leftParent.leftChild = this;
        }

        if (rightParent != null) {
            NetworkNode halfSibling = rightParent.rightChild;
            rightParent.leftChild = null;
            if (halfSibling != null) halfSibling.reorientParents();
            rightParent.rightChild = this;
        }
    }

    public void addLabel(String label) {
        labels.add(label);
    }

    public void removeLabel(String label) {
        labels.remove(label);
    }

    public boolean hasLabel(String label) {
        return labels.contains(label);
    }


    /**
     * many other methods below
     */

    @Override
    public String toString() {
        // System.out.println(getID());
        return buildNewick(Double.POSITIVE_INFINITY, true);
    }

    public String buildNewick(Double parentHeight, boolean isLeft) {
        StringBuilder subtreeString = new StringBuilder();
        if (leftParent == null || rightParent == null || isLeft) { // only add children to left-attached reticulation nodes
            String leftSubtreeString = null;
            String rightSubtreeString = null;
            if (leftChild != null) leftSubtreeString = leftChild.buildNewick(height, true);
            if (rightChild != null) rightSubtreeString = rightChild.buildNewick(height, false);

            if (leftSubtreeString != null || rightSubtreeString != null) {
                subtreeString.append("(");
                if (leftSubtreeString != null) subtreeString.append(leftSubtreeString);
                if (leftSubtreeString != null && rightSubtreeString != null) subtreeString.append(",");
                if (rightSubtreeString != null) subtreeString.append(rightSubtreeString);
                subtreeString.append(")");
            }
        }
        subtreeString.append(getID());
        if (parentHeight < Double.POSITIVE_INFINITY) {
            String branchLengthSuffix = String.format(":%.8g", parentHeight - height);
            while (branchLengthSuffix.charAt(branchLengthSuffix.length() - 1) == '0')
                branchLengthSuffix = branchLengthSuffix.substring(0, branchLengthSuffix.length() - 1);
            subtreeString.append(branchLengthSuffix);
        }
        return subtreeString.toString();
    }

    /**
     * number the network leafs as 0, 1, ..., leafNodeCount-1
     * number the speciation nodes (except root) as leafNodeCount, ..., leafNodeCount+speciationNodeCount-1
     * number the reticulation nodes as leafNodeCount+speciationNodeCount, ..., nodeCount-2
     * number the root as nodeCount-1
     * @param traversalDirection 0 -> left, 1 -> right
     * @return branch number
     */
    public int getBranchNumber(final int traversalDirection) {
        // if this is a reticulation node
        if (leftParent != null && rightParent != null) {
            // this is the index of the first reticulation node
            final int reticulationNodeOffset = network.getNodeCount() - network.getReticulationNodeCount() - 1;
            final int reticulationNumber = labelNr - reticulationNodeOffset;  // 0, 1, 2, ...
            return reticulationNodeOffset + (reticulationNumber * 2) + traversalDirection;
        } else if (leftParent == null && rightParent == null) {
            // this is a root node
            return network.getBranchCount() - 1;
        } else {
            // leaf and non-root speciation nodes have equivalent branch and node numbers
            return labelNr;
        }
    }

    public int getLeftBranchNumber() {
        return getBranchNumber(0);
    }

    public int getRightBranchNumber() {
        return getBranchNumber(1);
    }
}

