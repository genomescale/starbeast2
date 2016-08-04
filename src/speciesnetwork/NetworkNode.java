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
    protected double inheritProb = 0.5;

    /**
     * height of this node
     */
    protected double height = Double.MAX_VALUE;

    /**
     * children and branch numbers of the children of this node
     */
    Set<NetworkNode> children;
    int leftChildBranchNumber = -1;
    int rightChildBranchNumber = -1;

    /**
     * branch numbers of the parents of this node
     */
    Set<NetworkNode> parents;
    int gammaParentBranchNumber = -1;
    int otherParentBranchNumber = -1;

    protected int nParents;
    protected int nChildren;

    // a set of labels
    Set<String> labels = new TreeSet<>();

    protected NetworkNode clone;

    /**
     * status of this node after an operation is performed on the state
     */
    int isDirty = Network.IS_CLEAN;

    protected void updateRelationships() {
        parents = getParents();
        children = getChildren();
        nParents = parents.size();
        nChildren = children.size();
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
        if (nParents < 2) // not a reticulation node
            throw new RuntimeException();
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
        for (NetworkNode c: children) {
            c.isDirty |= Network.IS_DIRTY;
        }
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

    public Set<NetworkNode> getParents() {
        final Set<NetworkNode> parents = new HashSet<>();

        if (gammaParentBranchNumber >= 0) {
            final NetworkNode gammaParent = network.getNode(gammaParentBranchNumber / 2);
            parents.add(gammaParent);
        }

        if (otherParentBranchNumber >= 0) {
            final NetworkNode otherParent = network.getNode(otherParentBranchNumber / 2);
            parents.add(otherParent);
        }

        return parents;
    }

    public Set<NetworkNode> getChildren() {
        final Set<NetworkNode> children = new HashSet<>();

        if (leftChildBranchNumber >= 0) {
            final NetworkNode leftChild = network.getNode(leftChildBranchNumber / 2);
            children.add(leftChild);
        }

        if (rightChildBranchNumber >= 0) {
            final NetworkNode rightChild = network.getNode(rightChildBranchNumber / 2);
            children.add(rightChild);
        }

        return children;
    }

    public void removeParent(final NetworkNode pNode) {
        final int pNodeNumber = pNode.labelNr;

        if ((gammaParentBranchNumber / 2) == pNodeNumber) {
            gammaParentBranchNumber = -1;
        } else if ((otherParentBranchNumber / 2) == pNodeNumber) {
            otherParentBranchNumber = -1;
        } else {
            throw new RuntimeException("Node is not a parent that can be removed.");
        }

        isDirty = Network.IS_FILTHY;
        updateRelationships();
    }

    public void addParent(final NetworkNode pNode) {
        final int pNodeNumber = pNode.labelNr;

        if (gammaParentBranchNumber == -1) {
            gammaParentBranchNumber = pNodeNumber * 2;
        } else if (otherParentBranchNumber == -1) {
            otherParentBranchNumber = pNodeNumber * 2;
        } else {
            throw new RuntimeException("This node already has two parents.");
        }

        isDirty = Network.IS_FILTHY;
        updateRelationships();
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

    // whether the node has been visited, say by a recursive method
    protected boolean visited = false; // this should be used in a public method

    protected boolean touched = false; // this is only used inside this class

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

    // returns total node count (leaf, internal including root) of subtree defined by this node
    public int getNodeCount() {
        network.resetAllTouched();
        return recurseNodeCount();
    }

    private int recurseNodeCount() {
        if (touched) return 0;

        int nodeCount = 1;
        for (NetworkNode c: children) {
            nodeCount += c.recurseNodeCount();
        }

        touched = true;
        return nodeCount;
    }

    public int getLeafNodeCount() {
        network.resetAllTouched();
        return recurseLeafNodeCount();
    }

    private int recurseLeafNodeCount() {
        if (touched)
            return 0;
        else if (nChildren == 0)
            return 1;

        int nodeCount = 0;
        for (NetworkNode c: children) {
            nodeCount += c.recurseLeafNodeCount();
        }

        touched = true;
        return nodeCount;
    }

    public int getSpeciationNodeCount() {
        network.resetAllTouched();
        return recurseSpeciationNodeCount();
    }

    private int recurseSpeciationNodeCount() {
        if (touched) return 0;

        // don't count reticulation nodes
        int nodeCount = (nChildren == 2) ? 1 : 0;
        for (NetworkNode c: children) {
            nodeCount += c.recurseSpeciationNodeCount();
        }

        touched = true;
        return nodeCount;
    }

    public int getReticulationNodeCount() {
        network.resetAllTouched();
        return recurseReticulationNodeCount();
    }

    private int recurseReticulationNodeCount() {
        if (touched) return 0;

        // only count reticulation nodes
        int nodeCount = (nParents == 2) ? 1 : 0;
        for (NetworkNode c: children) {
            nodeCount += c.recurseReticulationNodeCount();
        }

        touched = true;
        return nodeCount;
    }

    /**
     * get all child node under this node, if this node is leaf then list.size() = 0.
     */
    public List<NetworkNode> getAllChildNodes() {
        final List<NetworkNode> childNodes = new ArrayList<>();
        network.resetAllTouched();
        getAllChildNodes(childNodes);
        return childNodes;
    }
    // recursive
    private void getAllChildNodes(final List<NetworkNode> childNodes) {
        if (touched) return;

        childNodes.add(this);
        for (NetworkNode c: children) {
            c.getAllChildNodes();
        }

        touched = true;
    }

    /**
     * scale height of this node and all its descendants
     * @param scale scale factor
     */
    public void scale (final double scale) {
        network.resetAllTouched();
        recursiveScale(scale);
    }

    private void recursiveScale(final double scale) {
        if (touched) return;

        startEditing();
        height *= scale;
        isDirty |= Network.IS_DIRTY;
        for (NetworkNode c: children) c.recursiveScale(scale);
        touched = true;
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

    @Override
    public String toString() {
        return buildNewick(Double.POSITIVE_INFINITY, -1);
    }

    public String buildNewick(Double parentHeight, int parentBranchNumber) {
        StringBuilder subtreeString = new StringBuilder();
        // only add children to the gamma parent attached hybrid node
        // this is an arbitrary choice
        if (nParents != 2 || parentBranchNumber == gammaParentBranchNumber) {
            String leftSubtreeString = null;
            String rightSubtreeString = null;
            if (leftChildBranchNumber >= 0) {
                final NetworkNode leftChild = network.networkNodes[leftChildBranchNumber / 2];
                leftSubtreeString = leftChild.buildNewick(height, leftChildBranchNumber);
            }
            if (rightChildBranchNumber >= 0) {
                final NetworkNode rightChild = network.networkNodes[rightChildBranchNumber / 2];
                rightSubtreeString = rightChild.buildNewick(height, rightChildBranchNumber);
            }

            if (leftSubtreeString != null || rightSubtreeString != null) {
                subtreeString.append("(");
                if (leftSubtreeString != null) subtreeString.append(leftSubtreeString);
                if (leftSubtreeString != null && rightSubtreeString != null) subtreeString.append(",");
                if (rightSubtreeString != null) subtreeString.append(rightSubtreeString);
                subtreeString.append(")");
            }

            subtreeString.append(getID());
        } else if (nParents == 2) {
            subtreeString.append(getID());

            if (parentBranchNumber == gammaParentBranchNumber) {
                final String gammaStr = String.format("[&gamma=%.8g]", inheritProb);
                subtreeString.append(gammaStr);
            } else if (parentBranchNumber == otherParentBranchNumber) {
                final String gammaStr = String.format("[&gamma=%.8g]", 1.0 - inheritProb);
                subtreeString.append(gammaStr);
            } else {
                throw new RuntimeException("blah");
            }
        } else {
            subtreeString.append(getID());
        }

        if (parentHeight < Double.POSITIVE_INFINITY) {
            String branchLengthSuffix = String.format(":%.8g", parentHeight - height);
            while (branchLengthSuffix.charAt(branchLengthSuffix.length() - 1) == '0')
                branchLengthSuffix = branchLengthSuffix.substring(0, branchLengthSuffix.length() - 1);
            subtreeString.append(branchLengthSuffix);
        }
        return subtreeString.toString();
    }

    public Double getGamma() {
        return inheritProb;
    }

    public void setGamma(final Double newGamma) {
        startEditing();
        inheritProb = newGamma;
        isDirty |= Network.IS_DIRTY;
    }
}
