package speciesnetwork;

import java.util.*;

import beast.core.Description;
import beast.evolution.tree.Node;

/**
 * NetworkNode is equivalent to Node but can have 2 parents or 2 children.
 * @author Chi Zhang
 * @author Huw Ogilvie
 */

@Description("Network node for binary rooted network")
public class NetworkNode {
    /**
     * the taxonomic name of this node
     */
    protected String label;

    /**
     * inheritance probability associated with the left parent branch
     */
    protected double inheritProb;

    /**
     * height of this node
     */
    protected double height;

    /**
     * children and parents of this node
     */
    protected Set<Integer> childBranchNumbers;
    protected Set<NetworkNode> children;
    protected Set<NetworkNode> parents;

    /**
     * counts of children and parents of this node
     */
    protected int nParents;
    protected int nChildren;

    /**
     * status of this node after an operation is performed on the state
     */
    int isDirty;

    protected void updateRelationships() {
        parents = getParents();
        children = getChildren();
        nParents = parents.size();
        nChildren = children.size();
    }

    /**
     * meta-data contained in square brackets in Newick
     */
    protected String metaDataString;

    /**
     * arbitrarily labeled metadata on this node
     */
    protected Map<String, Object> metaData = new TreeMap<>();

    /**
     * the network that this node is a part of
     */
    protected Network network;

    public NetworkNode() {
        label = "";
        inheritProb = 0.5;
        height = 0.0;
        childBranchNumbers = new HashSet<>();
        children = new HashSet<>();
        parents = new HashSet<>();
        nParents = 0;
        nChildren = 0;
        isDirty = Network.IS_DIRTY;
    }

    // instantiate a new network node with the same height, labels and metadata as a tree node
    // this does not copy the parents or children
    public NetworkNode(Node treeNode) {
        height = treeNode.getHeight();
        metaDataString = treeNode.metaDataString;
        for (String metaDataKey: treeNode.getMetaDataNames()) {
            Object metaDataValue = treeNode.getMetaData(metaDataKey);
            metaData.put(metaDataKey, metaDataValue);
        }
    }
    
    public void copyTo(NetworkNode dst) {
        dst.label = label;
        dst.inheritProb = inheritProb;
        dst.height = height;
        dst.childBranchNumbers.clear();
        dst.childBranchNumbers.addAll(childBranchNumbers);
        dst.children.clear();
        dst.children.addAll(children);
        dst.parents.clear();
        dst.parents.addAll(parents);
        dst.nParents = nParents;
        dst.nChildren = nChildren;
        dst.isDirty = isDirty;
    }

    public void copyFrom(NetworkNode src) {
        label = src.label;
        inheritProb = src.inheritProb;
        height = src.height;
        childBranchNumbers.clear();
        childBranchNumbers.addAll(src.childBranchNumbers);
        children.clear();
        children.addAll(src.children);
        parents.clear();
        parents.addAll(src.parents);
        nParents = src.nParents;
        nChildren = src.nChildren;
        isDirty = src.isDirty;
    }

    public Network getNetwork() {
        return network;
    }

    public int getNr() {
        for (int i = 0; i < network.networkNodes.length; i++) {
            if (network.networkNodes[i] == this) return i; 
        }

        throw new RuntimeException("Node is not attached to the network!");
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

        for (NetworkNode n: network.networkNodes) {
            if (n.children.contains(this)) {
                parents.add(n);
            }
        }

        return parents;
    }

    public Set<NetworkNode> getChildren() {
        final Set<NetworkNode> children = new HashSet<>();

        for (Integer i: childBranchNumbers) {
            final NetworkNode childNode = network.networkNodes[i];
            children.add(childNode);
        }

        return children;
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

    @Override
    public String toString() {
        return buildNewick(Double.POSITIVE_INFINITY, -1);
    }

    public String buildNewick(Double parentHeight, Integer parentBranchNumber) {
        final StringBuilder subtreeString = new StringBuilder();
        // only add children to the gamma parent attached reticulation node
        if ((nChildren > 0 && nParents < 2) || parentBranchNumber % 2 == 0) {
            subtreeString.append("(");
            int i = 0;
            for (Integer childBranchNumber: childBranchNumbers) {
                if (i > 0) subtreeString.append(",");
                NetworkNode childNode = network.networkNodes[childBranchNumber / 2];
                subtreeString.append(childNode.buildNewick(height, childBranchNumber));
                i++;
            }
            subtreeString.append(")");
        }

        subtreeString.append(label);

        // add inheritance probabilities to reticulation nodes
        if (nParents == 2) {
            subtreeString.append("[&gamma=");

            if (parentBranchNumber % 2 == 0) subtreeString.append(inheritProb);
            else subtreeString.append(1.0 - inheritProb);

            subtreeString.append("]");
        }

        if (parentHeight < Double.POSITIVE_INFINITY) {
            final double branchLength = parentHeight - height;
            subtreeString.append(":");
            subtreeString.append(branchLength);
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
