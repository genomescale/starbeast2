package speciesnetwork;

import java.io.PrintStream;
import java.util.List;
import java.util.Set;
import java.util.HashSet;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

/**
 * Network class to replace Tree class
 * It includes tip (in-degree 1/out-degree 0), bifurcation (1/2), and reticulation (2/1) nodes.
 * @author Chi Zhang
 * @author Huw Ogilvie
 */

@Description("Network representing reticulate evolution of species")
public class Network extends StateNode {
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
    protected int storedNodeCount = -1;
    protected int speciationNodeCount = -1;
    protected int storedSpeciationNodeCount = -1;
    protected int leafNodeCount = -1;
    protected int storedLeafNodeCount = -1;
    protected int reticulationNodeCount = -1;
    protected int storedReticulationNodeCount = -1;

    /**
     * array of all nodes in the network
     */
    protected NetworkNode[] nodes = null;
    protected NetworkNode[] storedNodes = null;

    @Override
    public void initAndValidate() {
        if (nodeCount < 0) {
            if (taxonSetInput.get() != null) {
                makeCaterpillar(0, 1);
                updateRelationships();
            } else {
                // make dummy network with a single root node
                nodes = new NetworkNode[1];
                nodes[1] = new NetworkNode(this);
                nodeCount = leafNodeCount = 1;
                speciationNodeCount = reticulationNodeCount = 0;
            }
        }
    }

    public void updateRelationships() {
        for (NetworkNode n: nodes) {
            n.updateRelationships();
        }
    }

    @Override
    public int scale(final double scale) {
        for (NetworkNode n: nodes) {
            n.height = n.height * scale;
        }

        return speciationNodeCount + reticulationNodeCount;
    }

    private void makeCaterpillar(final double minInternalHeight, final double step) {
        // make a caterpillar
        final List<String> taxa = taxonSetInput.get().asStringList();
        leafNodeCount = taxa.size();
        speciationNodeCount = leafNodeCount - 1;
        reticulationNodeCount = 0;
        nodeCount = leafNodeCount + speciationNodeCount;
        nodes = new NetworkNode[nodeCount];

        int leftNr = 0;
        nodes[leftNr] = new NetworkNode(this);
        NetworkNode left = nodes[leftNr];
        left.label = taxa.get(leftNr);
        for (int rightNr = 1; rightNr < leafNodeCount; rightNr++) {
            nodes[rightNr] = new NetworkNode(this);
            final NetworkNode right = nodes[rightNr];
            right.height = 0.0;
            right.label = taxa.get(rightNr);
            final int parentNr = leafNodeCount + (rightNr - 1);
            nodes[parentNr] = new NetworkNode(this);
            final NetworkNode parent = nodes[parentNr];
            parent.height = minInternalHeight + rightNr * step;
            parent.childBranchNumbers.add(rightNr);
            parent.childBranchNumbers.add(leftNr);
            left = parent;
            leftNr = parentNr;
        }
    }

    public NetworkNode getRoot() {
        return nodes[nodeCount - 1];
    }

    public void swapRoot(final int replacementNodeNumber) {
        final int rootNodeNumber = nodeCount - 1;
        swapNodes(replacementNodeNumber, rootNodeNumber);
    }

    public void swapNodes(final int nodeI, final int nodeJ) {
        final NetworkNode tmp = nodes[nodeI];
        nodes[nodeI] = nodes[nodeJ];
        nodes[nodeJ] = tmp;
    }

    /**
     * @return the number of nodes
     */
    public int getNodeCount() {
        return nodeCount;
    }

    public int getSpeciationNodeCount() {
        return speciationNodeCount;
    }

    public int getLeafNodeCount() {
        return leafNodeCount;
    }

    public int getReticulationNodeCount() {
        return reticulationNodeCount;
    }

    /**
     * @return the index of the first reticulation node
     */
    public int getReticulationOffset() {
        return (leafNodeCount + speciationNodeCount) - 1;
    }

    /**
     * @return get the total number of branches in the tree
     */
    public int getBranchCount() {
        return nodeCount + reticulationNodeCount;
    }

    /**
     * @return the number of branches at the given time
     */
    public int getBranchCount(double time) {
        int nB = 1;
        for (NetworkNode node : nodes) {
            if (node.getHeight() > time) {
                if (node.isReticulation()) nB--;
                else nB++;
            }
        }
        return nB;
    }

    /**
     * @return a list of leaf nodes contained in this network
     */
    public Set<NetworkNode> getLeafNodes() {
        final Set<NetworkNode> lNodes = new HashSet<>();
        for (int i = 0; i < leafNodeCount; i++) {
            lNodes.add(nodes[i]);
        }
        return lNodes;
    }

    /**
     * @return a list of internal nodes contained in this network
     */
    public Set<NetworkNode> getInternalNodes() {
        final Set<NetworkNode> iNodes = new HashSet<>();
        for (int i = leafNodeCount; i < nodeCount; i++) {
            iNodes.add(nodes[i]);
        }
        return iNodes;
    }

    /**
     * @return a list of reticulation nodes contained in this network
     */
    public Set<NetworkNode> getReticulationNodes() {
        final int reticulationOffset = getReticulationOffset();
        final Set<NetworkNode> rNodes = new HashSet<>();
        for (int i = 0; i < reticulationNodeCount; i++) {
            rNodes.add(nodes[i + reticulationOffset]);
        }
        return rNodes;
    }

    public NetworkNode getNode(final int nodeI) {
        return nodes[nodeI];
    }

    public int getNodeNumber(final String query) {
        for (int i = 0; i < nodeCount; i++) {
            final NetworkNode n = nodes[i];
            if (n != null && n.label != null) {
                if (n.label.equals(query)) return i;
            }
        }

        return -1; // no match
    }

    /**
     * @return an array of all the nodes in this network
     */
    public NetworkNode[] getNodes() {
        final NetworkNode[] nodesCopy = new NetworkNode[nodeCount];
        System.arraycopy(nodes, 0, nodesCopy, 0, nodeCount);
        return nodesCopy;
    }

    public double getNetworkLength() {
        return getRoot().getSubnetworkLength();
    }

    public String toString() {
        return getRoot().toString();
    }

    @Override
    public void setEverythingDirty(final boolean isDirty) {
        setSomethingIsDirty(isDirty);
        if (!isDirty) {
            for(NetworkNode node : nodes) {
                node.isDirty = IS_CLEAN;
            }
        } else {
            for(NetworkNode node : nodes) {
                node.isDirty = IS_FILTHY;
            }
        }
    }

    /**
     * @return a deep copy of this network
     */
    @Override
    public Network copy() {
        Network copy = new Network();
        copy.index = index;
        copyNetwork(this, copy);
        return copy;
    }
    
    /**
     * copy of all values into existing network
     */
    @Override
    public void assignTo(final StateNode other) {
        final Network dst = (Network) other;
        copyNetwork(this, dst);
    }

    /**
     * copy of all values from existing network
     */
    @Override
    public void assignFrom(final StateNode other) {
        final Network src = (Network) other;
        copyNetwork(src, this);
    }

    protected static void copyNetwork(Network src, Network dst) {
        final int copyNodeCount = src.nodeCount;

        dst.setID(src.getID());
        dst.nodeCount = copyNodeCount;
        dst.speciationNodeCount = src.speciationNodeCount;
        dst.leafNodeCount = src.leafNodeCount;
        dst.reticulationNodeCount = src.reticulationNodeCount;

        dst.nodes = new NetworkNode[copyNodeCount];
        for (int i = 0; i < copyNodeCount; i++) {
            dst.nodes[i] = new NetworkNode();
            dst.nodes[i].copyFrom(src.nodes[i]);
        }

        dst.updateRelationships();
    }

    /**
     * as assignFrom, but assumes this network has been initialized
     * with the same dimensions as the source network
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        final Network src = (Network) other;

        for (int i = 0; i < nodeCount; i++) {
            nodes[i].copyFrom(src.nodes[i]);
        }

        updateRelationships();
    }

    /**
     * reconstruct tree from XML fragment in the form of a DOM node
     */
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        final String newick = node.getTextContent();
        final TreeParser parser = new TreeParser();
        try {
            parser.thresholdInput.setValue(1e-10, parser);
        } catch (Exception e1) {
            e1.printStackTrace();
        }
        try {
            parser.offsetInput.setValue(0, parser);
            parser.parseNewick(newick); // TODO covert to network
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    protected void store() {
        storedNodeCount = nodeCount;
        storedSpeciationNodeCount = speciationNodeCount;
        storedLeafNodeCount = leafNodeCount;
        storedReticulationNodeCount = reticulationNodeCount;
        storedNodes = nodes;

        nodes = new NetworkNode[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            nodes[i] = new NetworkNode();
            nodes[i].copyFrom(storedNodes[i]);
        }

        updateRelationships();
    }

    @Override
    public void restore() {
        int tmpNodeCount = nodeCount;
        nodeCount = storedNodeCount;
        storedNodeCount = tmpNodeCount;

        int tmpSpeciationNodeCount = speciationNodeCount;
        speciationNodeCount = storedSpeciationNodeCount;
        storedSpeciationNodeCount = tmpSpeciationNodeCount;

        int tmpLeafNodeCount = leafNodeCount;
        leafNodeCount = storedLeafNodeCount;
        storedLeafNodeCount = tmpLeafNodeCount;

        int tmpReticulationNodeCount = reticulationNodeCount;
        reticulationNodeCount = storedReticulationNodeCount;
        storedReticulationNodeCount = tmpReticulationNodeCount;

        NetworkNode[] tmpNodes = nodes;
        nodes = storedNodes;
        storedNodes = tmpNodes;
    }

    /** Loggable interface implementation follows **/

    @Override
    public void init(PrintStream out) {
        out.println("#NEXUS\n");
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + leafNodeCount + ";");
        out.println("\t\tTaxlabels");
        for (int i = 0; i < leafNodeCount; i++) {
            out.println("\t\t\t" + nodes[i].label);
        }
        out.println("\t\t\t;");
        out.println("End;");
        out.println("Begin trees;");
    }

    @Override
    public void log(int sample, PrintStream out) {
        Network network = (Network) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        final String newick = network.getRoot().toString();
        out.print(newick);
        out.print(";");
    }

    /**
     * @see beast.core.Loggable *
     */
    @Override
    public void close(PrintStream out) {
        out.print("End;");
    }

    @Override
    public int getDimension() {
        return nodeCount;
    }

    @Override
    public double getArrayValue() {
        return getRoot().height;
    }

    @Override
    public double getArrayValue(final int nodeI) {
        return nodes[nodeI].height;
    }

    public boolean isDirty() {
        for (NetworkNode n: nodes) {
            if (n.isDirty != IS_CLEAN) return true;
        }
        return false;
    }

    protected void resetAllTouched() {
        for (NetworkNode n: nodes) {
            n.touched = false;
        }
    }

    protected void resetAllVisited() {
        for (NetworkNode n: nodes) {
            n.visited = false;
        }
    }

    /**
     * Returns branch number that corresponds to a node number
     */
    public int getBranchNumber(final int nodeNumber) {
        final int reticulationOffset = getReticulationOffset();
        if (nodeNumber < reticulationOffset) {
            return nodeNumber;
        } else {
            return (nodeNumber * 2) - reticulationOffset;
        }
    }

    /**
     * Vice versa
     */
    public int getNodeNumber(final int branchNumber) {
        final int reticulationOffset = getReticulationOffset();
        if (branchNumber < reticulationOffset) {
            return branchNumber;
        } else {
            return ((branchNumber - reticulationOffset) / 2) + reticulationOffset;
        }
    }
}
