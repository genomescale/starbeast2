package speciesnetwork;

import java.io.PrintStream;
import java.util.List;
import java.util.Set;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

/**
 * Network class to replace Tree class
 * It includes both bifurcation node (in-degree 1, out-degree 2 or 0(tip)) and reticulation node (in degree 2, out-degree 1).
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

    /**
     * array of all nodes in the network
     */
    protected NetworkNode[] networkNodes = null;
    protected NetworkNode[] storedNetworkNodes = null;

    @Override
    public void initAndValidate() {
        if (nodeCount < 0) {
            if (taxonSetInput.get() != null) {
                makeCaterpillar(0, 1);
                updateRelationships();
            } else {
                // make dummy network with a single root node
                networkNodes = new NetworkNode[1];
                networkNodes[1] = new NetworkNode(this);
                nodeCount = leafNodeCount = 1;
                speciationNodeCount = 0;
            }
        }
    }

    public void updateRelationships() {
        for (NetworkNode n: networkNodes) {
            n.updateRelationships();
        }
    }

    private void makeCaterpillar(final double minInternalHeight, final double step) {
        // make a caterpillar
        final List<String> taxa = taxonSetInput.get().asStringList();
        leafNodeCount = taxa.size();
        speciationNodeCount = leafNodeCount - 1;
        nodeCount = leafNodeCount + speciationNodeCount;
        networkNodes = new NetworkNode[nodeCount];

        int leftNr = 0;
        networkNodes[leftNr] = new NetworkNode(this);
        NetworkNode left = networkNodes[leftNr];
        left.label = taxa.get(leftNr);
        for (int rightNr = 1; rightNr < leafNodeCount; rightNr++) {
            networkNodes[rightNr] = new NetworkNode(this);
            final NetworkNode right = networkNodes[rightNr];
            right.height = 0.0;
            right.label = taxa.get(rightNr);
            final int parentNr = leafNodeCount + (rightNr - 1);
            networkNodes[parentNr] = new NetworkNode(this);
            final NetworkNode parent = networkNodes[parentNr];
            parent.height = minInternalHeight + rightNr * step;
            parent.childBranchNumbers.add(rightNr * 2);
            parent.childBranchNumbers.add(leftNr * 2);
            left = parent;
            leftNr = parentNr;
        }
    }

    public NetworkNode getRoot() {
        return networkNodes[nodeCount - 1];
    }

    public void swapRoot(final int replacementNodeNumber) {
        final int rootNodeNumber = nodeCount -1;
        swapNodes(replacementNodeNumber, rootNodeNumber);
    }

    public void swapNodes(final int nodeI, final int nodeJ) {
        final NetworkNode tmp = networkNodes[nodeI];
        networkNodes[nodeI] = networkNodes[nodeJ];
        networkNodes[nodeJ] = tmp;
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
        return nodeCount - (leafNodeCount + speciationNodeCount);
    }

    // the index of the first reticulation node
    public int getReticulationNodeOffset() {
        return leafNodeCount + speciationNodeCount - 1;
    }

    /**
     * @return get the total number of branches in the tree
     */
    public int getBranchCount() {
        return (nodeCount * 2) - (leafNodeCount + speciationNodeCount);
    }

    /**
     * @return the number of branches at the given time
     */
    public int getBranchCount(double time) {
        int nB = 1;
        for (NetworkNode node : networkNodes) {
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
            lNodes.add(networkNodes[i]);
        }
        return lNodes;
    }

    /**
     * @return a list of internal nodes contained in this network
     */
    public Set<NetworkNode> getInternalNodes() {
        final Set<NetworkNode> iNodes = new HashSet<>();
        for (int i = leafNodeCount; i < nodeCount; i++) {
            iNodes.add(networkNodes[i]);
        }
        return iNodes;
    }

    /**
     * @return a list of reticulation nodes contained in this network
     */
    public Set<NetworkNode> getReticulationNodes() {
        final Set<NetworkNode> rNodes = new HashSet<>();
        for (int i = 0; i < nodeCount; i++) {
            final int reticulationNodeNumber = getReticulationNodeOffset() + i;
            rNodes.add(networkNodes[reticulationNodeNumber]);
        }
        return rNodes;
    }

    public NetworkNode getNode(final int nodeI) {
        return networkNodes[nodeI];
    }

    public int getNodeNr(final String query) {
        int matchingNodeNr = -1;
        for (int i = 0; i < nodeCount; i++) {
            final NetworkNode n = networkNodes[i];
            if (n != null && n.label != null) {
                if (n.label.equals(query)) return i;
            }
        }
        return matchingNodeNr;
    }

    /**
     * @return an array of all the nodes in this network
     */
    public NetworkNode[] getAllNodesAsArray() {
        final NetworkNode[] nodesCopy = new NetworkNode[nodeCount];
        System.arraycopy(networkNodes, 0, nodesCopy, 0, nodeCount);
        return nodesCopy;
    }

    public double getNetworkLength() {
        double length = 0.0;
        for (NetworkNode n: networkNodes) {
            for (NetworkNode c: n.children) {
                length += n.height - c.height;
            }
        }

        return length;
    }

    /**
     * @return an array of taxon names in order of their node numbers
     */
    public String[] getTaxaNames() {
        final String[] taxaNames = new String[leafNodeCount];
        for (int i = 0; i < leafNodeCount; i++) {
            taxaNames[i] = networkNodes[i].label;
        }

        return taxaNames;
    }

    public String toString() {
        return getRoot().toString();
    }

    @Override
    public void setEverythingDirty(final boolean isDirty) {
        setSomethingIsDirty(isDirty);
        if (!isDirty) {
            for(NetworkNode node : networkNodes) {
                node.isDirty = IS_CLEAN;
            }
        } else {
            for(NetworkNode node : networkNodes) {
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

        dst.networkNodes = new NetworkNode[copyNodeCount];
        for (int i = 0; i < copyNodeCount; i++) {
            dst.networkNodes[i] = new NetworkNode();
            dst.networkNodes[i].copyFrom(src.networkNodes[i]);
        }

        for (int i = 0; i < copyNodeCount; i++) {
            dst.networkNodes[i].updateRelationships();
        }
    }

    /**
     * as assignFrom, but only copy network structure
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        final Network src = (Network) other;
        for (int i = 0; i < nodeCount; i++) {
            networkNodes[i].copyFrom(src.networkNodes[i]);
        }
    }

    @Override
    public int scale(final double scale) {
        getRoot().scale(scale);
        return getSpeciationNodeCount() + getReticulationNodeCount();
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

        if (storedNetworkNodes == null || storedNodeCount != nodeCount) { // rebuild array
            storedNetworkNodes = new NetworkNode[nodeCount];
            for (int i = 0; i < nodeCount; i++) {
                storedNetworkNodes[i] = new NetworkNode();
                storedNetworkNodes[i].copyFrom(networkNodes[i]);
            }
        } else {
            for (int i = 0; i < nodeCount; i++) {
                storedNetworkNodes[i].copyFrom(networkNodes[i]);
            }
        }
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

        NetworkNode[] tmpNetworkNodes = networkNodes;
        networkNodes = storedNetworkNodes;
        storedNetworkNodes = tmpNetworkNodes;

        for (int i = 0; i < nodeCount; i++) {
            networkNodes[i].updateRelationships();
        }
    }

    /** Loggable interface implementation follows **/

    /**
     * print translate block for NEXUS beast.tree file
     */
    public static void printTranslate(final NetworkNode node, final PrintStream out, final int nodeCount) {
        final List<String> translateLines = new ArrayList<>();
        printTranslate(node, translateLines, nodeCount);
        Collections.sort(translateLines);
        for (final String line : translateLines) {
            out.println(line);
        }
    }

    static public int taxaTranslationOffset = 1;

    /**
     * need this helper so that we can sort list of entries *
     */
    static void printTranslate(NetworkNode node, List<String> translateLines, int nodeCount) {
        if (node.isLeaf()) {
            final String nr = (node.getNr() + taxaTranslationOffset) + "";
            String line = "\t\t" + "    ".substring(nr.length()) + nr + " " + node.label;
            if (node.getNr() < nodeCount) {
                line += ",";
            }
            translateLines.add(line);
        } else {
            for (NetworkNode c: node.children) {
                printTranslate(c, translateLines, nodeCount);
            }
        }
    }

    public static void printTaxa(final NetworkNode node, final PrintStream out, final int nodeCount) {
        final List<String> translateLines = new ArrayList<>();
        printTranslate(node, translateLines, nodeCount);
        Collections.sort(translateLines);
        for (String line : translateLines) {
            line = line.split("\\s+")[2];
            out.println("\t\t\t" + line.replace(',', ' '));
        }
    }

    @Override
    public void init(PrintStream out) {
        NetworkNode node = getRoot();
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
        Network network = (Network) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final String newick = network.getRoot().toNewick();
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
        return networkNodes[nodeI].height;
    }

    public void addReticulationNode(final int nodeI) {
        // TODO
    }

    public void removeReticulationNode(final int nodeI) {
        // TODO
    }

    public boolean isDirty() {
        for (NetworkNode n: networkNodes) {
            if (n.isDirty != IS_CLEAN) return true;
        }
        return false;
    }

    protected void resetAllTouched() {
        for (NetworkNode n: networkNodes) {
            n.touched = false;
        }
    }

    protected void resetAllVisited() {
        for (NetworkNode n: networkNodes) {
            n.visited = false;
        }
    }
}
