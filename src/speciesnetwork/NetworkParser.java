package speciesnetwork;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.core.Input.Validate;

/**
 * Parse the network of extended Newick format.
 * @author Huw Ogilvie
 */

@Description("Parse the network of extended Newick format.")
public class NetworkParser extends Network implements StateNodeInitialiser {
    final public Input<Network> networkInput = new Input<>("initial", "Network to initialize.");
    final public Input<Tree> treeInput = new Input<>("tree", "Tree initialized from extended newick string", Validate.REQUIRED);

    private int nextLeafNr;
    private int nextSpeciationNr;
    private int nextReticulationNr;

    @Override
    public void initAndValidate() {
        final Tree tree = treeInput.get();
        final Node treeRoot = tree.getRoot();

        // Step (1) is to initialize the node counts and array
        leafNodeCount = 0;
        speciationNodeCount = 0;
        int hybridNodeCount = 0;

        for (Node n: tree.getNodesAsArray()) {
            if (n.getID() != null && n.getID().startsWith("#H")) {
                hybridNodeCount++;
            } else {
                if (n.isLeaf()) {
                    leafNodeCount++;
                } else {
                    speciationNodeCount++;
                }
            }
        }

        assert hybridNodeCount % 2 == 0;
        reticulationNodeCount = hybridNodeCount / 2;
        nodeCount = leafNodeCount + speciationNodeCount + reticulationNodeCount;
        nodes = new NetworkNode[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            nodes[i] = new NetworkNode(this); 
        }

        nextLeafNr = 0;
        nextSpeciationNr = leafNodeCount;
        nextReticulationNr = (leafNodeCount + speciationNodeCount) - 1;

        // Step (2) is to recursively copy the tree to the network
        rebuildNetwork(treeRoot);

        // Update the cached parents and children for each node
        updateRelationships();

        super.initAndValidate();
    }

    private Integer rebuildNetwork(final Node treeNode) {
        int branchNumber;
        NetworkNode newNode;

        final String nodeLabel = treeNode.getID();
        final double nodeHeight = treeNode.getHeight();

        final int matchingNodeNr = getNodeNumber(nodeLabel);
        if (matchingNodeNr < 0) {
            int newNodeNumber;
            double inheritProb = 0.5;

            if (treeNode.isRoot()) {
                newNodeNumber = nodeCount - 1;
            } else if (nodeLabel.startsWith("#H")) {
                if (treeNode.getMetaDataNames().contains("gamma"))
                    inheritProb = (Double) treeNode.getMetaData("gamma");
                newNodeNumber = nextReticulationNr;
                nextReticulationNr++;
            } else if (treeNode.isLeaf()) {
                newNodeNumber = nextLeafNr;
                nextLeafNr++;
            } else {
                newNodeNumber = nextSpeciationNr;
                nextSpeciationNr++;
            }

            newNode = nodes[newNodeNumber];
            newNode.label = nodeLabel;
            newNode.height = nodeHeight;
            newNode.inheritProb = inheritProb;

            branchNumber = getBranchNumber(newNodeNumber);
        } else {
            newNode = nodes[matchingNodeNr];
            branchNumber = getBranchNumber(matchingNodeNr) + 1;
        }

        for (Node c: treeNode.getChildren()) {
            final int childBranchNumber = rebuildNetwork(c);
            newNode.childBranchNumbers.add(childBranchNumber);
        }

        return branchNumber;
    }

    @Override
    public void initStateNodes() {
        if (networkInput.get() != null) {
            networkInput.get().assignFrom(this);
        }
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
        if (networkInput.get() != null) {
            stateNodes.add(networkInput.get());
        }
    }
}
