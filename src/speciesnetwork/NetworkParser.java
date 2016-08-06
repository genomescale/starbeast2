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
        final int reticulationNodeCount = hybridNodeCount / 2;
        nodeCount = leafNodeCount + speciationNodeCount + reticulationNodeCount;
        // System.out.println(String.format("%d = %d + %d + %d", nodeCount, leafNodeCount, speciationNodeCount, reticulationNodeCount));
        networkNodes = new NetworkNode[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            networkNodes[i] = new NetworkNode(); 
        }

        nextLeafNr = 0;
        nextSpeciationNr = leafNodeCount;
        nextReticulationNr = (leafNodeCount + speciationNodeCount) - 1;
        rebuildNetwork(treeRoot);
        updateRelationships();

        super.initAndValidate();
    }

    private Integer rebuildNetwork(final Node treeNode) {
        int branchNumber;
        int nodeNumber;
        NetworkNode newNode;

        final String nodeLabel = treeNode.getID();
        final double nodeHeight = treeNode.getHeight();

        final int matchingNodeNr = getNodeNr(nodeLabel);
        if (matchingNodeNr < 0) {
            if (treeNode.isRoot()) {
                nodeNumber = nodeCount -1;
            } else if (nodeLabel.startsWith("#H")) {
                nodeNumber = nextReticulationNr;
                nextReticulationNr++;
            } else if (treeNode.isLeaf()) {
                nodeNumber = nextLeafNr;
                nextLeafNr++;
            } else {
                nodeNumber = nextSpeciationNr;
                nextSpeciationNr++;
            }

            networkNodes[nodeNumber] = new NetworkNode(this);
            branchNumber = nodeNumber * 2;
        } else {
            nodeNumber = matchingNodeNr;
            branchNumber = (matchingNodeNr * 2) + 1;
        }

        newNode = networkNodes[nodeNumber];
        newNode.label = nodeLabel;
        newNode.height = nodeHeight;

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
