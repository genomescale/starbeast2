package starbeast2;

import java.util.List;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.core.Input.Validate;

public class NetworkParser extends Network implements StateNodeInitialiser {
    public final Input<TreeInterface> treeInput = new Input<>("tree", "Tree initialized from extended newick string", Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
        Node treeRoot = treeInput.get().getRoot();
        NetworkNode newRoot = new NetworkNode(treeRoot);
        setRoot(newRoot);
        speciationNodeCount = 1;
        leafNodeCount = reticulationNodeCount = 0;
        super.initAndValidate();

        rebuildNetwork(treeRoot.getLeft(), root, true);
        rebuildNetwork(treeRoot.getRight(), root, false);
        root.updateSizes();
    }

    private void rebuildNetwork(final Node treeNode, final NetworkNode parentNode, final boolean isLeft) throws Exception {
        final Node leftChild = treeNode.getLeft();
        final Node rightChild = treeNode.getRight();
        final String nodeLabel = treeNode.getID();
        boolean reticulation;
        NetworkNode networkNode;
        if (nodeLabel == null) {
            networkNode = null;
            reticulation = false;
        } else {
            networkNode = getNode(nodeLabel);
            final int hStart = nodeLabel.indexOf('#') + 1;
            reticulation = (hStart > 0) && (nodeLabel.length() > hStart) && (nodeLabel.charAt(hStart) == 'H');
        }

        if (networkNode == null) {
            networkNode = new NetworkNode(treeNode);
            if (reticulation) addReticulationNode(networkNode);
            else if (treeNode.isLeaf()) addLeafNode(networkNode);
            else addSpeciationNode(networkNode);
        }

        // complete reticulation nodes are always left-attached
        if (reticulation && leftChild == null && rightChild == null) networkNode.setLeftParent(parentNode);
        // incomplete (no children) reticulation nodes are always right-attached
        else if (reticulation) networkNode.setRightParent(parentNode);
        else if (isLeft) networkNode.setLeftParent(parentNode);
        else networkNode.setRightParent(parentNode);

        if (leftChild != null) rebuildNetwork(leftChild, networkNode, true);
        if (rightChild != null) rebuildNetwork(rightChild, networkNode, false);

        networkNode.updateSizes();
    }

    @Override
    public void initStateNodes() throws Exception {
        // TODO Auto-generated method stub

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        // TODO Auto-generated method stub

    }
}
