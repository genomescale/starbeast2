package speciesnetwork;

import java.util.List;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.core.Input.Validate;

/**
 * parse the network of extended newick format
 * @author Huw Ogilvie
 */

public class NetworkParser extends Network implements StateNodeInitialiser {
    public final Input<TreeInterface> treeInput = new Input<>("tree", "Tree initialized from extended newick string", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        Node treeRoot = treeInput.get().getRoot();
        NetworkNode newRoot = new NetworkNode(treeRoot);
        setRoot(newRoot);
        speciationNodeCount = 1;
        leafNodeCount = reticulationNodeCount = 0;
        super.initAndValidate();

        // System.out.println(String.format("Root labels = %s -> %s <- %s", treeRoot.getLeft().getID(), treeRoot.getID(), treeRoot.getRight().getID()));
        // System.out.println(String.format("Root numbers = %d -> %d <- %d", treeRoot.getLeft().getNr(), treeRoot.getNr(), treeRoot.getRight().getNr()));
        rebuildNetwork(treeRoot.getLeft(), root, true);
        rebuildNetwork(treeRoot.getRight(), root, false);
        root.updateSizes();
    }

    private void rebuildNetwork(final Node treeNode, final NetworkNode parentNode, final boolean isLeft) {
        // System.out.println(String.format("Current branch: %s - %s a.k.a. %d - %d", parentNode.getID(), treeNode.getID(), parentNode.getNr(), treeNode.getNr()));
        final Node leftChild = treeNode.getLeft();
        final Node rightChild = treeNode.getRight();
        final String nodeLabel = treeNode.getID();
        final int hStart = (nodeLabel == null) ? 0 : nodeLabel.indexOf('#') + 1;
        boolean reticulation = false;
        NetworkNode networkNode = null;
        reticulation = (hStart > 0) && (nodeLabel.length() > hStart) && (nodeLabel.charAt(hStart) == 'H');
        if (reticulation) networkNode = getNode(nodeLabel);

        if (networkNode == null) {
            networkNode = new NetworkNode(treeNode);
            if (reticulation) addReticulationNode(networkNode);
            else if (treeNode.isLeaf()) addLeafNode(networkNode);
            else addSpeciationNode(networkNode);
        }

        // partial reticulation nodes are always right-attached
        if (reticulation && leftChild == null && rightChild == null) networkNode.setRightParent(parentNode);
        // complete reticulation nodes are always left-attached
        else if (reticulation) networkNode.setLeftParent(parentNode);
        else if (isLeft) networkNode.setLeftParent(parentNode);
        else networkNode.setRightParent(parentNode);

        // System.out.println(String.format("%d = %d + %d + %d", nodeCount, leafNodeCount, speciationNodeCount, reticulationNodeCount));
        if (leftChild != null) rebuildNetwork(leftChild, networkNode, true);
        if (rightChild != null) rebuildNetwork(rightChild, networkNode, false);

        networkNode.updateSizes();
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
