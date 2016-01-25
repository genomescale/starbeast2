package starbeast2;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
        super.initAndValidate();

        TreeInterface inputTree = treeInput.get(); 
        for (Node treeLeaf: inputTree.getExternalNodes()) {
            NetworkNode netLeaf = new NetworkNode();
            netLeaf.setID(treeLeaf.getID());
            netLeaf.height = treeLeaf.getHeight();
            addLeafNode(netLeaf);
        }

        Map<String, NetworkNode> networkNodeNames = new HashMap<>();
        Node treeRoot = inputTree.getRoot();
        root.setID(treeRoot.getID());
        root.height = treeRoot.getHeight();
        recursiveCreateNetwork(treeRoot, root, networkNodeNames);
        root.updateSizes();
    }

    private void recursiveCreateNetwork(Node treeNode, NetworkNode netNode, Map<String, NetworkNode> networkNodeNames) {
        for (Node treeChild: treeNode.getChildren()) {
            String childName = treeChild.getID();
            final int hStart = childName.indexOf('#');
            if (hStart != -1) childName = childName.substring(0, hStart);
            NetworkNode netChild = networkNodeNames.get(childName);

            if (netChild == null) {
                netChild = new NetworkNode();
                networkNodeNames.put(childName, netChild);
                netChild.setID(childName);
                netChild.height = treeChild.getHeight();
                addInternalNode(netChild);
                recursiveCreateNetwork(treeChild, netChild, networkNodeNames);
            }

            if (netNode.leftChild == null) {
                if (netChild.getLeftParent() != null) netChild.reorientParents();
                netChild.setLeftParent(netNode);
            } else {
                if (netChild.getRightParent() != null) netChild.reorientParents();
                netChild.setRightParent(netNode);
            }

            netChild.updateSizes();
        }
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
