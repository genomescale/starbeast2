package speciesnetwork;

import java.util.List;

import beast.evolution.tree.Node;

public final class SanityChecks {
    protected boolean checkTreeSanity(Node node) {
        final List<Node> children = node.getChildren();
        final int nChildren = children.size();

        for (Node childNode: children) {
            assert childNode.getParent() == node;
            assert childNode.getHeight() <= node.getHeight();
            if (!node.isLeaf()) {
                checkTreeSanity(childNode);
            }
        }

        if (node.isLeaf()) {
            assert nChildren == 0;
        } else {
            assert nChildren == 2;
        }

        return true;
    }

    public boolean checkNetworkSanity(NetworkNode node) {
        final List<NetworkNode> children = node.getChildren();
        final List<NetworkNode> parents = node.getParents();

        final int nChildren = children.size();
        final int nParents = parents.size();

        for (NetworkNode child: children) {
            if (node.getLeftChild() == child)
                assert child.getLeftParent() == node;
            if (node.getRightChild() == child)
                assert child.getRightParent() == node;
            assert child.getHeight() <= node.getHeight();
            if (!node.isLeaf()) {
                checkNetworkSanity(child);
            }
        }

        if (node.isLeaf()) {
            assert nChildren == 0 && nParents == 1;
        } else if (node.isReticulation()){
            assert nChildren == 1 && nParents == 2;
        } else if (node.isRoot()){
            assert nChildren == 2 && nParents == 0;
        } else {
            assert nChildren == 2 && nParents == 1;
        }

        return true;
    }
}
