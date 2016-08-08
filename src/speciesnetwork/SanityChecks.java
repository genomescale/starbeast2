package speciesnetwork;

import java.util.List;

import com.google.common.collect.Multiset;

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

    public static void checkNetworkSanity(NetworkNode node) {
        final Multiset<NetworkNode> children = node.getChildren();
        final Multiset<NetworkNode> parents = node.getParents();

        final int nChildren = children.size();
        final int nParents = parents.size();

        assert nChildren <= 2;
        assert nParents <= 2;

        if (node.getHeight() == 0.0) {
            assert nChildren == 0;
            assert nParents == 1;
        } else {
            if (nParents == 2) {
                assert nChildren == 1;
            } else {
                assert nChildren == 2;
            }
        }

        for (NetworkNode child: children) {
            assert child.getParents().contains(node);
            assert child.getHeight() <= node.getHeight();
            checkNetworkSanity(child);
        }
    }
}
