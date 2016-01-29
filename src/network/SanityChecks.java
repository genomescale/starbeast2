package network;

import java.util.List;

import beast.evolution.tree.Node;

final class SanityChecks {
    protected boolean computeCoalescentTimes(List<GeneTree> geneTrees) {
        for (GeneTree geneTree: geneTrees) {
            if (!geneTree.computeCoalescentTimes()) {
                // this gene tree IS NOT compatible with the species tree
                return false;
            }
        }

        return true;
    }

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
}
