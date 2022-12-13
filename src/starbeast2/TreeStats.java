package starbeast2;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

final public class TreeStats {
    protected static double getLength(TreeInterface tree) {
        final Node treeRoot = tree.getRoot();
        final double treeHeight = treeRoot.getHeight();
        final double treeLength = recurseLength(treeRoot, treeHeight);
        return treeLength;
    }

    private static double recurseLength(final Node treeNode, final double parentHeight) {
        if (treeNode.isLeaf()) {
            return parentHeight;
        } else {
            double subtreeLength = 0.0;

            final double nodeHeight = treeNode.getHeight();
            subtreeLength += parentHeight - nodeHeight;

            final double leftChildLength = recurseLength(treeNode.getLeft(), nodeHeight);
            final double rightChildLength = recurseLength(treeNode.getRight(), nodeHeight);
            subtreeLength += leftChildLength;
            subtreeLength += rightChildLength;

            return subtreeLength;
        }
    }
}
