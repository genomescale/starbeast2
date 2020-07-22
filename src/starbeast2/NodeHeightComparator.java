package starbeast2;

import java.util.Comparator;

import beast.evolution.tree.Node;

// unambiguous, so only identical nodes are considered equal
// if height is equal, consider distance from root (depth)
// if depth is equal, consider assigned node number
public final class NodeHeightComparator implements Comparator<Node> {
    final int lessThan = -1;
    final int greaterThan = 1;

    @Override
    public int compare(Node nodeA, Node nodeB) {
        final double heightA = nodeA.getHeight();
        final double heightB = nodeB.getHeight();
        if (heightA == heightB) {
            final int depthA = calculateNodeDepth(nodeA);
            final int depthB = calculateNodeDepth(nodeB);
            if (depthA == depthB) {
                final int nodeNumberA = nodeA.getNr();
                final int nodeNumberB = nodeB.getNr();
                if (nodeNumberA == nodeNumberB) {
                    return 0;
                }
                return nodeNumberA > nodeNumberB ? greaterThan : lessThan;
            }
            return depthA > depthB ? greaterThan : lessThan;
        }
        return heightA > heightB ? greaterThan : lessThan;
    }

    private int calculateNodeDepth(Node node) {
        if (node.isRoot()) {
            return 0;
        }

        return calculateNodeDepth(node.getParent()) - 1;
    }
}

