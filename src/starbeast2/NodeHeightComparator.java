package starbeast2;

import java.util.Comparator;

import beast.evolution.tree.Node;

class NodeHeightComparator implements Comparator<Node> {

    @Override
    public int compare(Node nodeA, Node nodeB) {
        final double heightA = nodeA.getHeight();
        final double heightB = nodeB.getHeight();
        if (heightA == heightB) {
            final int depthA = calculateNodeDepth(nodeA);
            final int depthB = calculateNodeDepth(nodeB);
            if (depthA == depthB) {
                return 0;
            }
            return depthA > depthB ? 1 : -1;
        }
        return heightA < heightB ? 1 : -1;
    }

    private int calculateNodeDepth(Node node) {
        if (node.isRoot()) {
            return 1;
        }

        return calculateNodeDepth(node.getParent()) + 1;
    }
}
