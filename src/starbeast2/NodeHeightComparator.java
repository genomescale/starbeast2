package starbeast2;

import java.util.Comparator;

import beast.evolution.tree.Node;

public class NodeHeightComparator implements Comparator<Node> {

    @Override
    public int compare(Node nodeA, Node nodeB) {
        final double heightA = nodeA.getHeight();
        final double heightB = nodeB.getHeight();
        return heightA < heightB ? 1 : heightA == heightB ? 0 : -1;
    }

}
