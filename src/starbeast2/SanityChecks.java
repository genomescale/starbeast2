package starbeast2;

import java.util.List;

import beast.evolution.tree.Node;

final class SanityChecks {
    protected static boolean checkTreeSanity(Node node) {
        final double nodeHeight = node.getHeight();
        final List<Node> children = node.getChildren();
        assert children.size() == 2;

        for (Node childNode: children) {
            final double childHeight = childNode.getHeight();
            assert childNode.getParent() == node;

            if (childNode.isLeaf()) {
                // direct ancestor branches have zero height
                // so equal height allowed in this case
                assert childHeight <= nodeHeight;
            } else{
                if (childHeight >= nodeHeight) {
                    System.out.println(node.toNewick());
                }
                assert childHeight < nodeHeight;
                checkTreeSanity(childNode);
            }
        }

        return true;
    }
}
