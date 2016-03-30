package speciesnetwork.operators;

import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * @author Alexei Drummond
 * @author Chi Zhang
 */

@Description("Flip a random gene tree lineage with all its descendants in one side of the loop to another side in the network.")
public class FlipNetworkLoop extends Operator {
    public Input<Tree> geneTreeInput =
            new Input<>("geneTree", "The gene tree.", Input.Validate.REQUIRED);
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    public Input<IntegerParameter> embeddingInput =
            new Input<>("embedding", "The matrix to embed the gene tree within the species network.", Input.Validate.REQUIRED);

    private Multimap<NetworkNode, String> pathDirections = HashMultimap.create(); // all directions in the loop
    private Map<Node, String> lineagePathDir = new HashMap<>();  // lineage traversing direction
    private int speciesLeafCount;
    private IntegerParameter embedding;
    private String traverseDirection;

    @Override
    public void initAndValidate() {
        Network speciesNetwork = speciesNetworkInput.get();
        speciesLeafCount = speciesNetwork.getLeafNodes().size();
        embedding = embeddingInput.get();
    }

    @Override
    public double proposal() {
        Tree geneTree = geneTreeInput.get();
        Network speciesNetwork = speciesNetworkInput.get();

        List<NetworkNode> hybridNodes = speciesNetwork.getReticulationNodes();
        // if there is no reticulation node, this operator doesn't apply
        if (hybridNodes.isEmpty()) return Double.NEGATIVE_INFINITY;

        pathDirections.clear(); // clear before filling in
        lineagePathDir.clear();
        // pick a hybrid node randomly from the network
        int rnd = Randomizer.nextInt(hybridNodes.size());
        NetworkNode hybridNode = hybridNodes.get(rnd);
        // find the top node of the minimal loop with the hybrid node at the bottom
        NetworkNode topNode = findLoop(hybridNode, true);  // pathDirections is also set
        final int topNodeNr = topNode.getNr();

        final int geneNodeCount = geneTree.getNodeCount();
        for (int j = 0; j < geneNodeCount; j++) {
            // find the gene lineages traversing the loop-top network node
            if (embedding.getMatrixValue(topNodeNr-speciesLeafCount, j) > -1) {
                final Node geneNodeTop = geneTree.getNode(j);
                Set<NetworkNode> traversedNodes = new HashSet<>();
                traverseDirection = "";
                // check if all the descendants of geneNodeTop traverse one side of the loop
                if (allLineagesInLoop(geneNodeTop, topNode, hybridNode, traversedNodes))
                    lineagePathDir.put(geneNodeTop, traverseDirection);
            }
        }
        // if there is no lineage traversing, this operator doesn't apply
        if (lineagePathDir.isEmpty()) return Double.NEGATIVE_INFINITY;

        // pick a lineage randomly, with its traversing direction
        List<Node> keys = new ArrayList<>(lineagePathDir.keySet());
        rnd = Randomizer.nextInt(keys.size());
        final Node geneNodeTop = keys.get(rnd);
        final String currentPathDir = lineagePathDir.get(geneNodeTop);
        // delete the current direction from the list
        List<String> loopPathDirections = new ArrayList<>(pathDirections.get(hybridNode));
        loopPathDirections.remove(currentPathDir);

        // pick a new direction randomly
        assert (loopPathDirections.size() > 0);
        rnd = Randomizer.nextInt(loopPathDirections.size());
        final String proposedPathDir = loopPathDirections.get(rnd);

        // make the flip
        // first reset the embedding of the current path, then build the new embedding of the proposed path
        if (!setEmbedding(geneNodeTop, topNode, hybridNode, currentPathDir, true) ||
            !setEmbedding(geneNodeTop, topNode, hybridNode, proposedPathDir, false))
            return Double.NEGATIVE_INFINITY;

        return 0.0;
    }

    /**
     * @param hybridNode the hybrid node forming the bottom of the loop
     * @return the top network node on the minimal loop from the given hybridization node
     */
    private NetworkNode findLoop(NetworkNode hybridNode, boolean cleanup) {
        // check if hybridNode is actually hybrid
        if (!hybridNode.isReticulation()) throw new RuntimeException();

        NetworkNode[] returnNode = new NetworkNode[1];
        // traverse left, label A; traverse right, label B
        label(hybridNode.getLeftParent(), "A", null, null);
        label(hybridNode.getRightParent(), "B", "A", returnNode);

        NetworkNode topNode = returnNode[0];
        pathDirections.put(topNode, "");  // initialize as empty string
        // find all the path directions from top node to hybrid node
        getPathDirections(topNode, hybridNode, "A");
        getPathDirections(topNode, hybridNode, "B");

        if (cleanup) {
            unlabel(hybridNode.getLeftParent(), "A");
            unlabel(hybridNode.getRightParent(), "B");
        }

        return returnNode[0];
    }

    private void unlabel(NetworkNode node, String label) {
        node.removeLabel(label);

        if (node.getLeftParent() != null) unlabel(node.getLeftParent(), label);
        if (node.getRightParent() != null) unlabel(node.getRightParent(), label);
    }

    private void label(NetworkNode node, String label, String checkLabel, NetworkNode[] returnNode) {
        node.addLabel(label);

        if (checkLabel != null && node.hasLabel(checkLabel)) {
            if (returnNode[0] == null || node.getHeight() < returnNode[0].getHeight())
                returnNode[0] = node;
        }

        if (node.getLeftParent() != null) label(node.getLeftParent(), label, checkLabel, returnNode);
        if (node.getRightParent() != null) label(node.getRightParent(), label, checkLabel, returnNode);
    }

    /**
     * get (forward in time) directions of the loop, from top node of the loop to
     * @param bottomNode bottom node of the loop (hybrid node)
     * @param checkLabel label
     */
    private void getPathDirections(NetworkNode node, NetworkNode bottomNode, String checkLabel) {
        if (node == bottomNode)
            return;

        NetworkNode leftNode = node.getLeftChild();
        if (leftNode != null && (leftNode.hasLabel(checkLabel) || leftNode == bottomNode)) {
            for (final String s : pathDirections.get(node)) {
                pathDirections.put(leftNode, s + "0");  // traversing left
            }
            getPathDirections(leftNode, bottomNode, checkLabel);
        }

        NetworkNode rightNode = node.getRightChild();
        if (rightNode != null && (rightNode.hasLabel(checkLabel) || rightNode == bottomNode)) {
            for (final String s : pathDirections.get(node)) {
                pathDirections.put(rightNode, s + "1"); // traversing right
            }
            getPathDirections(rightNode, bottomNode, checkLabel);
        }
    }

    /**
     * Are geneNode and its descendants traversing one side of the loop
     */
    private boolean allLineagesInLoop (Node geneNode, NetworkNode netNode, NetworkNode bottomNode, Set<NetworkNode> traversedNodes) {
        if (netNode.isLeaf()) return false;
        traversedNodes.add(netNode); // add the node to the path set
        if (netNode == bottomNode) return true;

        if (geneNode.getHeight() < netNode.getHeight()) {
            final int traversalNodeNr = netNode.getNr()-speciesLeafCount;
            final int geneNodeNr = geneNode.getNr();
            final NetworkNode leftNode = netNode.getLeftChild();
            final NetworkNode rightNode = netNode.getRightChild();

            if (embedding.getMatrixValue(traversalNodeNr, geneNodeNr) == 0 && leftNode != null) {
                if (!traversedNodes.contains(leftNode))
                    traverseDirection += "0";
                return allLineagesInLoop(geneNode, leftNode, bottomNode, traversedNodes);
            }
            else if (embedding.getMatrixValue(traversalNodeNr, geneNodeNr) == 1 && rightNode != null) {
                if (!traversedNodes.contains(rightNode))
                    traverseDirection += "1";
                return allLineagesInLoop(geneNode, rightNode, bottomNode, traversedNodes);
            } else {
                return false; // something is wrong
            }
        } else {
            return !geneNode.isLeaf() &&
                    allLineagesInLoop(geneNode.getLeft(), netNode, bottomNode, traversedNodes) &&
                    allLineagesInLoop(geneNode.getRight(), netNode, bottomNode, traversedNodes);
        }
    }

    /**
     * set or reset the embedding according to the traversal direction of one side of the loop
     */
    private boolean setEmbedding (Node geneNode, NetworkNode netNode, NetworkNode bottomNode,
                                  String traverseDirection, boolean reset) {
        if (netNode.isLeaf()) return false;
        if (netNode == bottomNode) return true;

        if (geneNode.getHeight() < netNode.getHeight()) {
            final int traversalNodeNr = netNode.getNr()-speciesLeafCount;
            final int geneNodeNr = geneNode.getNr();
            final NetworkNode leftNode = netNode.getLeftChild();
            final NetworkNode rightNode = netNode.getRightChild();

            if (leftNode != null && getDirection(leftNode, bottomNode, traverseDirection) == '0') {
                if (reset)
                    embedding.setMatrixValue(traversalNodeNr, geneNodeNr, -1);
                else
                    embedding.setMatrixValue(traversalNodeNr, geneNodeNr, 0);
                return setEmbedding(geneNode, leftNode, bottomNode, traverseDirection, reset);
            }
            else if (rightNode != null && getDirection(rightNode, bottomNode, traverseDirection) == '1') {
                if (reset)
                    embedding.setMatrixValue(traversalNodeNr, geneNodeNr, -1);
                else
                    embedding.setMatrixValue(traversalNodeNr, geneNodeNr, 1);
                return setEmbedding(geneNode, rightNode, bottomNode, traverseDirection, reset);
            } else {
                return false; // something is wrong
            }
        } else {
            return !geneNode.isLeaf() &&
                    setEmbedding(geneNode.getLeft(), netNode, bottomNode, traverseDirection, reset) &&
                    setEmbedding(geneNode.getRight(), netNode, bottomNode, traverseDirection, reset);
        }
    }

    /** get traversal direction to node given the end node and the traversal direction
     * @return 0 -> left, 1 -> right, -1 wrong
     */
    private int getDirection(NetworkNode node, NetworkNode end, String direction) {
        // can go backward as child and parent directions are corresponded the same
        int i = direction.length() - 1;
        while (end != node && i >= 0) {
            NetworkNode left = end.getLeftParent();
            NetworkNode right = end.getRightParent();
            if (left != null && direction.charAt(i) == '0') end = left;
            else if (right != null && direction.charAt(i) == '1') end = right;
            else return -1;  // something is wrong
            i--;
        }
        return end == node ? direction.charAt(i) : -1;
    }
}
