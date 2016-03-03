package speciesnetwork;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;

/**
 * Flip a random gene tree lineage with all its descendants in one side of the loop to the other side in the network.
 * @author Alexei Drummond
 * @author Chi Zhang
 */
public class FlipNetworkLoop extends Operator {
    public Input<Tree> geneTreeInput =
            new Input<>("geneTree", "The gene tree.", Input.Validate.REQUIRED);
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    public Input<IntegerParameter> embeddingInput =
            new Input<>("embedding", "The matrix to embed the gene tree within the species network.", Input.Validate.REQUIRED);
    final private Multimap<NetworkNode, String> pathDirections = HashMultimap.create();

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Tree geneTree = geneTreeInput.get();
        Network speciesNetwork = speciesNetworkInput.get();
        IntegerParameter embedding = embeddingInput.get();

        List<NetworkNode> hybridNodes = speciesNetwork.getReticulationNodes();
        // if there is no reticulation node, this operator doesn't apply
        if (hybridNodes.isEmpty()) return Double.NEGATIVE_INFINITY;

        // pick a hybrid node randomly from the network
        int rnd = Randomizer.nextInt(hybridNodes.size());
        NetworkNode hybridNode = hybridNodes.get(rnd);
        // find the top node of the minimal loop with the hybrid node at the bottom
        NetworkNode topNode = findLoop(hybridNode, true);
        final int topNodeNr = topNode.getNr();
        final int bottomNodeNr = hybridNode.getNr();

        // find the gene lineages traversing the loop-top network node
        List<Node> topGeneNodes = new ArrayList<>();
        final int geneNodeCount = geneTree.getNodeCount();
        for (int j = 0; j < geneNodeCount; j++) {
            if (embedding.getMatrixValue(topNodeNr, j) > -1) {
                final Node geneNodeTop = geneTree.getNode(j);

                // check if all the descendants of geneNodeTop traverse one side of the loop


                topGeneNodes.add(geneNodeTop);
            }
        }

        // if there is no lineage traversing, this operator doesn't apply
        if (topGeneNodes.isEmpty()) return Double.NEGATIVE_INFINITY;

        // pick a lineage randomly, flip it (and all its descendant lineages) to the other side of the loop
        rnd = Randomizer.nextInt(topGeneNodes.size());
        final Node geneNodeTop = topGeneNodes.get(rnd);



        return 0.0;
    }

    /**
     * @param hybridNode the hybrid node forming the bottom of the loop
     * @return the top network node on the minimal loop from the given hybridization node
     */
    private NetworkNode findLoop(NetworkNode hybridNode, boolean cleanup) {
        // check if hybridNode is actually hybrid
        if (!hybridNode.isReticulation()) throw new RuntimeException();

        NetworkNode topNode = new NetworkNode();

        // traverse left, label A; traverse right, label B
        label(hybridNode.getLeftParent(), "A", null, null);
        label(hybridNode.getRightParent(), "B", "A", topNode);

        // find all the paths connecting top node and hybrid node
        getPathDirections(topNode, topNode, hybridNode, "A");
        getPathDirections(topNode, topNode, hybridNode, "B");

        if (cleanup) {
            unlabel(hybridNode.getLeftParent(), "A");
            unlabel(hybridNode.getRightParent(), "B");
        }

        return topNode;
    }

    private void unlabel(NetworkNode node, String label) {
        node.removeLabel(label);

        if (node.getLeftParent() != null) unlabel(node.getLeftParent(), label);
        if (node.getRightParent() != null) unlabel(node.getRightParent(), label);
    }

    private void label(NetworkNode node, String label, String checkLabel, NetworkNode returnNode) {
        node.addLabel(label);

        if (checkLabel != null && node.hasLabel(checkLabel)) {
            if (returnNode == null || node.getHeight() < returnNode.getHeight())
                returnNode = node;
        }

        if (node.getLeftParent() != null) label(node.getLeftParent(), label, checkLabel, returnNode);
        if (node.getRightParent() != null) label(node.getRightParent(), label, checkLabel, returnNode);
    }

    /**
     * get (forward in time) directions of the loop
     * @param topNode    top node of the loop (speciation node)
     * @param bottomNode bottom node of the loop (hybrid node)
     * @param checkLabel label
     */
    private void getPathDirections(NetworkNode node, NetworkNode topNode, NetworkNode bottomNode, String checkLabel) {
        if (node == topNode) {
            pathDirections.put(node, "");  // initialize as empty string
        }
        if (node == bottomNode)
            return;

        NetworkNode leftNode = node.getLeftChild();
        if (leftNode != null && (leftNode.hasLabel(checkLabel) || leftNode == bottomNode)) {
            for (final String s : pathDirections.get(node)) {
                pathDirections.put(leftNode, s + "0");  // traversing left
            }
            getPathDirections(leftNode, topNode, bottomNode, checkLabel);
        }

        NetworkNode rightNode = node.getRightChild();
        if (rightNode != null && (rightNode.hasLabel(checkLabel) || rightNode == bottomNode)) {
            for (final String s : pathDirections.get(node)) {
                pathDirections.put(rightNode, s + "1"); // traversing right
            }
            getPathDirections(rightNode, topNode, bottomNode, checkLabel);
        }
    }

    // is tree node ancNode the ancestor of childNode?
    boolean isAncestor(Node ancNode, Node childNode) {
        Node nextNode = childNode;
        while (nextNode != ancNode && nextNode != null)
            nextNode = nextNode.getParent();
        return nextNode != null;
    }
}
