package speciesnetwork;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;

import java.util.ArrayList;
import java.util.List;

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

    private Tree geneTree;
    private Network speciesNetwork;
    private IntegerParameter embedding;

    @Override
    public void initAndValidate() {
    }

    /**
     *
     */
    @Override
    public double proposal() {
        geneTree = geneTreeInput.get();
        speciesNetwork = speciesNetworkInput.get();
        embedding = embeddingInput.get();

        List<NetworkNode> hybridNodes = speciesNetwork.getReticulationNodes();
        // if there is no reticulation node, this operator doesn't apply
        if (hybridNodes.isEmpty()) return Double.NEGATIVE_INFINITY;

        // pick a hybrid node randomly from the network
        int rnd = Randomizer.nextInt(hybridNodes.size());
        NetworkNode hybridNode = hybridNodes.get(rnd);
        // find the top node of the loop with the hybrid node
        NetworkNode topNode = findLoop(hybridNode, true);
        final int topNodeNum = topNode.getNr();

        // find the gene lineages traversing the top network node
        final int geneNodeCount = geneTree.getNodeCount();
        List<Integer> geneNodeNums = new ArrayList<>();
        for (int j = 0; j < geneNodeCount; j++) {
            if (embedding.getMatrixValue(topNodeNum, j) > -1)
                geneNodeNums.add(j);
        }
        // if there is no lineage traversing, this operator doesn't apply
        if (geneNodeNums.isEmpty()) return Double.NEGATIVE_INFINITY;

        /* pick a lineage randomly from all the lineages traversing the top node
         * if this lineage (and all its descendant lineages) traverses the hybrid node,
         * flip it (them) to the other side of the loop; otherwise, abort
         */
        rnd = Randomizer.nextInt(geneNodeNums.size());
        int geneNodeNr = geneNodeNums.get(rnd);



        return 0.0;
    }

    /**
     * @param hybridNode the hybrid node forming the bottom of the loop
     * @return the top network node on the minimal loop from the given hybrization node
     */
    private NetworkNode findLoop(NetworkNode hybridNode, boolean cleanup) {
        // check if hybridNode is actually hybrid
        if (!hybridNode.isReticulation()) throw new RuntimeException();

        NetworkNode[] returnNode = new NetworkNode[1];

        // traverse left, label A; traverse right, label B
        label(hybridNode.getLeftParent(), "A", null, null);
        label(hybridNode.getRightParent(), "B", "A", returnNode);

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


    /* public static void main(String[] args) throws Exception {
        final String testNetwork =
            "((((A:0.1)#H1:0.1)#H2:0.3,#H4:0.1)#S2:0.1,((((#H1:0.1)#H3:0.1,#H2:0.1)#S1:0.1)#H4:0.1,#H3:0.3)#S3:0.1)#R";

        TreeParser treeParser = new TreeParser();
        NetworkParser networkParser = new NetworkParser();
        treeParser.initByName("newick", testNetwork, "IsLabelledNewick", true, "adjustTipHeights", false);
        networkParser.initByName("tree", treeParser);

        NetworkNode hybridNode = networkParser.getNode("#H3");
        NetworkNode topNode = findLoop(hybridNode, true);

        System.out.println(testNetwork);
        System.out.println("Hybrid Node =" + hybridNode.getID());
        System.out.println("Top Node =" + topNode.getID());
    } */
}
