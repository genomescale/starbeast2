package speciesnetwork;

import beast.core.Operator;
import beast.util.TreeParser;

/**
 * @author Alexei Drummond
 * @author Chi Zhang
 */
public class FlipNetworkLoop extends Operator {

    @Override
    public void initAndValidate() {
    }

    /**
     *
     */
    @Override
    public double proposal() {

        return 0.0;
    }

    /**
     * @param network
     * @param hybridNode the hybrid node forming the bottom of the loop
     * @return the top network node on the minimal loop from the given hybrization node
     */
    static NetworkNode findLoop(Network network, NetworkNode hybridNode, boolean cleanup) {

        // check if hybridNode is actually hybrid
        if (!hybridNode.isReticulation()) throw new RuntimeException();

        NetworkNode[] returnNode = new NetworkNode[1];

        // traverse left, label A
        label(hybridNode.getLeftParent(), "A", null, null);
        label(hybridNode.getRightParent(), "B", "A", returnNode);

        if (cleanup) {
            unlabel(hybridNode.getLeftParent(), "A");
            unlabel(hybridNode.getRightParent(), "B");
        }

        return returnNode[0];

    }

    static void unlabel(NetworkNode node, String label) {
        node.removeLabel(label);

        if (node.getLeftParent() != null) unlabel(node.getLeftParent(), label);
        if (node.getRightParent() != null) unlabel(node.getRightParent(), label);
    }

    static void label(NetworkNode node, String label, String checkLabel, NetworkNode[] returnNode) {

        node.addLabel(label);

        if (checkLabel != null && node.hasLabel(checkLabel)) {
            if (returnNode[0] == null) {
                returnNode[0] = node;
            } else {
                double currentHeight = returnNode[0].getHeight();
                if (node.getHeight() < currentHeight) {
                    returnNode[0] = node;
                }
            }
        }

        if (node.getLeftParent() != null) label(node.getLeftParent(), label, checkLabel, returnNode);
        if (node.getRightParent() != null) label(node.getRightParent(), label, checkLabel, returnNode);
    }

    /*
    public static void main(String[] args) throws Exception {

        final String testNetwork = "((((A:0.1)#H1:0.1)#H2:0.3,#H4:0.1)#S2:0.1,((((#H1:0.1)#H3:0.1,#H2:0.1)#S1:0.1)#H4:0.1,#H3:0.3)#S3:0.1)#R";
        TreeParser treeParser = new TreeParser();
        NetworkParser networkParser = new NetworkParser();

        treeParser.initByName("newick", testNetwork, "IsLabelledNewick", true, "adjustTipHeights", false);
        networkParser.initByName("tree", treeParser);

        NetworkNode hybridNode = networkParser.getNode("#H4");

        NetworkNode topNode = findLoop(networkParser, hybridNode, true);

        System.out.println(testNetwork);

        System.out.println("Hybrid Node =" + hybridNode.getID());
        System.out.println("Top Node =" + topNode.getID());

    } */
}
