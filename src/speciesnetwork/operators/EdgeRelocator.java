package speciesnetwork.operators;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;

import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

/**
 * @author Chi Zhang
 */

@Description("Relocate the source of an edge starting with speciation node, " +
        "or the destination of an edge ending with hybridization node.")
public class EdgeRelocator extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);

    private enum Direction {LEFT, RIGHT}

    // empty constructor to facilitate construction by XML + initAndValidate
    public EdgeRelocator() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();

        // pick an internal node randomly
        List<NetworkNode> intNodes = speciesNetwork.getInternalNodes();
        NetworkNode snNode = intNodes.get(Randomizer.nextInt(intNodes.size()));

        // start moving
        Direction edgeDirection;
        double lower, upper;
        if (snNode.isReticulation()) {
            // move either the left parent or the right parent branch
            if (Randomizer.nextBoolean()) {
                edgeDirection = Direction.LEFT;
                upper = snNode.getLeftParent().getHeight();
            } else {
                edgeDirection = Direction.RIGHT;
                upper = snNode.getRightParent().getHeight();
            }

            // look for all the candidates branches to attach to
            for (NetworkNode node: speciesNetwork.getAllNodesAsArray()) {
                if (node.getHeight() < upper && node != snNode) {

                }
            }


        } else {
            // move either the left child or the right child branch
            if (Randomizer.nextBoolean()) edgeDirection = Direction.LEFT;
            else edgeDirection = Direction.RIGHT;



        }


        return 0.0;
    }

}
