package speciesnetwork;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;

import java.util.List;

/**
 * @author Chi Zhang
 */

@Description("Randomly selects an internal network node and move its height using an uniform sliding window.")
public class NodeSlider extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    final public Input<Double> windowSizeInput =
            new Input<>("windowSize", "The size of the sliding window, default 0.1.", 0.1);

    // empty constructor to facilitate construction by XML + initAndValidate
    public NodeSlider() {
    }

    @Override
    public void initAndValidate() {
    }

    /**
     * Propose a new network-node height from a uniform distribution.
     * If the new value is outside the boundary, the excess is reflected back into the interval.
     * The proposal ratio is 1.0.
     */
    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();
        final double windowSize = windowSizeInput.get();

        // pick an internal node randomly
        List<NetworkNode> intNodes = speciesNetwork.getInternalNodes();
        NetworkNode snNode = intNodes.get(Randomizer.nextInt(intNodes.size()));

        // determine the lower and upper bounds
        NetworkNode leftParent = snNode.getLeftParent();
        NetworkNode rightParent = snNode.getRightParent();
        double upper = Double.MAX_VALUE;
        if (leftParent != null)
            upper = leftParent.getHeight();
        if (rightParent != null && upper > rightParent.getHeight())
            upper = rightParent.getHeight();
        NetworkNode leftChild = snNode.getLeftChild();
        NetworkNode rightChild = snNode.getRightChild();
        double lower = Double.MIN_NORMAL;
        if (leftChild != null)
            lower = leftChild.getHeight();
        if (rightChild != null && lower < rightChild.getHeight())
            lower = rightChild.getHeight();
        if (lower >= upper) return Double.NEGATIVE_INFINITY;  // something is wrong

        // propose a new height, reflect it back if it's outside the boundary
        final double oldHeight = snNode.getHeight();
        double newHeight = oldHeight + (Randomizer.nextDouble() - 0.5) * windowSize;
        while (newHeight < lower || newHeight > upper) {
            if (newHeight < lower)
                newHeight = 2.0 * lower - newHeight;
            if (newHeight > upper)
                newHeight = 2.0 * upper - newHeight;
        }

        // update the new node height
        snNode.setHeight(newHeight);

        return 0.0;
    }
}
