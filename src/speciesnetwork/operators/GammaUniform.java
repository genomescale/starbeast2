package speciesnetwork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

@Description("Changes the value of gamma by applying a random walk to the logit of gamma.")
public class GammaUniform extends Operator {
    public Input<Network> speciesNetworkInput = new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        // TODO Auto-generated method stub

    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();

        final int nReticulations = speciesNetwork.getReticulationNodeCount();
        final int randomNodeIndex = Randomizer.nextInt(nReticulations) + speciesNetwork.getReticulationOffset();
        final NetworkNode randomNode = speciesNetwork.getNode(randomNodeIndex);

        final Double newGamma = Randomizer.nextDouble();
        randomNode.setGamma(newGamma);

        return 0.0;
    }
}
