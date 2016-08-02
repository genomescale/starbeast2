package speciesnetwork.operators;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

@Description("Changes the value of gamma by applying a random walk to the logit of gamma.")
public class GammaRandomWalk extends Operator {
    public Input<Network> speciesNetworkInput = new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    public Input<Double> windowSizeInput = new Input<>("windowSize", "The size of the sliding window, default 1.0.", 1.0);

    @Override
    public void initAndValidate() {
        // TODO Auto-generated method stub

    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final Double windowSize = windowSizeInput.get();

        final List<NetworkNode> reticulationNodes = speciesNetwork.getReticulationNodes();
        final int nReticulations = reticulationNodes.size();
        final int randomNodeIndex = Randomizer.nextInt(nReticulations);
        final NetworkNode randomNode = reticulationNodes.get(randomNodeIndex);

        final Double logOddsShift = (Randomizer.nextDouble() * windowSize * 2) - windowSize;
        final Double currentGamma = randomNode.getGamma();
        final Double currentLogOdds = Math.log(currentGamma / (1.0 - currentGamma));
        final Double newLogOdds = currentLogOdds + logOddsShift;
        final Double newGamma = 1.0 / (1.0 + Math.exp(-newLogOdds));
        randomNode.setGamma(newGamma);

        return 0.0;
    }
}
