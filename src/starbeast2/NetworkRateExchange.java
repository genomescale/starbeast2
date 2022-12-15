package starbeast2;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

/**
 * Changes value of any tree branch rate by magnitude delta,
 * and balances the sum of rates by changing a different branch rate
 * by negative delta.
 *
 * @author Huw A. Ogilvie
 * 
 */

public class NetworkRateExchange extends AdaptiveOperator {
    final public Input<RealParameter> treeRatesInput = new Input<>("treeRates", "The branch rates.", Validate.REQUIRED);
    final public Input<Double> deltaInput = new Input<>("delta", "Magnitude of change for two randomly picked values.", 1.0);

    private int nNodes;
    private double deltaScaleFactor;
    private double lowerBound;
    private double upperBound;

    @Override
    public void initAndValidate() {
        final RealParameter treeRates = treeRatesInput.get();
        nNodes = treeRates.getDimension();
        lowerBound = treeRates.getLower();
        upperBound = treeRates.getUpper();

        deltaScaleFactor = 2.0 * deltaInput.get() / nNodes;

        setLimits(2, nNodes);
        super.initAndValidate();
    }

    // symmetric proposal distribution
    @Override
    public double proposal() {
        final RealParameter treeRates = treeRatesInput.get();
        final double[] treeRatesArray = treeRates.getDoubleValues();
        final int[] network = chooseK(nNodes);

        // exchange a new delta between all pairs of nodes in 'network'
        for (int i = 0; i < (discreteK - 1); i++) {
            for (int j = i + 1; j < discreteK; j++) {
                final double delta = (Randomizer.nextDouble() - 0.5) * deltaScaleFactor;
                treeRatesArray[network[i]] += delta;
                treeRatesArray[network[j]] -= delta;
            }
        }

        for (int i = 0; i < discreteK; i++) {
            final int nodeNumber = network[i];
            final double newRate = treeRatesArray[nodeNumber];
            if (newRate < lowerBound || newRate > upperBound) {
                return Double.NEGATIVE_INFINITY;
            } else {
                treeRates.setValue(nodeNumber, newRate);
            }
        }

        return 0.0;
    }
}
