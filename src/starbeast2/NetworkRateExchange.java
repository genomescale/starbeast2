package starbeast2;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;
import beast.util.Randomizer;

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

        deltaScaleFactor = deltaInput.get() / nNodes;

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
                final double delta = (Randomizer.nextDouble() - 0.5) * deltaScaleFactor * (4.0 / discreteK);
                treeRatesArray[network[i]] += delta;
                treeRatesArray[network[j]] -= delta;
            }
        }

        for (int i = 0; i < discreteK; i++) {
            final double newRateI = treeRatesArray[network[i]];
            if (newRateI < lowerBound || newRateI > upperBound) {
                return Double.NEGATIVE_INFINITY;
            } else {
                treeRates.setValue(network[i], newRateI);
            }
        }

        return 0.0;
    }
}
