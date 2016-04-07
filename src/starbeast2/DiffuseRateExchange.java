package starbeast2;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;
import beast.util.Randomizer;

/**
 * Changes value of any tree branch rate by magnitude delta,
 * and balances the sum of rates by changing the remaining rates
 * by negative delta / (number of branches - 1)
 *
 * @author Huw A. Ogilvie
 * 
 */

public class DiffuseRateExchange extends AdaptiveOperator {
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

        setLimits(1, nNodes / 2);
        super.initAndValidate();
    }

    // symmetric proposal distribution
    @Override
    public double proposal() {
        final RealParameter treeRates = treeRatesInput.get();
        final double delta = (Randomizer.nextDouble() - 0.5) * deltaScaleFactor;
        final double reciprocalDelta = delta * discreteK / (nNodes - discreteK);
        final double[] treeRatesArray = treeRates.getDoubleValues();

        final int[] permutation = chooseK(nNodes);
        for (int i = 0; i < nNodes; i++) {
            final int nodeNumber = permutation[i];
            final double originalBranchRate = treeRatesArray[nodeNumber];
            final double newBranchRate = (i < discreteK) ? originalBranchRate + delta : originalBranchRate - reciprocalDelta;
            if (newBranchRate < lowerBound || newBranchRate > upperBound) {
                return Double.NEGATIVE_INFINITY;
            } else {
                treeRatesArray[nodeNumber] = newBranchRate;
            }
        }

        for (int i = 0; i < nNodes; i++) {
            final double newRate = treeRatesArray[i];
            treeRates.setValue(i, newRate);
        }

        return 0.0;
    }
}
