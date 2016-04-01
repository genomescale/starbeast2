package starbeast2;

import beast.core.Input;
import beast.core.Operator;
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

public class PairRateExchange extends Operator {
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
    }

    @Override
    public double proposal() {
        // symmetric proposal distribution
        final double delta = (Randomizer.nextDouble() - 0.5) * deltaScaleFactor * 2.0;

        final RealParameter treeRates = treeRatesInput.get();
        final int nodeNumber = Randomizer.nextInt(nNodes);
        final double originalBranchRate = treeRates.getValue(nodeNumber);
        final double newBranchRate = originalBranchRate + delta;
        if (newBranchRate < lowerBound || newBranchRate > upperBound) {
            return Double.NEGATIVE_INFINITY;
        } else {
            treeRates.setValue(nodeNumber, newBranchRate);
        }

        int altNodeNumber = Randomizer.nextInt(nNodes);
        while (altNodeNumber == nodeNumber) {
            altNodeNumber = Randomizer.nextInt(nNodes);
        }

        final double originalAltBranchRate = treeRates.getValue(altNodeNumber);
        final double newAltBranchRate = originalAltBranchRate - delta;
        if (newAltBranchRate < lowerBound || newAltBranchRate > upperBound) {
            return Double.NEGATIVE_INFINITY;
        } else {
            treeRates.setValue(altNodeNumber, newAltBranchRate);
        }

        return 0.0;
    }
}
