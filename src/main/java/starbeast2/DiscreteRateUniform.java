package starbeast2;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.Int;
import beast.base.spec.inference.parameter.IntVectorParam;
import beast.base.util.Randomizer;

/**
 *
 * @author Huw A. Ogilvie
 *
 */

public class DiscreteRateUniform extends AdaptiveOperator {
    final public Input<IntVectorParam<Int>> treeRatesInput = new Input<>("treeRates", "The branch rates.", Validate.REQUIRED);

    private int nNodes;
    private int lowerBound;
    private int integerRange;

    @Override
    public void initAndValidate() {
        final IntVectorParam<Int> treeRates = treeRatesInput.get();
        nNodes = treeRates.size();
        lowerBound = treeRates.getLower();
        integerRange = 1 + treeRates.getUpper() - lowerBound;

        setLimits(1, nNodes);
        super.initAndValidate();
    }

    // symmetric proposal distribution
    @Override
    public double proposal() {
        final IntVectorParam<Int> treeRates = treeRatesInput.get();
        final int[] cycle = chooseK(nNodes);

        for (int i = 0; i < discreteK; i++) {
            final int nodeNumber = cycle[i];
            final int newRate = lowerBound + Randomizer.nextInt(integerRange);
            treeRates.set(nodeNumber, newRate);
        }

        return 0.0;
    }
}
