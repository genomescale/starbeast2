package starbeast2;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.Int;
import beast.base.spec.inference.parameter.IntVectorParam;

/**
 *
 * @author Huw A. Ogilvie
 *
 */

public class DiscreteRateCycle extends AdaptiveOperator {
    final public Input<IntVectorParam<Int>> treeRatesInput = new Input<>("treeRates", "The branch rates.", Validate.REQUIRED);

    private int nNodes;

    @Override
    public void initAndValidate() {
        final IntVectorParam<Int> treeRates = treeRatesInput.get();
        nNodes = treeRates.size();
        setLimits(2, nNodes);
        super.initAndValidate();
    }

    // symmetric proposal distribution
    @Override
    public double proposal() {
        final IntVectorParam<Int> treeRates = treeRatesInput.get();
        final int[] treeRatesArray = treeRates.getValues();
        final int[] cycle = chooseK(nNodes);

        final int lastNodeNumber = cycle[discreteK - 1];
        int previousRate = treeRatesArray[lastNodeNumber];
        for (int i = 0; i < discreteK; i++) {
            final int nodeNumber = cycle[i];
            treeRates.set(nodeNumber, previousRate);
            previousRate = treeRatesArray[nodeNumber];
        }

        return 0.0;
    }
}
