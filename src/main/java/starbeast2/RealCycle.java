package starbeast2;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;

/**
 *
 * @author Huw A. Ogilvie
 *
 */

public class RealCycle extends AdaptiveOperator {
    final public Input<RealVectorParam<Real>> parameterInput = new Input<>("parameter", "The branch rates.", Validate.REQUIRED);

    private int nNodes;

    @Override
    public void initAndValidate() {
        final RealVectorParam<Real> parameter = parameterInput.get();
        nNodes = parameter.size();
        setLimits(2, nNodes);
        super.initAndValidate();
    }

    // symmetric proposal distribution
    @Override
    public double proposal() {
        final RealVectorParam<Real> parameter = parameterInput.get();
        final double[] parameterArray = parameter.getValues();
        final int[] cycle = chooseK(nNodes);

        final int lastNodeNumber = cycle[discreteK - 1];
        double previousRate = parameterArray[lastNodeNumber];
        for (int i = 0; i < discreteK; i++) {
            final int nodeNumber = cycle[i];
            parameter.set(nodeNumber, previousRate);
            previousRate = parameterArray[nodeNumber];
        }

        return 0.0;
    }
}
