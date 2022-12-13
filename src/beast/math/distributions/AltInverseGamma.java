package beast.math.distributions;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.distribution.InverseGamma;
import beast.base.inference.parameter.RealParameter;

@Description("Inverse Gamma distribution, used as prior. Parameterized by its mean instead of scale.")

public class AltInverseGamma extends InverseGamma {
    final public Input<RealParameter> meanInput = new Input<>("mean", "mean of the distribution, defaults to 1");

    void refresh() {
        double alpha;
        double mean;

        if (alphaInput.get() == null) {
            alpha = 2;
        } else {
            alpha = alphaInput.get().getArrayValue();
        }

        if (meanInput.get() == null) {
            mean = 1;
        } else {
            mean = meanInput.get().getValue();
        }

        final double beta = mean * (alpha - 1);

        dist.setAlphaBeta(alpha, beta);
    }
}
