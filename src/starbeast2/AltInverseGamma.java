package starbeast2;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;

@Description("Inverse Gamma distribution, used as prior. Parameterized using shape (alpha) and mean instead of shape and scale (beta).")
public class AltInverseGamma extends ParametricDistribution {
    final public Input<RealParameter> alphaInput = new Input<>("alpha", "shape parameter, defaults to 2.0");
    final public Input<RealParameter> meanInput = new Input<>("mean", "mean of the distribution, defaults to 1.0");

    InverseGammaImpl dist = new InverseGammaImpl(2.0, 1.0);

    @Override
    public void initAndValidate() {
        refresh();
    }

    /**
     * ensure internal state is up to date *
     */
    void refresh() {
        final double alpha = (alphaInput.get() == null) ? 2.0 : alphaInput.get().getValue();
        final double mean = (meanInput.get() == null) ? 1.0 : meanInput.get().getValue();
        final double beta = mean * (alpha - 1.0);
        dist.setAlphaBeta(alpha, beta);
    }

    @Override
    public Distribution getDistribution() {
        refresh();
        return dist;
    }

    class InverseGammaImpl implements ContinuousDistribution {
        double m_fAlpha;
        double m_fBeta;
        // log of the constant beta^alpha/Gamma(alpha)
        double C;

        InverseGammaImpl(double alpha, double beta) {
            setAlphaBeta(alpha, beta);
        }

        void setAlphaBeta(double alpha, double beta) {
            m_fAlpha = alpha;
            m_fBeta = beta;
            C = m_fAlpha * Math.log(m_fBeta) - org.apache.commons.math.special.Gamma.logGamma(m_fAlpha);
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
            throw new MathException("Not implemented yet");
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            throw new MathException("Not implemented yet");
        }

        @Override
        public double inverseCumulativeProbability(double p) throws MathException {
            throw new MathException("Not implemented yet");
        }

        @Override
        public double density(double x) {
            double logP = logDensity(x);
            return Math.exp(logP);
        }

        @Override
        public double logDensity(double x) {
            double logP = -(m_fAlpha + 1.0) * Math.log(x) - (m_fBeta / x) + C;
            return logP;
        }
    } // class OneOnXImpl


    @Override
    public double getMean() {
    	return meanInput.get().getValue();
    }
    
} // class InverseGamma
