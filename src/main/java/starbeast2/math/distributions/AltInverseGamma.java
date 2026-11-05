package starbeast2.math.distributions;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.RealScalar;
import org.apache.commons.numbers.gamma.LogGamma;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.GammaDistribution;

import java.util.List;

@Description("Inverse Gamma distribution, used as prior. Parameterized by its mean instead of scale.")
public class AltInverseGamma extends ScalarDistribution<RealScalar<PositiveReal>, Double> {
    final public Input<RealScalar<PositiveReal>> alphaInput = new Input<>("alpha", "shape parameter, defaults to 2");
    final public Input<RealScalar<PositiveReal>> meanInput = new Input<>("mean", "mean of the distribution, defaults to 1");

    private GammaDistribution dist = GammaDistribution.of(2, 1);
    private ContinuousDistribution.Sampler sampler;

    private double alpha;
    private double beta;
    // log of the constant beta^alpha/Gamma(alpha)
    private double C;

    @Override
    public void initAndValidate() {
        refresh();
        super.initAndValidate();
    }

    /**
     * ensure internal state is up to date, re-deriving beta from alpha and the mean *
     */
    @Override
    public void refresh() {
        alpha = (alphaInput.get() != null) ? alphaInput.get().get() : 2.0;
        final double mean = (meanInput.get() != null) ? meanInput.get().get() : 1.0;
        beta = mean * (alpha - 1);
        C = alpha * Math.log(beta) - LogGamma.value(alpha);

        // Floating point comparison
        if (isNotEqual(dist.getShape(), alpha) || isNotEqual(dist.getScale(), 1.0 / beta)) {
            dist = GammaDistribution.of(alpha, 1.0 / beta);
        }
    }

    @Override
    public double calculateLogP() {
        logP = logDensity(param.get()); // no unboxing needed, faster
        return logP;
    }

    @Override
    protected double calcLogP(Double value) {
        return logDensity(value); // scalar
    }

    // handle offset in one place
    public double logDensity(double x) {
        refresh(); // this make sure distribution parameters are updated if they are sampled during MCMC
        return -(alpha + 1.0) * Math.log(x) - (beta / x) + C;
    }

    @Override
    public double density(double x) {
        return Math.exp(logDensity(x));
    }

    @Override
    public List<Double> sample() {
        if (sampler == null) {
            // Ensure sampler exists
            sampler = dist.createSampler(rng);
        }
        final double y = sampler.sample(); // sample from Gamma
        final double x = 1.0 / y; // sample from inverse Gamma
        return List.of(x);
    }

    @Override
    protected Object getApacheDistribution() {
        refresh(); // this make sure distribution parameters are updated if they are sampled during MCMC
        return dist;
    }
}
