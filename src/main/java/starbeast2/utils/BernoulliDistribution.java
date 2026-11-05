package starbeast2.utils;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.distribution.TensorDistribution;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;

import java.util.Arrays;
import java.util.List;

// Diff vs beast-base's beast.base.spec.inference.distribution.Bernoulli:
// - beast-base's Bernoulli extends ScalarDistribution<BoolScalar, Boolean> and models a single
//   true Boolean outcome; a vector of trials is obtained externally by wrapping it in an IID
//   distribution over a BoolVectorParam.
// - This class instead extends TensorDistribution<RealVector<UnitInterval>, Double> directly,
//   preserving the original BEAST2 behaviour of looping over a real-valued (0.0/1.0) vector
//   param in one class, since starbeast2's indicator parameters are real-valued rather than
//   true booleans. It cannot be replaced by beast-base's Bernoulli + IID without also migrating
//   those callers' data from RealVectorParam to BoolVectorParam, which would be a breaking change.
@Description("Bernoulli distribution over a vector of independent 0/1-valued components, " +
        "each with the same success probability.")
public class BernoulliDistribution extends TensorDistribution<RealVector<UnitInterval>, Double> {
    final public Input<RealScalar<UnitInterval>> pInput = new Input<>("p", "success probability ", Validate.REQUIRED);

    private double p;

    @Override
    public void initAndValidate() {
        refresh();
        super.initAndValidate();
    }

    /**
     * ensure internal state is up to date *
     */
    @Override
    public void refresh() {
        p = pInput.get().get();
    }

    @Override
    public double calculateLogP() {
        logP = calcLogP(param.getElements());
        return logP;
    }

    @Override
    protected double calcLogP(Double... value) {
        return calcLogP(Arrays.asList(value));
    }

    private double calcLogP(List<Double> value) {
        refresh(); // this make sure distribution parameters are updated if they are sampled during MCMC

        double logP = 0;
        for (double x : value) {
            if (x == 1) {
                logP += Math.log(p);
            } else if (x == 0) {
                logP += Math.log(1 - p);
            } else {
                throw new IllegalArgumentException("value to calculate binomial is not 0 nor 1");
            }
        }
        return logP;
    }

    @Override
    public List<Double> sample() {
        return null;
    }

    @Override
    public Double getLowerBoundOfParameter() {
        return 0.0;
    }

    @Override
    public Double getUpperBoundOfParameter() {
        return 1.0;
    }
}
