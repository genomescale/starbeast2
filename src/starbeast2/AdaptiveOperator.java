package starbeast2;



import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;

/**
 * Common methods for executing the same operator multiple times in one step.
 * Useful for cases where acceptance ratios are way too high.
 *
 * @author Huw A. Ogilvie
 */
@Description("A generic operator for use with a sum-constrained (possibly weighted) vector parameter.")
public abstract class AdaptiveOperator extends Operator {
    public final Input<Integer> kInput = new Input<>("k", "Number of operations to perform per step.", 1);
    public final Input<Boolean> autoOptimizeInput = new Input<>("autoOptimize", "Adjust 'k' during the MCMC run to improve mixing.", true);

    protected boolean autoOptimize;
    protected int discreteK;
    protected double continuousK;

    protected double lower;
    protected double upper;

    @Override
    public void initAndValidate() {
        discreteK = kInput.get();
        continuousK = discreteK;
        autoOptimize = autoOptimizeInput.get();
    }

    @Override
    public double getCoercableParameterValue() {
        return continuousK;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        continuousK = Math.min(Math.max(lower, value), upper);
        discreteK = (int) Math.round(continuousK);
    }

    @Override
    public void optimize(final double logAlpha) {
        if (autoOptimize) {
            double _delta = calcDelta(logAlpha);
            _delta += Math.log(continuousK);
            setCoercableParameterValue(Math.exp(_delta));
            continuousK = Math.max(0.5000000001, continuousK);
            discreteK = (int) Math.round(continuousK);
        }
    }

    @Override
    public final String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double newContinuousK = continuousK * ratio;
        final int newDiscreteK = (newContinuousK < lower) ? 1 : (int) Math.round(newContinuousK); 

        if (newDiscreteK != discreteK && (prob < 0.10 || prob > 0.40)) {
            return String.format("Try setting 'k' to about %d", newDiscreteK);
        } else {
            return "";
        }
    }

    protected void setLimits(final int lower, final int upper) {
        this.lower = (double) lower;
        this.upper = (double) upper;
        this.lower -= 0.4999999999;
        this.upper += 0.4999999999;
    }

    // chooses 'k' numbers from between 0 and n - 1
    // i.e., random sampling without replacement
    // this is a Durstenfeld shuffle which terminates after k loops
    protected int[] chooseK(final int n) {
        final int[] sequential = new int[n];

        for (int i = 0; i < n; i++) sequential[i] = i;

        for (int i = 0; i < discreteK; i++) {
            final int j = Randomizer.nextInt(n - i);
            if (j > 0) { // swap [i] with [i + j]
                final int i_temp = sequential[i];
                final int iplusj = i + j;
                sequential[i] = sequential[iplusj];
                sequential[iplusj] = i_temp;
            }
        }

        final int[] sample = Arrays.copyOf(sequential, discreteK);

        return sample;
    }
}
