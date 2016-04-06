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
    protected int n; // as in binomial coefficient
    protected int k; // as in binomial coefficient
    private double continuousK;

    @Override
    public void initAndValidate() {
        k = kInput.get();
        continuousK = k;
        autoOptimize = autoOptimizeInput.get();
    }

    @Override
    public double getCoercableParameterValue() {
        return continuousK;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        continuousK = value;
    }

    @Override
    public void optimize(final double logAlpha) {
        if (autoOptimize) {
            double _delta = calcDelta(logAlpha);
            _delta += Math.log(continuousK);
            continuousK = Math.exp(_delta);
            continuousK = Math.max(0.5000000001, continuousK);
            k = (int) Math.round(continuousK);
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
        final int newK = (newContinuousK < 0.5000000001) ? 1 : (int) Math.round(newContinuousK); 

        if (prob < 0.10 || prob > 0.40) {
            return String.format("Try setting 'k' to about %d", newK);
        } else {
            return "";
        }
    }

    // chooses 'k' numbers from between 0 and n - 1
    // i.e., random sampling without replacement
    // this is a Durstenfeld shuffle which terminates after k loops
    protected int[] chooseK() {
        final int[] sequential = new int[n];

        for (int i = 0; i < n; i++) sequential[i] = i;

        for (int i = 0; i < k; i++) {
            final int j = Randomizer.nextInt(n - i);
            if (j > 0) { // swap [i] with [i + j]
                final int i_temp = sequential[i];
                final int iplusj = i + j;
                sequential[i] = sequential[iplusj];
                sequential[iplusj] = i_temp;
            }
        }

        final int[] sample = Arrays.copyOf(sequential, k);

        return sample;
    }
}
