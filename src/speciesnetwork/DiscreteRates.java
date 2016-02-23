package speciesnetwork;

import org.apache.commons.math.MathException;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.math.distributions.LogNormalDistributionModel;

public class DiscreteRates extends BranchRateModel.Base implements SpeciesNetworkRates {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
    final public Input<Boolean> estimateRootInput = new Input<>("estimateRoot", "Estimate rate of the root branch.", false);
    final public Input<IntegerParameter> branchRatesInput = new Input<>("rates", "Discrete per-branch rates.", Input.Validate.REQUIRED);
    final public Input<LogNormalDistributionModel> rateDistributionInput = new Input<>("distr", "The distribution governing the rates among branches. Must have mean of 1.", Input.Validate.REQUIRED);

    private int nBins;
    private double currentLogNormalStdev;
    private double storedLogNormalStdev;
    private double[] binRates;
    private double[] storedBinRates;
    private double[] ratesArray;
    private double[] storedRatesArray;

    private int nEstimatedRates;
    private int rootNodeNumber;
    private boolean estimateRoot;
    private boolean needsUpdate;
    private boolean binRatesNeedsUpdate;

    @Override
    public boolean requiresRecalculation() {
        // This check is not reliable
        /* if (rateDistributionInput.get().isDirtyCalculation()) {
            binRatesNeedsUpdate = true;
        } else {
            binRatesNeedsUpdate = false;
        } */

        final double proposedLogNormalStdev = rateDistributionInput.get().SParameterInput.get().getValue();
        if (proposedLogNormalStdev != currentLogNormalStdev) {
            binRatesNeedsUpdate = true;
            currentLogNormalStdev = proposedLogNormalStdev;
        } else {
            binRatesNeedsUpdate = false;
        }

        needsUpdate = binRatesNeedsUpdate || branchRatesInput.isDirty() || meanRateInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
        storedLogNormalStdev = currentLogNormalStdev;
        System.arraycopy(binRates, 0, storedBinRates, 0, binRates.length);
        System.arraycopy(ratesArray, 0, storedRatesArray, 0, ratesArray.length);
        super.store();
    }

    @Override
    public void restore() {
        double tmpLogNormalStdev = currentLogNormalStdev;
        double[] tmpBinRates = binRates;
        double[] tmpRatesArray = ratesArray;

        currentLogNormalStdev = storedLogNormalStdev;
        binRates = storedBinRates;
        ratesArray = storedRatesArray;
        
        storedLogNormalStdev = tmpLogNormalStdev;
        storedBinRates = tmpBinRates;
        storedRatesArray = tmpRatesArray;

        super.restore();
    }

    @Override
    public void initAndValidate() {
        final IntegerParameter branchRates = branchRatesInput.get();
        final TreeInterface speciesTree = treeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();
        estimateRoot = estimateRootInput.get().booleanValue();
        rootNodeNumber = speciesTree.getRoot().getNr();
        ratesArray = new double[speciesNodes.length];
        storedRatesArray = new double[speciesNodes.length];

        if (estimateRoot) {
            nEstimatedRates = speciesNodes.length;
        } else {
            nEstimatedRates = speciesNodes.length - 1;
        }

        nBins = nEstimatedRates;

        branchRates.setDimension(nEstimatedRates);
        branchRates.setLower(0);
        branchRates.setUpper(nBins - 1);

        currentLogNormalStdev = -1.0;
        storedLogNormalStdev = -1.0;

        binRates = new double[nBins];
        storedBinRates = new double[nBins];

        binRatesNeedsUpdate = true;
        needsUpdate = true;
    }

    private void update() {
        final LogNormalDistributionModel rateDistribution = rateDistributionInput.get();

        if (binRatesNeedsUpdate) {
            try {
                for (int i = 0; i < nBins; i++) {
                    binRates[i] = rateDistribution.inverseCumulativeProbability((i + 0.5) / nBins);
                }
            } catch (MathException e) {
                throw new RuntimeException("Failed to compute inverse cumulative probability!");
            }
        }

        Double estimatedMean;
        final RealParameter estimatedMeanParameter = meanRateInput.get();
        if (estimatedMeanParameter == null) {
            estimatedMean = 1.0;
        } else {
            estimatedMean = estimatedMeanParameter.getValue();
        }

        final Integer[] branchRatePointers = branchRatesInput.get().getValues();
        for (int i = 0; i < nEstimatedRates; i++) {
            int b = branchRatePointers[i];
            ratesArray[i] = binRates[b] * estimatedMean;
        }

        if (!estimateRoot) ratesArray[rootNodeNumber] = estimatedMean;

        /* StringBuffer x = new StringBuffer();
        x.append(treeInput.get().getID());
        for (int i = 0; i < ratesArray.length; i++) {
            x.append(" ");
            x.append(ratesArray[i]);
        }
        System.out.println(x); */
    }

    @Override
    public double[] getRatesArray() {
        if (needsUpdate) {
            synchronized (this) {
                update();
                needsUpdate = false;
            }
        }

        return ratesArray;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate) {
            synchronized (this) {
                update();
                needsUpdate = false;
            }
        }

        assert ratesArray[node.getNr()] > 0.0;
        return ratesArray[node.getNr()];
    }
}
