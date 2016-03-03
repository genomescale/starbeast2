package starbeast2;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

@Description("Uses continuous rates stored in clock space, which should have a standard normal prior distribution." +
"These are scaled to have a spread of 'stdev' and a mean in real space of 1.")
public class RealRates extends BranchRateModel.Base implements SpeciesTreeRates {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
    final public Input<Boolean> estimateRootInput = new Input<>("estimateRoot", "Estimate rate of the root branch.", false);
    final public Input<RealParameter> treeRatesInput = new Input<>("treeRates", "Per-branch rates. Must have a log standard normal prior distribution.", Input.Validate.REQUIRED);

    private double[] realRatesArray;
    private double[] storedRatesArray;
    private int nEstimatedRates;
    private int rootNodeNumber;
    private boolean estimateRoot;
    private boolean needsUpdate;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = meanRateInput.isDirty() || treeRatesInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
        System.arraycopy(realRatesArray, 0, storedRatesArray, 0, realRatesArray.length);
        super.store();
    }

    @Override
    public void restore() {
        double[] tmpRatesArray = realRatesArray;
        realRatesArray = storedRatesArray;
        storedRatesArray = tmpRatesArray;
        super.restore();
    }

    @Override
    public void initAndValidate() {
        final RealParameter branchRates = treeRatesInput.get();
        final TreeInterface speciesTree = treeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();
        estimateRoot = estimateRootInput.get().booleanValue();
        rootNodeNumber = speciesTree.getRoot().getNr();
        realRatesArray = new double[speciesNodes.length];
        storedRatesArray = new double[speciesNodes.length];

        if (estimateRoot) {
            nEstimatedRates = speciesNodes.length;
        } else {
            nEstimatedRates = speciesNodes.length - 1;
        }

        branchRates.setDimension(nEstimatedRates);

        needsUpdate = true;
    }

    private void update() {
        final RealParameter estimatedMeanParameter = meanRateInput.get();
        Double estimatedMean;
        if (estimatedMeanParameter == null) {
            estimatedMean = 1.0;
        } else {
            estimatedMean = estimatedMeanParameter.getValue();
        }

        final Double[] treeRatesArray = treeRatesInput.get().getValues();
        for (int i = 0; i < nEstimatedRates; i++) {
            realRatesArray[i] = treeRatesArray[i] * estimatedMean;
        }

        if (!estimateRoot) realRatesArray[rootNodeNumber] = estimatedMean;
    }

    @Override
    public double[] getRatesArray() {
        if (needsUpdate) {
            synchronized (this) {
                update();
                needsUpdate = false;
            }
        }

        return realRatesArray;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate) {
            synchronized (this) {
                update();
                needsUpdate = false;
            }
        }

        return realRatesArray[node.getNr()];
    }
}
