package starbeast2;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

/**
 * Changes value of an internal tree branch rate by magnitude delta,
 * and balances the sum of rates by changing one of the child branch rates
 * by negative delta.
 *
 * @author Huw A. Ogilvie
 * 
 */

public class ChildRateExchange extends AdaptiveOperator {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "The tree with branch rates.", Validate.REQUIRED);
    final public Input<RealParameter> treeRatesInput = new Input<>("treeRates", "The branch rates.", Validate.REQUIRED);
    final public Input<Double> deltaInput = new Input<>("delta", "Magnitude of change for two randomly picked values.", 1.0);

    private TreeInterface tree;

    private int nNodes;
    private int nLeafNodes;
    private int nInternalNodes;
    private double deltaScaleFactor;
    private double lowerBound;
    private double upperBound;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        nLeafNodes = tree.getLeafNodeCount();
        deltaScaleFactor = deltaInput.get();

        final RealParameter treeRates = treeRatesInput.get();
        nNodes = treeRates.getDimension();
        nInternalNodes = nNodes - nLeafNodes;
        lowerBound = treeRates.getLower();
        upperBound = treeRates.getUpper();

        deltaScaleFactor = 2.0 * deltaInput.get() / nNodes;

        setLimits(1, nInternalNodes);
        super.initAndValidate();
    }

    // symmetric proposal distribution
    @Override
    public double proposal() {
        tree = treeInput.get();
        final RealParameter treeRates = treeRatesInput.get();
        final double delta = (Randomizer.nextDouble() - 0.5) * deltaScaleFactor;
        final double[] treeRatesArray = treeRates.getDoubleValues();
        final int[] permutation = chooseK(nInternalNodes);

        for (int i = 0; i < discreteK; i++ ) {
            final int parentNodeNumber = nLeafNodes + permutation[i];
            final Node parentNode = tree.getNode(parentNodeNumber);
            final Node childNode = Randomizer.nextBoolean() ? parentNode.getLeft() : parentNode.getRight();
            final int childNodeNumber = childNode.getNr();
            treeRatesArray[parentNodeNumber] += delta;
            treeRatesArray[childNodeNumber] -= delta;
        }

        for (int i = 0; i < nNodes; i++) {
            final double newRate = treeRatesArray[i];
            treeRates.setValue(i, newRate);
            if (newRate < lowerBound || newRate > upperBound) {
                return Double.NEGATIVE_INFINITY;
            } else {
                treeRates.setValue(i, newRate);
            }
        }

        return 0.0;
    }
}
