package starbeast2;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

/**
 * Changes value of an internal tree branch rate by magnitude delta,
 * and balances the sum of rates by changing both child branch rates
 * by negative half delta.
 *
 * @author Huw A. Ogilvie
 * 
 */

public class ChildrenRateExchange extends Operator {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "The tree with branch rates.", Validate.REQUIRED);
    final public Input<RealParameter> treeRatesInput = new Input<>("treeRates", "The branch rates.", Validate.REQUIRED);
    final public Input<Double> deltaInput = new Input<>("delta", "Magnitude of change for two randomly picked values.", 1.0);

    private TreeInterface tree;

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
        nInternalNodes = treeRates.getDimension() - nLeafNodes;
        lowerBound = treeRates.getLower();
        upperBound = treeRates.getUpper();
    }

    @Override
    public double proposal() {
        tree = treeInput.get();

        // symmetric proposal distribution
        final double delta = (Randomizer.nextDouble() - 0.5) * deltaScaleFactor * 2.0;

        final RealParameter treeRates = treeRatesInput.get();
        final int parentNodeNumber = nLeafNodes + Randomizer.nextInt(nInternalNodes);
        final double originalParentBranchRate = treeRates.getValue(parentNodeNumber);
        final double newParentBranchRate = originalParentBranchRate + delta;
        if (newParentBranchRate < lowerBound || newParentBranchRate > upperBound) {
            return Double.NEGATIVE_INFINITY;
        } else {
            treeRates.setValue(parentNodeNumber, newParentBranchRate);
        }

        final Node parentNode = tree.getNode(parentNodeNumber);
        for (Node childNode: parentNode.getChildren()) {
            final int childNodeNumber = childNode.getNr();
            final double originalChildBranchRate = treeRates.getValue(childNodeNumber);
            final double newChildBranchRate = originalChildBranchRate - (delta * 0.5);
            if (newChildBranchRate < lowerBound || newChildBranchRate > upperBound) {
                return Double.NEGATIVE_INFINITY;
            } else {
                treeRates.setValue(childNodeNumber, newChildBranchRate);
            }
        }

        return 0.0;
    }
}
