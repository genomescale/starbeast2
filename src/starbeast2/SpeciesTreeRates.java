package starbeast2;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public abstract class SpeciesTreeRates extends CalculationNode implements BranchRateModel {
    public Input<TreeInterface> speciesTreeInput = new Input<>("tree", "Species tree to apply per-branch rates to.");

    abstract Double[] getRatesArray();

    @Override
    public double getRateForBranch(Node node) {
        final Double[] ratesArray = getRatesArray();

        return ratesArray[node.getNr()];
    }
}
