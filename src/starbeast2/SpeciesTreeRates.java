package starbeast2;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.tree.TreeInterface;

public abstract class SpeciesTreeRates extends CalculationNode {
    public Input<TreeInterface> speciesTreeInput = new Input<>("tree", "Species tree to apply per-branch rates to.");
    abstract Double[] getRatesArray();
}
