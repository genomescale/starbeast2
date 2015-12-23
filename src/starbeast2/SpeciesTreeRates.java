package starbeast2;

import beast.core.CalculationNode;

public abstract class SpeciesTreeRates extends CalculationNode {
    abstract Double[] getRatesArray();
}
