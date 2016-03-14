package speciesnetwork;

import java.text.DecimalFormat;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;

/**
 * @author Huw Ogilvie
 * @author Joseph Heled
 */

@Description("Calculates probability of coalescence events within a branch based on a demographic function.")
public abstract class PopulationSizeModel extends CalculationNode {
    abstract public double branchLogP(int speciesNetworkPopNumber, double[] perGenePloidy,
                                      List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts);

    // Sets the appropriate dimension size of each population size state node
    // To successfully resume from a saved state, this must be called via an initAndValidate method
    abstract public void initPopSizes(final int nPopulation);

    // Sets model-compatible default population sizes
    // To successfully begin a run, this must be called from a StateNodeInitializer
    abstract public void initPopSizes(final double initialPopSizes);

    // Per-branch population size information which will be added to a Newick string.
    // If no information is available, do not override the superclass method
    abstract public void serialize(NetworkNode speciesNetworkNode, StringBuffer buf, DecimalFormat df);
}
