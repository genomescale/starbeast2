package starbeast2;


import java.text.DecimalFormat;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

@Description("Calculates probability of coalescence events within a branch based on a demographic function.")
public abstract class MultispeciesPopulationModel extends CalculationNode {
    public abstract double branchLogP(int geneTreeNumber, int speciesTreeNodeNumber, Node speciesTreeNode, double perGenePloidy, double[] branchCoalescentTimes, int branchLineageCounts, int branchEventCounts);

    /* return unique gene tree number, to be used as argument to branchLogP */
    int geneTreeCount = 0;
    int speciesNodeCount = 0;

    public int getGeneTreeNumber(int speciesNodeCount) {
    	this.speciesNodeCount = speciesNodeCount;
    	geneTreeCount++;
    	return geneTreeCount - 1;
    }

    // Sets the appropriate dimension size of each population size state node
    // To successfully resume from a saved state, this must be called via an initAndValidate method
    public void initPopSizes(final int nSpeciesBranches) {
    }

    // Sets model-compatible default population sizes
    // To successfully begin a run, this must be called from a StateNodeInitializer
    public void initPopSizes(final double initialPopSizes) {
    }

    // Per-branch population size information which will be added to a Newick string.
    // If no information is available, do not override the superclass method
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
    }

	public void doStore() {
		// override if something needs to be stored for commonContributionLogP()
	}

	public void doRetore() {
		// override if something needs to be stored for commonContributionLogP()
	}

    protected static double constantLogP(double popSize, double perGenePloidy, double[] branchCoalescentTimes, int branchLineageCounts, int branchEventCounts) {
        int branchQ = 0;
        double branchLogR = 0.0;
        double branchGamma = 0.0;

        final int geneN = branchLineageCounts;
        final double[] geneCoalescentTimes = branchCoalescentTimes;
        final int geneK = branchEventCounts;
        final double genePloidy = perGenePloidy; 
        branchLogR -= geneK * Math.log(genePloidy);
        branchQ += geneK;

        double partialGamma = 0.0;
        for (int i = 0; i < geneK; i++) {
            partialGamma += (geneCoalescentTimes[i + 1] - geneCoalescentTimes[i]) * (geneN - i) * (geneN - (i + 1.0)) / 2.0;
        }
        
        if (geneN - geneK > 1) {
            partialGamma += (geneCoalescentTimes[geneK + 1] - geneCoalescentTimes[geneK]) * (geneN - geneK) * (geneN - (geneK + 1.0)) / 2.0;
        }

        branchGamma += partialGamma / genePloidy;

        final double logP = branchLogR - (branchQ * Math.log(popSize)) - (branchGamma / popSize);
        return logP;
    }
}
