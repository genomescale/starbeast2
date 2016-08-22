package starbeast2;


import java.text.DecimalFormat;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
* @author Joseph Heled
 */

@Description("Calculates probability of coalescence events within a branch based on a demographic function.")
public abstract class MultispeciesPopulationModel extends CalculationNode {
    public abstract double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double[] perGenePloidy, double[][] branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts);

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

    protected static double constantLogP(double popSize, double[] perGenePloidy, double[][] branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final int nGenes = perGenePloidy.length;

        int branchQ = 0;
        double branchLogR = 0.0;
        double branchGamma = 0.0;

        for (int j = 0; j < nGenes; j++) {
            final int geneN = branchLineageCounts[j];
            final double[] geneCoalescentTimes = branchCoalescentTimes[j];
            final int geneK = branchEventCounts[j];
            final double genePloidy = perGenePloidy[j]; 
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
        }

        final double logP = branchLogR - (branchQ * Math.log(popSize)) - (branchGamma / popSize);
        return logP;
    }

    // copied from *BEAST v2.3, with small modifications to work with starbeast2
    protected static double linearLogP(double topPopSize, double tipPopSize, double[] perGenePloidy, double[][] branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final int nGenes = perGenePloidy.length;

        double logP = 0.0;
        for (int j = 0; j < nGenes; j++) {
            final double fPopSizeTop = topPopSize * perGenePloidy[j];
            final double fPopSizeBottom = tipPopSize * perGenePloidy[j];
            final double[] fTimes = branchCoalescentTimes[j];
            final int k = branchEventCounts[j];
            final int nLineagesBottom = branchLineageCounts[j];

            final double d5 = fPopSizeTop - fPopSizeBottom;
            final double fTime0 = fTimes[0];
            final double a = d5 / (fTimes[k + 1] - fTime0);
            final double b = fPopSizeBottom;

            if (Math.abs(d5) < 1e-10) {
                // use approximation for small values to bypass numerical instability
                for (int i = 0; i <= k; i++) {
                    final double fTimeip1 = fTimes[i + 1];
                    final double fPopSize = a * (fTimeip1 - fTime0) + b;
                    if( i < k ) {
                        logP += -Math.log(fPopSize);
                    }
                    // slope = 0, so population function is constant

                    final int i1 = nLineagesBottom - i;
                    logP -= (i1 * (i1 - 1.0) / 2.0) * (fTimeip1 - fTimes[i]) / fPopSize;
                }
            } else {
                final double vv = b - a * fTime0;
                for (int i = 0; i <= k; i++) {
                    final double fPopSize = a * fTimes[i + 1] + vv;
                    if( i < k ) {
                        logP += -Math.log(fPopSize);
                    }
                    final double f = fPopSize / (a * fTimes[i] + vv);
    
                    final int i1 = nLineagesBottom - i;
                    logP += -(i1 * (i1 - 1.0) / 2.0) * Math.log(f) / a;
                }
            }
        }

        return logP;
    }
}
