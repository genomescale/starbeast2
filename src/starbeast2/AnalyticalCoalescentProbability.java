package starbeast2;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.core.State;
import beast.evolution.tree.Node;

/**
* @author Remco Bouckaert
* @author Joseph Heled
* @author Huw Ogilvie
 */

@Description("Calculates probability of gene trees conditioned on a species tree (the multi-species coalescent).")
public class AnalyticalCoalescentProbability extends Distribution {
    final public Input<List<GeneTree>> geneTreesInput = new Input<>("geneTree", "Gene tree within the species tree.", new ArrayList<>());
    final public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree.", Validate.REQUIRED);
    final public Input<RealParameter> populationShapeInput = new Input<>("populationShape", "Shape of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);
    final public Input<RealParameter> populationMeanInput = new Input<>("populationMean", "Mean of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);
    final public Input<BooleanParameter> disableInput = new Input<>("disable", "Disable calculating analytical probability.", Validate.REQUIRED);

    private SpeciesTree speciesTree;
    private List<GeneTree> geneTrees;
    private RealParameter invGammaShape;
    private RealParameter invGammaMean;

    private int nGeneTrees;
    private int speciesNodeCount;
    private double[] perGenePloidy;

    private double alpha;
    private double beta;

    private double storedAlpha;
    private double storedBeta;

    private int[] allLineageCounts;
    private int[] allEventCounts;
    private double[][] allCoalescentTimes;

    private int[] storedLineageCounts;
    private int[] storedEventCounts;
    private double[][] storedCoalescentTimes;

    private double[] perBranchLogP;
    private double[] storedPerBranchLogP;

    private boolean dontCalculate;

    @Override
    public void store() {
        if (dontCalculate) return;

        storedAlpha = alpha;
        storedBeta = beta;

        System.arraycopy(allLineageCounts, 0, storedLineageCounts, 0, allLineageCounts.length);
        System.arraycopy(allEventCounts, 0, storedEventCounts, 0, allEventCounts.length);
        for (int i = 0; i < allCoalescentTimes.length; i++)
            System.arraycopy(allCoalescentTimes, 0, storedCoalescentTimes, 0, allCoalescentTimes.length);
        System.arraycopy(perBranchLogP, 0, storedPerBranchLogP, 0, perBranchLogP.length);

        super.store();
    }

    @Override
    public void restore() {
        if (dontCalculate) return;

        double tmpAlpha = alpha;
        double tmpBeta = beta;
        int[] tmpLineageCounts = allLineageCounts;
        int[] tmpEventCounts = allEventCounts;
        double[][] tmpCoalescentTimes = allCoalescentTimes;
        double[] tmpPerBranchLogP = perBranchLogP;

        alpha = storedAlpha;
        beta = storedBeta;
        allLineageCounts = storedLineageCounts;
        allEventCounts = storedEventCounts;
        allCoalescentTimes = storedCoalescentTimes;
        perBranchLogP = storedPerBranchLogP;

        storedAlpha = tmpAlpha;
        storedBeta = tmpBeta;
        storedLineageCounts = tmpLineageCounts;
        storedEventCounts = tmpEventCounts;
        storedCoalescentTimes = tmpCoalescentTimes;
        storedPerBranchLogP = tmpPerBranchLogP;

        super.restore();
    }

    @Override
    public void initAndValidate() {
        dontCalculate = disableInput.get().getValue();
        if (dontCalculate) return;

        speciesTree = speciesTreeInput.get();
        speciesNodeCount = speciesTree.getNodeCount();

        geneTrees = geneTreesInput.get();
        nGeneTrees = geneTrees.size();

        allLineageCounts = new int[speciesNodeCount*nGeneTrees];
        allEventCounts = new int[speciesNodeCount*nGeneTrees];
        allCoalescentTimes = new double[speciesNodeCount*nGeneTrees][];
        perBranchLogP = new double[speciesNodeCount];

        perGenePloidy = new double[nGeneTrees];
        for (int geneI = 0; geneI < nGeneTrees; geneI++) {
            final GeneTree geneTreeI = geneTrees.get(geneI);
            perGenePloidy[geneI] = geneTreeI.ploidy;
        }
    }

    @Override
	public double calculateLogP() {
        logP = 0.0;
        if (dontCalculate) return logP;

        // recompute coalescent counts and times
        geneTrees = geneTreesInput.get();
        for (int geneI = 0; geneI < nGeneTrees; geneI++) {
            final GeneTree geneTree = geneTrees.get(geneI);
            if (!geneTree.computeCoalescentTimes()) {
                assert SanityChecks.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane
                logP = Double.NEGATIVE_INFINITY;
                return logP;
            }
        }

        // need to recompute all branches if the parameters of the prior distribution have changed
        invGammaShape = populationShapeInput.get();
        invGammaMean = populationMeanInput.get();

        boolean updatedPrior = false;
        if (invGammaShape.isDirty(0) || alpha == 0.0) {
            alpha = invGammaShape.getValue();
            updatedPrior = true;
        }
        if (invGammaMean.isDirty(0) || beta == 0.0) {
            beta = invGammaMean.getValue() * (alpha - 1.0);
            updatedPrior = true;
        }

        final int[] branchLineageCounts = new int[nGeneTrees];
        final int[] branchEventCounts = new int[nGeneTrees];
        final double[][] branchCoalescentTimes = new double[nGeneTrees][];

        // rebuild species start and end times (if necessary)
        final Node[] speciesTreeNodes = speciesTree.getNodesAsArray();
        int nodeGeneI = 0;
        for (int nodeI = 0; nodeI < speciesNodeCount; nodeI++) {
            final Node speciesNode = speciesTreeNodes[nodeI];
            final Node parentNode = speciesNode.getParent();

            final double speciesEndTime = speciesNode.getHeight();
            final double speciesStartTime = (parentNode == null) ? Double.POSITIVE_INFINITY : parentNode.getHeight();

            boolean dirtyBranch = false;
            for (int geneI = 0; geneI < nGeneTrees; geneI++) {
                final GeneTree geneTree = geneTrees.get(geneI);

                if (geneTree.isDirtyBranch(geneI)) {
                    dirtyBranch = true;

                    final double [] tmpCoalescentTimes = geneTree.getCoalescentTimes(nodeI);
                    final int geneBranchEventCount = tmpCoalescentTimes.length;
                    final double [] geneBranchCoalescentTimes = new double[geneBranchEventCount + 2];
                    geneBranchCoalescentTimes[0] = speciesEndTime;
                    if (geneBranchEventCount > 0)
                        System.arraycopy(tmpCoalescentTimes, 0, geneBranchCoalescentTimes, 1, geneBranchEventCount);
                    geneBranchCoalescentTimes[geneBranchEventCount + 1] = speciesStartTime;

                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts[nodeI];
                    allLineageCounts[nodeGeneI] = geneBranchLineageCount;
                    allEventCounts[nodeGeneI] = geneBranchEventCount;
                    allCoalescentTimes[nodeGeneI] = geneBranchCoalescentTimes;
                }

                branchLineageCounts[geneI] = allLineageCounts[nodeGeneI];
                branchEventCounts[geneI] = allEventCounts[nodeGeneI];
                branchCoalescentTimes[geneI] = allCoalescentTimes[nodeGeneI];

                nodeGeneI++;
            }

            if (updatedPrior || dirtyBranch)
                perBranchLogP[nodeI] = branchLogP(alpha, beta, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);

            logP += perBranchLogP[nodeI];
        }

        return logP;
    }

    static private double branchLogP(double alpha, double beta, double[] perGenePloidy, double[][] branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
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

        double logGammaRatio = 0.0;
        for (int i = 0; i < branchQ; i++) {
            logGammaRatio += Math.log(alpha + i);
        }

        final double logP = branchLogR + (alpha * Math.log(beta)) - ((alpha + branchQ) * Math.log(beta + branchGamma)) + logGammaRatio;

        return logP;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
}
