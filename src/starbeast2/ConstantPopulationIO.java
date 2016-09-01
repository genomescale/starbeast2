package starbeast2;


import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulationIO extends MultispeciesPopulationModel {
    public Input<RealParameter> populationShapeInput = new Input<RealParameter>("populationShape", "Shape of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);
    public Input<RealParameter> populationMeanInput = new Input<RealParameter>("populationMean", "Mean of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);

    RealParameter invgammaShape;
    RealParameter invgammaMean;
    
    // the following arrays store results per genetree per species branch
    // use flattened matrix for performance
    double [] branchGamma;
    double [] storedBranchGamma;
    int [] branchQ;
    int [] storedBranchQ;
    

    @Override
    public void initAndValidate() {
        invgammaShape = populationShapeInput.get();
        invgammaMean = populationMeanInput.get();
    }
    
    @Override
    public int getGeneTreeNumber(int speciesNodeCount) {
    	this.speciesNodeCount = speciesNodeCount;
    	geneTreeCount++;
    	branchGamma = new double[speciesNodeCount * geneTreeCount];
    	storedBranchGamma = new double[speciesNodeCount * geneTreeCount];
    	branchQ = new int[speciesNodeCount * geneTreeCount];
    	storedBranchQ = new int[speciesNodeCount * geneTreeCount];
    	return geneTreeCount - 1;
    }


    @Override
    public double branchLogP(int geneTreeNumber, int speciesTreeNodeNumber, Node speciesTreeNode, double perGenePloidy, double[] branchCoalescentTimes, int branchLineageCounts, int branchEventCounts) {

        //int branchQ = 0;
        double branchLogR = 0.0;
        //double branchGamma = 0.0;

        final int geneN = branchLineageCounts;
        final double[] geneCoalescentTimes = branchCoalescentTimes;
        final int geneK = branchEventCounts;
        final double genePloidy = perGenePloidy; 
        branchLogR -= geneK * Math.log(genePloidy);
        int k = speciesTreeNodeNumber * geneTreeCount + geneTreeNumber;
        branchQ [k]= geneK;

        double partialGamma = 0.0;
        for (int i = 0; i < geneK; i++) {
            partialGamma += (geneCoalescentTimes[i + 1] - geneCoalescentTimes[i]) * (geneN - i) * (geneN - (i + 1.0)) / 2.0;
        }
        
        if (geneN - geneK > 1) {
            partialGamma += (geneCoalescentTimes[geneK + 1] - geneCoalescentTimes[geneK]) * (geneN - geneK) * (geneN - (geneK + 1.0)) / 2.0;
        }

        branchGamma[k] = partialGamma / genePloidy;

        final double logP = branchLogR;

        return logP;
    }
    
    @Override
    public double commonContributionLogP() {
        final double alpha = invgammaShape.getValue();
        final double beta = invgammaMean.getValue() * (alpha - 1.0);
        double logP = 0;
        
        for (int j = 0; j < speciesNodeCount; j++) {
        	int branchQTotal = 0;
        	int k = j * geneTreeCount;
        	for (int i = 0; i < geneTreeCount; i++) {
        		branchQTotal += branchQ[k++];
        	}
        	double branchGammaTotal = 0;
        	k = j * geneTreeCount;
        	for (int i = 0; i < geneTreeCount; i++) {
        		branchGammaTotal += branchGamma[k++];
        	}
        	
	        double logGammaRatio = 0.0;
	        for (int i = 0; i < branchQTotal; i++) {
	            logGammaRatio += Math.log(alpha + i);
	        }
	        
	        logP += (alpha * Math.log(beta)) - ((alpha + branchQTotal) * Math.log(beta + branchGammaTotal)) + logGammaRatio;
        }
    	return logP;
    }

    
    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double[] perGenePloidy, double[][] branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter invgammaShape = populationShapeInput.get();
        final RealParameter invgammaMean = populationMeanInput.get();
        final double alpha = invgammaShape.getValue();
        final double beta = invgammaMean.getValue() * (alpha - 1.0);

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
    protected boolean requiresRecalculation() {
    	return super.requiresRecalculation();
    }

    
    @Override
	public void doStore() {
    	System.arraycopy(branchQ, 0, storedBranchQ, 0, branchQ.length);
    	System.arraycopy(branchGamma, 0, storedBranchGamma, 0, branchGamma.length);
	}

    @Override
	public void doRetore() {
    	int [] tmp = branchQ;
    	branchQ = storedBranchQ;
    	storedBranchQ = tmp;
    	
    	double [] tmp2 = branchGamma;
    	branchGamma = storedBranchGamma;
    	storedBranchGamma = tmp2;
	}

}
