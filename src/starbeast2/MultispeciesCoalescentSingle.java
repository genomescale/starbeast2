package starbeast2;

import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

/**
* @author Remco Bouckaert
* @author Joseph Heled
* @author Huw Ogilvie
 */

@Description("Calculates probability of a single gene tree conditioned on a species tree (the multi-species coalescent).")
public class MultispeciesCoalescentSingle extends Distribution {
    public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree.", Validate.REQUIRED);
    public Input<GeneTree> geneTreeInput = new Input<>("geneTree", "Gene tree within the species tree.", Validate.REQUIRED);
    public Input<MultispeciesPopulationModel> populationModelInput = new Input<>("populationModel", "The species tree population model.", Validate.REQUIRED);

    private int speciesTreeNodeCount;
    private double perGenePloidy;

    TreeInterface speciesTree;
    
    double [] logPBranchContribution;
    double [] storedLogPBranchContribution;
    
    double[] speciesStartTimes;
    double[] speciesEndTimes;
    GeneTree geneTrees;

    // used for optimising pop model calculation
    int geneTreeNumber = 0;
    
    @Override
    public void store() {
    	System.arraycopy(logPBranchContribution, 0, storedLogPBranchContribution, 0, logPBranchContribution.length);
    	super.store();
    }
    
    @Override
    public void restore() {
    	double [] tmp = logPBranchContribution;
    	logPBranchContribution = storedLogPBranchContribution;
    	storedLogPBranchContribution = tmp;
    	
    	super.restore();
    }
    
    @Override
    public void initAndValidate() {
        geneTrees = geneTreeInput.get();

        // initialize gene trees and store ploidy
        perGenePloidy = geneTrees.ploidy;

        final MultispeciesPopulationModel populationModel = populationModelInput.get();
        speciesTree = speciesTreeInput.get().getTree();
        speciesTreeNodeCount = speciesTree.getNodeCount();
        populationModel.initPopSizes(speciesTreeNodeCount);

        logPBranchContribution = new double[speciesTreeNodeCount];
        storedLogPBranchContribution = new double[speciesTreeNodeCount];

        geneTreeNumber = populationModel.getGeneTreeNumber(speciesTreeNodeCount);
    }
    
    @Override
	public double calculateLogP() {
        final MultispeciesPopulationModel populationModel = populationModelInput.get();

        speciesTreeNodeCount = speciesTree.getNodeCount();
        speciesStartTimes = new double[speciesTreeNodeCount]; // the earlier date (rootward end)
        speciesEndTimes = new double[speciesTreeNodeCount]; // the later date (tipward end)

        geneTrees = geneTreeInput.get();

        assert SanityChecks.checkTreeSanity(speciesTree.getRoot()); // species tree should not be insane

        for (int i = 0; i < speciesTreeNodeCount; i++) {
            final Node speciesNode = speciesTree.getNode(i);
            final Node parentNode = speciesNode.getParent();

            speciesEndTimes[i] = speciesNode.getHeight();
            
            if (parentNode == null) {
                speciesStartTimes[i] = Double.POSITIVE_INFINITY;
            } else {
                speciesStartTimes[i] = parentNode.getHeight();
            }
        }
        
        
        final GeneTree geneTree = geneTrees;
        assert SanityChecks.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either
        logP = 0.0;
        boolean needsUpdate = populationModel.isDirtyCalculation();
        if (geneTree.computeCoalescentTimes()) {
            for (int i = 0; i < speciesTreeNodeCount; i++) { // for each species tree node/branch "i"
            	if (needsUpdate || geneTree.isDirtyBranch(i)) {
	                final double [] geneBranchCoalescentTimes = geneTree.getCoalescentTimes(i);
	                final int geneBranchEventCount = geneBranchCoalescentTimes.length;
	                final double [] coalescentTimesIJ = new double[geneBranchEventCount + 2];
	                coalescentTimesIJ[0] = speciesEndTimes[i];
	                if (geneBranchEventCount > 0) {
	                	System.arraycopy(geneBranchCoalescentTimes, 0, coalescentTimesIJ, 1, geneBranchEventCount);
	                }
	                coalescentTimesIJ[geneBranchEventCount + 1] = speciesStartTimes[i];
	
	                final int geneBranchLineageCount = geneTree.coalescentLineageCounts[i];

		            final Node speciesTreeNode = speciesTree.getNode(i); 
		            final double branchLogP = populationModel.branchLogP(geneTreeNumber, i, speciesTreeNode, 
		            		perGenePloidy, coalescentTimesIJ, geneBranchLineageCount, geneBranchEventCount);
	        		logPBranchContribution[i] = branchLogP;
		            logP += branchLogP;            	
            	} else {
            		logP += logPBranchContribution[i];
            	}
            }
        } else { // this gene tree IS NOT compatible with the species tree
            logP = Double.NEGATIVE_INFINITY;
        }
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
