package starbeast2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

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

@Description("Calculates probability of gene trees conditioned on a species tree (the multi-species coalescent).")
public class MultispeciesCoalescent extends Distribution {
    public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree.", Validate.REQUIRED);
    public Input<List<GeneTree>> geneTreeInput = new Input<>("geneTree", "Gene tree within the species tree.", new ArrayList<>());
    public Input<MultispeciesPopulationModel> populationModelInput = new Input<>("populationModel", "The species tree population model.", Validate.REQUIRED);

    private int nGeneTrees;
    private int speciesTreeNodeCount;
    private double[] perGenePloidy;

    private int[] allLineageCounts;// = new ArrayList<>();
    private int[] allEventCounts;// = new ArrayList<>();
    private double[][][] allCoalescentTimes;// = new ArrayList<>();
    
    public static ExecutorService pool;
    private final List<Callable<Double>> likelihoodCallers = new ArrayList<Callable<Double>>();
    int threadCount = 1;
    double [] threadResult;
    TreeInterface speciesTree;
    
    double [] logPBranchContribution;
    double [] storedLogPBranchContribution;
    
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
        final List<GeneTree> geneTrees = geneTreeInput.get();
        nGeneTrees = geneTrees.size();
        if (nGeneTrees > 1) {
        	throw new IllegalArgumentException("Expected only 1 gene tree");
        }

        // initialize gene trees and store ploidy
        perGenePloidy = new double[nGeneTrees];
        for (int i = 0; i < nGeneTrees; i++) {
            final GeneTree geneTreeI = geneTrees.get(i);
            perGenePloidy[i] = geneTreeI.ploidy;
        }

        final MultispeciesPopulationModel populationModel = populationModelInput.get();
        speciesTree = speciesTreeInput.get().getTree();
        speciesTreeNodeCount = speciesTree.getNodeCount();
        populationModel.initPopSizes(speciesTreeNodeCount);

        allLineageCounts = new int[speciesTreeNodeCount*nGeneTrees];
        allEventCounts = new int[speciesTreeNodeCount*nGeneTrees];
        allCoalescentTimes  = new double[speciesTreeNodeCount][nGeneTrees][];
        
        if (threadCount> 1) {
	        pool = Executors.newFixedThreadPool(threadCount);
        }
        int d = geneTrees.size();
        for (int i = 0; i < threadCount; i++) {
        	likelihoodCallers.add(new LikelihoodCaller(i, i * d / threadCount, (i+1) * d / threadCount));
        }
        threadResult = new double[threadCount];
        
        logPBranchContribution = new double[speciesTreeNodeCount];
        storedLogPBranchContribution = new double[speciesTreeNodeCount];
        
        
        geneTreeNumber = populationModel.getGeneTreeNumber(speciesTreeNodeCount);

    }
    
    
    class LikelihoodCaller implements Callable<Double> {
        private final int start;
        private final int end;
        private final int thread;

        public LikelihoodCaller(int thread, int start, int end) {
        	this.thread = thread;
            this.start = start;
            this.end = end;
        }

        public Double call() {
  		  	try {
  		        for (int j = start; j < end; j++) { // for each gene "j"
  		            final GeneTree geneTree = geneTrees.get(j);
  		            if (geneTree.treeInput.isDirty()) {
	  		            assert SanityChecks.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either
	  		            if (geneTree.computeCoalescentTimes()) {
	  		                for (int i = 0; i < speciesTreeNodeCount; i++) { // for each species tree node/branch "i"
	  		                    final double [] geneBranchCoalescentTimes = geneTree.getCoalescentTimes(i);
	  		                    final int geneBranchEventCount = geneBranchCoalescentTimes.length;
	  		                    final double [] coalescentTimesIJ = new double[geneBranchEventCount + 2];
	  		                    coalescentTimesIJ[0] = speciesEndTimes[i];
	  		                    if (geneBranchEventCount > 0) {
	  		                    	System.arraycopy(geneBranchCoalescentTimes, 0, coalescentTimesIJ, 1, geneBranchEventCount);
	  		                    }
	  		                    coalescentTimesIJ[geneBranchEventCount + 1] = speciesStartTimes[i];
	
	  		                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts[i];
	  		                    final int k = i * nGeneTrees + j;
	  		                    allLineageCounts[k] = geneBranchLineageCount;
	  		                    allEventCounts[k] = geneBranchEventCount;
	  		                    allCoalescentTimes[i][j] = coalescentTimesIJ;
	  		                }
	  		            } else { // this gene tree IS NOT compatible with the species tree
	  		            	threadResult[thread] = Double.NEGATIVE_INFINITY;
	  		                return Double.NEGATIVE_INFINITY;
	  		            }
  		            }
  		        }
  		        for (int j = start; j < end; j++) { // for each gene "j"
  		            final GeneTree geneTree = geneTrees.get(j);
  		            if (!geneTree.treeInput.isDirty()) {
	  		            assert SanityChecks.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either
	  		            if (geneTree.computeCoalescentTimes()) {
	  		                for (int i = 0; i < speciesTreeNodeCount; i++) { // for each species tree node/branch "i"
	  		                    final double [] geneBranchCoalescentTimes = geneTree.getCoalescentTimes(i);
	  		                    final int geneBranchEventCount = geneBranchCoalescentTimes.length;
	  		                    final double [] coalescentTimesIJ = new double[geneBranchEventCount + 2];
	  		                    coalescentTimesIJ[0] = speciesEndTimes[i];
	  		                    if (geneBranchEventCount > 0) {
	  		                    	System.arraycopy(geneBranchCoalescentTimes, 0, coalescentTimesIJ, 1, geneBranchEventCount);
	  		                    }
	  		                    coalescentTimesIJ[geneBranchEventCount + 1] = speciesStartTimes[i];
	
	  		                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts[i];
	  		                    final int k = i * nGeneTrees + j;
	  		                    allLineageCounts[k] = geneBranchLineageCount;
	  		                    allEventCounts[k] = geneBranchEventCount;
	  		                    allCoalescentTimes[i][j] = coalescentTimesIJ;
	  		                }
	  		            } else { // this gene tree IS NOT compatible with the species tree
	  		            	threadResult[thread] = Double.NEGATIVE_INFINITY;
	  		                return Double.NEGATIVE_INFINITY;
	  		            }
  		            }
  		        }
  		  	} catch (Exception e) {
  		  		System.err.println("Something went wrong in thread " + start + " " + end);
				e.printStackTrace();
				System.exit(0);
			}
  		  	threadResult[thread] = 0.0;
            return 0.0;
        }

    }


    double[] speciesStartTimes;
    double[] speciesEndTimes;
    List<GeneTree> geneTrees;

    
	public double calculateLogP0() {
        final TreeInterface speciesTree = speciesTreeInput.get().getTree();
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
        
        
		if (threadCount > 1) {
            try {
				pool.invokeAll(likelihoodCallers);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} else {
				try {
					likelihoodCallers.get(0).call();
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		}

    	for (double f : threadResult) {
    		if (Double.isInfinite(f)) {
    			logP = Double.NEGATIVE_INFINITY;
    			return logP;
    		}
    	}
        
//        // transpose gene-branch list of lists to branch-gene list of lists
//        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
//            final GeneTree geneTree = geneTrees.get(j);
//            assert SanityChecks.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either
//            if (geneTree.computeCoalescentTimes()) {
//                for (int i = 0; i < speciesTreeNodeCount; i++) { // for each species tree node/branch "i"
////                    final List<Double> timesView = geneTree.getCoalescentTimes(i);
////                    final int geneBranchEventCount = timesView.size();
////                    final Double[] geneBranchCoalescentTimes = new Double[geneBranchEventCount];
////                    timesView.toArray(geneBranchCoalescentTimes);
////                    Arrays.sort(geneBranchCoalescentTimes);
//                    final double [] geneBranchCoalescentTimes = geneTree.getCoalescentTimes(i);
//                    final int geneBranchEventCount = geneBranchCoalescentTimes.length;
//
//                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts[i];
//
//                    final double[] coalescentTimesIJ = new double[geneBranchEventCount + 2];
//                    coalescentTimesIJ[0] = speciesEndTimes[i];
//                    if (geneBranchEventCount > 0) {
//                    	System.arraycopy(geneBranchCoalescentTimes, 0, coalescentTimesIJ, 1, geneBranchEventCount);
//                    }
////                    for (int k = 0; k < geneBranchEventCount; k++) {
////                        coalescentTimesIJ[k + 1] = geneBranchCoalescentTimes[k];
////                    }
//                    coalescentTimesIJ[geneBranchEventCount + 1] = speciesStartTimes[i];
//
//                    final int k = i * nGeneTrees + j;
//                    allLineageCounts[k] = geneBranchLineageCount;
//                    allEventCounts[k] = geneBranchEventCount;
//                    allCoalescentTimes[i][j]=coalescentTimesIJ;
//                }
//            } else { // this gene tree IS NOT compatible with the species tree
//                logP = Double.NEGATIVE_INFINITY;
//                return logP;
//            }
//        }
//
//		}
		
        logP = 0.0;
        int[] branchLineageCounts = new int[nGeneTrees];
        int[] branchEventCounts = new int[nGeneTrees];
        
        
        
        final GeneTree geneTree = geneTrees.get(0);
        for (int i = 0; i < speciesTreeNodeCount; i++) {
	            final Node speciesTreeNode = speciesTree.getNode(i); 
	            final double[][] branchCoalescentTimes = allCoalescentTimes[i];
	            final int k = i * nGeneTrees;
	            System.arraycopy(allLineageCounts, k, branchLineageCounts, 0, nGeneTrees);
	            System.arraycopy(allEventCounts, k, branchEventCounts, 0, nGeneTrees);
	            final double branchLogP = populationModel.branchLogP(i, speciesTreeNode, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
	            if (getID().equals("speciescoalescent12") && i == 0) {
//	            	System.err.println("\n" + Arrays.toString(allCoalescentTimes[i][0]) + " " + 
//	            			Arrays.toString(branchLineageCounts) + " " + 
//	            					Arrays.toString(branchEventCounts));
	            }
	        	if (geneTree.isDirtyBranch(i) || populationModel.isDirtyCalculation()) {
	        		logPBranchContribution[i] = branchLogP;
		            logP += branchLogP;
	        	} else {
	        		if (logPBranchContribution[i] != branchLogP ) {
	        			int h = 3;
	        			h++;
	        		}
	        		logPBranchContribution[i] = branchLogP;
	        		logP += logPBranchContribution[i];
	        	}

            /* for (int j = 0; j < branchCoalescentTimes.size(); j++) {
                Double[] geneTimes = branchCoalescentTimes.get(j);
                for (int k = 0; k < geneTimes.length; k++) {
                    System.out.println(String.format("%d/%d/%d: %f", i, j, k, geneTimes[k]));
                }
                System.out.println(String.format("%d/%d: %d, %d", i, j, branchLineageCounts[j], branchEventCounts[j]));
            }
            System.out.println(String.format("%d: %f", i, branchLogP)); */
        }

        return logP;
    }

    
    
    //@Override
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
        
        
        final GeneTree geneTree = geneTrees.get(0);
        assert SanityChecks.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either
        logP = 0.0;
        int[] branchLineageCounts = new int[nGeneTrees];
        int[] branchEventCounts = new int[nGeneTrees];
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
	                final int k = i * nGeneTrees;
	                allLineageCounts[k] = geneBranchLineageCount;
	                allEventCounts[k] = geneBranchEventCount;
	                allCoalescentTimes[i][0] = coalescentTimesIJ;

		            final Node speciesTreeNode = speciesTree.getNode(i); 
		            final double[][] branchCoalescentTimes = allCoalescentTimes[i];
		            branchLineageCounts[0] = allLineageCounts[i];
		            branchEventCounts[0] = allEventCounts[i];
		            //final double branchLogP = populationModel.branchLogP(i, speciesTreeNode, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
		            final double branchLogP = populationModel.branchLogP(geneTreeNumber, i, speciesTreeNode, 
		            		perGenePloidy[0], branchCoalescentTimes[0], branchLineageCounts[0], branchEventCounts[0]);
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
