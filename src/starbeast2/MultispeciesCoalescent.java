package starbeast2;

import java.util.ArrayList;
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

    @Override
    public void initAndValidate() {
        final List<GeneTree> geneTrees = geneTreeInput.get();
        nGeneTrees = geneTrees.size();

        // initialize gene trees and store ploidy
        perGenePloidy = new double[nGeneTrees];
        for (int i = 0; i < nGeneTrees; i++) {
            final GeneTree geneTreeI = geneTrees.get(i);
            perGenePloidy[i] = geneTreeI.ploidy;
        }

        final MultispeciesPopulationModel populationModel = populationModelInput.get();
        final TreeInterface speciesTree = speciesTreeInput.get().getTree();
        speciesTreeNodeCount = speciesTree.getNodeCount();
        populationModel.initPopSizes(speciesTreeNodeCount);

        allLineageCounts = new int[speciesTreeNodeCount*nGeneTrees];
        allEventCounts = new int[speciesTreeNodeCount*nGeneTrees];
        allCoalescentTimes  = new double[speciesTreeNodeCount][nGeneTrees][];
    }

    @SuppressWarnings("unchecked")
	public double calculateLogP() {
        final TreeInterface speciesTree = speciesTreeInput.get().getTree();
        final MultispeciesPopulationModel populationModel = populationModelInput.get();

        speciesTreeNodeCount = speciesTree.getNodeCount();
        double[] speciesStartTimes = new double[speciesTreeNodeCount]; // the earlier date (rootward end)
        double[] speciesEndTimes = new double[speciesTreeNodeCount]; // the later date (tipward end)

        final List<GeneTree> geneTrees = geneTreeInput.get();

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
        
        
        // transpose gene-branch list of lists to branch-gene list of lists
        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
            final GeneTree geneTree = geneTrees.get(j);
            assert SanityChecks.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either
            if (geneTree.computeCoalescentTimes()) {
                for (int i = 0; i < speciesTreeNodeCount; i++) { // for each species tree node/branch "i"
//                    final List<Double> timesView = geneTree.getCoalescentTimes(i);
//                    final int geneBranchEventCount = timesView.size();
//                    final Double[] geneBranchCoalescentTimes = new Double[geneBranchEventCount];
//                    timesView.toArray(geneBranchCoalescentTimes);
//                    Arrays.sort(geneBranchCoalescentTimes);
                    final double [] geneBranchCoalescentTimes = geneTree.getCoalescentTimes(i);
                    final int geneBranchEventCount = geneBranchCoalescentTimes.length;

                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts[i];

                    final double[] coalescentTimesIJ = new double[geneBranchEventCount + 2];
                    coalescentTimesIJ[0] = speciesEndTimes[i];
                    if (geneBranchEventCount > 0) {
                    	System.arraycopy(geneBranchCoalescentTimes, 0, coalescentTimesIJ, 1, geneBranchEventCount);
                    }
//                    for (int k = 0; k < geneBranchEventCount; k++) {
//                        coalescentTimesIJ[k + 1] = geneBranchCoalescentTimes[k];
//                    }
                    coalescentTimesIJ[geneBranchEventCount + 1] = speciesStartTimes[i];

                    final int k = i * nGeneTrees + j;
                    allLineageCounts[k] = geneBranchLineageCount;
                    allEventCounts[k] = geneBranchEventCount;
                    allCoalescentTimes[i][j]=coalescentTimesIJ;
                }
            } else { // this gene tree IS NOT compatible with the species tree
                logP = Double.NEGATIVE_INFINITY;
                return logP;
            }
        }

        logP = 0.0;
        int[] branchLineageCounts = new int[nGeneTrees];
        int[] branchEventCounts = new int[nGeneTrees];
        for (int i = 0; i < speciesTreeNodeCount; i++) {
            final Node speciesTreeNode = speciesTree.getNode(i); 
            final double[][] branchCoalescentTimes = allCoalescentTimes[i];
            final int k = i * nGeneTrees;
            System.arraycopy(allLineageCounts, k, branchLineageCounts, 0, nGeneTrees);
            System.arraycopy(allEventCounts, k, branchEventCounts, 0, nGeneTrees);
            final double branchLogP = populationModel.branchLogP(i, speciesTreeNode, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
            logP += branchLogP;

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
