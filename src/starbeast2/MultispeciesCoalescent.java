package starbeast2;

import java.util.ArrayList;
import java.util.Arrays;
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

    final private List<int[]> allLineageCounts = new ArrayList<>();
    final private List<int[]> allEventCounts = new ArrayList<>();
    final private List<List<Double[]>> allCoalescentTimes = new ArrayList<>();

    final static SanityChecks sc = new SanityChecks();

    @Override
    public void initAndValidate() throws Exception {
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
    }

    public double calculateLogP() throws Exception {
        final TreeInterface speciesTree = speciesTreeInput.get().getTree();
        final MultispeciesPopulationModel populationModel = populationModelInput.get();

        speciesTreeNodeCount = speciesTree.getNodeCount();
        double[] speciesStartTimes = new double[speciesTreeNodeCount]; // the earlier date (rootward end)
        double[] speciesEndTimes = new double[speciesTreeNodeCount]; // the later date (tipward end)

        final List<GeneTree> geneTrees = geneTreeInput.get();

        assert sc.checkTreeSanity(speciesTree.getRoot()); // species tree should not be insane

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

        allLineageCounts.clear();
        allEventCounts.clear();
        allCoalescentTimes.clear();

        for (int i = 0; i < speciesTreeNodeCount; i++) {
            allLineageCounts.add(new int[nGeneTrees]);
            allEventCounts.add(new int[nGeneTrees]);
            allCoalescentTimes.add(new ArrayList<>());
        }

        // transpose gene-branch list of lists to branch-gene list of lists
        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
            final GeneTree geneTree = geneTrees.get(j);
            assert sc.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either
            if (geneTree.computeCoalescentTimes()) {
                for (int i = 0; i < speciesTreeNodeCount; i++) { // for each species tree node/branch "i"
                    final List<Double> timesView = geneTree.coalescentTimes.get(i);
                    final int geneBranchEventCount = timesView.size();
                    final Double[] geneBranchCoalescentTimes = new Double[geneBranchEventCount];
                    timesView.toArray(geneBranchCoalescentTimes);
                    Arrays.sort(geneBranchCoalescentTimes);

                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts.count(i);

                    final Double[] coalescentTimesIJ = new Double[geneBranchEventCount + 2];
                    coalescentTimesIJ[0] = speciesEndTimes[i];
                    for (int k = 0; k < geneBranchEventCount; k++) {
                        coalescentTimesIJ[k + 1] = geneBranchCoalescentTimes[k];
                    }
                    coalescentTimesIJ[geneBranchEventCount + 1] = speciesStartTimes[i];

                    allLineageCounts.get(i)[j] = geneBranchLineageCount;
                    allEventCounts.get(i)[j] = geneBranchEventCount;
                    allCoalescentTimes.get(i).add(coalescentTimesIJ);
                }
            } else { // this gene tree IS NOT compatible with the species tree
                logP = Double.NEGATIVE_INFINITY;
                return logP;
            }
        }

        logP = 0.0;
        for (int i = 0; i < speciesTreeNodeCount; i++) {
            final Node speciesTreeNode = speciesTree.getNode(i); 
            final List<Double[]> branchCoalescentTimes = allCoalescentTimes.get(i);
            final int[] branchLineageCounts = allLineageCounts.get(i);
            final int[] branchEventCounts = allEventCounts.get(i);
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
