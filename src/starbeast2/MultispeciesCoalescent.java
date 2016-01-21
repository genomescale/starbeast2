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

/**
* @author Remco Bouckaert
* @author Joseph Heled
* @author Huw Ogilvie
 */

@Description("Calculates probability of gene trees conditioned on a species tree (the multi-species coalescent).")
public class MultispeciesCoalescent extends Distribution {
    public Input<SpeciesNetwork> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<List<GeneTree>> geneTreeInput =
            new Input<>("geneTrees", "Gene trees within the species network.", new ArrayList<>());
    public Input<PopulationSizeModel> populationModelInput =
            new Input<>("populationModel", "The species network population model.", Validate.REQUIRED);

    private int nGeneTrees;
    private int speciesNetworkNodeCount;
    private double[] perGenePloidy;

    final private List<int[]> allLineageCounts = new ArrayList<>();
    final private List<int[]> allEventCounts = new ArrayList<>();
    final private List<List<Double[]>> allCoalescentTimes = new ArrayList<>();

    final static SanityChecks sc = new SanityChecks();

    @Override
    public void initAndValidate() throws Exception {
        final List<GeneTree> geneTrees = geneTreeInput.get();
        nGeneTrees = geneTrees.size();

        perGenePloidy = new double[nGeneTrees];
        for (int i = 0; i < nGeneTrees; i++) {
            final GeneTree geneTreeI = geneTrees.get(i);
            perGenePloidy[i] = geneTreeI.ploidy;
        }

        final PopulationSizeModel populationModel = populationModelInput.get();
        final Network speciesNetwork = speciesNetworkInput.get().getNetwork();
        speciesNetworkNodeCount = speciesNetwork.getNodeCount();
        populationModel.initPopSizes(speciesNetworkNodeCount * 2);
    }

    @Override
    public double calculateLogP() throws Exception {
        final Network speciesNetwork = speciesNetworkInput.get().getNetwork();
        final PopulationSizeModel populationModel = populationModelInput.get();
        speciesNetworkNodeCount = speciesNetwork.getNodeCount();
        double[] speciesStartTimes = new double[2*speciesNetworkNodeCount]; // the earlier date (rootward end)
        double[] speciesEndTimes = new double[2*speciesNetworkNodeCount]; // the later date (tipward end)

        final List<GeneTree> geneTrees = geneTreeInput.get();
        // assert sc.checkNetworkSanity(speciesNetwork.getRoot()); // species network should not be insane

        // [2*i] for the left branch and [2i+1] for the right branch (if right branch exists)
        for (int i = 0; i < speciesNetworkNodeCount; i++) {
            final NetworkNode speciesNode = speciesNetwork.getNode(i);
            final NetworkNode leftParent = speciesNode.getLeftParent();
            final NetworkNode rightParent = speciesNode.getLeftParent();

            speciesEndTimes[i*2+1] = speciesEndTimes[i*2] = speciesNode.getHeight();
            if (speciesNode.isRoot()) {
                speciesStartTimes[i*2+1] = speciesStartTimes[i*2] = Double.POSITIVE_INFINITY;
            } else {
                speciesStartTimes[i*2+1] = speciesStartTimes[i*2] = leftParent.getHeight();
                if (rightParent != null)
                    speciesStartTimes[i*2+1] = rightParent.getHeight();
            }
        }

        allLineageCounts.clear();
        allEventCounts.clear();
        allCoalescentTimes.clear();
        for (int i = 0; i < 2 * speciesNetworkNodeCount; i++) {
            allLineageCounts.add(new int[nGeneTrees]);
            allEventCounts.add(new int[nGeneTrees]);
            allCoalescentTimes.add(new ArrayList<>());
        }

        // transpose gene-branch list of lists to branch-gene list of lists
        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
            final GeneTree geneTree = geneTrees.get(j);
            // assert sc.checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either

            if (geneTree.computeCoalescentTimes()) {
                for (int i = 0; i < 2 * speciesNetworkNodeCount; i++) { // for each species network branch "i"
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
        for (int i = 0; i < 2 * speciesNetworkNodeCount; i++) {
            final NetworkNode speciesNetworkNode = speciesNetwork.getNode(i/2);
            final List<Double[]> branchCoalescentTimes = allCoalescentTimes.get(i);
            final int[] branchLineageCounts = allLineageCounts.get(i);
            final int[] branchEventCounts = allEventCounts.get(i);

            // linearPopulation uses i and speciesTreeNode, which needs double check later???
            final double branchLogP = populationModel.branchLogP(i, speciesNetworkNode, perGenePloidy,
                                                branchCoalescentTimes, branchLineageCounts, branchEventCounts);
            logP += branchLogP;
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
