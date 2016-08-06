package speciesnetwork;

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
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<List<GeneTreeInSpeciesNetwork>> geneTreeWrapperInput =
            new Input<>("geneTreeWithin", "Gene trees within the species network.", new ArrayList<>());
    public Input<PopulationSizeModel> populationModelInput =
            new Input<>("populationModel", "The species network population model.", Validate.REQUIRED);

    private int nGeneTrees;
    private double[] perGenePloidy;

    final private List<int[]> allLineageCounts = new ArrayList<>();
    final private List<int[]> allEventCounts = new ArrayList<>();
    final private List<List<Double[]>> allCoalescentTimes = new ArrayList<>();

    final static SanityChecks sc = new SanityChecks();

    @Override
    public void initAndValidate() {
        final List<GeneTreeInSpeciesNetwork> geneTrees = geneTreeWrapperInput.get();
        nGeneTrees = geneTrees.size();

        perGenePloidy = new double[nGeneTrees];
        for (int i = 0; i < nGeneTrees; i++) {
            final GeneTreeInSpeciesNetwork geneTreeI = geneTrees.get(i);
            perGenePloidy[i] = geneTreeI.ploidy;
        }

        final Network speciesNetwork = speciesNetworkInput.get();
        final int speciesBranchCount = speciesNetwork.getBranchCount();
        final PopulationSizeModel populationModel = populationModelInput.get();
        populationModel.initPopSizes(speciesBranchCount);
    }

    private void buildTimes(NetworkNode node, int branchNumber, double parentHeight, double[] startTimes, double[] endTimes) {
        final double nodeHeight = node.height;
        endTimes[branchNumber] = nodeHeight;
        startTimes[branchNumber] = nodeHeight;
        for (int childBranchNr: node.childBranchNumbers) {
            final NetworkNode childNode = node.getChild(childBranchNr);
            buildTimes(childNode, childBranchNr, nodeHeight, startTimes, endTimes);
        }
    }

    @Override
    public double calculateLogP() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final NetworkNode speciesRoot = speciesNetwork.getRoot();
        assert sc.checkNetworkSanity(speciesRoot); // species network should not be insane
        final int speciesNodeCount = speciesNetwork.getNodeCount();
        final int rootBranchNr = (speciesNodeCount - 1) * 2;

        double[] speciesStartTimes = new double[rootBranchNr + 1]; // the earlier date (rootward end)
        double[] speciesEndTimes = new double[rootBranchNr + 1]; // the later date (tipward end)

        buildTimes(speciesRoot, rootBranchNr, Double.POSITIVE_INFINITY, speciesStartTimes, speciesEndTimes);

        allLineageCounts.clear();
        allEventCounts.clear();
        allCoalescentTimes.clear();
        for (int i = 0; i <= rootBranchNr; i++) {
            allLineageCounts.add(new int[nGeneTrees]);
            allEventCounts.add(new int[nGeneTrees]);
            allCoalescentTimes.add(new ArrayList<>());
        }

        final List<GeneTreeInSpeciesNetwork> geneTrees = geneTreeWrapperInput.get();
        // transpose gene-branch list of lists to branch-gene list of lists
        logP = 0.0;
        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
            final GeneTreeInSpeciesNetwork geneTree = geneTrees.get(j);
            geneTree.computeCoalescentTimes();
            logP += geneTree.logGammaSum;
            for (int i = 0; i <= rootBranchNr; i++) { // for each species network branch "i"
                final List<Double> timesView = geneTree.coalescentTimes.get(i);
                final int geneBranchEventCount = timesView.size();
                final Double[] geneBranchCoalescentTimes = new Double[geneBranchEventCount];
                timesView.toArray(geneBranchCoalescentTimes);
                Arrays.sort(geneBranchCoalescentTimes);

                final int geneBranchLineageCount = geneTree.coalescentLineageCounts.count(i);

                final Double[] coalescentTimesIJ = new Double[geneBranchEventCount + 2];
                coalescentTimesIJ[0] = speciesEndTimes[i];
                System.arraycopy(geneBranchCoalescentTimes, 0, coalescentTimesIJ, 1, geneBranchEventCount);
                coalescentTimesIJ[geneBranchEventCount + 1] = speciesStartTimes[i];

                allLineageCounts.get(i)[j] = geneBranchLineageCount;
                allEventCounts.get(i)[j] = geneBranchEventCount;
                allCoalescentTimes.get(i).add(coalescentTimesIJ);
            }
        }

        final PopulationSizeModel populationModel = populationModelInput.get();
        for (int i = 0; i <= rootBranchNr; i++) {
            final List<Double[]> branchCoalescentTimes = allCoalescentTimes.get(i);
            final int[] branchLineageCounts = allLineageCounts.get(i);
            final int[] branchEventCounts = allEventCounts.get(i);

            logP += populationModel.branchLogP(i, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
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
