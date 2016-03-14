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
    public Input<List<GeneTreeInSpeciesNetwork>> geneTreeInput =
            new Input<>("geneTree", "Gene trees within the species network.", new ArrayList<>());
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
    public void initAndValidate() {
        final List<GeneTreeInSpeciesNetwork> geneTrees = geneTreeInput.get();
        nGeneTrees = geneTrees.size();

        perGenePloidy = new double[nGeneTrees];
        for (int i = 0; i < nGeneTrees; i++) {
            final GeneTreeInSpeciesNetwork geneTreeI = geneTrees.get(i);
            perGenePloidy[i] = geneTreeI.ploidy;
        }

        final PopulationSizeModel populationModel = populationModelInput.get();
        final Network speciesNetwork = speciesNetworkInput.get();
        final int speciesBranchCount = speciesNetwork.getBranchCount();

        speciesNetworkNodeCount = speciesNetwork.getNodeCount();
        populationModel.initPopSizes(speciesBranchCount);
    }

    @Override
    public double calculateLogP() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final int speciesNodeCount = speciesNetwork.getNodeCount();
        final int reticulationNodeCount = speciesNetwork.getReticulationNodeCount();
        // each reticulation node has two branches
        final int speciesBranchCount = speciesNodeCount + reticulationNodeCount;

        final PopulationSizeModel populationModel = populationModelInput.get();
        speciesNetworkNodeCount = speciesNetwork.getNodeCount();
        double[] speciesStartTimes = new double[speciesBranchCount]; // the earlier date (rootward end)
        double[] speciesEndTimes = new double[speciesBranchCount]; // the later date (tipward end)

        final List<GeneTreeInSpeciesNetwork> geneTrees = geneTreeInput.get();
        // assert sc.checkNetworkSanity(speciesNetwork.getRoot()); // species network should not be insane

        for (int i = 0; i < speciesNetworkNodeCount; i++) {
            final NetworkNode speciesNetworkNode = speciesNetwork.getNode(i);
            final NetworkNode leftParent = speciesNetworkNode.getLeftParent();
            final NetworkNode rightParent = speciesNetworkNode.getRightParent();

            final int speciesBranchNumber = speciesNetworkNode.getLeftBranchNumber();
            speciesEndTimes[speciesBranchNumber] = speciesNetworkNode.getHeight();

            if (leftParent != null && rightParent != null) { // this is a reticulation node
                speciesEndTimes[speciesBranchNumber + 1] = speciesNetworkNode.getHeight();
                speciesStartTimes[speciesBranchNumber] = leftParent.getHeight();
                speciesStartTimes[speciesBranchNumber + 1] = rightParent.getHeight();
            } else if (leftParent != null) {
                speciesStartTimes[speciesBranchNumber] = leftParent.getHeight();
            } else if (rightParent != null) {
                speciesStartTimes[speciesBranchNumber] = rightParent.getHeight();
            } else { // this is the root node
                speciesStartTimes[speciesBranchNumber] = Double.POSITIVE_INFINITY;
            }
        }

        allLineageCounts.clear();
        allEventCounts.clear();
        allCoalescentTimes.clear();
        for (int i = 0; i < speciesBranchCount; i++) {
            allLineageCounts.add(new int[nGeneTrees]);
            allEventCounts.add(new int[nGeneTrees]);
            allCoalescentTimes.add(new ArrayList<>());
        }

        // transpose gene-branch list of lists to branch-gene list of lists
        logP = 0.0;
        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
            final GeneTreeInSpeciesNetwork geneTree = geneTrees.get(j);
            geneTree.computeCoalescentTimes();
            logP += geneTree.logGammaSum;
            for (int i = 0; i < speciesBranchCount; i++) { // for each species network branch "i"
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

        for (int i = 0; i < speciesBranchCount; i++) {
            final List<Double[]> branchCoalescentTimes = allCoalescentTimes.get(i);
            final int[] branchLineageCounts = allLineageCounts.get(i);
            final int[] branchEventCounts = allEventCounts.get(i);

            final double branchLogP = populationModel.branchLogP(i, perGenePloidy,
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
