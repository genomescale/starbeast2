package speciesnetwork;

import beast.core.Input;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;

public class StarBeastClock extends BranchRateModel.Base {
    public Input<GeneTreeInSpeciesNetwork> geneTreeInput =
            new Input<>("geneTree", "The gene tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    public Input<SpeciesNetworkRates> speciesTreeRatesInput =
            new Input<>("speciesTreeRates", "The per-branch rates for the species tree", Input.Validate.REQUIRED);

    private int geneNodeCount;
    private double[] branchRates;
    private double[] storedBranchRates;
    private boolean needsUpdate;

    @Override
    public void initAndValidate() {
        geneNodeCount = geneTreeInput.get().getTree().getNodeCount();
        branchRates = new double[geneNodeCount];
        storedBranchRates = new double[geneNodeCount];
        needsUpdate = true;
    }

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = geneTreeInput.isDirty() || speciesTreeRatesInput.isDirty() || meanRateInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
        System.arraycopy(branchRates, 0, storedBranchRates, 0, branchRates.length);
        super.store();
    }

    @Override
    public void restore() {
        double[] tmpRatesArray = branchRates;
        branchRates = storedBranchRates;
        storedBranchRates = tmpRatesArray;
        super.restore();
    }

    private void update() throws Exception {
        final double geneTreeRate = meanRateInput.get().getValue();
        final double[] speciesTreeRates = speciesTreeRatesInput.get().getRatesArray();
        final double[][] speciesOccupancy = geneTreeInput.get().getSpeciesOccupancy();

        final int speciesNodeCount = speciesTreeRates.length;
        for (int i = 0; i < geneNodeCount - 1; i++) {
            double weightedSum = 0.0;
            double branchLength = 0.0;
            for (int j = 0; j < speciesNodeCount; j++) {
                // System.out.println(String.format("%d, %d: %f, %f", i, j, speciesTreeRates[j], speciesOccupancy[i][j]));
                weightedSum += speciesTreeRates[j] * speciesOccupancy[i][j];
                branchLength += speciesOccupancy[i][j];
            }

            branchRates[i] = geneTreeRate * weightedSum / branchLength;
        }
        // set the rate for the root branch of this gene to equal the input mean rate
        branchRates[geneNodeCount - 1] = geneTreeRate;
        
        needsUpdate = false;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate) {
            try {
                update();
            } catch (Exception e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        return branchRates[node.getNr()];
    }
}

