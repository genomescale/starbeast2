package starbeast2;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import beast.core.Input;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public class StarBeastClock extends BranchRateModel.Base {
    public Input<GeneTree> geneTreeInput = new Input<>("geneTree", "The gene tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    public Input<SpeciesTreeRates> speciesTreeRatesInput = new Input<>("speciesTreeRates", "The per-branch rates for the species tree", Input.Validate.REQUIRED);

    private int geneNodeCount;
    private double[] branchRates;
    private double[] storedBranchRates;
    private boolean needsUpdate;

    @Override
    public void initAndValidate() throws Exception {
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

    private void update() {
        final double geneTreeRate = meanRateInput.get().getValue();
        final double[] speciesTreeRates = speciesTreeRatesInput.get().getRatesArray();
        final double[][] speciesOccupancy = geneTreeInput.get().getSpeciesOccupancy();

        final int speciesNodeCount = speciesTreeRates.length;
        for (int i = 0; i < geneNodeCount; i++) {
            double weightedRate = 0.0;
            double branchLength = 0.0;

            for (int j = 0; j < speciesNodeCount; j++) {
                weightedRate += speciesTreeRates[j] * speciesOccupancy[i][j];
                branchLength += speciesOccupancy[i][j];
            }

            branchRates[i] = geneTreeRate * weightedRate / branchLength;
        }
        
        needsUpdate = false;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate) {
            update();
        }

        return branchRates[node.getNr()];
    }

    // for testing purposes
    public double getRate(TreeInterface geneTree, String[] targetNames) {
        final Set<String> s = new HashSet<>(Arrays.asList(targetNames));

        for (Node n: geneTree.getNodesAsArray()) {
            final HashSet<String> comparison = new HashSet<>();
            if (n.isLeaf()) {
                final String leafName = n.getID();
                comparison.add(leafName);
            } else {
                for (Node l: n.getAllLeafNodes()) {
                    final String leafName = l.getID();
                    comparison.add(leafName);
                }
            }

            if (s.equals(comparison)) {
                return getRateForBranch(n);
            }
        }

        return Double.NEGATIVE_INFINITY;
    }
}

