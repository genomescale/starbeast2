package starbeast2;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import beast.core.Input;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public class StarBeastClock extends BranchRateModel.Base {
    public Input<MultispeciesCoalescent> multispeciesCoalescentInput = new Input<MultispeciesCoalescent>("multispeciesCoalescent", "The multispecies coalescent calculation node.", Input.Validate.REQUIRED);
    public Input<SpeciesTreeRates> speciesTreeRatesInput = new Input<>("speciesTreeRates", "The per-branch rates for the species tree", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double getRateForBranch(Node node) {
        final double geneTreeRate = meanRateInput.get().getValue();
        if (node.isRoot()) {
            return geneTreeRate;
        }

        final Double[] speciesTreeRates = speciesTreeRatesInput.get().getRatesArray();
        final double[] speciesTreeOccupancy = multispeciesCoalescentInput.get().getOccupancy(node);

        final int nRates = speciesTreeRates.length;
        
        assert speciesTreeOccupancy != null;

        double relaxedRate = 0.0;
        double totalOccupancy = 0.0;
        for (int i = 0; i < nRates; i++) {
            relaxedRate += speciesTreeRates[i] * speciesTreeOccupancy[i];
            totalOccupancy += speciesTreeOccupancy[i];
        }

        final double geneTreeBranchRate = relaxedRate * geneTreeRate / totalOccupancy;

        return geneTreeBranchRate;
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

