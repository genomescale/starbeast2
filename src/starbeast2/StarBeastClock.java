package starbeast2;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;

public class StarBeastClock extends CalculationNode implements BranchRateModel {
    public Input<MultispeciesCoalescent> multispeciesCoalescentInput = new Input<MultispeciesCoalescent>("multispeciesCoalescent", "The multispecies coalescent calculation node.", Input.Validate.REQUIRED);
    public Input<SpeciesTreeRates> speciesTreeRatesInput = new Input<>("speciesTreeRates", "The per-branch rates for the species tree", Input.Validate.REQUIRED);
    public Input<RealParameter> geneTreeRateInput = new Input<>("geneTreeRate", "Clock rate multiplier for this gene tree.", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double getRateForBranch(Node node) {
        final double geneTreeRate = geneTreeRateInput.get().getValue();
        if (node.isRoot()) {
            return geneTreeRate;
        }

        final double[] speciesTreeRates = speciesTreeRatesInput.get().getRatesArray();
        final double[] speciesTreeOccupancy = multispeciesCoalescentInput.get().getOccupancy(node);

        final int nRates = speciesTreeRates.length;

        double relaxedRate = 0.0;
        double totalOccupancy = 0.0;
        for (int i = 0; i < nRates; i++) {
            relaxedRate += speciesTreeRates[i] * speciesTreeOccupancy[i];
            totalOccupancy += speciesTreeOccupancy[i];
        }

        final double geneTreeBranchRate = relaxedRate * geneTreeRate / totalOccupancy;

        return geneTreeBranchRate;
    }
}

