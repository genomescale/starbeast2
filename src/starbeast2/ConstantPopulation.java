package starbeast2;


import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulation extends MultispeciesPopulationModel {
    public Input<RealParameter> popSizesInput = new Input<RealParameter>("popSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double[] perGenePloidy, List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter popSizes = popSizesInput.get();
        final double popSize = popSizes.getValue(speciesTreeNodeNumber);
        double logP = constantLogP(popSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);

        return logP;
    }

    @Override
    public void initPopSizes(double[] popInitial) {
        final RealParameter popSizes = popSizesInput.get();
        final int nBranches = popInitial.length;

        popSizes.setDimension(nBranches);

        for (int i = 0; i < popSizes.getDimension(); i++) {
            popSizes.setValue(i, popInitial[i]);
        }
    }

    @Override
    public String serialize(Node speciesTreeNode) {
        final RealParameter popSizes = popSizesInput.get();
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();
        final double branchPopSize = popSizes.getValue(speciesTreeNodeNumber);
        final String dmv = String.format("dmv=%f", branchPopSize);

        return dmv;
    }
}
