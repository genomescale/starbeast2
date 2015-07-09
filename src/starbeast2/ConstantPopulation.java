package starbeast2;


import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulation extends MultispeciesPopulationModel {
    public Input<RealParameter> popSizesInput = new Input<RealParameter>("popSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    private RealParameter popSizes;

    @Override
    public void initAndValidate() throws Exception {
        popSizes = popSizesInput.get();
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, double[] perGenePloidy, List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final double popSize = popSizes.getValue(speciesTreeNodeNumber);
        double logP = constantLogP(popSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);

        return logP;
    }

    @Override
    public List<StateNode> initializePopSizes(TreeInterface speciesTree, double popInitial) {
        final int nBranches = speciesTree.getNodeCount();
        final List<StateNode> popSizeVectors = new ArrayList<StateNode>();

        popSizes.setDimension(nBranches);

        popSizeVectors.add(popSizes);

        for (int i = 0; i < popSizes.getDimension(); i++) {
            popSizes.setValue(i, popInitial);
        }

        return popSizeVectors;
    }

    @Override
    public String serialize(Node speciesTreeNode) {
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();
        final double branchPopSize = popSizes.getValue(speciesTreeNodeNumber);
        final String dmv = String.format("dmv=%f", branchPopSize);

        return dmv;
    }
}
