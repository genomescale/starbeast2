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

public class LinearWithConstantRoot extends MultispeciesPopulationModel {
    public Input<RealParameter> topPopSizesInput = new Input<RealParameter>("topPopSizes", "Population sizes at the top (rootward) end of each branch.", Validate.REQUIRED);
    public Input<RealParameter> tipPopSizesInput = new Input<RealParameter>("tipPopSizes", "Population sizes at the tips of leaf branches.", Validate.REQUIRED);

    private int nBranches;
    private int nSpecies;

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double[] perGenePloidy, List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final double branchTopPopSize = topPopSizes.getValue(speciesTreeNodeNumber);

        double logP;
        if (speciesTreeNode.isRoot()) {
            logP = constantLogP(branchTopPopSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
        } else {
            double branchTipPopSize;
            if (speciesTreeNode.isLeaf()) {
                branchTipPopSize = tipPopSizes.getValue(speciesTreeNodeNumber);
            } else {
                final int leftChildNodeNumber = speciesTreeNode.getLeft().getNr();
                final int rightChildNodeNumber = speciesTreeNode.getRight().getNr();
                branchTipPopSize = topPopSizes.getValue(leftChildNodeNumber) + topPopSizes.getValue(rightChildNodeNumber);
            }
            logP = linearLogP(branchTopPopSize, branchTipPopSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
        }

        return logP;
    }

    @Override
    public List<StateNode> initializePopSizes(TreeInterface speciesTree, double popInitial) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        nBranches = speciesTree.getNodeCount();
        nSpecies = speciesTree.getLeafNodeCount();
        final List<StateNode> popSizeVectors = new ArrayList<StateNode>();

        topPopSizes.setDimension(nBranches);
        tipPopSizes.setDimension(nSpecies);

        popSizeVectors.add(topPopSizes);
        popSizeVectors.add(tipPopSizes);

        for (int i = 0; i < topPopSizes.getDimension(); i++) {
            topPopSizes.setValue(i, popInitial / 2.0);
        }

        for (int i = 0; i < tipPopSizes.getDimension(); i++) {
            tipPopSizes.setValue(i, popInitial);
        }

        return popSizeVectors;
    }

    @Override
    public String serialize(Node speciesTreeNode) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();

        final double branchTopPopSize = topPopSizes.getValue(speciesTreeNodeNumber);

        double branchTipPopSize;
        
        if (speciesTreeNode.isRoot()) {
            branchTipPopSize = branchTopPopSize;
        } else if (speciesTreeNode.isLeaf()) {
            branchTipPopSize = tipPopSizes.getValue(speciesTreeNodeNumber);
        } else {
            final int leftChildNodeNumber = speciesTreeNode.getLeft().getNr();
            final int rightChildNodeNumber = speciesTreeNode.getRight().getNr();
            branchTipPopSize = topPopSizes.getValue(leftChildNodeNumber) + topPopSizes.getValue(rightChildNodeNumber);
        }

        final String dmv = String.format("dmv={%f,%f}", branchTopPopSize, branchTipPopSize);

        return dmv;
    }
}
