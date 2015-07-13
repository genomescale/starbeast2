package starbeast2;


import java.util.List;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class LinearWithConstantRoot extends LinearPopulation {
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
