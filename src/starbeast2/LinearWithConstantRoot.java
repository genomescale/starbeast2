package starbeast2;


import java.text.DecimalFormat;
import java.util.List;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class LinearWithConstantRoot extends LinearPopulation {
    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double[] perGenePloidy, double[][] branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
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
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
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

        buf.append("dmv={");
        if (df == null) {
            buf.append(branchTopPopSize);
            buf.append(",");
            buf.append(branchTipPopSize);
        } else {
            buf.append(df.format(branchTopPopSize));
            buf.append(",");
            buf.append(df.format(branchTipPopSize));
        }
        buf.append("}");    }
}
