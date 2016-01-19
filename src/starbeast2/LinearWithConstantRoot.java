package starbeast2;


import java.text.DecimalFormat;
import java.util.List;

import beast.core.parameter.RealParameter;

/**
* @author Huw Ogilvie
 */

public class LinearWithConstantRoot extends LinearPopulationSize {
    @Override
    public double branchLogP(int speciesNetworkNodeNumber, NetworkNode speciesNetworkNode, double[] perGenePloidy,
                             List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final double branchTopPopSize = topPopSizes.getValue(speciesNetworkNodeNumber);

        double logP;
        if (speciesNetworkNode.isRoot()) {
            logP = constantLogP(branchTopPopSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
        } else {
            double branchTipPopSize;
            if (speciesNetworkNode.isLeaf()) {
                branchTipPopSize = tipPopSizes.getValue(speciesNetworkNodeNumber);
            } else {
                final int leftChildNodeNumber = speciesNetworkNode.getLeftChild().getNr();
                final int rightChildNodeNumber = speciesNetworkNode.getRightChild().getNr();
                branchTipPopSize = topPopSizes.getValue(leftChildNodeNumber) + topPopSizes.getValue(rightChildNodeNumber);
            }
            logP = linearLogP(branchTopPopSize, branchTipPopSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
        }

        return logP;
    }

    @Override
    public void serialize(NetworkNode speciesNetworkNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final int speciesNetworkNodeNumber = speciesNetworkNode.getNr();

        final double branchTopPopSize = topPopSizes.getValue(speciesNetworkNodeNumber);

        double branchTipPopSize;
        
        if (speciesNetworkNode.isRoot()) {
            branchTipPopSize = branchTopPopSize;
        } else if (speciesNetworkNode.isLeaf()) {
            branchTipPopSize = tipPopSizes.getValue(speciesNetworkNodeNumber);
        } else {
            final int leftChildNodeNumber = speciesNetworkNode.getLeftChild().getNr();
            final int rightChildNodeNumber = speciesNetworkNode.getRightChild().getNr();
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
