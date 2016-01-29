package network;

import java.text.DecimalFormat;
import java.util.List;

import beast.core.parameter.RealParameter;

/**
* @author Huw Ogilvie
 */

public class LinearWithConstantRoot extends LinearPopulation {
    @Override
    public double branchLogP(int speciesNetworkPopNumber, NetworkNode speciesNetworkNode, double[] perGenePloidy,
                             List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();

        final double branchTopPopSize = topPopSizes.getValue(speciesNetworkPopNumber);

        final double logP;
        if (speciesNetworkNode.isRoot()) {
            logP = ConstantPopulation.constantLogP(branchTopPopSize, perGenePloidy,
                              branchCoalescentTimes, branchLineageCounts, branchEventCounts);
        } else {
            final double branchTipPopSize;
            if (speciesNetworkNode.isLeaf()) {
                branchTipPopSize = tipPopSizes.getValue(speciesNetworkPopNumber);
            } else if (speciesNetworkNode.isReticulation()) {
                final int childNodeNumber;
                if (speciesNetworkNode.getLeftChild() != null)
                    childNodeNumber = speciesNetworkNode.getLeftChild().getNr();
                else
                    childNodeNumber = speciesNetworkNode.getRightChild().getNr();
                branchTipPopSize = topPopSizes.getValue(childNodeNumber*2);
            } else {
                final int leftChildNodeNumber = speciesNetworkNode.getLeftChild().getNr();
                final int rightChildNodeNumber = speciesNetworkNode.getRightChild().getNr();
                branchTipPopSize = topPopSizes.getValue(leftChildNodeNumber*2) +
                                   topPopSizes.getValue(rightChildNodeNumber*2+1);
            }
            logP = linearLogP(branchTopPopSize, branchTipPopSize, perGenePloidy,
                              branchCoalescentTimes, branchLineageCounts, branchEventCounts);
        }

        return logP;
    }

    @Override
    public void serialize(NetworkNode speciesNetworkNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final int speciesNetworkNodeNumber = speciesNetworkNode.getNr();

        // TODO ???
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
