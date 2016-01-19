package starbeast2;


import java.text.DecimalFormat;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulationSize extends PopulationSizeModel {
    public Input<RealParameter> popSizesInput = new Input<RealParameter>("popSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double branchLogP(int speciesNetworkNodeNumber, NetworkNode speciesNetworkNode, double[] perGenePloidy,
                             List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter popSizes = popSizesInput.get();
        final double popSize = popSizes.getValue(speciesNetworkNodeNumber);
        double logP = constantLogP(popSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);

        return logP;
    }

    @Override
    public void initPopSizes(int nBranches) {
        final RealParameter popSizes = popSizesInput.get();
        popSizes.setDimension(nBranches);
    }

    @Override
    public void initPopSizes(double popInitial) {
        final RealParameter popSizes = popSizesInput.get();

        for (int i = 0; i < popSizes.getDimension(); i++) {
            popSizes.setValue(i, popInitial);
        }
    }

    @Override
    public void serialize(NetworkNode speciesNetworkNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter popSizes = popSizesInput.get();
        final int speciesNetworkNodeNumber = speciesNetworkNode.getNr();
        final double branchPopSize = popSizes.getValue(speciesNetworkNodeNumber);

        buf.append("dmv={");
        if (df == null) {
            buf.append(branchPopSize);
        } else {
            buf.append(df.format(branchPopSize));
        }
        buf.append("}");
    }
}
