package starbeast2;


import java.text.DecimalFormat;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulation extends MultispeciesPopulationModel {
    public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double[] perGenePloidy, double[][] branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter popSizes = popSizesInput.get();
        final double popSize = popSizes.getValue(speciesTreeNodeNumber);
        double logP = constantLogP(popSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);

        return logP;
    }

    @Override
    public double branchLogP(int geneTreeNode, int speciesTreeNodeNumber, Node speciesTreeNode, double perGenePloidy, double[] branchCoalescentTimes, int branchLineageCounts, int branchEventCounts) {
        final RealParameter popSizes = popSizesInput.get();
        final double popSize = popSizes.getValue(speciesTreeNodeNumber);
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
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter popSizes = popSizesInput.get();
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();
        final double branchPopSize = popSizes.getValue(speciesTreeNodeNumber);

        buf.append("dmv={");
        if (df == null) {
            buf.append(branchPopSize);
        } else {
            buf.append(df.format(branchPopSize));
        }
        buf.append("}");
    }
}
