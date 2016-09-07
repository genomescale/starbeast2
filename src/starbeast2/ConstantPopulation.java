package starbeast2;


import java.text.DecimalFormat;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulation extends PopulationModel {
    public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        final RealParameter popSizes = popSizesInput.get();
        final double popSize = popSizes.getValue(speciesTreeNodeNumber);
        double logP = constantLogP(popSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);

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

    @Override
    public boolean isDirtyBranch(Node speciesNode) {
        final RealParameter popSizes = popSizesInput.get();
        return popSizes.isDirty(speciesNode.getNr());
    }

    protected static double constantLogP(double popSize, double genePloidy, double[] geneCoalescentTimes, int geneN, int geneK) {
        double partialGamma = 0.0;
        for (int i = 0; i < geneK; i++) {
            partialGamma += (geneCoalescentTimes[i + 1] - geneCoalescentTimes[i]) * (geneN - i) * (geneN - (i + 1.0)) / 2.0;
        }
        
        if (geneN - geneK > 1) {
            partialGamma += (geneCoalescentTimes[geneK + 1] - geneCoalescentTimes[geneK]) * (geneN - geneK) * (geneN - (geneK + 1.0)) / 2.0;
        }

        final double branchGamma = partialGamma / genePloidy;
        final double branchLogR = -geneK * Math.log(genePloidy);
        final double logP = branchLogR - (geneK * Math.log(popSize)) - (branchGamma / popSize);
        return logP;
    }
}
