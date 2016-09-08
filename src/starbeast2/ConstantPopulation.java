package starbeast2;


import java.text.DecimalFormat;
import java.util.Arrays;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulation extends PopulationModel {
    public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    private boolean needsUpdate;
    private boolean[] speciesBranchStatus;
    private int speciesNodeCount;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = true;
        return needsUpdate;
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        speciesNodeCount = speciesTree.getNodeCount();
        speciesBranchStatus = new boolean[speciesNodeCount];
        needsUpdate = true;
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        final RealParameter popSizes = popSizesInput.get();
        final double popSize = popSizes.getValue(speciesTreeNodeNumber);
        double logP = constantLogP(popSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);

        return logP;
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
        if (needsUpdate) {
            final RealParameter popSizes = popSizesInput.get();

            Arrays.fill(speciesBranchStatus, false);

            for (int nodeI = 0; nodeI < speciesNodeCount; nodeI++)
                speciesBranchStatus[nodeI] |= popSizes.isDirty(nodeI);

            needsUpdate = false;
        }

        return speciesBranchStatus[speciesNode.getNr()];
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
