package starbeast2;


import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;

import java.text.DecimalFormat;
import java.util.Arrays;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulations extends CalculationNode implements PopulationModel {
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.", Validate.REQUIRED);
    public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    private SpeciesTreeInterface speciesTree;

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
        speciesTree = speciesTreeInput.get();
        speciesNodeCount = speciesTree.getNodeCount();
        popSizesInput.get().setDimension(speciesNodeCount);
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
        final double lower = popSizes.getLower();
        final double upper = popSizes.getUpper();

        if (popSizes.isEstimatedInput.get() && popInitial > lower && popInitial < upper) {
	        for (int i = 0; i < popSizes.getDimension(); i++) {
	            popSizes.setValue(i, popInitial);
	        }
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
                speciesBranchStatus[nodeI] = popSizes.isDirty(nodeI);

            needsUpdate = false;
        }

        return speciesBranchStatus[speciesNode.getNr()];
    }

    protected static double constantLogP(double popSize, double ploidy, double[] geneTimes, int geneN, int geneK) {
        double partialGamma = 0.0;
        for (int i = 0; i < geneK; i++) {
            partialGamma += (geneTimes[i + 1] - geneTimes[i]) * (geneN - i) * (geneN - (i + 1.0)) / 2.0;
        }
        
        if (geneN - geneK > 1) {
            partialGamma += (geneTimes[geneK + 1] - geneTimes[geneK]) * (geneN - geneK) * (geneN - (geneK + 1.0)) / 2.0;
        }

        final double branchGamma = partialGamma / ploidy;
        final double branchLogR = -geneK * Math.log(ploidy);
        final double logP = branchLogR - (geneK * Math.log(popSize)) - (branchGamma / popSize);

        /* StringBuffer sb = new StringBuffer();
        sb.append(popSize + "/" + ploidy + "/" + geneN + "/" + geneK + ": ");
        for (int i = 0; i < geneTimes.length; i++) {
            sb.append(geneTimes[i] + ", ");
        }
        sb.append(logP);
        System.out.println(sb); */

        return logP;
    }

    @Override
    public PopulationModel getBaseModel() {
        return this;
    }
}
