package starbeast2;


import java.text.DecimalFormat;
import java.util.Arrays;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class UniformPopulations extends CalculationNode implements PopulationModel {
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.", Validate.REQUIRED);
    public Input<RealParameter> universalSizeInput = new Input<RealParameter>("universalSize", "Universal constant population size.", Validate.REQUIRED);

    private SpeciesTreeInterface speciesTree;

    private boolean needsUpdate;
    private int speciesNodeCount;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = true;
        return needsUpdate;
    }

    @Override
    public void initAndValidate() {
        speciesTree = speciesTreeInput.get();
        needsUpdate = true;
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        final double popSize = universalSizeInput.get().getValue();
        double logP = uniformLogP(popSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);

        return logP;
    }

    @Override
    public void initPopSizes(double popInitial) {
        final RealParameter universalSize = universalSizeInput.get();
        final double lower = universalSize.getLower();
        final double upper = universalSize.getUpper();

        if (universalSize.isEstimatedInput.get() && popInitial > lower && popInitial < upper) {
            universalSize.setValue(popInitial);
        }
    }

    @Override
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        final double popSize = universalSizeInput.get().getValue();

        buf.append("dmv={");
        if (df == null) {
            buf.append(popSize);
        } else {
            buf.append(df.format(popSize));
        }
        buf.append("}");
    }

    @Override
    public boolean isDirtyBranch(Node speciesNode) {
        final RealParameter universalSize = universalSizeInput.get();

        return universalSize.isDirty(0);
    }

    protected static double uniformLogP(double popSize, double ploidy, double[] geneTimes, int geneN, int geneK) {
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
