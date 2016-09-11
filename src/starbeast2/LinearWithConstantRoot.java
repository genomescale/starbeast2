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
* @author Joseph Heled
 */

public class LinearWithConstantRoot extends CalculationNode implements PopulationModel {
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.");
    public Input<RealParameter> topPopSizesInput = new Input<>("topPopSizes", "Population sizes at the top (rootward) end of each branch.", Validate.REQUIRED);
    public Input<RealParameter> tipPopSizesInput = new Input<>("tipPopSizes", "Population sizes at the tips of leaf branches.", Validate.REQUIRED);

    private SpeciesTreeInterface speciesTree;

    private boolean needsUpdate;
    private boolean[] speciesBranchStatus;
    private int rootNodeNumber;
    private int leafNodeCount;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = true;
        return needsUpdate;
    }

    @Override
    public void initAndValidate() {
        speciesTree = speciesTreeInput.get();
        final int speciesNodeCount = speciesTree.getNodeCount();
        rootNodeNumber = speciesNodeCount - 1;
        leafNodeCount = speciesTree.getLeafNodeCount();
        topPopSizesInput.get().setDimension(rootNodeNumber);
        tipPopSizesInput.get().setDimension(leafNodeCount);
        speciesBranchStatus = new boolean[speciesNodeCount];
        needsUpdate = true;
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();

        double branchTipPopSize;
        if (speciesTreeNode.isLeaf()) {
            branchTipPopSize = tipPopSizes.getValue(speciesTreeNodeNumber);
        } else {
            final int leftChildNodeNumber = speciesTreeNode.getLeft().getNr();
            final int rightChildNodeNumber = speciesTreeNode.getRight().getNr();
            branchTipPopSize = topPopSizes.getValue(leftChildNodeNumber) + topPopSizes.getValue(rightChildNodeNumber);
        }

        if (speciesTreeNode.isRoot()) {
            return ConstantPopulation.constantLogP(branchTipPopSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
        } else {
            final double branchTopPopSize = topPopSizes.getValue(speciesTreeNodeNumber);
            return linearLogP(branchTopPopSize, branchTipPopSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
        }
    }

    @Override
    public void initPopSizes(double popInitial) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();

        for (int i = 0; i < topPopSizes.getDimension(); i++) {
            topPopSizes.setValue(i, popInitial / 2.0);
        }

        for (int i = 0; i < tipPopSizes.getDimension(); i++) {
            tipPopSizes.setValue(i, popInitial);
        }
    }

    @Override
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();

        final double branchTopPopSize = topPopSizes.getValue(speciesTreeNodeNumber);

        double branchTipPopSize;
        if (speciesTreeNode.isLeaf()) {
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
        buf.append("}");
    }

    @Override
    public boolean isDirtyBranch(Node speciesNode) {
        if (needsUpdate) {
            final RealParameter topPopSizes = topPopSizesInput.get();
            final RealParameter tipPopSizes = tipPopSizesInput.get();

            Arrays.fill(speciesBranchStatus, false);
            Node[] speciesNodes = speciesTree.getNodesAsArray();
            
            // non-root nodes (linear population sizes)
            for (int nodeI = 0; nodeI < rootNodeNumber; nodeI++) {
                speciesBranchStatus[nodeI] |= topPopSizes.isDirty(nodeI);

                if (nodeI < leafNodeCount) { // a leaf node
                    speciesBranchStatus[nodeI] |= tipPopSizes.isDirty(nodeI);
                } else {
                    final int leftChildNumber = speciesNodes[nodeI].getLeft().getNr();
                    final int rightChildNumber = speciesNodes[nodeI].getRight().getNr();
                    speciesBranchStatus[nodeI] |= topPopSizes.isDirty(leftChildNumber);
                    speciesBranchStatus[nodeI] |= topPopSizes.isDirty(rightChildNumber);
                }
            }

            // the root node (constant population size)
            final int leftChildNumber = speciesNodes[rootNodeNumber].getLeft().getNr();
            final int rightChildNumber = speciesNodes[rootNodeNumber].getRight().getNr();
            speciesBranchStatus[rootNodeNumber] |= topPopSizes.isDirty(leftChildNumber);
            speciesBranchStatus[rootNodeNumber] |= topPopSizes.isDirty(rightChildNumber);

            needsUpdate = false;
        }

        return speciesBranchStatus[speciesNode.getNr()];
    }

    static double linearLogP(double topPopSize, double tipPopSize, double ploidy, double[] fTimes, int nLineagesBottom, int k) {
        final double fPopSizeTop = topPopSize * ploidy;
        final double fPopSizeBottom = tipPopSize * ploidy;

        final double d5 = fPopSizeTop - fPopSizeBottom;
        final double fTime0 = fTimes[0];
        final double a = d5 / (fTimes[k + 1] - fTime0);
        final double b = fPopSizeBottom;

        double logP = 0.0;
        if (Math.abs(d5) < 1e-10) {
            // use approximation for small values to bypass numerical instability
            for (int i = 0; i <= k; i++) {
                final double fTimeip1 = fTimes[i + 1];
                final double fPopSize = a * (fTimeip1 - fTime0) + b;
                if( i < k ) {
                    logP += -Math.log(fPopSize);
                }
                // slope = 0, so population function is constant

                final int i1 = nLineagesBottom - i;
                logP -= (i1 * (i1 - 1.0) / 2.0) * (fTimeip1 - fTimes[i]) / fPopSize;
            }
        } else {
            final double vv = b - a * fTime0;
            for (int i = 0; i <= k; i++) {
                final double fPopSize = a * fTimes[i + 1] + vv;
                if( i < k ) {
                    logP += -Math.log(fPopSize);
                }
                final double f = fPopSize / (a * fTimes[i] + vv);

                final int i1 = nLineagesBottom - i;
                logP += -(i1 * (i1 - 1.0) / 2.0) * Math.log(f) / a;
            }
        }

        /* StringBuffer sb = new StringBuffer();
        sb.append(topPopSize + "/" + tipPopSize + "/" + ploidy + "/" + nLineagesBottom + "/" + k + ": ");
        for (int i = 0; i < fTimes.length; i++) {
            sb.append(fTimes[i] + ", ");
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
