package starbeast2;


import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;

import java.text.DecimalFormat;
import java.util.Arrays;

/**
* @author Huw Ogilvie
* @author Joseph Heled
 */

public class LinearWithConstantRoot extends CalculationNode implements PopulationModel {
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.", Validate.REQUIRED);
    public Input<RealParameter> tipPopSizesInput = new Input<>("tipPopulationSizes", "Population sizes at the tips of leaf branches.", Validate.REQUIRED);
    public Input<RealParameter> topPopSizesInput = new Input<>("topPopulationSizes", "Population sizes at the top of non-root branches.", Validate.REQUIRED);

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
        leafNodeCount = speciesTree.getLeafNodeCount(); // also the number of "tip" population sizes
        rootNodeNumber = speciesNodeCount - 1; // also the number of "top" population sizes
        tipPopSizesInput.get().setDimension(leafNodeCount);
        topPopSizesInput.get().setDimension(rootNodeNumber);
        speciesBranchStatus = new boolean[speciesNodeCount];
        needsUpdate = true;
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final RealParameter topPopSizes = topPopSizesInput.get();

        double branchTipPopSize;
        if (speciesTreeNode.isLeaf()) {
            branchTipPopSize = tipPopSizes.getValue(speciesTreeNodeNumber);
        } else {
            final int leftChildTopI = speciesTreeNode.getLeft().getNr();
            final int rightChildTopI = speciesTreeNode.getRight().getNr();
            branchTipPopSize = topPopSizes.getValue(leftChildTopI) + topPopSizes.getValue(rightChildTopI);
        }

        if (speciesTreeNode.isRoot()) {
            return ConstantPopulations.constantLogP(branchTipPopSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
        } else {
            final int speciesTopI = speciesTreeNodeNumber;
            final double branchTopPopSize = topPopSizes.getValue(speciesTopI);
            return linearLogP(branchTopPopSize, branchTipPopSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
        }
    }

    @Override
    public void initPopSizes(double tipInitial) {
        final double topInitial = tipInitial * 0.5;

        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final RealParameter topPopSizes = topPopSizesInput.get();

        final double tipLower = tipPopSizes.getLower();
        final double tipUpper = tipPopSizes.getUpper();

        final double topLower = topPopSizes.getLower();
        final double topUpper = topPopSizes.getUpper();

        if (tipPopSizes.isEstimatedInput.get() && tipInitial > tipLower && tipInitial < tipUpper) {
            for (int i = 0; i < tipPopSizes.getDimension(); i++)
                tipPopSizes.setValue(i, tipInitial);
        }

        if (topPopSizes.isEstimatedInput.get() && topInitial > topLower && topInitial < topUpper) {
	        for (int i = 0; i < topPopSizes.getDimension(); i++)
	            topPopSizes.setValue(i, topInitial);
        }
    }

    @Override
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final RealParameter topPopSizes = topPopSizesInput.get();
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();

        double branchTipPopSize;
        if (speciesTreeNode.isLeaf()) {
            branchTipPopSize = tipPopSizes.getValue(speciesTreeNodeNumber);
        } else {
            final int leftChildTopI = speciesTreeNode.getLeft().getNr();
            final int rightChildTopI = speciesTreeNode.getRight().getNr();
            branchTipPopSize = topPopSizes.getValue(leftChildTopI) + topPopSizes.getValue(rightChildTopI);
        }

        final double branchTopPopSize = (speciesTreeNode.isRoot()) ? branchTipPopSize : topPopSizes.getValue(speciesTreeNode.getNr());

        if (df == null) buf.append("dmv={" + branchTopPopSize + "," + branchTipPopSize + "}");
        else buf.append("dmv={" + df.format(branchTopPopSize) + "," + df.format(branchTipPopSize) + "}");
    }

    @Override
    public boolean isDirtyBranch(Node speciesNode) {
        if (needsUpdate) {
            final RealParameter tipPopSizes = tipPopSizesInput.get();
            final RealParameter topPopSizes = topPopSizesInput.get();

            Arrays.fill(speciesBranchStatus, false);
            Node[] speciesNodes = speciesTree.getNodesAsArray();

            // non-root nodes (linear population sizes)
            for (int nodeI = 0; nodeI < speciesNodes.length; nodeI++) {
                // if the "top" population is dirty, no need to check the tip
                if (nodeI < rootNodeNumber && topPopSizes.isDirty(nodeI)) { // not the root node
                    speciesBranchStatus[nodeI] = true;
                    continue;
                }

                if (nodeI < leafNodeCount) { // is a leaf node
                    speciesBranchStatus[nodeI] = tipPopSizes.isDirty(nodeI);
                } else { // is an internal node
                    final Node leftChild = speciesNodes[nodeI].getLeft();
                    final Node rightChild = speciesNodes[nodeI].getRight();
                    final int leftChildTopI = leftChild.getNr();
                    final int rightChildTopI = rightChild.getNr();
                    speciesBranchStatus[nodeI] = topPopSizes.isDirty(leftChildTopI) ||
                            topPopSizes.isDirty(rightChildTopI) ||
                            leftChild.isDirty() != Tree.IS_CLEAN ||
                            rightChild.isDirty() != Tree.IS_CLEAN;
                }
            }

            //Arrays.fill(speciesBranchStatus, true);
            needsUpdate = false;
        }

        return speciesBranchStatus[speciesNode.getNr()];
    }

    static double linearLogP(double topPopSize, double lwcrPopSize, double ploidy, double[] fTimes, int nLineagesBottom, int k) {
        final double fPopSizeTop = topPopSize * ploidy;
        final double fPopSizeBottom = lwcrPopSize * ploidy;

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
        sb.append(topPopSize + "/" + lwcrPopSize + "/" + ploidy + "/" + nLineagesBottom + "/" + k + ": ");
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
