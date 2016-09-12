package starbeast2;


import java.text.DecimalFormat;
import java.util.Arrays;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
* @author Huw Ogilvie
* @author Joseph Heled
 */

public class LinearWithConstantRoot extends CalculationNode implements PopulationModel {
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.", Validate.REQUIRED);
    public Input<RealParameter> lwcrPopSizesInput = new Input<>("populationSizes", "Population sizes at the tips of leaf branches.", Validate.REQUIRED);

    private SpeciesTreeInterface speciesTree;

    private boolean needsUpdate;
    private boolean[] speciesBranchStatus;
    private int rootNodeNumber;
    private int leafNodeCount;

    // scale top population sizes by half, so that tip sizes
    // have the same expectation as top sizes for a given prior
    public static final double TOP_SCALE_FACTOR = 0.5;

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
        lwcrPopSizesInput.get().setDimension(leafNodeCount + rootNodeNumber);
        speciesBranchStatus = new boolean[speciesNodeCount];
        needsUpdate = true;
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        final RealParameter lwcrPopSizes = lwcrPopSizesInput.get();

        double branchTipPopSize;
        if (speciesTreeNode.isLeaf()) {
            branchTipPopSize = lwcrPopSizes.getValue(speciesTreeNodeNumber);
        } else {
            final int leftChildTopI = leafNodeCount + speciesTreeNode.getLeft().getNr();
            final int rightChildTopI = leafNodeCount + speciesTreeNode.getRight().getNr();
            branchTipPopSize = (lwcrPopSizes.getValue(leftChildTopI) + lwcrPopSizes.getValue(rightChildTopI)) * TOP_SCALE_FACTOR;
        }

        if (speciesTreeNode.isRoot()) {
            return ConstantPopulations.constantLogP(branchTipPopSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
        } else {
            final int speciesTopI = leafNodeCount + speciesTreeNodeNumber;
            final double branchTopPopSize = lwcrPopSizes.getValue(speciesTopI) * TOP_SCALE_FACTOR;
            return linearLogP(branchTopPopSize, branchTipPopSize, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
        }
    }

    @Override
    public void initPopSizes(double popInitial) {
        final RealParameter lwcrPopSizes = lwcrPopSizesInput.get();

        for (int i = 0; i < lwcrPopSizes.getDimension(); i++)
            lwcrPopSizes.setValue(i, popInitial);
    }

    @Override
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter lwcrPopSizes = lwcrPopSizesInput.get();
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();

        double branchTipPopSize;
        if (speciesTreeNode.isLeaf()) {
            branchTipPopSize = lwcrPopSizes.getValue(speciesTreeNodeNumber);
        } else {
            final int leftChildTopI = leafNodeCount + speciesTreeNode.getLeft().getNr();
            final int rightChildTopI = leafNodeCount + speciesTreeNode.getRight().getNr();
            branchTipPopSize = (lwcrPopSizes.getValue(leftChildTopI) + lwcrPopSizes.getValue(rightChildTopI)) * TOP_SCALE_FACTOR;
        }

        if (speciesTreeNode.isRoot()) {
            buf.append("dmv=");
            if (df == null) buf.append(branchTipPopSize);
            else buf.append(df.format(branchTipPopSize));
            buf.append("}");
        } else {
            final int speciesTopI = leafNodeCount + speciesTreeNodeNumber;
            final double branchTopPopSize = lwcrPopSizes.getValue(speciesTopI) * TOP_SCALE_FACTOR;
            buf.append("dmv={");
            if (df == null) buf.append(branchTopPopSize + "," + branchTipPopSize);
            else buf.append(df.format(branchTopPopSize) + "," + df.format(branchTipPopSize));
            buf.append("}");
        }
    }

    @Override
    public boolean isDirtyBranch(Node speciesNode) {
        if (needsUpdate) {
            final RealParameter lwcrPopSizes = lwcrPopSizesInput.get();

            Arrays.fill(speciesBranchStatus, false);
            Node[] speciesNodes = speciesTree.getNodesAsArray();

            // non-root nodes (linear population sizes)
            for (int nodeI = 0; nodeI < speciesNodes.length; nodeI++) {
                final int speciesTopI = leafNodeCount + nodeI;
                // if the "top" population is dirty, no need to check the tip
                if (nodeI < rootNodeNumber && lwcrPopSizes.isDirty(speciesTopI)) { // not the root node
                    speciesBranchStatus[nodeI] = true;
                    continue;
                }

                if (nodeI < leafNodeCount) { // is a leaf node
                    speciesBranchStatus[nodeI] = lwcrPopSizes.isDirty(nodeI);
                } else { // is an internal node
                    final Node leftChild = speciesNodes[nodeI].getLeft();
                    final Node rightChild = speciesNodes[nodeI].getRight();
                    final int leftChildTopI = leafNodeCount + leftChild.getNr();
                    final int rightChildTopI = leafNodeCount + rightChild.getNr();
                    speciesBranchStatus[nodeI] = lwcrPopSizes.isDirty(leftChildTopI) ||
                            lwcrPopSizes.isDirty(rightChildTopI) ||
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
