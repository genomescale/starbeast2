package starbeast2;


import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

public class LinearPopulation extends MultispeciesPopulationModel {
    public Input<RealParameter> topPopSizesInput = new Input<RealParameter>("topPopSizes", "Population sizes at the top (rootward) end of each branch.", Validate.REQUIRED);
    public Input<RealParameter> tipPopSizesInput = new Input<RealParameter>("tipPopSizes", "Population sizes at the tips of leaf branches.", Validate.REQUIRED);

    private int nBranches;
    private int nSpecies;

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double[] perGenePloidy, List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final int nGenes = perGenePloidy.length;

        final double branchTopPopSize = topPopSizes.getValue(speciesTreeNodeNumber);

        double branchTipPopSize;
        if (speciesTreeNode.isLeaf()) {
            branchTipPopSize = tipPopSizes.getValue(speciesTreeNodeNumber);
        } else {
            final int leftChildNodeNumber = speciesTreeNode.getLeft().getNr();
            final int rightChildNodeNumber = speciesTreeNode.getRight().getNr();
            branchTipPopSize = topPopSizes.getValue(leftChildNodeNumber) + topPopSizes.getValue(rightChildNodeNumber);
        }

        // set the root branch heights for each gene tree (to equal the highest gene tree root node)
        double tallestGeneTreeHeight = 0.0;
        if (speciesTreeNode.isRoot()) {
            for (int j = 0; j < nGenes; j++) {
                tallestGeneTreeHeight = Math.max(tallestGeneTreeHeight, branchCoalescentTimes.get(j)[branchEventCounts[j]]);
            }
            for (int j = 0; j < nGenes; j++) {
                branchCoalescentTimes.get(j)[branchEventCounts[j] + 1] = tallestGeneTreeHeight;
            }
        }

        final double logP = linearLogP(branchTopPopSize, branchTipPopSize, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);

        /*// for debugging
        if (speciesTreeNode.isRoot()) {
            System.out.println(String.format("Tallest gene tree height = %f", tallestGeneTreeHeight));
            System.out.println(String.format("Root node %d logP = %f", speciesTreeNodeNumber, logP));
        } else if (speciesTreeNode.isLeaf()) {
            System.out.println(String.format("Leaf node %d logP = %f", speciesTreeNodeNumber, logP));
        } else { 
            System.out.println(String.format("Internal node %d logP = %f", speciesTreeNodeNumber, logP));
        }
        */

        return logP;
    }

    @Override
    public void initPopSizes(int nBranches) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final int nSpecies = (nBranches + 1) / 2;

        topPopSizes.setDimension(nBranches);
        tipPopSizes.setDimension(nSpecies);
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
    public String serialize(Node speciesTreeNode) {
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

        final String dmv = String.format("dmv={%f,%f}", branchTopPopSize, branchTipPopSize);

        return dmv;
    }
}
