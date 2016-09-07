package starbeast2;


import java.text.DecimalFormat;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

@Description("Calculates probability of coalescence events within a branch for a single gene based on a demographic function.")
public class PopulationModel extends CalculationNode {
    public Input<PopulationModel> childModelInput = new Input<>("childModel", "Pass calculations onwards to another model");
    public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.");

    PopulationModel childModel;
    SpeciesTree speciesTree;

    @Override
    public void initAndValidate() {
        childModel = childModelInput.get();
        speciesTree = speciesTreeInput.get();
    }

    // Calculate the truncated coalescent probability for a single species tree branch and gene
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        if (childModel != null) return childModel.branchLogP(speciesTreeNodeNumber, speciesTreeNode, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
        else return 0.0;
    }

    // Sets the appropriate dimension size of each population size state node
    // To successfully resume from a saved state, this must be called via an initAndValidate method
    public void initPopSizes(final int nSpeciesBranches) {
        if (childModel != null) childModel.initPopSizes(nSpeciesBranches);
    }

    // Sets model-compatible default population sizes
    // To successfully begin a run, this must be called from a StateNodeInitializer
    public void initPopSizes(final double initialPopSizes) {
        if (childModel != null) childModel.initPopSizes(initialPopSizes);
    }

    // Per-branch population size information which will be added to a Newick string.
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        if (childModel != null) childModel.serialize(speciesTreeNode, buf, df);
    }

    public boolean isDirtyBranch(Node speciesTreeNode) {
        if (childModel != null) return childModel.isDirtyBranch(speciesTreeNode);
        else return false;
    }
}
