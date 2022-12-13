package starbeast2;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;

import java.text.DecimalFormat;

/**
* @author Huw Ogilvie
 */

@Description("Pass calculations onwards to another model.")
public class PassthroughModel extends CalculationNode implements PopulationModel {
    public Input<PopulationModel> childModelInput = new Input<>("childModel", "Pass calculations onwards to another model", Validate.REQUIRED);

    PopulationModel childModel;

    @Override
    public void initAndValidate() {
        childModel = childModelInput.get();
    }

    // Calculate the truncated coalescent probability for a single species tree branch and gene
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        return childModel.branchLogP(speciesTreeNodeNumber, speciesTreeNode, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
    }

    // Sets the appropriate dimension size of each population size state node
    // To successfully resume from a saved state, this must be called via an initAndValidate method
    public void initPopSizes(final int nSpeciesBranches) {
        childModel.initPopSizes(nSpeciesBranches);
    }

    // Sets model-compatible default population sizes
    // To successfully begin a run, this must be called from a StateNodeInitializer
    public void initPopSizes(final double initialPopSizes) {
        childModel.initPopSizes(initialPopSizes);
    }

    // Per-branch population size information which will be added to a Newick string.
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        childModel.serialize(speciesTreeNode, buf, df);
    }

    public boolean isDirtyBranch(Node speciesTreeNode) {
        return childModel.isDirtyBranch(speciesTreeNode);
    }

    @Override
    public PopulationModel getBaseModel() {
        return childModel;
    }
}
