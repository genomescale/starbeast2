package starbeast2;


import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;

import java.text.DecimalFormat;

/**
* @author Huw Ogilvie
 */

// Used to help BEAUTi with Analytical population sizes
public class DummyModel extends CalculationNode implements PopulationModel {
    @Override
    public void initAndValidate() {}

    @Override
    public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy, double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
        throw new RuntimeException("Cannot calculate the log probability of a branch for a single gene tree");
    }

    @Override
    public void initPopSizes(double initialPopSizes) {}

    @Override
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {}

    @Override
    public boolean isDirtyBranch(Node speciesTreeNode) {
        return false;
    }

    @Override
    public PopulationModel getBaseModel() {
        return this;
    }
}
