package starbeast2;

import beast.evolution.tree.Node;

/**
* @author Huw Ogilvie
 */

// for use with analytical, not MCMC, integration of population sizes
public class DummyPopulation extends MultispeciesPopulationModel {
    @Override
    public void initAndValidate() {
    }

    @Override
    public double branchLogP(int geneTreeNumber, int speciesTreeNodeNumber, Node speciesTreeNode, double perGenePloidy,
            double[] branchCoalescentTimes, int branchLineageCounts, int branchEventCounts) {
        return 0;
    }
}
