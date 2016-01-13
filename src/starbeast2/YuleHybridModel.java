package starbeast2;


import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.speciation.*;


/**
 * pure birth model for the species network
 */
@Description("Pure birth model (i.e. no deaths) with hybridization")  // speciation times
public class YuleHybridModel extends SpeciesTreeDistribution {
    public Input<RealParameter> originInput =
            new Input<RealParameter>("origin", "the time when the process started");
    public Input<RealParameter> diversificationInput =
            new Input<RealParameter>("diversificationRate", "net diversification rate, lambda - mu", Validate.REQUIRED);
    public Input<RealParameter> hybridizationInput =
            new Input<RealParameter>("hybridizationRate", "hybridization rate, nu");
    public Input<RealParameter> rhoProbInput =
            new Input<RealParameter>("rho", "sampling prob. of extant species, rho");

    public Input<Boolean> conditionalOnRootInput =
            new Input<Boolean>("conditionalOnRoot", "condition on the root (default false: on the time of origin)", false);

    protected boolean conditionalOnRoot;
    protected boolean conditionalOnOrigin;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        conditionalOnRoot = conditionalOnRootInput.get();
        conditionalOnOrigin = originInput.get() != null;

        if (conditionalOnRoot && conditionalOnOrigin) {
            throw new RuntimeException("ERROR: Cannot condition on both root and origin.");
        }

        // make sure that all tips are at the same height,
        // otherwise this Yule Model is not appropriate
        TreeInterface tree = treeInput.get();
        if (tree == null) {
            tree = treeIntervalsInput.get().treeInput.get();
        }
        List<Node> leafs = tree.getExternalNodes();
        double height = leafs.get(0).getHeight();
        for (Node leaf : leafs) {
            if (Math.abs(leaf.getHeight() - height) > 1e-8) {
                System.err.println("WARNING: Yule Model cannot handle dated tips. Use for example a coalescent prior instead.");
                break;
            }
        }
    }

    @Override
    public double calculateTreeLogLikelihood(final TreeInterface tree) {
        if (conditionalOnOrigin && tree.getRoot().getHeight() > originInput.get().getValue())
            return Double.NEGATIVE_INFINITY;

        final double rho = rhoProbInput.get() == null ? 1.0 : rhoProbInput.get().getValue();
        final double lambda = diversificationInput.get().getValue();
        final double nu = hybridizationInput.get().getValue();

        double logP= 0;

        // get nodes sorted by their ages (from extant to origin?), each node has its age/height
        // each node should also have a label indicating either speciation or hybridization
        final Node[] nodes = tree.getNodesAsArray();

        // need another array to hold number of species in each interval?


        // calculate probability of the species network


        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {  // ?????
        return super.requiresRecalculation()
                || diversificationInput.get().somethingIsDirty()
                || (conditionalOnOrigin && originInput.get().somethingIsDirty());
    }

    @Override
    public boolean canHandleTipDates() {
        return false;
    }

}
