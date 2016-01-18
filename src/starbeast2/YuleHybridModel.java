package starbeast2;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.Distribution;

/**
 * pure birth model for the species network
 */
@Description("Pure birth model (i.e. no deaths) with hybridization")  // speciation times
public class YuleHybridModel extends Distribution {
    public Input<Network> networkInput =
            new Input<>("network", "The species network.", Validate.REQUIRED);
    final public Input<RealParameter> diversificationInput =
            new Input<>("diversificationRate", "Net diversification rate, lambda - mu.", Validate.REQUIRED);
    final public Input<RealParameter> hybridizationInput =
            new Input<>("hybridizationRate", "Hybridization rate, nu.", Validate.REQUIRED);
    final public Input<RealParameter> rhoProbInput =
            new Input<>("rho", "Sampling prob. of extant species, rho.");
    // final public Input<RealParameter> originHeightInput =
    //      new Input<>("originHeight", "the height of the point of origin of the process");
    // public Input<Boolean> conditionalOnRootInput =
    //      new Input<>("conditionalOnRoot", "condition on the root (otherwise: on the time of origin)", true);

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        // make sure that all tips are at the same height,
        // otherwise this Yule Model is not appropriate
        final Network network = networkInput.get();
        List<NetworkNode> leafs = network.getLeafNodes();
        double height = leafs.get(0).getHeight();
        for (NetworkNode leaf : leafs) {
            if (Math.abs(leaf.getHeight() - height) > 1e-8) {
                System.err.println("WARNING: Yule Model cannot handle dated tips.");
                break;
            }
        }
    }

    @Override
    public double calculateLogP() throws Exception {
        final Network network = networkInput.get();
        final double lambda = diversificationInput.get().getValue();
        final double nu = hybridizationInput.get().getValue();
        // final double rho = rhoProbInput.get() == null ? 1.0 : rhoProbInput.get().getValue();

        // sort the internal nodes according to their heights
        final List<NetworkNode> nodes = network.getInternalNodes();
        Collections.sort(nodes, new heightComparator());

        double logP= 0;
        // calculate probability of the network
        for (int i = 0; i < nodes.size(); i++) {
            final double nodeHeight = nodes.get(i).getHeight();
            final double nextHeight;
            if (i == 0)  // the youngest internal node
                nextHeight = 0.0;
            else
                nextHeight = nodes.get(i-1).getHeight();

            final int nBranch = network.getBranchCount((nodeHeight + nextHeight) /2.0);
            logP += (nBranch * lambda + nu * nBranch * (nBranch -1) /2) * (nextHeight - nodeHeight);
            if (i > 0 && nodes.get(i-1).isReticulation())
                logP += Math.log(nu);
            else if (i > 0)
                logP += Math.log(lambda);
        }
        return logP;
    }

    private class heightComparator implements Comparator<NetworkNode> {
        @Override
        public int compare(NetworkNode n1, NetworkNode n2) {
            return Double.compare(n1.getHeight(), n2.getHeight());
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation() || diversificationInput.get().somethingIsDirty();
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) { }
}
