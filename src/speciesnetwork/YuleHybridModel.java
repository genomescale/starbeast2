package speciesnetwork;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.distribution.BetaDistributionImpl;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.Distribution;

/**
 * Pure birth model for the species network.
 * @author Chi Zhang
 */

@Description("Pure birth model (i.e. no deaths) with hybridization")  // speciation times
public class YuleHybridModel extends Distribution {
    public Input<Network> networkInput =
            new Input<>("network", "The species network.", Validate.REQUIRED);
    final public Input<RealParameter> diversificationInput =
            new Input<>("diversificationRate", "Speciation rate, lambda.", Validate.REQUIRED);
    final public Input<RealParameter> hybridizationInput =
            new Input<>("hybridizationRate", "Hybridization rate, nu.", Validate.REQUIRED);
    final public Input<RealParameter> rhoProbInput =
            new Input<>("rho", "Sampling prob. of extant species, rho.");
    final public Input<RealParameter> betaShapeInput =
            new Input<>("betaShape", "Shape of the symmetric beta prior distribution on gammas.", Validate.REQUIRED);
    // final public Input<RealParameter> originHeightInput =
    //      new Input<>("originHeight", "the height of the point of origin of the process");
    // public Input<Boolean> conditionalOnRootInput =
    //      new Input<>("conditionalOnRoot", "condition on the root (otherwise: on the time of origin)", true);

    final private static Comparator<NetworkNode> hc = new NodeHeightComparator();

    private BetaDistributionImpl gammaPrior;
    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // make sure that all tips are at the same height,
        // otherwise this Yule Model is not appropriate
        final Network network = networkInput.get();
        final double firstHeight = network.nodes[0].height;
        for (int i = 1; i < network.leafNodeCount; i++) {
            final double height = network.nodes[i].height;
            if (Math.abs(firstHeight - height) > 1e-8) {
                System.err.println("WARNING: Yule Model cannot handle dated tips.");
                break;
            }
        }

        final double betaShape = betaShapeInput.get().getValue();
        gammaPrior = new BetaDistributionImpl(betaShape, betaShape);
    }

    @Override
    public double calculateLogP() {
        final Network network = networkInput.get();
        final double lambda = diversificationInput.get().getValue();
        final double nu = hybridizationInput.get().getValue();
        // final double rho = rhoProbInput.get() == null ? 1.0 : rhoProbInput.get().getValue();

        // sort the internal nodes according to their heights
        List<NetworkNode> nodes = new ArrayList<>();
        for (NetworkNode n: network.getInternalNodes()) {
            nodes.add(n);
        }
        nodes.sort(hc);

        double logP= 0;
        // calculate probability of the network
        for (int i = 0; i < nodes.size(); i++) {
            final NetworkNode node = nodes.get(i);
            final double nodeHeight = node.getHeight();
            final double nextHeight;
            if (i == 0)  // the youngest internal node
                nextHeight = 0.0;  // the tip
            else
                nextHeight = nodes.get(i-1).getHeight();

            final int nBranch = network.getBranchCount((nodeHeight + nextHeight) /2.0);
            logP += (nBranch * lambda + nu * nBranch * (nBranch -1) /2) * (nextHeight - nodeHeight);
            if (i > 0 && nodes.get(i-1).isReticulation())
                logP += Math.log(nu);
            else if (i > 0)
                logP += Math.log(lambda);

            if (node.isReticulation()) logP += gammaPrior.logDensity(node.inheritProb);
            if (!(logP > -Double.MAX_VALUE && logP < Double.MAX_VALUE)) System.out.println("???");
        }
        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        //return super.requiresRecalculation() || diversificationInput.get().somethingIsDirty() || networkInput.get().isDirty();
        return true;
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
