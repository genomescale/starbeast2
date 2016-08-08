package speciesnetwork.operators;

import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;
import beast.core.Operator;
import beast.core.StateNode;

/**
 * @author Chi Zhang
 */

@Description("Combine a gene tree operator with RebuildEmbedding.")
public class NodeSlider extends Operator {
    public Input<List<RebuildEmbedding>> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene tree within species network.", new ArrayList<>());
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    public Input<Double> windowSizeInput =
            new Input<>("windowSize", "The size of the sliding window, default 0.01.", 0.01);

    int nGenes;

    @Override
    public void initAndValidate() {
        nGenes = rebuildEmbeddingInput.get().size();
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final List<RebuildEmbedding> reembedOps = rebuildEmbeddingInput.get();
        final double windowSize = windowSizeInput.get();

        // count the number of alternative traversing choices for the current state (n0)
        int oldChoices = 0;
        for (int i = 0; i < nGenes; i++) {
            final int nChoices = reembedOps.get(i).getNumberOfChoices();
            if (nChoices < 0)
                throw new RuntimeException("Developer ERROR: current embedding invalid! geneTree " + i);
            oldChoices += nChoices;
        }

        // pick an internal node randomly
        final int leafNodeCount = speciesNetwork.getLeafNodeCount();
        final int traversalNodeCount = speciesNetwork.getNodeCount() - leafNodeCount;
        final int randomNodeNumber = leafNodeCount + Randomizer.nextInt(traversalNodeCount);
        NetworkNode snNode = speciesNetwork.getNode(randomNodeNumber);

        // determine the lower and upper bounds
        double upper = Double.MAX_VALUE;
        for (NetworkNode p: snNode.getParents()) {
            upper = Math.min(upper, p.getHeight());
        }

        double lower = 0.0;
        for (NetworkNode c: snNode.getChildren()) {
            lower = Math.max(lower, c.getHeight());
        }

        // propose a new height, reflect it back if it's outside the boundary
        final double oldHeight = snNode.getHeight();
        double newHeight = oldHeight + (Randomizer.nextDouble() - 0.5) * windowSize;
        while (newHeight < lower || newHeight > upper) {
            if (newHeight < lower)
                newHeight = 2.0 * lower - newHeight;
            if (newHeight > upper)
                newHeight = 2.0 * upper - newHeight;
        }

        // update the new node height
        speciesNetwork.startEditing(this);
        snNode.setHeight(newHeight);
        SanityChecks.checkNetworkSanity(speciesNetwork.getRoot());

        // System.out.println(String.format("%d: %f -> %f [%f-%f]", randomNodeNumber, oldHeight, snNode.getHeight(), Math.min(100, upper), lower));
        // update the embedding in the new species network
        int newChoices = 0;
        for (int i = 0; i < nGenes; i++) {
            final RebuildEmbedding reembedOp = reembedOps.get(i);
            final int nChoices = reembedOp.initializeEmbedding();
            // System.out.println(String.format("Gene tree %d: %d choices", i, nChoices));
            if (!reembedOp.listStateNodes().isEmpty()) // copied from JointOperator
                reembedOp.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();
            if (nChoices < 0) return Double.NEGATIVE_INFINITY;
            newChoices += nChoices;
        }

        final double logHr = (newChoices - oldChoices) * Math.log(2);
        return logHr;
    }

    @Override
    public List<StateNode> listStateNodes() {
        final List<RebuildEmbedding> reembedOps = rebuildEmbeddingInput.get();
        List<StateNode> stateNodeList = new ArrayList<>();

        stateNodeList.addAll(super.listStateNodes());
        for (int i = 0; i < nGenes; i++) {
            stateNodeList.addAll(reembedOps.get(i).listStateNodes());
        }

        return stateNodeList;
    }
}
