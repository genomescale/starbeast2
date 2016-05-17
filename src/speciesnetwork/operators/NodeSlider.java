package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

/**
 * @author Chi Zhang
 */

@Description("Randomly selects an internal network node and move its height using an uniform sliding window.")
public class NodeSlider extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    public Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "list of gene trees embedded in species network", new ArrayList<>());
    public Input<List<IntegerParameter>> embeddingsInput =
            new Input<>("embedding", "The matrices to embed the gene trees in the species network.", new ArrayList<>());
    public Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonSuperset", "Super-set of taxon sets mapping lineages to species.", Input.Validate.REQUIRED);
    public Input<Double> windowSizeInput =
            new Input<>("windowSize", "The size of the sliding window, default 0.1.", 0.1);

    // empty constructor to facilitate construction by XML + initAndValidate
    public NodeSlider() {
    }

    @Override
    public void initAndValidate() {
    }

    /**
     * Propose a new network-node height from a uniform distribution.
     * If the new value is outside the boundary, the excess is reflected back into the interval.
     * The proposal ratio is 1.0.
     */
    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();
        List<Tree> geneTrees = geneTreesInput.get();
        List<IntegerParameter> embeddings = embeddingsInput.get();
        final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
        final double windowSize = windowSizeInput.get();

        // StringBuffer sb = new StringBuffer();
        // sb.append(speciesNetwork.toString());
        // sb.append("\n");
        // check the embedding in the current species network
        for (int ig = 0; ig < geneTrees.size(); ig++) {
            IntegerParameter embedding = embeddings.get(ig);
            Tree geneTree = geneTrees.get(ig);

            RebuildEmbedding rebuildOperator = new RebuildEmbedding();
            rebuildOperator.initByName("speciesNetwork", speciesNetwork, "taxonSuperset", taxonSuperSet,
                    "geneTree", geneTree, "embedding", embedding);
            if(rebuildOperator.getNumberOfChoices() < 0)
                throw new RuntimeException("Developer ERROR: current embedding invalid! geneTree " + ig);
        }

        // pick an internal node randomly
        List<NetworkNode> intNodes = speciesNetwork.getInternalNodes();
        NetworkNode snNode = intNodes.get(Randomizer.nextInt(intNodes.size()));

        // determine the lower and upper bounds
        NetworkNode leftParent = snNode.getLeftParent();
        NetworkNode rightParent = snNode.getRightParent();
        double upper = Double.MAX_VALUE;
        if (leftParent != null)
            upper = leftParent.getHeight();
        if (rightParent != null && upper > rightParent.getHeight())
            upper = rightParent.getHeight();
        NetworkNode leftChild = snNode.getLeftChild();
        NetworkNode rightChild = snNode.getRightChild();
        double lower = Double.MIN_NORMAL;
        if (leftChild != null)
            lower = leftChild.getHeight();
        if (rightChild != null && lower < rightChild.getHeight())
            lower = rightChild.getHeight();
        if (lower >= upper)
            throw new RuntimeException("Developer ERROR: upper bound must be larger than lower bound!");

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
        snNode.setHeight(newHeight);

        // sb.append(speciesNetwork.toString());
        // sb.append("\n");
        // System.out.println(sb);
        // check the embedding in the new species network
        for (int ig = 0; ig < geneTrees.size(); ig++) {
            IntegerParameter embedding = embeddings.get(ig);
            Tree geneTree = geneTrees.get(ig);

            RebuildEmbedding rebuildOperator = new RebuildEmbedding();
            rebuildOperator.initByName("speciesNetwork", speciesNetwork, "taxonSuperset", taxonSuperSet,
                    "geneTree", geneTree, "embedding", embedding);
            if(rebuildOperator.getNumberOfChoices() < 0)
                return Double.NEGATIVE_INFINITY;  // not a valid embedding
        }

        return 0.0;
    }
}
