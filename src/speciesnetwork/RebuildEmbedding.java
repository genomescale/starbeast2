package speciesnetwork;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * Rebuild the embedding of a gene tree in the species network.
 * @author Huw Ogilvie
 */

public class RebuildEmbedding extends Operator {
    public Input<Tree> geneTreeInput = new Input<>("geneTree", "The gene tree.", Validate.REQUIRED);
    public Input<Network> speciesNetworkInput = new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonSuperset", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);
    public Input<IntegerParameter> embeddingInput =
            new Input<>("embedding", "The matrix to embed the gene tree within the species network.", Validate.REQUIRED);

    final private Map<String, NetworkNode> tipMap = new HashMap<>();
    // heirs are the gene tree leaf tip numbers below each gene tree node or species network node
    final private Multimap<Node, Integer> geneNodeHeirs = HashMultimap.create();
    final private Multimap<NetworkNode, Integer> speciesNodeHeirs = HashMultimap.create();
    
    private int speciesLeafCount;
    private Tree geneTree;
    private Network speciesNetwork;
    private IntegerParameter embedding;

    @Override
    public void initAndValidate() {
        // generate map of species network tip names to species network tip nodes
        final Network speciesNetwork = speciesNetworkInput.get();
        final Map<String, NetworkNode> speciesNodeMap = new HashMap<>();

        final List<NetworkNode> speciesLeafNodes = speciesNetwork.getLeafNodes();
        speciesLeafCount = speciesLeafNodes.size();
        for (NetworkNode leafNode: speciesLeafNodes) {
            final String speciesName = leafNode.getID();
            speciesNodeMap.put(speciesName, leafNode);
        }

        // generate map of gene tree tip names to species network tip nodes
        final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
        final Set<Taxon> speciesSet = new HashSet<>(taxonSuperSet.taxonsetInput.get());

        for (Taxon species: speciesSet) {
            final String speciesName = species.getID();
            final NetworkNode speciesNode = speciesNodeMap.get(speciesName);
            final TaxonSet speciesTaxonSet = (TaxonSet) species;
            final Set<Taxon> tipSet = new HashSet<>(speciesTaxonSet.taxonsetInput.get());

            for (Taxon tip: tipSet) {
                final String tipName = tip.getID();
                tipMap.put(tipName, speciesNode);
            }
        }
    }

    /**
     * The proposal distribution depends only on the topology and times of the gene
     * tree and on the topology and times of the species network. Because this move
     * does not change topologies or times, the forward and reverse proposal
     * distributions are identical and the move is symmetric.
     */
    @Override
    public double proposal() {
        geneTree = geneTreeInput.get();
        speciesNetwork = speciesNetworkInput.get();
        embedding = embeddingInput.get();

        try {
            reinitializeEmbedding();
        } catch (Exception e) {
            e.printStackTrace();
        }

        for (final Node geneLeaf: geneTree.getExternalNodes()) {
            setLeafNodeHeirs(geneLeaf);
            recurseGeneHeirs(geneLeaf);
        }

        for (final NetworkNode speciesLeaf: speciesNetwork.getLeafNodes()) {
            recurseSpeciesHeirs(speciesLeaf);
        }

        // if no embedding of the gene trees is compatible with the species network
        if (!recurseRebuild(geneTree.getRoot(), speciesNetwork.getRoot())) {
            return Double.NEGATIVE_INFINITY;
        }

        // print matrix for debugging
        /* StringBuffer sb = new StringBuffer();
        for (int i = -1; i < embedding.getMinorDimension2(); i++) {
            for (int j = 0; j < embedding.getMinorDimension1(); j++) {
                if (i == -1) sb.append(j);
                else sb.append(embedding.getMatrixValue(i, j));
                sb.append(" ");
            }
            sb.append("\n");
        }
        System.out.println(sb); */

        return 0.0;
    }

    private void reinitializeEmbedding() throws Exception {
        // nothing traverses through species leaf nodes
        final int traversalNodeCount = speciesNetwork.getNodeCount() - speciesNetwork.getLeafNodeCount();
        final int geneNodeCount = geneTree.getNodeCount();

        embedding.setDimension(traversalNodeCount * geneNodeCount);
        embedding.setMinorDimension(geneNodeCount);

        // initialize embedding matrix to -1 (no traversal)
        for (int i = 0; i < traversalNodeCount; i++) {
            for (int j = 0; j < geneNodeCount; j++) {
                embedding.setMatrixValue(i, j, -1);
            }
        }
    }

    // the heir for each gene leaf node is itself
    // the heirs for each species leaf node is the associated gene leaf nodes
    private void setLeafNodeHeirs(final Node geneTreeNode) {
        final String tipName = geneTreeNode.getID();
        final NetworkNode speciesNetworkNode = tipMap.get(tipName);
        // System.out.println(String.format("%s - %d", speciesNetworkNode.getID(), geneTreeNode.getNr()));
        geneNodeHeirs.put(geneTreeNode, geneTreeNode.getNr());
        speciesNodeHeirs.put(speciesNetworkNode, geneTreeNode.getNr());
    }

    private void recurseGeneHeirs(final Node geneTreeNode) {
        if (!geneTreeNode.isLeaf()) {
            final Node leftChild = geneTreeNode.getLeft();
            final Node rightChild = geneTreeNode.getRight();
            geneNodeHeirs.putAll(geneTreeNode, geneNodeHeirs.get(leftChild));
            geneNodeHeirs.putAll(geneTreeNode, geneNodeHeirs.get(rightChild));
        }

        final Node parentNode = geneTreeNode.getParent();
        if (parentNode != null) recurseGeneHeirs(parentNode);
    }

    private void recurseSpeciesHeirs(final NetworkNode speciesNetworkNode) {
        if (!speciesNetworkNode.isLeaf()) {
            final NetworkNode leftChild = speciesNetworkNode.getLeftChild();
            final NetworkNode rightChild = speciesNetworkNode.getRightChild();
            if (leftChild != null) speciesNodeHeirs.putAll(speciesNetworkNode, speciesNodeHeirs.get(leftChild));
            if (rightChild != null) speciesNodeHeirs.putAll(speciesNetworkNode, speciesNodeHeirs.get(rightChild));
        }

        final NetworkNode leftParent = speciesNetworkNode.getLeftParent();
        final NetworkNode rightParent = speciesNetworkNode.getRightParent();

        if (leftParent != null) recurseSpeciesHeirs(leftParent);
        if (rightParent != null) recurseSpeciesHeirs(rightParent);
    }

    private boolean recurseRebuild(final Node geneTreeNode, final NetworkNode speciesNetworkNode) {
        // System.out.println(String.format("%s/%d:%d", speciesNetworkNode.getID(), speciesNetworkNode.getNr(), geneTreeNode.getNr()));
        final double geneTreeNodeHeight = geneTreeNode.getHeight();
        final double speciesNetworkNodeHeight = speciesNetworkNode.getHeight();
        // this coalescence node must be embedded in a descendant species network branch
        if (geneTreeNodeHeight < speciesNetworkNodeHeight) {
            final int traversalNodeNumber = speciesNetworkNode.getNr() - speciesLeafCount;
            final int geneTreeNodeNumber = geneTreeNode.getNr();

            final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);

            final NetworkNode leftSpecies = speciesNetworkNode.getLeftChild();
            final NetworkNode rightSpecies = speciesNetworkNode.getRightChild();

            /* StringBuffer sb = new StringBuffer();
            for (Integer i: requiredHeirs) {
                sb.append(i);
                sb.append(" ");
            }
            sb.append("- ");
            for (Integer i: speciesNodeHeirs.get(leftSpecies)) {
                sb.append(i);
                sb.append(" ");
            }
            sb.append("- ");
            for (Integer i: speciesNodeHeirs.get(rightSpecies)) {
                sb.append(i);
                sb.append(" ");
            }
            System.out.println(sb.toString()); */

            if (leftSpecies != null && speciesNodeHeirs.get(leftSpecies).containsAll(requiredHeirs)) {
                if (rightSpecies != null && speciesNodeHeirs.get(rightSpecies).containsAll(requiredHeirs)) {
                    // both species children are compatible with the gene tree, in which case the embedding becomes stochastic
                    if (Randomizer.nextBoolean()) {
                        embedding.setMatrixValue(traversalNodeNumber, geneTreeNodeNumber, 0);
                        return recurseRebuild(geneTreeNode, leftSpecies);
                    } else {
                        embedding.setMatrixValue(traversalNodeNumber, geneTreeNodeNumber, 1);
                        return recurseRebuild(geneTreeNode, rightSpecies);
                    }
                } else {
                    // only the left branch is compatible with the gene tree
                    embedding.setMatrixValue(traversalNodeNumber, geneTreeNodeNumber, 0);
                    return recurseRebuild(geneTreeNode, leftSpecies);
                }
            } else if (rightSpecies != null && speciesNodeHeirs.get(rightSpecies).containsAll(requiredHeirs)) {
                // only the right branch is compatible with the gene tree
                embedding.setMatrixValue(traversalNodeNumber, geneTreeNodeNumber, 1);
                return recurseRebuild(geneTreeNode, rightSpecies);
            } else {
                // neither child branch is compatible with the gene tree
                return false;
            }
        } else if (geneTreeNode.isLeaf()) {
            return true;
        } else {
            // embed both gene tree children
            return recurseRebuild(geneTreeNode.getLeft(), speciesNetworkNode) &&
                   recurseRebuild(geneTreeNode.getRight(), speciesNetworkNode);
        }
    }
}
