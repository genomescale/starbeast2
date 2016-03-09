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
 * @author Chi Zhang
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
    
    private int speciesLeafCount, nChoices;
    private Tree geneTree;
    private Network speciesNetwork;
    private IntegerParameter embedding;

    @Override
    public void initAndValidate() {
        // generate map of species network tip names to species network tip nodes
        geneTree = geneTreeInput.get();
        speciesNetwork = speciesNetworkInput.get();
        embedding = embeddingInput.get();
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
     */
    @Override
    public double proposal() {
        // count the number of alternative traversing choices for the current state (n0)
        nChoices = 0;
        countNumberOfChoices(geneTree.getRoot(), speciesNetwork.getRoot());
        final int oldChoices = nChoices;

        try {
            reinitializeEmbedding();
        } catch (Exception e) {
            e.printStackTrace();
        }

        for (final Node geneLeaf : geneTree.getExternalNodes()) {
            setLeafNodeHeirs(geneLeaf);
            recurseGeneHeirs(geneLeaf);
        }
        for (final NetworkNode speciesLeaf : speciesNetwork.getLeafNodes()) {
            recurseSpeciesHeirs(speciesLeaf);
        }

        // rebuild the embedding, and
        // count the number of alternative traversing choices for the proposed state (n1)
        nChoices = 0;
        if (!recurseRebuild(geneTree.getRoot(), speciesNetwork.getRoot())) {
            // if no embedding of the gene trees is compatible with the species network
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

        // the proposal ratio is (2^n1)/(2^n0)
        return (nChoices - oldChoices) * Math.log(2);
    }

    public boolean initializeEmbedding() {
        reinitializeEmbedding();

        for (final Node geneLeaf: geneTree.getExternalNodes()) {
            setLeafNodeHeirs(geneLeaf);
            recurseGeneHeirs(geneLeaf);
        }
        for (final NetworkNode speciesLeaf: speciesNetwork.getLeafNodes()) {
            recurseSpeciesHeirs(speciesLeaf);
        }
        // rebuild the embedding
        return recurseRebuild(geneTree.getRoot(), speciesNetwork.getRoot());
    }

    private void reinitializeEmbedding() {
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
                    nChoices++;
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

    private boolean countNumberOfChoices(final Node geneTreeNode, final NetworkNode speciesNetworkNode) {
        final double geneTreeNodeHeight = geneTreeNode.getHeight();
        final double speciesNetworkNodeHeight = speciesNetworkNode.getHeight();
        // this coalescence node must be embedded in a descendant species network branch
        if (geneTreeNodeHeight < speciesNetworkNodeHeight) {
            final int traversalNodeNumber = speciesNetworkNode.getNr() - speciesLeafCount;
            final int geneTreeNodeNumber = geneTreeNode.getNr();
            final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
            final NetworkNode leftSpecies = speciesNetworkNode.getLeftChild();
            final NetworkNode rightSpecies = speciesNetworkNode.getRightChild();

            if (leftSpecies != null && speciesNodeHeirs.get(leftSpecies).containsAll(requiredHeirs)) {
                if (rightSpecies != null && speciesNodeHeirs.get(rightSpecies).containsAll(requiredHeirs)) {
                    // both species children are compatible with the gene tree
                    nChoices++;
                    if (embedding.getMatrixValue(traversalNodeNumber, geneTreeNodeNumber) == 0)
                        return countNumberOfChoices(geneTreeNode, leftSpecies);
                    else // embedding.getMatrixValue(traversalNodeNumber, geneTreeNodeNumber) == 1
                        return countNumberOfChoices(geneTreeNode, rightSpecies);
                } else {
                    // only the left branch is compatible with the gene tree
                    return countNumberOfChoices(geneTreeNode, leftSpecies);
                }
            } else if (rightSpecies != null && speciesNodeHeirs.get(rightSpecies).containsAll(requiredHeirs)) {
                // only the right branch is compatible with the gene tree
                return countNumberOfChoices(geneTreeNode, rightSpecies);
            } else {
                // neither child branch is compatible with the gene tree
                return false;  // for a valid embedding, should never go here
            }
        } else if (geneTreeNode.isLeaf()) {
            return true;
        } else {
            return countNumberOfChoices(geneTreeNode.getLeft(), speciesNetworkNode) &&
                   countNumberOfChoices(geneTreeNode.getRight(), speciesNetworkNode);
        }
    }
}
