package speciesnetwork.simulator;

import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.Operator;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import speciesnetwork.*;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * @author Chi Zhang
 */

@Description("Simulate gene trees given a species network (multispecies coalescent).")
public class CoalescentSimulator extends Operator {
    public Input<GeneTreeInSpeciesNetwork> geneTreeWrapperInput =
            new Input<>("geneTreeWithin", "Gene trees within the species network.", Validate.REQUIRED);
    public Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonSuperset", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);
    public Input<RealParameter> popSizesInput =
            new Input<>("popSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    private Network speciesNetwork;
    private Tree geneTree;
    private IntegerParameter embedding;
    private RealParameter popSizes;
    private double ploidy;

    private Multimap<NetworkNode, Node> networkNodeGeneLineagesMap = HashMultimap.create();

    @Override
    public void initAndValidate() {
        speciesNetwork = geneTreeWrapperInput.get().speciesNetworkInput.get();
        geneTree = geneTreeWrapperInput.get().geneTreeInput.get();
        embedding = geneTreeWrapperInput.get().embeddingInput.get();
        ploidy = geneTreeWrapperInput.get().ploidyInput.get();
        popSizes = popSizesInput.get();
        // assert (popSizes.getDimension() == speciesNetwork.getBranchCount());

        networkNodeGeneLineagesMap.clear();

        // generate map of tip names to tip nodes
        final Map<String, NetworkNode> speciesNodeMap = new HashMap<>();
        for (NetworkNode leafNode: speciesNetwork.getLeafNodes()) {
            final String speciesName = leafNode.getID();
            speciesNodeMap.put(speciesName, leafNode);
        }
        final Map<String, Node> geneNodeMap = new HashMap<>();
        for (Node leafNode: geneTree.getExternalNodes()) {
            final String geneName = leafNode.getID();
            geneNodeMap.put(geneName, leafNode);
        }
        // generate multimap of species network tip node to gene tree tip nodes
        final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
        for (Taxon speciesTip: taxonSuperSet.taxonsetInput.get()) {
            final NetworkNode speciesNode = speciesNodeMap.get(speciesTip.getID());
            final TaxonSet speciesTaxonSet = (TaxonSet) speciesTip;
            for (Taxon geneTip: speciesTaxonSet.taxonsetInput.get()) {
                final Node geneNode = geneNodeMap.get(geneTip.getID());
                networkNodeGeneLineagesMap.put(speciesNode, geneNode);
            }
        }
    }

    @Override
    public double proposal() {
        // initialize embedding matrix to -1 (no traversal)
        final int traversalNodeCount = speciesNetwork.getNodeCount() - speciesNetwork.getLeafNodeCount();
        final int geneNodeCount = geneTree.getNodeCount();
        for (int i = 0; i < traversalNodeCount; i++) {
            for (int j = 0; j < geneNodeCount; j++) {
                embedding.setMatrixValue(i, j, -1);
            }
        }

        NetworkNode networkRoot = speciesNetwork.getRoot();
        // reset visited indicator
        networkRoot.recursiveResetVisited();

        // simulate the gene tree
        simulate(networkRoot);

        return 0.0;
    }

    // recursively simulate lineages coalescent in each population
    private void simulate(NetworkNode snNode) {
        if (snNode.isVisited())
            return;
        if (snNode.getLeftChild() != null) {
            simulate(snNode.getLeftChild());
        }
        if (snNode.getRightChild() != null) {
            simulate(snNode.getRightChild());
        }

        snNode.setVisited();  // set visited indicator

        final Collection<Node> lineagesAtBottom = networkNodeGeneLineagesMap.get(snNode);

        final NetworkNode leftParent = snNode.getLeftParent();
        final NetworkNode rightParent = snNode.getLeftParent();

        if (snNode.isReticulation()) {



        } else {
            final int speciesBranchNumber = snNode.getNr();
            final double popSize = popSizes.getArrayValue(speciesBranchNumber);
            final double bottomHeight = snNode.getHeight();
            final double topHeight;
            if (leftParent == null && rightParent == null)  // network root
                topHeight = Double.POSITIVE_INFINITY;
            else if (leftParent != null)
                topHeight = snNode.getLeftParent().getHeight();
            else
                topHeight = snNode.getRightParent().getHeight();

            List<Node> lineagesAtTop = simulateCoalescentEvents(lineagesAtBottom, bottomHeight, topHeight, ploidy*popSize);
            if (leftParent != null)
                networkNodeGeneLineagesMap.putAll(leftParent, lineagesAtTop);
            else if (rightParent != null)
                networkNodeGeneLineagesMap.putAll(rightParent, lineagesAtTop);
            else
                geneTree.setRoot(lineagesAtTop.get(0));
        }
    }

    private List<Node> simulateCoalescentEvents(Collection<Node> lineagesAtBottom, double bottomHeight, double topHeight, double pNu) {
        // start from the lineages at the bottom
        List<Node> currentLineages = new ArrayList<>(lineagesAtBottom);
        double currentHeight = bottomHeight;

        // then go up backward in time
        while (currentLineages.size() > 1 && currentHeight < topHeight) {
            // generate a coalescent waiting time
            final int nLineage = currentLineages.size();
            final double coalescentRate = nLineage * (nLineage - 1) / (2 * pNu);
            final double waitingTime = Randomizer.nextExponential(coalescentRate);
            currentHeight += waitingTime;

            if (currentHeight < topHeight) {
                // randomly pick two lineages to coalescence
                int rnd = Randomizer.nextInt(nLineage);
                final Node left = currentLineages.get(rnd);
                currentLineages.remove(left);
                rnd = Randomizer.nextInt(nLineage - 1);
                final Node right = currentLineages.get(rnd);
                currentLineages.remove(right);
                // create a new node as the parent of the two picked nodes
                final Node node = new Node();
                node.setLeft(left);
                node.setRight(right);
                left.setParent(node);
                right.setParent(node);
                node.setHeight(currentHeight);
                //node.setNr();
                currentLineages.add(node);
            }
        }

        return currentLineages;
    }
}
