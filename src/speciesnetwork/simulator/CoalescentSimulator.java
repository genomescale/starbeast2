package speciesnetwork.simulator;

import java.io.IOException;
import java.util.*;

import beast.app.seqgen.*;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.XMLProducer;
import speciesnetwork.*;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * @author Chi Zhang
 */

@Description("Simulate gene trees given a species network (multispecies coalescent).")
public class CoalescentSimulator extends Runnable {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network for embedding the gene tree.", Validate.REQUIRED);
    public Input<RealParameter> gammaInput =
            new Input<>("gamma", "Inheritance probabilities (traversing left backward in time).", Validate.REQUIRED);
    public Input<RealParameter> popSizesInput =
            new Input<>("popSizes", "Constant per-branch population sizes.", Validate.REQUIRED);
    public Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonSuperset", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);

    public Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", new ArrayList<>());
    public Input<List<IntegerParameter>> embeddingsInput =
            new Input<>("embedding", "Map of gene tree traversal within the species network.", new ArrayList<>());
    public Input<RealParameter> ploidiesInput =
            new Input<>("ploidy", "Ploidy (copy number) for each gene (default is 2).");
    public Input<State> startStateInput =
            new Input<>("state", "elements of the state space", Validate.REQUIRED);

    public Input<List<SequenceSimulator>> seqSimulatorsInput =
            new Input<>("sequenceSimulator", "Sequence simulator.", new ArrayList<>());

    private Network speciesNetwork;
    private RealParameter gammaP;
    private RealParameter popSizes;
    private List<Tree> geneTrees;
    private List<IntegerParameter> embeddings;
    private RealParameter ploidies;

    private Multimap<NetworkNode, Node> networkNodeGeneLineagesMap = HashMultimap.create();
    private int nodeIndex;  //node index number
    private int speciesLeafNodeCount;

    private List<SequenceSimulator> seqSimulators;

    @Override
    public void initAndValidate() {
        speciesNetwork = speciesNetworkInput.get();
        gammaP = gammaInput.get();
        popSizes = popSizesInput.get();
        geneTrees = geneTreesInput.get();
        embeddings = embeddingsInput.get();
        ploidies = ploidiesInput.get();
        seqSimulators = seqSimulatorsInput.get();

        // sanity check
        if (geneTrees == null || embeddings == null || geneTrees.size() != embeddings.size())
            throw new RuntimeException("Check gene tree and embedding input!");

        // initialize state nodes, essential
        State state = startStateInput.get();
        state.initialise();
    }

    @Override
    public void run() throws IOException {
        // check correctness of parameter dimensions
        final int speciesBranchCount = speciesNetwork.getBranchCount();
        final int speciesReticulationNodeCount = speciesNetwork.getReticulationNodeCount();
        if (popSizes.getDimension() != speciesBranchCount)
            popSizes.setDimension(speciesBranchCount);
        if (gammaP.getDimension() != speciesReticulationNodeCount)
            gammaP.setDimension(speciesReticulationNodeCount);

        final int traversalNodeCount = speciesNetwork.getNodeCount() - speciesNetwork.getLeafNodeCount();
        // simulate each gene tree and alignment
        for (int ig = 0; ig < geneTrees.size(); ig++) {
            Tree geneTree = geneTrees.get(ig);
            IntegerParameter embedding = embeddings.get(ig);

            // initialize embedding matrix to -1 (no traversal)
            final int geneNodeCount = geneTree.getNodeCount();
            embedding.setDimension(traversalNodeCount * geneNodeCount);
            embedding.setMinorDimension(geneNodeCount);
            for (int i = 0; i < traversalNodeCount; i++)
                for (int j = 0; j < geneNodeCount; j++)
                    embedding.setMatrixValue(i, j, -1);

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
            // multimap of species network tip node to gene tree tip nodes
            final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
            for (Taxon speciesTip: taxonSuperSet.taxonsetInput.get()) {
                final NetworkNode speciesNode = speciesNodeMap.get(speciesTip.getID());
                final TaxonSet speciesTaxonSet = (TaxonSet) speciesTip;
                for (Taxon geneTip: speciesTaxonSet.taxonsetInput.get()) {
                    final Node geneNode = geneNodeMap.get(geneTip.getID());
                    networkNodeGeneLineagesMap.put(speciesNode, geneNode);
                }
            }

            // reset visited indicator
            NetworkNode networkRoot = speciesNetwork.getRoot();
            networkRoot.resetAllVisited();

            // simulate the gene tree
            nodeIndex = 0;
            speciesLeafNodeCount = speciesNetwork.getLeafNodeCount();
            double ploidy = 2.0;
            if (ploidies != null && ploidies.getDimension() > ig)
                ploidy = ploidies.getArrayValue(ig);
            simulateGeneTree(networkRoot, geneTree, embedding, ploidy);

            // simulate alignment on the gene tree
            if (seqSimulators != null)
                simulateSequences(ig);
        }
    }

    // recursively simulate lineages coalescent in each population
    private void simulateGeneTree(NetworkNode snNode, Tree geneTree, IntegerParameter embedding, double ploidy) {
        if (snNode.isVisited())
            return;
        if (snNode.getLeftChild() != null) {
            simulateGeneTree(snNode.getLeftChild(), geneTree, embedding, ploidy);
        }
        if (snNode.getRightChild() != null) {
            simulateGeneTree(snNode.getRightChild(), geneTree, embedding, ploidy);
        }

        snNode.setVisited();  // set visited indicator

        final Collection<Node> lineagesAtBottom = networkNodeGeneLineagesMap.get(snNode);
        final NetworkNode leftParent = snNode.getLeftParent();
        final NetworkNode rightParent = snNode.getRightParent();

        if (snNode.isReticulation()) {
            // assign lineages at the bottom to the left and right populations
            final int reticulationNumber = snNode.getReticulationNumber();
            final double leftP = gammaP.getArrayValue(reticulationNumber);
            final Collection<Node> lineagesAtLeft = new HashSet<>();
            final Collection<Node> lineagesAtRight = new HashSet<>();
            for (Node lineage : lineagesAtBottom) {
                if (Randomizer.nextDouble() < leftP)
                    lineagesAtLeft.add(lineage);
                else
                    lineagesAtRight.add(lineage);
            }

            final double bottomHeight = snNode.getHeight();
            final int leftBranchNumber = snNode.getLeftBranchNumber();
            final double lPopSize = popSizes.getArrayValue(leftBranchNumber);
            final double lTopHeight = snNode.getLeftParent().getHeight();
            List<Node> lineagesAtLeftTop =
                    simulateCoalescentEvents(lineagesAtLeft, bottomHeight, lTopHeight, ploidy*lPopSize, geneTree);
            final int rightBranchNumber = snNode.getRightBranchNumber();
            final double rPopSize = popSizes.getArrayValue(rightBranchNumber);
            final double rTopHeight = snNode.getRightParent().getHeight();
            List<Node> lineagesAtRightTop =
                    simulateCoalescentEvents(lineagesAtRight, bottomHeight, rTopHeight, ploidy*rPopSize, geneTree);

            networkNodeGeneLineagesMap.putAll(leftParent, lineagesAtLeftTop);
            networkNodeGeneLineagesMap.putAll(rightParent, lineagesAtRightTop);
            // update embedding
            final int traversalLeftNumber = leftParent.getNr() - speciesLeafNodeCount;
            for (final Node geneNode : lineagesAtLeftTop)
                embedding.setMatrixValue(traversalLeftNumber, geneNode.getNr(), 0);
            final int traversalRightNumber = rightParent.getNr() - speciesLeafNodeCount;
            for (final Node geneNode : lineagesAtRightTop)
                embedding.setMatrixValue(traversalRightNumber, geneNode.getNr(), 1);
        }
        else {
            final int speciesBranchNumber = snNode.getBranchNumber(0);
            final double popSize = popSizes.getArrayValue(speciesBranchNumber);
            final double bottomHeight = snNode.getHeight();
            final double topHeight;
            if (leftParent == null && rightParent == null)  // network root
                topHeight = Double.POSITIVE_INFINITY;
            else if (leftParent != null)
                topHeight = snNode.getLeftParent().getHeight();
            else
                topHeight = snNode.getRightParent().getHeight();

            List<Node> lineagesAtTop =
                    simulateCoalescentEvents(lineagesAtBottom, bottomHeight, topHeight, ploidy*popSize, geneTree);
            if (leftParent != null) {
                networkNodeGeneLineagesMap.putAll(leftParent, lineagesAtTop);
                // update embedding
                final int traversalLeftNumber = leftParent.getNr() - speciesLeafNodeCount;
                for (final Node geneNode : lineagesAtTop)
                    embedding.setMatrixValue(traversalLeftNumber, geneNode.getNr(), 0);
            } else if (rightParent != null) {
                networkNodeGeneLineagesMap.putAll(rightParent, lineagesAtTop);
                // update embedding
                final int traversalRightNumber = rightParent.getNr() - speciesLeafNodeCount;
                for (final Node geneNode : lineagesAtTop)
                    embedding.setMatrixValue(traversalRightNumber, geneNode.getNr(), 1);
            } else {
                geneTree.setRoot(lineagesAtTop.get(0));  // bad idea? no
            }
        }
    }

    private List<Node> simulateCoalescentEvents(Collection<Node> lineages, double bottomHeight,
                                                double topHeight, double pNu, Tree geneTree) {
        // start from the lineages at the bottom
        List<Node> currentLineages = new ArrayList<>(lineages);
        double currentHeight = bottomHeight;

        List<Node> internalNodes = geneTree.getInternalNodes();

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
                // deal with the parent of the two picked nodes
                final Node node = internalNodes.get(nodeIndex++);
                node.setLeft(left);   node.setRight(right);
                left.setParent(node); right.setParent(node);
                node.setHeight(currentHeight);
                currentLineages.add(node);
            }
        }

        return currentLineages;
    }

    private void simulateSequences(int locus) {
        SequenceSimulator seqSimulator = seqSimulators.get(locus);

        Alignment alignment = seqSimulator.simulate();

        // print the alignment to screen TODO: print to a XML file
        System.out.println(new XMLProducer().toRawXML(alignment));
    }
}
