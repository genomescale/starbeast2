package speciesnetwork.simulator;

import java.io.*;
import java.util.*;

import beast.app.seqgen.*;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.*;
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
public class CoalescentSimulator extends Runnable {
    final public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network for embedding the gene tree.", Validate.REQUIRED);
    final public Input<RealParameter> gammaInput =
            new Input<>("gamma", "Inheritance probabilities (traversing left backward in time).", Validate.REQUIRED);
    final public Input<RealParameter> popSizesInput =
            new Input<>("popSizes", "Constant per-branch population sizes.", Validate.REQUIRED);
    final public Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonSuperset", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);

    final public Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", new ArrayList<>());
    final public Input<List<IntegerParameter>> embeddingsInput =
            new Input<>("embedding", "Map of gene tree traversal within the species network.", new ArrayList<>());
    final public Input<RealParameter> ploidiesInput =
            new Input<>("ploidy", "Ploidy (copy number) for each gene (default is 2).");
    final public Input<State> startStateInput =
            new Input<>("state", "elements of the state space", Validate.REQUIRED);

    final public Input<List<SequenceSimulator>> seqSimulatorsInput =
            new Input<>("sequenceSimulator", "Sequence simulator.", new ArrayList<>());
    final public Input<String> outputFileNameInput =
            new Input<>("outputFileName", "If provided, write to this file rather than to standard out.");

    final public Input<Integer> iterationsInput = new Input<>("iterations","Number of iterations to simulate.");

    private Network speciesNetwork;
    private RealParameter gammaP;
    private RealParameter popSizes;
    private List<Tree> geneTrees;
    private List<IntegerParameter> embeddings;
    private RealParameter ploidies;

    private int nrOfGeneTrees;
    private int speciesLeafNodeCount;
    private Multimap<NetworkNode, Node> networkNodeGeneLineagesMap = HashMultimap.create();
    private int nodeIndex;  //node index number

    private List<SequenceSimulator> seqSimulators;
    private List<Alignment> alignments = new ArrayList<>();

    @Override
    public void initAndValidate() {
        speciesNetwork = speciesNetworkInput.get();
        gammaP = gammaInput.get();
        popSizes = popSizesInput.get();
        geneTrees = geneTreesInput.get();
        embeddings = embeddingsInput.get();
        ploidies = ploidiesInput.get();
        seqSimulators = seqSimulatorsInput.get();
        speciesLeafNodeCount = speciesNetwork.getLeafNodeCount();

        // sanity check
        if (geneTrees == null || embeddings == null || geneTrees.size() != embeddings.size())
            throw new RuntimeException("Check gene tree and embedding input!");
        nrOfGeneTrees = geneTrees.size();

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
        if (ploidies == null) ploidies = new RealParameter("2.0");
        ploidies.setDimension(nrOfGeneTrees);

        final int traversalNodeCount = speciesNetwork.getNodeCount() - speciesNetwork.getLeafNodeCount();
        // simulate each gene tree and alignment
        for (int ig = 0; ig < nrOfGeneTrees; ig++) {
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
            for (NetworkNode leafNode : speciesNetwork.getLeafNodes()) {
                final String speciesName = leafNode.getID();
                speciesNodeMap.put(speciesName, leafNode);
            }
            final Map<String, Node> geneNodeMap = new HashMap<>();
            for (Node leafNode : geneTree.getExternalNodes()) {
                final String geneName = leafNode.getID();
                geneNodeMap.put(geneName, leafNode);
            }
            // multimap of species network tip node to gene tree tip nodes
            final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
            for (Taxon speciesTip : taxonSuperSet.taxonsetInput.get()) {
                final NetworkNode speciesNode = speciesNodeMap.get(speciesTip.getID());
                final TaxonSet speciesTaxonSet = (TaxonSet) speciesTip;
                for (Taxon geneTip : speciesTaxonSet.taxonsetInput.get()) {
                    final Node geneNode = geneNodeMap.get(geneTip.getID());
                    networkNodeGeneLineagesMap.put(speciesNode, geneNode);
                }
            }

            // reset visited indicator
            NetworkNode networkRoot = speciesNetwork.getRoot();
            networkRoot.resetAllVisited();

            // simulate the gene tree
            nodeIndex = 0;
            simulateGeneTree(networkRoot, geneTree, embedding, ploidies.getValue(ig));

            // simulate alignment on the gene tree
            if (seqSimulators.size() > ig) {
                alignments.add(seqSimulators.get(ig).simulate());
            }
        }

        // output
        writeXMLOutput();
    }

    private void writeXMLOutput() throws IOException {
        String outputFileName = outputFileNameInput.get();
        PrintStream out;  // where to print
        if (outputFileName == null) {
            out = System.out;
        } else {
            String msg = "Writing";
            if (new File(outputFileName).exists())
                msg = "Warning: Overwriting";
            System.err.println(msg + " file " + outputFileName);
            out = new PrintStream(outputFileName);
        }

        // print header
        out.println("<?xml version='1.0' encoding='UTF-8'?>");
        out.println("<beast namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:" +
                "beast.core.util:beast.evolution.operators:beast.evolution.sitemodel:" +
                "beast.evolution.substitutionmodel:beast.evolution.likelihood\" version=\"2.0\">");
        // print sequence data
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("    <data id=\"gene" + (i+1) + "\" name=\"alignment\">");
            if (seqSimulators.size() > i) {  // have simulated alignments
                Alignment alignment = alignments.get(i);
                List<Sequence> sequences = alignment.sequenceInput.get();
                for (Sequence seq : sequences)
                    out.println("        <sequence taxon=\"" + seq.getTaxon() + "\" value=\"" + seq.getData() + "\"/>");
            } else {
                Tree geneTree = geneTrees.get(i);
                for (Node leaf : geneTree.getExternalNodes())
                    out.println("        <sequence taxon=\"" + leaf.getID() + "\" totalcount=\"4\" value=\"-\"/>");
            }
            out.println("    </data>");
        }
        out.println("");  // mappings
        out.println("    <map name=\"Uniform\">beast.math.distributions.Uniform</map>\n" +
                    "    <map name=\"Exponential\">beast.math.distributions.Exponential</map>\n" +
                    "    <map name=\"LogNormal\">beast.math.distributions.LogNormalDistributionModel</map>\n" +
                    "    <map name=\"Normal\">beast.math.distributions.Normal</map>\n" +
                    "    <map name=\"Beta\">beast.math.distributions.Beta</map>\n" +
                    "    <map name=\"Gamma\">beast.math.distributions.Gamma</map>\n" +
                    "    <map name=\"LaplaceDistribution\">beast.math.distributions.LaplaceDistribution</map>\n" +
                    "    <map name=\"InverseGamma\">beast.math.distributions.InverseGamma</map>\n" +
                    "    <map name=\"OneOnX\">beast.math.distributions.OneOnX</map>\n" +
                    "    <map name=\"prior\">beast.math.distributions.Prior</map>\n");
        // print initial species network
        out.println("    <init spec=\"beast.util.TreeParser\" id=\"newick:species\" IsLabelledNewick=\"true\" " +
                            "adjustTipHeights=\"false\" newick=\"" + speciesNetwork.toString() + "\"/>\n");
        out.println("    <run chainLength=\"1000000\" id=\"mcmc\" spec=\"MCMC\">");  // MCMC block
        out.println("        <state id=\"state\" storeEvery=\"1000\">");  // states
        // print state nodes
        StringBuilder buf = new StringBuilder();
        for (int k = 0; k < gammaP.getDimension(); k++) {
            buf.append(gammaP.getValue(k));  buf.append(" ");
        }
        out.println("            <parameter id=\"gammaP\" lower=\"0.0\" upper=\"1.0\" name=\"stateNode\">" + buf + "</parameter>");
        out.println("            <stateNode id=\"network:species\" spec=\"speciesnetwork.NetworkParser\" tree=\"@newick:species\">");
        out.println("                <taxonset id=\"taxonSuperset\" spec=\"TaxonSet\">");
        final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
        for (Taxon speciesTip : taxonSuperSet.taxonsetInput.get()) {
            out.println("                    <taxon id=\"" + speciesTip.getID() + "\" spec=\"TaxonSet\">");
            final TaxonSet speciesTaxonSet = (TaxonSet) speciesTip;
            for (Taxon geneTip : speciesTaxonSet.taxonsetInput.get())
                out.println("                        <taxon id=\"" + geneTip.getID() + "\" spec=\"Taxon\"/>");
            out.println("                    </taxon>");
        }
        out.println("                </taxonset>");
        out.println("            </stateNode>");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("            <tree id=\"tree:gene" + (i+1) + "\" name=\"stateNode\">");
            out.println("                <taxonset alignment=\"@gene" + (i+1) + "\" " +
                                                 "id=\"taxonset:gene" + (i+1) + "\" spec=\"TaxonSet\"/>");
            out.println("            </tree>");
            // print true embedding (doesn't make sense as gene node number may change)
            IntegerParameter embedding = embeddings.get(i);
            out.println("            <stateNode id=\"embedding:gene" + (i+1) + "\" spec=\"parameter.IntegerParameter\" " +
                                        "dimension=\"" + embedding.getDimension() + "\" minordimension=\"" +
                                        embedding.getMinorDimension1() + "\">" + (-1) + "</stateNode>");
        }
        out.println("        </state>\n");  // end of states
        // print initial/true gene trees
        for (int i = 0; i < nrOfGeneTrees; i++) {
            Tree geneTree = geneTrees.get(i);
            out.println("        <init spec=\"beast.util.TreeParser\" id=\"newick:gene" + (i+1) + "\" " +
                    "initial=\"@tree:gene" + (i+1) + "\" taxa=\"@gene" + (i+1) + "\" IsLabelledNewick=\"true\" " +
                    "newick=\"" + geneTree.getRoot().toNewick() + "\"/>");
        }
        // starbeast initializer
        out.println("        <init estimate=\"false\" id=\"initializer\" method=\"random\" " +
                                "spec=\"speciesnetwork.StarBeastInitializer\" speciesNetwork=\"@network:species\">");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
            out.println("            <rebuildEmbedding id=\"rebuildEmbedding:gene" + (i+1) + "\" taxonSuperset=\"@taxonSuperset\" " +
                    "spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" " +
                    "geneTree=\"@tree:gene" + (i+1) + "\" embedding=\"@embedding:gene" + (i+1) + "\" weight=\"0.0\"/>");
        }
        out.println("        </init>\n");
        // print posterior, prior, and likelihood stuff
        out.println("        <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">");
        out.println("            <distribution id=\"prior\" spec=\"util.CompoundDistribution\">");  // prior
        out.println("                <distribution id=\"coalescent\" spec=\"speciesnetwork.MultispeciesCoalescent\" " +
                                                    "speciesNetwork=\"@network:species\">");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("                    <geneTreeWithin id=\"geneTree:gene" + (i+1) + "\" ploidy=\"2.0\" " +
                    "spec=\"speciesnetwork.GeneTreeInSpeciesNetwork\" speciesNetwork=\"@network:species\" " +
                    "geneTree=\"@tree:gene" + (i+1) + "\" embedding=\"@embedding:gene" + (i+1) + "\" gamma=\"@gammaP\"/>");
        }
        buf = new StringBuilder();
        for (int k = 0; k < popSizes.getDimension(); k++) {
            buf.append(popSizes.getValue(k));  buf.append(" ");
        }
        out.println("                    <populationModel id=\"popModel\" popSizes=\"" + buf + "\" " +
                                                    "spec=\"speciesnetwork.ConstantPopulation\"/>");
        out.println("                    <!-- populationModel alpha=\"2.0\" beta=\"1.0\" id=\"popModel\" " +
                                                    "spec=\"speciesnetwork.ConstantPopulationIO\"/ -->");
        out.println("                </distribution>");
        out.println("                <prior id=\"gamma.prior\" name=\"distribution\" x=\"@gammaP\">");
        out.println("                    <Beta id=\"beta.distr\" name=\"distr\" alpha=\"1\" beta=\"1\"/>");
        out.println("                </prior>");
        out.println("            </distribution>");
        out.println("            <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">");  // likelihood
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("                <distribution data=\"@gene" + (i+1) + "\" id=\"likelihood:gene" + (i+1) + "\" " +
                                                 "tree=\"@tree:gene" + (i+1) + "\" spec=\"TreeLikelihood\">");
            out.println("                    <siteModel id=\"siteModel:gene" + (i+1) + "\" mutationRate=\"1.0\" " +
                                                    "proportionInvariant=\"0.0\" spec=\"SiteModel\">");
            out.println("                        <substModel id=\"jc:gene" + (i+1) + "\" spec=\"JukesCantor\"/>");
            out.println("                    </siteModel>");
            out.println("                    <branchRateModel clock.rate=\"1.0\" id=\"clock:gene" + (i+1) + "\" " +
                                                "spec=\"beast.evolution.branchratemodel.StrictClockModel\"/>");
            out.println("                </distribution>");
        }
        out.println("            </distribution>");
        out.println("        </distribution>\n");
        // print operators
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("        <operator id=\"scaleAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators." +
                    "JointReembedding\" rebuildEmbedding=\"@rebuildEmbedding:gene" + (i+1) + "\" weight=\"3.0\">");
            out.println("            <operator id=\"scale:gene" + (i+1) + "\" spec=\"ScaleOperator\" " +
                                        "scaleFactor=\"0.5\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"scaleRootAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators." +
                    "JointReembedding\" rebuildEmbedding=\"@rebuildEmbedding:gene" + (i+1) + "\" weight=\"3.0\">");
            out.println("            <operator id=\"scaleRoot:gene" + (i+1) + "\" spec=\"ScaleOperator\" " +
                      "rootOnly=\"true\" scaleFactor=\"0.5\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"uniformAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators." +
                    "JointReembedding\" rebuildEmbedding=\"@rebuildEmbedding:gene" + (i+1) + "\" weight=\"30.0\">");
            out.println("            <operator id=\"uniform:gene" + (i+1) + "\" spec=\"Uniform\" " +
                                                            "tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"subSlideAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators." +
                    "JointReembedding\" rebuildEmbedding=\"@rebuildEmbedding:gene" + (i+1) + "\" weight=\"15.0\">");
            out.println("            <operator id=\"subSlide:gene" + (i+1) + "\" spec=\"SubtreeSlide\" " +
                                                            "tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"narrowAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators." +
                    "JointReembedding\" rebuildEmbedding=\"@rebuildEmbedding:gene" + (i+1) + "\" weight=\"15.0\">");
            out.println("            <operator id=\"narrow:gene" + (i+1) + "\" spec=\"Exchange\" " +
                                                            "tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"wideAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators." +
                    "JointReembedding\" rebuildEmbedding=\"@rebuildEmbedding:gene" + (i+1) + "\" weight=\"3.0\">");
            out.println("            <operator id=\"wide:gene" + (i+1) + "\" spec=\"Exchange\" " +
                                         "isNarrow=\"false\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"WBAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators." +
                    "JointReembedding\" rebuildEmbedding=\"@rebuildEmbedding:gene" + (i+1) + "\" weight=\"3.0\">");
            out.println("            <operator id=\"WilsonBalding:gene" + (i+1) + "\" spec=\"WilsonBalding\" " +
                                                            "tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>\n");
        }
        out.println("        <operator id=\"gammaScaler\" spec=\"ScaleOperator\" parameter=\"@gammaP\" " +
                                "scaleFactor=\"0.25\" weight=\"10.0\"/>");
        out.println("        <operator id=\"nodeSlider\" spec=\"speciesnetwork.operators.NodeSlider\" " +
                "speciesNetwork=\"@network:species\" taxonSuperset=\"@taxonSuperset\" weight=\"50.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("            <geneTree idref=\"tree:gene" + (i + 1) + "\"/>");
            out.println("            <embedding idref=\"embedding:gene" + (i + 1) + "\"/>");
        }
        out.println("        </operator>");

        // print loggers
        out.println("");
        out.println("        <logger id=\"screenlog\" logEvery=\"1000\" model=\"@posterior\">");
        out.println("            <log idref=\"posterior\"/>");
        out.println("            <log idref=\"likelihood\"/>");
        out.println("            <log idref=\"prior\"/>");
        out.println("            <log idref=\"coalescent\"/>");
        out.println("        </logger>");
        out.println("        <logger fileName=\"" + outputFileName + ".trace.log\" id=\"tracelog\" " +
                            "logEvery=\"200\" model=\"@posterior\" sort=\"smart\">");
        out.println("            <log idref=\"posterior\"/>");
        out.println("            <log idref=\"likelihood\"/>");
        out.println("            <log idref=\"prior\"/>");
        out.println("            <log idref=\"coalescent\"/>");
        out.println("            <log idref=\"gammaP\"/>");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("            <log id=\"TH:gene" + (i+1) + "\" tree=\"@tree:gene" + (i+1) + "\" " +
                                        "spec=\"beast.evolution.tree.TreeHeightLogger\"/>");
        }
        out.println("        </logger>");
        out.println("        <logger fileName=\"" + outputFileName + ".species.trees\" id=\"treelog:species\" " +
                             "logEvery=\"200\" mode=\"tree\">");
        out.println("            <log id=\"networkLogger:species\" spec=\"speciesnetwork.SpeciesNetworkLogger\" " +
                                    "speciesNetwork=\"@network:species\"/>");
        out.println("        </logger>");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("        <logger fileName=\"" + outputFileName + ".gene" + (i+1) + ".log\" " +
                                    "id=\"embedlog:gene" + (i+1) + "\" logEvery=\"200\" sort=\"smart\">");
            out.println("            <log idref=\"embedding:gene" + (i+1) + "\"/>");
            out.println("        </logger>");
            out.println("        <logger fileName=\"" + outputFileName + ".gene" + (i+1) + ".trees\" " +
                                    "id=\"treelog:gene" + (i+1) + "\" logEvery=\"200\" mode=\"tree\">");
            out.println("            <log id=\"treeLogger:gene" + (i+1) + "\" tree=\"@tree:gene" + (i+1) + "\" " +
                                        "spec=\"beast.evolution.tree.TreeWithMetaDataLogger\"/>");
            out.println("        </logger>");
        }
        out.println("    </run>");  // end of MCMC
        out.println("</beast>");
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
            final double leftP = snNode.getGamma();
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
            final double lPopSize = popSizes.getValue(leftBranchNumber);
            final double lTopHeight = snNode.getLeftParent().getHeight();
            List<Node> lineagesAtLeftTop =
                    simulateCoalescentEvents(lineagesAtLeft, bottomHeight, lTopHeight, ploidy*lPopSize, geneTree);
            final int rightBranchNumber = snNode.getRightBranchNumber();
            final double rPopSize = popSizes.getValue(rightBranchNumber);
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
            final double popSize = popSizes.getValue(speciesBranchNumber);
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
                left.setParent(node); right.setParent(node);
                node.setLeft(left);   node.setRight(right);
                node.setHeight(currentHeight);
                currentLineages.add(node);
            }
        }

        return currentLineages;
    }
}
