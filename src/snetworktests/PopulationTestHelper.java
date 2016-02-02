package snetworktests;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import speciesnetwork.NetworkParser;
import speciesnetwork.GeneTreeInSpeciesNetwork;
import speciesnetwork.MultispeciesCoalescent;
import speciesnetwork.PopulationSizeModel;
import speciesnetwork.RebuildEmbedding;

abstract class PopulationTestHelper {
    State state = null;
    String newickSpeciesNetwork;
    List<String> newickGeneTrees = new ArrayList<>();

    TaxonSet speciesSuperset;
    TreeParser speciesTree;
    NetworkParser speciesNetwork;
    List<TreeParser> geneTrees = new ArrayList<>();
    List<GeneTreeInSpeciesNetwork> geneTreeWrappers = new ArrayList<>();
    List<IntegerParameter> geneTreeEmbedding = new ArrayList<>();

    MultispeciesCoalescent msc;

    double ploidy;
    double popSize;
    double expectedLogP;
    int nSpecies;

    final double allowedError = 10e-6;

    abstract public TaxonSet generateSuperset() throws Exception;
    abstract public PopulationSizeModel generatePopulationModel() throws Exception;

    @Test
    public void testLogP() throws Exception {
        speciesSuperset = generateSuperset();
        initializeSpeciesNetwork();
        initializeStateNodes();
        initializeGeneTrees();

        final int nBranches = (nSpecies * 2) - 1;
        final PopulationSizeModel populationModel = generatePopulationModel();
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("speciesNetwork", speciesNetwork, "geneTrees", geneTreeWrappers, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }

    public void initializeSpeciesNetwork() throws Exception {
        speciesTree = new TreeParser();
        speciesTree.initByName("newick", newickSpeciesNetwork, "IsLabelledNewick", true);
        speciesNetwork = new NetworkParser();
        speciesNetwork.initByName("tree", speciesTree);
    }

    public void initializeStateNodes() throws Exception {
        if (state == null) state = new State();
        for (int i = 0; i < newickGeneTrees.size(); i++) {
            IntegerParameter embedding = new IntegerParameter();
            geneTreeEmbedding.add(embedding);
            embedding.initByName("value", "1");
            state.initByName("stateNode", embedding);
        }
        state.initialise();
    }

    public void initializeGeneTrees() throws Exception {
        for (int i = 0; i < newickGeneTrees.size(); i++) {
            final String geneTreeNewick = newickGeneTrees.get(i);
            TreeParser geneTree = new TreeParser();
            geneTree.initByName("newick", geneTreeNewick, "IsLabelledNewick", true);
            geneTrees.add(geneTree);
            IntegerParameter embedding = geneTreeEmbedding.get(i);
            RebuildEmbedding rebuildOperator = new RebuildEmbedding();
            rebuildOperator.initByName("geneTree", geneTree, "speciesNetwork", speciesNetwork, "taxonSuperset", speciesSuperset, "embedding", embedding);
            assertEquals(rebuildOperator.proposal(), 0.0, allowedError);
            GeneTreeInSpeciesNetwork geneTreeWrapper = new GeneTreeInSpeciesNetwork();
            geneTreeWrapper.initByName("geneTree", geneTree, "ploidy", ploidy, "speciesNetwork", speciesNetwork, "embedding", embedding);
            geneTreeWrappers.add(geneTreeWrapper);
        }
    }
}
