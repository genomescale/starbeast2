package snetworktests;

import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import speciesnetwork.NetworkParser;
import speciesnetwork.GeneTreeInSpeciesNetwork;
import speciesnetwork.MultispeciesCoalescent;
import speciesnetwork.PopulationSizeModel;
import speciesnetwork.operators.RebuildEmbedding;

abstract class PopulationTestHelper {
    String newickSpeciesNetwork;
    List<String> newickGeneTrees = new ArrayList<>();
    List<IntegerParameter> geneTreeEmbeddings = new ArrayList<>();
    RealParameter gammaParameter;

    TaxonSet speciesSuperset;
    TreeParser speciesTree;
    NetworkParser speciesNetwork;
    List<TreeParser> geneTrees = new ArrayList<>();
    List<GeneTreeInSpeciesNetwork> geneTreeWrappers = new ArrayList<>();

    State state = null;
    MultispeciesCoalescent msc;
    int nSpecies;
    int nBranches;
    double popSize;
    double ploidy;
    double gamma;
    double expectedLogP;

    final double allowedError = 1e-6;

    abstract public TaxonSet generateSuperset();
    abstract public PopulationSizeModel generatePopulationModel();

    @Test
    public void testLogP() {
        speciesSuperset = generateSuperset();
        initializeSpeciesNetwork();
        initializeStateNodes();
        initializeGeneTrees(false);

        final PopulationSizeModel populationModel = generatePopulationModel();
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("speciesNetwork", speciesNetwork, "geneTreeWithin", geneTreeWrappers, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }

    private void initializeSpeciesNetwork() {
        speciesTree = new TreeParser();
        speciesTree.initByName("newick", newickSpeciesNetwork, "IsLabelledNewick", true, "adjustTipHeights", false);
        speciesNetwork = new NetworkParser();
        speciesNetwork.initByName("tree", speciesTree);
    }

    private void initializeStateNodes() {
        if (state == null) state = new State();
        assertEquals(newickGeneTrees.size(), geneTreeEmbeddings.size());
        for (int i = 0; i < newickGeneTrees.size(); i++) {
            IntegerParameter embedding = geneTreeEmbeddings.get(i);
            state.initByName("stateNode", embedding);
        }
        gammaParameter = new RealParameter();
        gammaParameter.initByName("value", String.valueOf(gamma));
        state.initByName("stateNode", gammaParameter);
        state.initialise();
    }

    private void initializeGeneTrees(boolean reembed) {
        for (int i = 0; i < newickGeneTrees.size(); i++) {
            final String geneTreeNewick = newickGeneTrees.get(i);
            TreeParser geneTree = new TreeParser();
            geneTree.initByName("newick", geneTreeNewick, "IsLabelledNewick", true);
            geneTrees.add(geneTree);
            IntegerParameter embedding = geneTreeEmbeddings.get(i);
            if (reembed) { // rebuild the embedding
                RebuildEmbedding rebuildOperator = new RebuildEmbedding();
                rebuildOperator.initByName("geneTree", geneTree, "speciesNetwork", speciesNetwork,
                                           "taxonSuperset", speciesSuperset, "embedding", embedding);
                assertTrue(rebuildOperator.initializeEmbedding());
            }
            GeneTreeInSpeciesNetwork geneTreeWrapper = new GeneTreeInSpeciesNetwork();
            geneTreeWrapper.initByName("geneTree", geneTree, "ploidy", ploidy, "speciesNetwork", speciesNetwork,
                                       "embedding", embedding);
            geneTreeWrappers.add(geneTreeWrapper);
        }
    }
}
