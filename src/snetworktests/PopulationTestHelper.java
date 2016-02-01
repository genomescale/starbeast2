package snetworktests;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import speciesnetwork.NetworkParser;
import speciesnetwork.GeneTree;
import speciesnetwork.MultispeciesCoalescent;
import speciesnetwork.PopulationSizeModel;
import speciesnetwork.SpeciesNetwork;

abstract class PopulationTestHelper {
    String newickSpeciesNetwork;
    List<String> newickGeneTrees = new ArrayList<>();

    NetworkParser speciesNetwork;
    List<TreeParser> geneTrees = new ArrayList<>();
    
    SpeciesNetwork speciesNetworkWrapper;
    List<GeneTree> geneTreeWrappers = new ArrayList<>();

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
        TaxonSet speciesSuperset = generateSuperset();
        initializeSpeciesNetwork(speciesSuperset);
        initializeGeneTrees();

        final int nBranches = (nSpecies * 2) - 1;
        final PopulationSizeModel populationModel = generatePopulationModel();
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("speciesNetwork", speciesNetworkWrapper, "geneTrees", geneTreeWrappers, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }

    public void initializeSpeciesNetwork(TaxonSet speciesSuperset) throws Exception {
        speciesNetwork = new NetworkParser();
        speciesNetwork.initByName("newick", newickSpeciesNetwork, "IsLabelledNewick", true);
        speciesNetworkWrapper = new SpeciesNetwork();
        speciesNetworkWrapper.initByName("network", speciesNetwork, "taxonSuperSet", speciesSuperset);
    }

    public void initializeGeneTrees() throws Exception {
        for (String geneTreeNewick: newickGeneTrees) {
            TreeParser geneTree = new TreeParser();
            geneTree.initByName("newick", geneTreeNewick, "IsLabelledNewick", true);
            geneTrees.add(geneTree);

            GeneTree geneTreeWrapper = new GeneTree();
            geneTreeWrapper.initByName("tree", geneTree, "ploidy", ploidy, "speciesTree", speciesNetworkWrapper);
            geneTreeWrappers.add(geneTreeWrapper);
        }
    }
}
