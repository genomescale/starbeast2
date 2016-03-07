package sb2tests;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import starbeast2.GeneTree;
import starbeast2.MultispeciesCoalescent;
import starbeast2.MultispeciesPopulationModel;
import starbeast2.SpeciesTree;

abstract class PopulationTestHelper {
    String newickSpeciesTree;
    List<String> newickGeneTrees = new ArrayList<>();

    TreeParser speciesTree;
    List<TreeParser> geneTrees = new ArrayList<>();
    
    SpeciesTree speciesTreeWrapper;
    List<GeneTree> geneTreeWrappers = new ArrayList<>();

    MultispeciesCoalescent msc;

    double ploidy;
    double popSize;
    double expectedLogP;
    int nSpecies;

    final double allowedError = 10e-6;

    abstract public TaxonSet generateSuperset() throws Exception;
    abstract public MultispeciesPopulationModel generatePopulationModel() throws Exception;

    @Test
    public void testLogP() throws Exception {
        TaxonSet speciesSuperset = generateSuperset();
        initializeSpeciesTree(speciesSuperset);
        initializeGeneTrees();

        final int nBranches = (nSpecies * 2) - 1;
        final MultispeciesPopulationModel populationModel = generatePopulationModel();
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("speciesTree", speciesTreeWrapper, "geneTree", geneTreeWrappers, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }

    public void initializeSpeciesTree(TaxonSet speciesSuperset) throws Exception {
        speciesTree = new TreeParser();
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
        speciesTreeWrapper = new SpeciesTree();
        speciesTreeWrapper.initByName("tree", speciesTree, "taxonSuperSet", speciesSuperset);
    }

    public void initializeGeneTrees() throws Exception {
        for (String geneTreeNewick: newickGeneTrees) {
            TreeParser geneTree = new TreeParser();
            geneTree.initByName("newick", geneTreeNewick, "IsLabelledNewick", true);
            geneTrees.add(geneTree);

            GeneTree geneTreeWrapper = new GeneTree();
            geneTreeWrapper.initByName("tree", geneTree, "ploidy", ploidy, "speciesTree", speciesTreeWrapper);
            geneTreeWrappers.add(geneTreeWrapper);
        }
    }
}