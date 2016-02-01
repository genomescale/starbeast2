package snetworktests;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import speciesnetwork.*;
import speciesnetwork.GeneTreeInSpeciesNetwork;

abstract class ExchangeTestHelper {
    String newickSpeciesTree;
    List<String> newickGeneTrees = new ArrayList<>();

    TreeParser speciesTree;
    List<TreeParser> geneTrees = new ArrayList<>();
    
    SpeciesTree speciesTreeWrapper;
    List<GeneTreeInSpeciesNetwork> geneTreeWrappers = new ArrayList<>();

    RealParameter popSizesParameter;
    PopulationSizeModel populationModel;
    MultispeciesCoalescent msc;

    double ploidy;
    double popSize;
    double expectedLogHR;
    String targetNodeLabel;
    boolean targetParent;

    final double allowedError = 10e-6;

    abstract public TaxonSet generateSuperset() throws Exception;

    @Test
    public void testLogHR() throws Exception {
        TaxonSet speciesSuperSet = generateSuperset();
        initializeSpeciesTree(speciesSuperSet);
        initializeGeneTrees();

        popSizesParameter = new RealParameter();
        popSizesParameter.initByName("value", String.valueOf(popSize));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", popSizesParameter);
        state.initialise();

        populationModel = new ConstantPopulation();
        populationModel.initByName("popSizes", popSizesParameter);

        msc = new MultispeciesCoalescent();
        msc.initByName("speciesTree", speciesTreeWrapper, "geneTree", geneTreeWrappers, "populationModel", populationModel);
        
        int nBranches = speciesTree.getNodeCount();
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        Node brother = null;
        for (Node n: speciesTree.getRoot().getAllLeafNodes()) {
            if (n.getID().equals(targetNodeLabel)) {
                if (targetParent) {
                    brother = n.getParent();
                } else {
                    brother = n;
                }
            }
        }

        CoordinatedExchange coex = new CoordinatedExchange();
        coex.initByName("tree", speciesTree, "speciesTree", speciesTreeWrapper, "geneTree", geneTrees);
        coex.manipulateSpeciesTree(brother);
        final double calculatedLogHR = coex.rearrangeGeneTrees();

        assertEquals(expectedLogHR, calculatedLogHR, allowedError);
    }

    public void initializeSpeciesTree(TaxonSet speciesSuperSet) throws Exception {
        speciesTree = new TreeParser();
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
        speciesTreeWrapper = new SpeciesTree();
        speciesTreeWrapper.initByName("tree", speciesTree, "taxonSuperSet", speciesSuperSet);
    }

    public void initializeGeneTrees() throws Exception {
        for (String geneTreeNewick: newickGeneTrees) {
            TreeParser geneTree = new TreeParser();
            geneTree.initByName("newick", geneTreeNewick, "IsLabelledNewick", true);
            geneTrees.add(geneTree);

            GeneTreeInSpeciesNetwork geneTreeWrapper = new GeneTreeInSpeciesNetwork();
            geneTreeWrapper.initByName("tree", geneTree, "ploidy", ploidy, "speciesTree", speciesTreeWrapper);
            geneTreeWrappers.add(geneTreeWrapper);
        }
    }
}
