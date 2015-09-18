package sb2tests;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import starbeast2.ConstantPopulation;
import starbeast2.CoordinatedExchange;
import starbeast2.GeneTreeWithinSpeciesTree;
import starbeast2.MultispeciesCoalescent;
import starbeast2.MultispeciesPopulationModel;

public class SmallCoordinatedExchangeTest {
    private final String newickSpeciesTree = "((A:1.0,B:1.0):2.0,C:3.0)";
    private final String newickGeneTreeA = "(((A1:1.5,B1:1.5):1.0,(A2:2.0,B2:2.0):0.5):1.5,C1:4.0)";

    final int nSpecies = 7;
    final int nBranches = (nSpecies * 2) - 1;
    final double[] popSizes = new double[nBranches];

    GeneTreeWithinSpeciesTree geneTreeTest;
    TaxonSet speciesSuperSet;
    final double allowedError = 10e-6;

    final TreeParser speciesTree = new TreeParser();
    final TreeParser geneTree = new TreeParser();
    
    final List<GeneTreeWithinSpeciesTree> geneTreeList = new ArrayList<GeneTreeWithinSpeciesTree>();
    
    MultispeciesPopulationModel populationModel;
    MultispeciesCoalescent msc;

    final double ploidy = 2.0;

    final double popSize = 0.3;
    final double expectedLogHR = Math.log(0.5) * 2; // this should be the right answer (calculated by hand)
    RealParameter popSizesParameter;

    private static List<Taxon> superSetList() throws Exception {
        List<Taxon> superSetList = new ArrayList<>();

        List<Taxon> taxonListA = new ArrayList<>();
        List<Taxon> taxonListB = new ArrayList<>();
        List<Taxon> taxonListC = new ArrayList<>();

        taxonListA.add(new Taxon("A1"));
        taxonListA.add(new Taxon("A2"));
        taxonListB.add(new Taxon("B1"));
        taxonListB.add(new Taxon("B2"));
        taxonListC.add(new Taxon("C1"));

        superSetList.add(new TaxonSet("A", taxonListA));
        superSetList.add(new TaxonSet("B", taxonListB));
        superSetList.add(new TaxonSet("C", taxonListC));

        return superSetList;
    }

    private void initializeTrees() throws Exception {
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
        geneTree.initByName("newick", newickGeneTreeA, "IsLabelledNewick", true);

        speciesSuperSet = new TaxonSet(superSetList());

        geneTreeTest = new GeneTreeWithinSpeciesTree();

        geneTreeTest.initByName("tree", geneTree, "ploidy", ploidy);

        geneTreeList.add(geneTreeTest);
    }

    @Test
    public void testLogP() throws Exception {
        initializeTrees();

        popSizesParameter = new RealParameter();
        popSizesParameter.initByName("value", String.valueOf(popSize));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", popSizesParameter);
        state.initialise();

        populationModel = new ConstantPopulation();
        populationModel.initByName("popSizes", popSizesParameter);

        msc = new MultispeciesCoalescent();
        msc.initByName("tree", speciesTree, "geneTree", geneTreeList, "taxonSuperSet", speciesSuperSet, "populationModel", populationModel);

        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        Node brother = null;
        for (Node n: speciesTree.getRoot().getAllLeafNodes()) {
            if (n.getID().equals("B")) {
                brother = n;
            }
        }

        CoordinatedExchange coex = new CoordinatedExchange();
        coex.initByName("multispeciesCoalescent", msc);
        coex.manipulateSpeciesTree(brother);
        final double calculatedLogHR = coex.rearrangeGeneTrees(msc);

        assertEquals(expectedLogHR, calculatedLogHR, allowedError);
    }
}
