package sb2tests;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import starbeast2.ConstantPopulationIO;
import starbeast2.ContinuousRates;
import starbeast2.GeneTreeWithinSpeciesTree;
import starbeast2.MultispeciesCoalescent;
import starbeast2.StarBeastClock;

public class ContinuousRatesTest extends PopulationTestHelper {
    private final String newickSpeciesTree = "((a:1.5,b:1.5):0.5,c:2.0)";
    private final String newickGeneTreeA = "((((a1:0.3,a2:0.3):1.6,(b1:1.8,b2:1.8):0.1):0.5,c1:2.4):0.6,c2:3.0)";

    private final double popSize = 0.3;
    private final double alpha = 1.5;
    private final double beta = 2.5;
    private final double geneRate = 1.5;
    private final double initialSpeciesRate = 1.0;

    private RealParameter alphaParameter;
    private RealParameter betaParameter;
    private RealParameter geneRateParameter;
    private RealParameter speciesRateParameter;

    private StarBeastClock clockModel;
    private ContinuousRates speciesRates;

    MultispeciesCoalescent msc;

    @Test
    public void testLogP() throws Exception {
        nSpecies = 3;
        nBranches = (nSpecies * 2) - 1;
        popSizes = new double[nBranches];

        initializeTrees(newickSpeciesTree, newickGeneTreeA);

        alphaParameter = new RealParameter();
        betaParameter = new RealParameter();
        geneRateParameter = new RealParameter();
        speciesRateParameter = new RealParameter();

        alphaParameter.initByName("value", String.valueOf(alpha));
        betaParameter.initByName("value", String.valueOf(beta));
        geneRateParameter.initByName("value", String.valueOf(geneRate));
        speciesRateParameter.initByName("value", String.valueOf(initialSpeciesRate));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", alphaParameter);
        state.initByName("stateNode", betaParameter);
        state.initByName("stateNode", geneRateParameter);
        state.initByName("stateNode", speciesRateParameter);
        state.initialise();

        populationModel = new ConstantPopulationIO();
        populationModel.initByName("alpha", alphaParameter, "beta", betaParameter);
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("tree", speciesTree, "geneTree", geneTreeList, "taxonSuperSet", speciesSuperSet, "populationModel", populationModel);

        speciesRates = new ContinuousRates();
        speciesRates.initByName("tree", speciesTree, "rates", speciesRateParameter);
        initializeRates();

        clockModel = new StarBeastClock();
        clockModel.initByName("multispeciesCoalescent", msc, "speciesTreeRates", speciesRates, "geneTreeRate", geneRateParameter);
        checkRates();
    }
    
    private void initializeTrees(final String newickSpeciesTree, final String newickGeneTree) throws Exception {
        String speciesCodes = "abc";
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
        geneTreeA.initByName("newick", newickGeneTree, "IsLabelledNewick", true);

        List<Taxon> superSetList = new ArrayList<>();
        for (int i = 0; i < nSpecies; i++) {
            final String speciesName = speciesCodes.substring(i, i + 1);
            List<Taxon> taxonList = new ArrayList<>();
            for (int j = 0; j < individualsPerSpecies; j++) {
                final String taxonName = String.format("%s%d", speciesName, j + 1);
                taxonList.add(new Taxon(taxonName));
            }
            superSetList.add(new TaxonSet(speciesName, taxonList));
        }
        speciesSuperSet = new TaxonSet(superSetList);

        geneTreeTestA = new GeneTreeWithinSpeciesTree();

        geneTreeTestA.initByName("tree", geneTreeA, "ploidy", ploidy);
        
        geneTreeList.add(geneTreeTestA);
    }
    
    private void initializeRates() {
        String[] node0 = {"a"};
        String[] node1 = {"b"};
        String[] node2 = {"c"};
        String[] node3 = {"a", "b"};
        String[] node4 = {"a", "b", "c"};

        double rate0 = 0.8;
        double rate1 = 1.2;
        double rate2 = 1.0;
        double rate3 = 1.6;
        double rate4 = 0.5;

        assertTrue(speciesRates.setRate(node0, rate0));
        assertTrue(speciesRates.setRate(node1, rate1));
        assertTrue(speciesRates.setRate(node2, rate2));
        assertTrue(speciesRates.setRate(node3, rate3));
        assertTrue(speciesRates.setRate(node4, rate4));
    }

    private void checkRates() {
        String[] node00 = {"a1"};
        String[] node01 = {"a2"};
        String[] node02 = {"b1"};
        String[] node03 = {"b2"};
        String[] node04 = {"c1"};
        String[] node05 = {"c2"};
        String[] node06 = {"a1", "a2"};
        String[] node07 = {"b1", "b2"};
        String[] node08 = {"a1", "a2", "b1", "b2"};
        String[] node09 = {"a1", "a2", "b1", "b2", "c1"};
        String[] node10 = {"a1", "a2", "b1", "b2", "c1", "c2"};

        double rate00 = 1.2;
        double rate01 = 1.2;
        double rate02 = 1.9;
        double rate03 = 1.9;
        double rate04 = 1.375;
        double rate05 = 1.25;
        double rate06 = 1.5;
        double rate07 = 2.4;
        double rate08 = 1.08;
        double rate09 = 0.75;
        double rate10 = 1.5;

        assertEquals(rate00, clockModel.getRate(geneTreeA, node00), allowedError);
        assertEquals(rate01, clockModel.getRate(geneTreeA, node01), allowedError);
        assertEquals(rate02, clockModel.getRate(geneTreeA, node02), allowedError);
        assertEquals(rate03, clockModel.getRate(geneTreeA, node03), allowedError);
        assertEquals(rate04, clockModel.getRate(geneTreeA, node04), allowedError);
        assertEquals(rate05, clockModel.getRate(geneTreeA, node05), allowedError);
        assertEquals(rate06, clockModel.getRate(geneTreeA, node06), allowedError);
        assertEquals(rate07, clockModel.getRate(geneTreeA, node07), allowedError);
        assertEquals(rate08, clockModel.getRate(geneTreeA, node08), allowedError);
        assertEquals(rate09, clockModel.getRate(geneTreeA, node09), allowedError);
        assertEquals(rate10, clockModel.getRate(geneTreeA, node10), allowedError);
    }
}
