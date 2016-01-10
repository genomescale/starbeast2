package sb2tests;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.util.TreeParser;
import starbeast2.ConstantPopulationIO;
import starbeast2.DiscreteRates;
import starbeast2.GeneTree;
import starbeast2.MultispeciesCoalescent;
import starbeast2.MultispeciesPopulationModel;
import starbeast2.SpeciesTree;
import starbeast2.StarBeastClock;

public class DiscreteRatesTest {
    private String newickSpeciesTree;
    private List<String> newickGeneTrees = new ArrayList<>();

    private TreeParser speciesTree = new TreeParser();
    private List<TreeParser> geneTrees = new ArrayList<>();
    
    private SpeciesTree speciesTreeWrapper = new SpeciesTree();
    private List<GeneTree> geneTreeWrappers = new ArrayList<>();

    private final int nSpecies = 3;
    private final int individualsPerSpecies = 2;
    private final double popSize = 0.3;
    private final double alpha = 1.5;
    private final double beta = 2.5;
    private final double geneRate = 1.5;
    private final int initialSpeciesRate = 1;
    private final double ploidy = 2.0;

    final double allowedError = 10e-6;

    private RealParameter alphaParameter;
    private RealParameter betaParameter;
    private RealParameter geneRateParameter;
    private IntegerParameter speciesRateParameter;

    private StarBeastClock clockModel;
    private DiscreteRates speciesRates;
    private LogNormalDistributionModel speciesRateDistribution;

    MultispeciesCoalescent msc;

    public DiscreteRatesTest() {
        newickSpeciesTree = "((a:1.5,b:1.5):0.5,c:2.0)";
        newickGeneTrees.add("((((a1:0.3,a2:0.3):1.6,(b1:1.8,b2:1.8):0.1):0.5,c1:2.4):0.6,c2:3.0)");
    }

    @Test
    public void testRates() throws Exception {
        TaxonSet speciesSuperSet = generateSuperset();
        
        initializeSpeciesTree(speciesSuperSet);
        initializeGeneTrees();

        alphaParameter = new RealParameter();
        betaParameter = new RealParameter();
        geneRateParameter = new RealParameter();
        speciesRateParameter = new IntegerParameter();

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

        final int nBranches = (2 * nSpecies) - 1;
        MultispeciesPopulationModel populationModel = new ConstantPopulationIO();
        populationModel.initByName("alpha", alphaParameter, "beta", betaParameter);
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("speciesTree", speciesTreeWrapper, "geneTree", geneTreeWrappers, "populationModel", populationModel);

        speciesRateDistribution = new LogNormalDistributionModel();
        speciesRateDistribution.initByName("M", "1.0", "S", "1.0", "meanInRealSpace", true);

        speciesRates = new DiscreteRates();
        speciesRates.initByName("tree", speciesTree, "rates", speciesRateParameter, "distr", speciesRateDistribution);
        initializeRates();

        clockModel = new StarBeastClock();
        clockModel.initByName("geneTree", geneTreeWrappers.get(0), "speciesTreeRates", speciesRates, "clock.rate", geneRateParameter);

        checkRates();
    }
    
    private TaxonSet generateSuperset() throws Exception {
        String speciesCodes = "abc";
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

        TaxonSet speciesSuperSet = new TaxonSet(superSetList);
        
        return speciesSuperSet;
    }

    private void initializeRates() {
        String[] node0 = {"a"};
        String[] node1 = {"b"};
        String[] node2 = {"c"};
        String[] node3 = {"a", "b"};
        String[] node4 = {"a", "b", "c"};

        int rate0 = 4;
        int rate1 = 0;
        int rate2 = 3;
        int rate3 = 2;
        int rate4 = 1;

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

        // correct gene tree rates were calculated by hand
        double rate00 = 3.2772894;
        double rate01 = 3.2772894;
        double rate02 = 0.3621035;
        double rate03 = 0.3621035;
        double rate04 = 1.370629;
        double rate05 = 1.204206;
        double rate06 = 2.685416;
        double rate07 = 0.909796;
        double rate08 = 0.6127731;
        double rate09 = 0.5385174;
        double rate10 = 1.5;

        Tree geneTree = geneTrees.get(0);
        assertEquals(rate00, clockModel.getRate(geneTree, node00), allowedError);
        assertEquals(rate01, clockModel.getRate(geneTree, node01), allowedError);
        assertEquals(rate02, clockModel.getRate(geneTree, node02), allowedError);
        assertEquals(rate03, clockModel.getRate(geneTree, node03), allowedError);
        assertEquals(rate04, clockModel.getRate(geneTree, node04), allowedError);
        assertEquals(rate05, clockModel.getRate(geneTree, node05), allowedError);
        assertEquals(rate06, clockModel.getRate(geneTree, node06), allowedError);
        assertEquals(rate07, clockModel.getRate(geneTree, node07), allowedError);
        assertEquals(rate08, clockModel.getRate(geneTree, node08), allowedError);
        assertEquals(rate09, clockModel.getRate(geneTree, node09), allowedError);
        assertEquals(rate10, clockModel.getRate(geneTree, node10), allowedError);
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

            GeneTree geneTreeWrapper = new GeneTree();
            geneTreeWrapper.initByName("tree", geneTree, "ploidy", ploidy, "speciesTree", speciesTreeWrapper);
            geneTreeWrappers.add(geneTreeWrapper);
        }
    }
}
