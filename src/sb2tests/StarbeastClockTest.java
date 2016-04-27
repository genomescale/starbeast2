package sb2tests;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.math.distributions.LogNormalDistributionModel;
import beast.util.TreeParser;
import starbeast2.DiscreteRates;
import starbeast2.GeneTree;
import starbeast2.SpeciesTree;
import starbeast2.StarBeastClock;

public class StarbeastClockTest {
    private String newickSpeciesTree;
    private String newickGeneTree;

    private TreeParser speciesTree = new TreeParser();
    private TreeParser geneTree = new TreeParser();
    
    private SpeciesTree speciesTreeWrapper = new SpeciesTree();
    private GeneTree geneTreeWrapper = new GeneTree();


    private final int nSpecies = 3;
    private final int individualsPerSpecies = 2;
    private final int initialBranchRate = 1;
    private final double meanRate = 1.5;
    private final double ploidy = 2.0;

    final double allowedError = 10e-6;

    private RealParameter meanRateParameter;
    private IntegerParameter branchRatesParameter;

    private StarBeastClock geneTreeClock;
    private DiscreteRates speciesTreeClock;
    private LogNormalDistributionModel speciesRateDistribution;

    public StarbeastClockTest() {
        newickSpeciesTree = "((a:1.5,b:1.5):0.5,c:2.0)";
        newickGeneTree = "((((a1:0.3,a2:0.3):1.6,(b1:1.8,b2:1.8):0.1):0.5,c1:2.4):0.6,c2:3.0)";
    }

    @Test
    public void testRates() throws Exception {
        TaxonSet speciesSuperSet = generateSuperset();
        initializeTrees(speciesSuperSet);

        meanRateParameter = new RealParameter();
        branchRatesParameter = new IntegerParameter();

        meanRateParameter.initByName("value", String.valueOf(meanRate));
        branchRatesParameter.initByName("value", String.valueOf(initialBranchRate));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", meanRateParameter);
        state.initByName("stateNode", branchRatesParameter);
        state.initialise();

        speciesRateDistribution = new LogNormalDistributionModel();
        speciesRateDistribution.initByName("M", "1.0", "S", "1.0", "meanInRealSpace", true);

        speciesTreeClock = new DiscreteRates();
        speciesTreeClock.initByName("tree", speciesTree, "rates", branchRatesParameter, "distr", speciesRateDistribution, "estimateRoot", true);
        initializeRates();

        geneTreeClock = new StarBeastClock();
        geneTreeClock.initByName("geneTree", geneTreeWrapper, "speciesTreeRates", speciesTreeClock, "clock.rate", meanRateParameter);
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
        String[] node00 = {"a"};
        String[] node01 = {"b"};
        String[] node02 = {"c"};
        String[] node03 = {"a", "b"};
        String[] node04 = {"a", "b", "c"};

        int rate00 = 0; // corresponds to 0.1683767
        int rate01 = 4; // corresponds to 2.1848596
        int rate02 = 4; // corresponds to 2.1848596
        int rate03 = 3; // corresponds to 1.0247006
        int rate04 = 1; // corresponds to 0.3590116

        assertTrue(setRate(rate00, speciesTree, node00));
        assertTrue(setRate(rate01, speciesTree, node01));
        assertTrue(setRate(rate02, speciesTree, node02));
        assertTrue(setRate(rate03, speciesTree, node03));
        assertTrue(setRate(rate04, speciesTree, node04));
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

        // correct tree rates were calculated by hand
        double rate00 = 0.2525650;
        double rate01 = 0.2525650;
        double rate02 = 2.98725;
        double rate03 = 2.98725;
        double rate04 = 2.820827;
        double rate05 = 2.364365;
        double rate06 = 0.5736865;
        double rate07 = 1.537051;
        double rate08 = 0.7382241;
        double rate09 = 0.5385174;
        double rate10 = 1.5;

        assertEquals(rate00, getRate(geneTree, node00), allowedError);
        assertEquals(rate01, getRate(geneTree, node01), allowedError);
        assertEquals(rate02, getRate(geneTree, node02), allowedError);
        assertEquals(rate03, getRate(geneTree, node03), allowedError);
        assertEquals(rate04, getRate(geneTree, node04), allowedError);
        assertEquals(rate05, getRate(geneTree, node05), allowedError);
        assertEquals(rate06, getRate(geneTree, node06), allowedError);
        assertEquals(rate07, getRate(geneTree, node07), allowedError);
        assertEquals(rate08, getRate(geneTree, node08), allowedError);
        assertEquals(rate09, getRate(geneTree, node09), allowedError);
        assertEquals(rate10, getRate(geneTree, node10), allowedError);
    }

    private boolean setRate(int rate, TreeParser tree, String[] target) {
        final Node targetNode = findNode(tree, target);
        if (targetNode == null) {
            return false;
        } else {
            branchRatesParameter.setValue(targetNode.getNr(), rate);
            return true;
        }
    }

    private double getRate(TreeParser tree, String[] target) {
        final Node targetNode = findNode(tree, target);
        return geneTreeClock.getRateForBranch(targetNode);
    }

    private Node findNode(final TreeParser tree, final String[] targetArray) {
        final Node[] treeNodes = tree.getNodesAsArray();
        final Set<String> targetSet = new HashSet<>();
        for (int i = 0; i < targetArray.length; i++) {
            targetSet.add(targetArray[i]);
        }

        for (Node node: treeNodes) {
            Set<String> nodeSet = new HashSet<>();

            if (node.isLeaf()) {
                nodeSet.add(node.getID());
            } else {
                final List<Node> leafNodes = node.getAllLeafNodes();
                for (Node leaf: leafNodes) {
                    nodeSet.add(leaf.getID());
                }
            }

            if (targetSet.equals(nodeSet)) {
                return node;
            }
        }

        return null;
    }

    public void initializeTrees(TaxonSet speciesSuperSet) throws Exception {
        speciesTree = new TreeParser();
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true, "taxonset", speciesSuperSet);
        speciesTreeWrapper = new SpeciesTree();
        speciesTreeWrapper.initByName("tree", speciesTree);

        geneTree = new TreeParser();
        geneTree.initByName("newick", newickGeneTree, "IsLabelledNewick", true);

        geneTreeWrapper = new GeneTree();
        geneTreeWrapper.initByName("tree", geneTree, "ploidy", ploidy, "speciesTree", speciesTreeWrapper);
    }
}
