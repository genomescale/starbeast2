package sb2tests;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import starbeast2.ConstantPopulation;
import starbeast2.CoordinatedExchange;
import starbeast2.GeneTree;
import starbeast2.MultispeciesCoalescent;
import starbeast2.PopulationModel;
import starbeast2.SpeciesTreeParser;

abstract class ExchangeTestHelper {
    String newickSpeciesTree;
    List<String> newickGeneTrees = new ArrayList<>();

    List<TreeParser> geneTrees = new ArrayList<>();
    
    SpeciesTreeParser speciesTreeWrapper;
    List<GeneTree> geneTreeWrappers = new ArrayList<>();

    RealParameter popsizeParameter;
    PopulationModel populationModel;
    MultispeciesCoalescent msc;

    double ploidy;
    double popSize;
    double expectedLogHR;
    String bTipLabel;
    String cTipLabel;
    boolean bIsParent;
    boolean cIsParent;

    final double allowedError = 10e-6;

    abstract public TaxonSet generateSuperset() throws Exception;

    @Test
    public void testLogHR() throws Exception {
        TaxonSet speciesSuperSet = generateSuperset();
        initializeSpeciesTree(speciesSuperSet);
        initializeGeneTrees();

        popsizeParameter = new RealParameter();
        popsizeParameter.initByName("value", String.valueOf(popSize));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", popsizeParameter);
        state.initialise();

        populationModel = new ConstantPopulation();
        populationModel.initByName("populationSizes", popsizeParameter, "speciesTree", speciesTreeWrapper);

        int nBranches = speciesTreeWrapper.getNodeCount();
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        Node cNode = null;
        Node bNode = null;
        for (Node n: speciesTreeWrapper.getRoot().getAllLeafNodes()) {
            if (n.getID().equals(bTipLabel)) {
                if (bIsParent) bNode = n.getParent();
                else bNode = n;
            } else if (n.getID().equals(cTipLabel)) {
                if (cIsParent) cNode = n.getParent();
                else cNode = n;
            }
        }

        Node yNode = bNode.getParent();
        Node zNode = yNode.getParent();
        Node aNode = (bNode == yNode.getRight()) ? yNode.getLeft() : yNode.getRight();

        if (cNode == null) {
            cNode = (yNode == zNode.getRight()) ? zNode.getLeft() : zNode.getRight();
        }

        CoordinatedExchange coex = new CoordinatedExchange();
        coex.initByName("tree", speciesTreeWrapper, "geneTree", geneTrees, "testing", true);
        coex.aNode = aNode;
        coex.bNode = bNode;
        coex.cNode = cNode;
        coex.yNode = yNode;
        coex.zNode = zNode;
        final double calculatedLogHR = coex.proposal();

        assertEquals(expectedLogHR, calculatedLogHR, allowedError);
    }

    public void initializeSpeciesTree(TaxonSet speciesSuperSet) throws Exception {
        speciesTreeWrapper = new SpeciesTreeParser();
        speciesTreeWrapper.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true, "taxonset", speciesSuperSet);
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
