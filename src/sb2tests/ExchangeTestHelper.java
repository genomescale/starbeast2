package sb2tests;

import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import org.junit.Test;
import starbeast2.*;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

abstract class ExchangeTestHelper {
    String newickSpeciesTree;
    List<String> newickGeneTrees = new ArrayList<>();

    List<TreeParser> geneTrees = new ArrayList<>();
    
    SpeciesTreeParser speciesTreeWrapper;
    List<GeneTree> geneTreeWrappers = new ArrayList<>();

    RealParameter popsizeParameter;
    ConstantPopulations populationModel;
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
        for (Tree gt: geneTrees) {
            state.initByName("stateNode", gt);
        }
        state.initialise();

        populationModel = new ConstantPopulations();
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
        coex.initByName("speciesTree", speciesTreeWrapper, "geneTree", geneTrees, "testing", true, "weight", 1.0);
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
