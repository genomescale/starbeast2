package snetworktests;

import java.util.ArrayList;
import java.util.List;

import beast.evolution.alignment.Taxon;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import speciesnetwork.*;

public class FlipLoopTest {
    String newickSpeciesNetwork;
    List<String> newickGeneTrees = new ArrayList<>();

    TaxonSet speciesSuperset;
    TreeParser speciesTree;
    NetworkParser speciesNetwork;
    List<TreeParser> geneTrees = new ArrayList<>();
    List<IntegerParameter> geneTreeEmbeddings = new ArrayList<>();
    RealParameter gammaParameter;

    int nSpecies;
    int nBranches;
    double popSize;
    double ploidy;
    double gamma;
    State state = null;

    double expectedLogP;
    final double allowedError = 1e-6;

    public FlipLoopTest() throws Exception {
        nSpecies = 3;
        nBranches = 8;
        popSize = 0.1;
        ploidy = 2.0;
        gamma = 0.6; // for branch #H1(5)-#S2(3)

        newickSpeciesNetwork = "((A:0.2,#H1:0.1)S1:0.3,((B:0.1)#H1:0.2,C:0.3)S2:0.2)R";
        newickGeneTrees.add("(((a1:0.07,a2:0.07):0.48,(b1:0.25,b2:0.25):0.30):0.08,(b3:0.35,c1:0.35):0.28)");

        IntegerParameter embedding = new IntegerParameter();
        embedding.initByName("value", "-1 -1 -1 -1  0  1 -1 -1 -1 -1 -1 " +
                                      "-1 -1  1  1 -1 -1  0 -1 -1 -1 -1 " +
                                      "-1 -1  0  0  0 -1 -1 -1 -1 -1 -1 " +
                                      "-1 -1 -1 -1 -1 -1  0  0 -1  1 -1", "dimension", "44", "minordimension", "11");
        geneTreeEmbeddings.add(embedding);
    }

    @Test
    public void testOperator() throws Exception {
        speciesSuperset = generateSuperset();
        initializeSpeciesNetwork();
        initializeStateNodes();
        initializeGeneTrees(true);

        // run the fliploop operator
        for (int i = 0; i < newickGeneTrees.size(); i++) {
            TreeParser geneTree = geneTrees.get(i);
            IntegerParameter embedding = geneTreeEmbeddings.get(i);
            FlipNetworkLoop flipOperator = new FlipNetworkLoop();
            flipOperator.initByName("geneTree", geneTree, "speciesNetwork", speciesNetwork, "embedding", embedding);
            assertEquals(flipOperator.proposal(), 0.0, allowedError);
        }
    }

    public TaxonSet generateSuperset() throws Exception {
        List<Taxon> superSetList = new ArrayList<>();

        List<Taxon> taxonListA = new ArrayList<>();
        taxonListA.add(new Taxon("a1"));
        taxonListA.add(new Taxon("a2"));
        superSetList.add(new TaxonSet("A", taxonListA));

        List<Taxon> taxonListB = new ArrayList<>();
        taxonListB.add(new Taxon("b1"));
        taxonListB.add(new Taxon("b2"));
        taxonListB.add(new Taxon("b3"));
        superSetList.add(new TaxonSet("B", taxonListB));

        List<Taxon> taxonListC = new ArrayList<>();
        taxonListC.add(new Taxon("c1"));
        superSetList.add(new TaxonSet("C", taxonListC));

        return new TaxonSet(superSetList);
    }

    public void initializeSpeciesNetwork() throws Exception {
        speciesTree = new TreeParser();
        speciesTree.initByName("newick", newickSpeciesNetwork, "IsLabelledNewick", true, "adjustTipHeights", false);
        speciesNetwork = new NetworkParser();
        speciesNetwork.initByName("tree", speciesTree);
    }

    public void initializeGeneTrees(boolean reembed) throws Exception {
        for (int i = 0; i < newickGeneTrees.size(); i++) {
            final String geneTreeNewick = newickGeneTrees.get(i);
            TreeParser geneTree = new TreeParser();
            geneTree.initByName("newick", geneTreeNewick, "IsLabelledNewick", true);
            geneTrees.add(geneTree);
            if (reembed) { // rebuild the embedding
                IntegerParameter embedding = geneTreeEmbeddings.get(i);
                RebuildEmbedding rebuildOperator = new RebuildEmbedding();
                rebuildOperator.initByName("geneTree", geneTree, "speciesNetwork", speciesNetwork,
                                           "taxonSuperset", speciesSuperset, "embedding", embedding);
            }
            /*
            GeneTreeInSpeciesNetwork geneTreeWrapper = new GeneTreeInSpeciesNetwork();
            geneTreeWrapper.initByName("geneTree", geneTree, "ploidy", ploidy, "speciesNetwork", speciesNetwork,
                                       "embedding", embedding, "gamma", gammaParameter);
            geneTreeWrappers.add(geneTreeWrapper); */
        }
    }

    public void initializeStateNodes() throws Exception {
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
}
