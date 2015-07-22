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

public class CoordinatedExchangeTest {
    final String newickSpeciesTree = "(T1:1.02081927766,((T2:0.388840623317,T3:0.388840623317):0.477223011734,(T4:0.521093885203,(T5:0.145649902275,(T6:0.0680533968312,T7:0.0680533968312):0.0775965054433):0.375443982929):0.344969749848):0.154755642606):0.902290428254;";
    final String newickGeneTree = "((T1_02:0.285597385824,((T1_03:0.0213928057436,T1_01:0.0213928057436):0.186943459751,T1_04:0.208336265495):0.0772611203298):0.975826233271,(((T4_03:0.207826294361,((T4_04:0.0392482649542,T4_02:0.0392482649542):0.0180238818382,T4_01:0.0572721467924):0.150554147569):0.408902767916,((T6_03:0.166960276323,T6_02:0.166960276323):0.0455888246696,((((T6_04:0.128915161032,T7_02:0.128915161032):0.0195933968798,((T7_04:0.0103290108377,T7_03:0.0103290108377):0.00721837753724,T7_01:0.0175473883749):0.130961169537):0.009992185274,(((T5_03:0.00536372179893,T5_04:0.00536372179893):0.0317199690567,T5_02:0.0370836908556):0.00279977835762,T5_01:0.0398834692132):0.118617273973):0.034382279669,T6_01:0.192883022855):0.0196660781378):0.404179961286):0.295076999889,(((T2_04:0.159180677136,(T2_01:0.00306956690132,T2_03:0.00306956690132):0.156111110235):0.171690618442,T2_02:0.330871295578):0.18855004693,((T3_03:0.000577428309595,T3_04:0.000577428309595):0.28575162438,(T3_02:0.00643512791271,T3_01:0.00643512791271):0.279893924777):0.233092289818):0.39238471966):0.349617556925):0.0;";

    final int nSpecies = 7;
    final int nBranches = (nSpecies * 2) - 1;
    final double[] popSizes = new double[nBranches];

    final int individualsPerSpecies = 4;
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
    final double expectedLogHR = Math.log(2) - Math.log(4); // this should be the right answer (calculated by hand)
    RealParameter popSizesParameter;

    public void initializeTrees() throws Exception {
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
        geneTree.initByName("newick", newickGeneTree, "IsLabelledNewick", true);

        List<Taxon> superSetList = new ArrayList<>();
        for (int i = 1; i <= nSpecies; i++) {
            final String speciesName = String.format("T%d", i);
            List<Taxon> taxonList = new ArrayList<>();
            for (int j = 1; j <= individualsPerSpecies; j++) {
                final String taxonName = String.format("T%d_0%d", i, j);
                taxonList.add(new Taxon(taxonName));
            }
            superSetList.add(new TaxonSet(speciesName, taxonList));
        }
        speciesSuperSet = new TaxonSet(superSetList);

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
            if (n.getID().equals("T6")) {
                brother = n.getParent();
            }
        }

        CoordinatedExchange coex = new CoordinatedExchange();
        coex.initByName("multispeciesCoalescent", msc);
        msc.computeCoalescentTimes();
        coex.manipulateSpeciesTree(brother);
        final double calculatedLogHR = coex.rearrangeGeneTrees(msc);

        assertEquals(expectedLogHR, calculatedLogHR, allowedError);
    }
}
