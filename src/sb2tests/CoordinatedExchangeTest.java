package sb2tests;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import starbeast2.ConstantPopulationIO;
import starbeast2.GeneTreeWithinSpeciesTree;
import starbeast2.MultispeciesCoalescent;
import starbeast2.MultispeciesPopulationModel;

public class CoordinatedExchangeTest {
    private final int individualsPerSpecies = 2;
    private final double allowedError = 10e-6;
    private final double expectedLogHR = 0.0; // this should be the right answer (calculated by hand)
    private final String newickGeneTree = "(('T8_02':2.33805784117,'T8_01':2.33805784117):4.29401094668,((('T2_02':1.50982163708,('T1_02':1.46121467289,(('T3_02':0.671306265393,'T3_01':0.671306265393):0.313144060841,'T2_01':0.984450326234):0.476764346659):0.0486069641853):0.326314720888,(('T7_01':0.756964337325,(('T5_01':0.261481658197,'T5_02':0.261481658197):0.137745763082,('T4_02':0.154201972895,'T4_01':0.154201972895):0.245025448384):0.357736916046):0.837134282273,('T7_02':0.85642850243,('T6_01':0.311122292472,'T6_02':0.311122292472):0.545306209958):0.737670117167):0.242037738369):0.247686919048,'T1_01':2.08382327701):4.54824551085):0.0;";
    private final String newickSpeciesTree = "(((T1:1.39640717064,(T2:0.760791939917,T3:0.760791939917):0.635615230727):0.118698784729,((T4:0.22999140259,T5:0.22999140259):0.502476265791,(T6:0.11342184614,T7:0.11342184614):0.619045822241):0.782638286992):1.25763388932,T8:2.77273984469):0.607360813616;";

    GeneTreeWithinSpeciesTree geneTreeTest;
    TaxonSet speciesSuperSet;

    MultispeciesCoalescent msc;

    private final TreeParser speciesTree = new TreeParser();
    private final TreeParser geneTree = new TreeParser();
    
    private final List<GeneTreeWithinSpeciesTree> geneTreeList = new ArrayList<GeneTreeWithinSpeciesTree>();
    
    MultispeciesPopulationModel populationModel;

    @Test
    public void testLogHR() throws Exception {
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
        geneTree.initByName("newick", newickGeneTree, "IsLabelledNewick", true);

        final int nSpecies = speciesTree.getNodeCount();

        List<Taxon> superSetList = new ArrayList<>();
        for (int i = 0; i < nSpecies; i++) {
            final String speciesName = String.format("T%d", i);
            List<Taxon> taxonList = new ArrayList<>();
            for (int j = 0; j < individualsPerSpecies; j++) {
                final String taxonName = String.format("T%d_0%d", i, j);
                taxonList.add(new Taxon(taxonName));
            }
            superSetList.add(new TaxonSet(speciesName, taxonList));
        }

        speciesSuperSet = new TaxonSet(superSetList);

        geneTreeTest = new GeneTreeWithinSpeciesTree();
        geneTreeTest.initByName("tree", geneTree);
        geneTreeList.add(geneTreeTest);

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initialise();

        populationModel = new ConstantPopulationIO();
        populationModel.initByName("alpha", "2.0", "beta", "2.0");

        msc = new MultispeciesCoalescent();
        msc.initByName("tree", speciesTree, "geneTree", geneTreeList, "taxonSuperSet", speciesSuperSet, "populationModel", populationModel);

        final double calculatedLogHR = coex.calculateLogP();
        assertEquals(expectedLogHR, calculatedLogHR, allowedError);
    }
}
