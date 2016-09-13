package sb2tests;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import starbeast2.MultispeciesCoalescent;
import starbeast2.GeneTree;
import starbeast2.SpeciesTreeParser;

public class MissingDataConstantIO {
    private String newickSpeciesTree;
    private List<String> newickGeneTrees = new ArrayList<>();

    private SpeciesTreeParser speciesTree;
    private List<TreeParser> geneTrees = new ArrayList<>();
    private List<GeneTree> geneTreeWrappers = new ArrayList<>();

    private MultispeciesCoalescent msc;

    private State state;
    private RealParameter alphaParameter;
    private RealParameter meanParameter;

    private double ploidy;
    private double expectedLogP;
    private int nSpecies;

    private final double allowedError = 10e-6;
    private final double alpha = 1.5;
    private final double beta = 2.5;
    private final double mean = beta / (alpha - 1.0);
    private final int individualsPerSpecies = 2;

    public MissingDataConstantIO() throws Exception {
        ploidy = 2.0;
        nSpecies = 4;
        expectedLogP = -10.956285249389675; // haven't checked this is the right answer

        state = new State();
        alphaParameter = new RealParameter();
        meanParameter = new RealParameter();

        alphaParameter.initByName("value", String.valueOf(alpha));
        meanParameter.initByName("value", String.valueOf(mean));
        
        newickSpeciesTree = "((s0:0.32057156677143211,s3:0.32057156677143211):1.2653250035015629,(s1:0.56540722294658641,s2:0.56540722294658641):1.0204893473264085)";
        newickGeneTrees.add("((((s0_tip1:0.3416660303037105,s3_tip0:0.3416660303037105):0.024561190897159135,s0_tip0:0.36622722120086965):0.0643095990846464,s3_tip1:0.43053682028551604):1.4201019862262891,(s2_tip0:0.19897724687831703,s2_tip1:0.19897724687831703):1.651661559633488)");
        newickGeneTrees.add("((s3_tip0:0.09482581277282173,s3_tip1:0.09482581277282173):1.6017973588278296,((s1_tip0:0.33170960882423645,s1_tip1:0.33170960882423645):0.29497523293318856,(s2_tip0:0.2908611340994834,s2_tip1:0.2908611340994834):0.3358237076579416):1.0699383298432266)");
    }

    @Test
    public void testLogP() throws Exception {
        TaxonSet speciesSuperset = generateSuperset();
        initialize(speciesSuperset);

        final double calculatedLogP = msc.calculateLogP();
        // System.out.println(String.format("expected %f, calculated %f", expectedLogP, calculatedLogP));
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }

    public TaxonSet generateSuperset() {
        List<Taxon> superSetList = new ArrayList<>();
        for (int i = 0; i < nSpecies; i++) {
            final String speciesName = String.format("s%d", i);
            List<Taxon> taxonList = new ArrayList<>();
            for (int j = 0; j < individualsPerSpecies; j++) {
                final String taxonName = String.format("s%d_tip%d", i, j);
                taxonList.add(new Taxon(taxonName));
            }
            superSetList.add(new TaxonSet(speciesName, taxonList));
        }

        TaxonSet speciesSuperset = new TaxonSet(superSetList);
        
        return speciesSuperset;
    }

    private void initialize(TaxonSet speciesSuperset) throws Exception {
        speciesTree = new SpeciesTreeParser();
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true, "taxonset", speciesSuperset);

        for (String geneTreeNewick: newickGeneTrees) {
            TreeParser geneTree = new TreeParser();
            geneTree.initByName("newick", geneTreeNewick, "IsLabelledNewick", true);
            geneTrees.add(geneTree);

            GeneTree geneTreeWrapper = new GeneTree();
            geneTreeWrapper.initByName("tree", geneTree, "ploidy", ploidy, "speciesTree", speciesTree);
            geneTreeWrappers.add(geneTreeWrapper);
        }

        state.initByName("stateNode", alphaParameter);
        state.initByName("stateNode", meanParameter);

        msc = new MultispeciesCoalescent();
        msc.initByName("populationShape", alphaParameter, "populationMean", meanParameter, "distribution", geneTreeWrappers);

        for (GeneTree gt: geneTreeWrappers) gt.calculateLogP();
    }
}
