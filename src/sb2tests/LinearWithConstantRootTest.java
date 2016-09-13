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
import starbeast2.GeneTree;
import starbeast2.LinearWithConstantRoot;
import starbeast2.SpeciesTreeParser;

public class LinearWithConstantRootTest {
    private String newickSpeciesTree;
    private List<String> newickGeneTrees = new ArrayList<>();

    private SpeciesTreeParser speciesTree;

    private LinearWithConstantRoot popModel;

    private State state;
    private RealParameter tipPopSizesParameter;
    private RealParameter topPopSizesParameter;

    private double ploidy;
    private double popSize;
    private double expectedLogP;
    private int nSpecies;

    private final double allowedError = 10e-6;
    private final int individualsPerSpecies = 2;

    public LinearWithConstantRootTest() throws Exception {
        popSize = 0.3;
        ploidy = 2.0;
        nSpecies = 4;
        expectedLogP = 2.8879769759752225; // need to double check this at some point

        state = new State();
        tipPopSizesParameter = new RealParameter();
        topPopSizesParameter = new RealParameter();

        newickSpeciesTree = "((s0:0.32057156677143211,s3:0.32057156677143211):1.2653250035015629,(s1:0.56540722294658641,s2:0.56540722294658641):1.0204893473264085)";
        newickGeneTrees.add("((((s0_tip1:0.3416660303037105,s3_tip0:0.3416660303037105):0.024561190897159135,s0_tip0:0.36622722120086965):0.0643095990846464,s3_tip1:0.43053682028551604):1.4201019862262891,((s1_tip0:0.14473698225381706,s1_tip1:0.14473698225381706):0.5135479407233198,(s2_tip0:0.19897724687831703,s2_tip1:0.19897724687831703):0.4593076760988198):1.1923538835346683)");
        newickGeneTrees.add("(((s0_tip0:0.04173231934154758,s0_tip1:0.04173231934154758):0.7845256741090114,(s3_tip0:0.09482581277282173,s3_tip1:0.09482581277282173):0.7314321806777372):0.8703651781500925,((s1_tip0:0.33170960882423645,s1_tip1:0.33170960882423645):0.29497523293318856,(s2_tip0:0.2908611340994834,s2_tip1:0.2908611340994834):0.3358237076579416):1.0699383298432266)");
    }

    @Test
    public void testLogP() throws Exception {
        TaxonSet speciesSuperset = generateSuperset();
        speciesTree = new SpeciesTreeParser();
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true, "taxonset", speciesSuperset);
        state.initByName("stateNode", speciesTree);

        tipPopSizesParameter.initByName("dimension", String.valueOf(nSpecies), "value", 1.0);
        topPopSizesParameter.initByName("dimension", String.valueOf(nSpecies * 2 - 2), "value", 1.0);
        state.initByName("stateNode", tipPopSizesParameter);
        state.initByName("stateNode", topPopSizesParameter);
        state.initialise();

        popModel = new LinearWithConstantRoot();
        popModel.initByName("tipPopulationSizes", tipPopSizesParameter, "topPopulationSizes", topPopSizesParameter, "speciesTree", speciesTree);
        popModel.initPopSizes(popSize);

        double calculatedLogP = 0.0;
        for (String geneTreeNewick: newickGeneTrees) {
            TreeParser geneTree = new TreeParser();
            geneTree.initByName("newick", geneTreeNewick, "IsLabelledNewick", true);

            GeneTree geneTreeWrapper = new GeneTree();
            geneTreeWrapper.initByName("tree", geneTree, "ploidy", ploidy, "speciesTree", speciesTree, "populationModel", popModel);
            calculatedLogP += geneTreeWrapper.calculateLogP();
        }

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
}
