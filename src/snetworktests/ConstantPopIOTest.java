package snetworktests;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import speciesnetwork.ConstantPopulationIO;
import speciesnetwork.PopulationSizeModel;

public class ConstantPopIOTest extends PopulationTestHelper {

    public ConstantPopIOTest() throws Exception {
        expectedLogP = -2.010226796; // this should be the right answer (calculated by hand)

        nSpecies = 3;
        nBranches = 8;
        popSize = 0.1;
        ploidy = 2.0;
        gamma = 0.6; // for branch #H1(5)-#S2(3)

        newickSpeciesNetwork = "((A:0.2,#H1:0.1)S1:0.3,((B:0.1)#H1:0.2,C:0.3)S2:0.2)R";
        newickGeneTrees.add("(((a1:0.07,a2:0.07):0.48,(b1:0.25,b2:0.25):0.30):0.08,(b3:0.35,c1:0.35):0.28)");
        newickGeneTrees.add("((((a1:0.10,a2:0.10):0.50,(b1:0.05,b2:0.05):0.55):0.05,b3:0.65):0.05,c1:0.70)");

        IntegerParameter embedding1 = new IntegerParameter();
        IntegerParameter embedding2 = new IntegerParameter();
        embedding1.initByName("value", "-1 -1 -1 -1  0  1 -1 -1 -1 -1 -1 " +
                                       "-1 -1  1  1 -1 -1  0 -1 -1 -1 -1 " +
                                       "-1 -1  0  0  0 -1 -1 -1 -1 -1 -1 " +
                                       "-1 -1 -1 -1 -1 -1  0  0 -1  1 -1", "dimension", "44", "minordimension", "11");
        embedding2.initByName("value", "-1 -1 -1 -1  0  1 -1  0 -1 -1 -1 " +
                                       "-1 -1 -1 -1 -1 -1  0 -1 -1 -1 -1 " +
                                       "-1 -1 -1 -1  0 -1 -1  0 -1 -1 -1 " +
                                       "-1 -1 -1 -1  1  1  0  1 -1 -1 -1", "dimension", "44", "minordimension", "11");
        geneTreeEmbeddings.add(embedding1);
        geneTreeEmbeddings.add(embedding2);
    }

    @Test
    public void testLogP() throws Exception {
        super.testLogP();
    }

    @Override
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

    @Override
    public PopulationSizeModel generatePopulationModel() throws Exception {
        final double alpha = 5.0, beta = 1.0;
        RealParameter alphaParameter = new RealParameter();
        RealParameter betaParameter = new RealParameter();
        alphaParameter.initByName("value", String.valueOf(alpha));
        betaParameter.initByName("value", String.valueOf(beta));

        state.initByName("stateNode", alphaParameter);
        state.initByName("stateNode", betaParameter);
        state.initialise();

        PopulationSizeModel populationModel = new ConstantPopulationIO();
        populationModel.initByName("alpha", alphaParameter, "beta", betaParameter);
        
        return populationModel;
    }
}
