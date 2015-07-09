package test;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import msc.ConstantPopulationIO;
import msc.MultispeciesCoalescent;
import test.beast.BEASTTestCase;

public class ConstantIOTest extends PopulationTestHelper {
    private final double popSize = 0.3;
    private final double alpha = 1.5;
    private final double beta = 2.5;
    private final double expectedLogP = -14.5233984762; // this should be the right answer (calculated by hand)

    private RealParameter alphaParameter;
    private RealParameter betaParameter;

    MultispeciesCoalescent msc;

    @Test
    public void testLogP() throws Exception {
        initializeSpeciesTree();
        initializeGeneTrees();

        alphaParameter = new RealParameter();
        betaParameter = new RealParameter();

        alphaParameter.initByName("value", String.valueOf(alpha));
        betaParameter.initByName("value", String.valueOf(beta));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", alphaParameter);
        state.initByName("stateNode", betaParameter);
        state.initialise();

        populationModel = new ConstantPopulationIO();
        populationModel.initByName("alpha", alphaParameter, "beta", betaParameter);
        populationModel.initializePopSizes(speciesTree, popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("tree", speciesTree, "geneTree", geneTreeList, "taxonSuperSet", speciesSuperSet, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, BEASTTestCase.PRECISION);
    }
}
