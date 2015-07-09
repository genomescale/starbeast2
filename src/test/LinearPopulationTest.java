package test;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import msc.LinearPopulation;
import msc.MultispeciesCoalescent;
import test.beast.BEASTTestCase;

public class LinearPopulationTest extends PopulationTestHelper {
    private final double popSize = 0.3;
    private final double expectedLogP = 3.622793963076785; // I have no idea if this is the right answer
    private RealParameter topPopSizesParameter;
    private RealParameter tipPopSizesParameter;

    MultispeciesCoalescent msc;

    @Test
    public void testLogP() throws Exception {
        initializeSpeciesTree();
        initializeGeneTrees();

        topPopSizesParameter = new RealParameter();
        tipPopSizesParameter = new RealParameter();
        topPopSizesParameter.initByName("value", String.valueOf(popSize));
        tipPopSizesParameter.initByName("value", String.valueOf(popSize));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", topPopSizesParameter);
        state.initByName("stateNode", tipPopSizesParameter);
        state.initialise();

        populationModel = new LinearPopulation();
        populationModel.initByName("topPopSizes", topPopSizesParameter, "tipPopSizes", tipPopSizesParameter);
        populationModel.initializePopSizes(speciesTree, popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("tree", speciesTree, "geneTree", geneTreeList, "taxonSuperSet", speciesSuperSet, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, BEASTTestCase.PRECISION);
    }
}
