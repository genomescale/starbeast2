package test;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import msc.ConstantPopulation;
import msc.MultispeciesCoalescent;
import test.beast.BEASTTestCase;

public class ConstantPopulationTest extends PopulationTestHelper {
    private final double popSize = 0.3;
    private final double expectedLogP = 2.0784641098; // this should be the right answer (calculated by hand)
    private RealParameter popSizesParameter;

    MultispeciesCoalescent msc;

    @Test
    public void testLogP() throws Exception {
        initializeSpeciesTree();
        initializeGeneTrees();

        popSizesParameter = new RealParameter();
        popSizesParameter.initByName("value", String.valueOf(popSize));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", popSizesParameter);
        state.initialise();

        populationModel = new ConstantPopulation();
        populationModel.initByName("popSizes", popSizesParameter);
        populationModel.initializePopSizes(speciesTree, popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("tree", speciesTree, "geneTree", geneTreeList, "taxonSuperSet", speciesSuperSet, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, BEASTTestCase.PRECISION);
    }
}
