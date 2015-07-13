package sb2tests;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import starbeast2.LinearWithConstantRoot;
import starbeast2.MultispeciesCoalescent;

public class LinearWithConstantRootTest extends PopulationTestHelper {
    private final double popSize = 0.3;
    private final double expectedLogP = 3.6484899411510012; // I have no idea if this is the right answer
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

        populationModel = new LinearWithConstantRoot();
        populationModel.initByName("topPopSizes", topPopSizesParameter, "tipPopSizes", tipPopSizesParameter);

        msc = new MultispeciesCoalescent();
        msc.initByName("tree", speciesTree, "geneTree", geneTreeList, "taxonSuperSet", speciesSuperSet, "populationModel", populationModel);

        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }
}
