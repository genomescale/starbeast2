package sb2tests;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.core.State;
import beast.util.TreeParser;
import starbeast2.TreeLengthLogger;

public class TreeLengthLoggerTest {
    private String newickTree = "((((a1:0.3,a2:0.3):1.6,(b1:1.8,b2:1.8):0.1):0.5,c1:2.4):0.6,c2:3.0)";
    private TreeParser testTree;

    final double expectedLength = 12.4;
    final double allowedError = 10e-6;
    private TreeLengthLogger treeLengthLogger;

    @Test
    public void testRates() throws Exception {
        testTree = new TreeParser();
        testTree.initByName("newick", newickTree, "IsLabelledNewick", true);

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initialise();

        treeLengthLogger = new TreeLengthLogger();
        treeLengthLogger.initByName("tree", testTree);
        final double computedLength = treeLengthLogger.getArrayValue();

        assertEquals(expectedLength, computedLength, allowedError);
    }
}
