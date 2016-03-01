package snetworktests;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.util.TreeParser;
import speciesnetwork.NetworkParser;

public class NetworkParserTest {
    TreeParser treeParser;
    NetworkParser networkParser;

    // "Partial" reticulation nodes (that is, #H nodes without any directly attached children in the extended newick string)
    // should be on the right. "Complete" reticulation nodes (with directly attached children) should be on the left.
    // cannot parse if there is :: inheritProb in the string
    final String testNetwork2 = "((A:0.2,#H1:0.1)#S1:0.3,((B:0.1)#H1:0.2,C:0.3)#S2:0.2)#R";
    // final String testNetwork2 = "(((B:0.1)#H1:0.1,A:0.2)#S1:0.3,(C:0.3,#H1:0.2)#S2:0.2)#R";

    final String testNetwork3 = "((((A:0.1)#H1:0.1)#H2:0.3,#H4:0.1)#S2:0.1,((((#H1:0.1)#H3:0.1,#H2:0.1)#S1:0.1)#H4:0.1,#H3:0.3)#S3:0.1)#R";

    public NetworkParserTest() {
        treeParser = new TreeParser();
        networkParser = new NetworkParser();

        try {
            treeParser.initByName("newick", testNetwork2, "IsLabelledNewick", true, "adjustTipHeights", false);
            networkParser.initByName("tree", treeParser);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testLogHR() throws Exception {
        // System.out.println(networkParser.toString());
        assertEquals(testNetwork2, networkParser.toString());
    }
}
