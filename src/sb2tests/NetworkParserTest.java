package sb2tests;

import org.junit.Test;

import beast.util.TreeParser;
import starbeast2.NetworkParser;

public class NetworkParserTest {
    TreeParser treeParser;
    NetworkParser networkParser;
    
    final String testNetwork = "((1:1.4,((2:0.7,(3:0.2,(4:0.1)Y#H1:0.1)g:0.5)e:0.3,(((Y#H1:0.2,5:0.3)h:0.5,6:0.8)f:0.1)X#H2:0.1)c:0.4)a:0.5,((X#H2:0.1,7:1.0)d:0.6,8:1.6)b:0.3)r";
    public NetworkParserTest() {
        treeParser = new TreeParser();
        networkParser = new NetworkParser();

        try {
            treeParser.initByName("newick", testNetwork, "IsLabelledNewick", true, "adjustTipHeights", false);
            networkParser.initByName("tree", treeParser);
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
    }

    @Test
    public void testLogHR() throws Exception {
        System.out.println(networkParser.toString());
    }
}
