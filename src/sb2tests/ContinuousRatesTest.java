package sb2tests;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import starbeast2.ContinuousRates;

public class ContinuousRatesTest {
    private String newickTree = "((((a1:0.3,a2:0.3):1.6,(b1:1.8,b2:1.8):0.1):0.5,c1:2.4):0.6,c2:3.0)";
    private TreeParser testTree;

    private final double meanRate = 1.5;
    private final double stdev = 0.2;
    private final double initialBranchRate = 1.0;

    final double allowedError = 10e-6;

    private RealParameter meanRateParameter;
    private RealParameter branchRatesParameter;
    private RealParameter branchRateSpread;

    private ContinuousRates clockModel;

    @Test
    public void testRates() throws Exception {
        initializeTree();

        meanRateParameter = new RealParameter();
        branchRatesParameter = new RealParameter();
        branchRateSpread = new RealParameter();

        meanRateParameter.initByName("value", String.valueOf(meanRate));
        branchRatesParameter.initByName("value", String.valueOf(initialBranchRate));
        branchRateSpread.initByName("value", String.valueOf(stdev));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", meanRateParameter);
        state.initByName("stateNode", branchRatesParameter);
        state.initByName("stateNode", branchRateSpread);
        state.initialise();

        clockModel = new ContinuousRates();
        clockModel.initByName("tree", testTree, "logRates", branchRatesParameter, "stdev", branchRateSpread, "estimateRoot", false, "clock.rate", meanRateParameter);

        initializeRates();
        checkRates();
    }

    private void initializeRates() {
        String[] node00 = {"a1"};
        String[] node01 = {"a2"};
        String[] node02 = {"b1"};
        String[] node03 = {"b2"};
        String[] node04 = {"c1"};
        String[] node05 = {"c2"};
        String[] node06 = {"a1", "a2"};
        String[] node07 = {"b1", "b2"};
        String[] node08 = {"a1", "a2", "b1", "b2"};
        String[] node09 = {"a1", "a2", "b1", "b2", "c1"};
        
        double rate00 = -2.1309141;
        double rate01 =  0.4201902;
        double rate02 =  1.1237991;
        double rate03 = -1.8328272;
        double rate04 =  0.8646902;
        double rate05 =  0.2305576;
        double rate06 =  0.9259551;
        double rate07 = -0.5233576;
        double rate08 = -1.0421507;
        double rate09 =  0.1124448;

        assertTrue(setRate(rate00, testTree, node00));
        assertTrue(setRate(rate01, testTree, node01));
        assertTrue(setRate(rate02, testTree, node02));
        assertTrue(setRate(rate03, testTree, node03));
        assertTrue(setRate(rate04, testTree, node04));
        assertTrue(setRate(rate05, testTree, node05));
        assertTrue(setRate(rate06, testTree, node06));
        assertTrue(setRate(rate07, testTree, node07));
        assertTrue(setRate(rate08, testTree, node08));
        assertTrue(setRate(rate09, testTree, node09));
    }

    private void checkRates() {
        String[] node00 = {"a1"};
        String[] node01 = {"a2"};
        String[] node02 = {"b1"};
        String[] node03 = {"b2"};
        String[] node04 = {"c1"};
        String[] node05 = {"c2"};
        String[] node06 = {"a1", "a2"};
        String[] node07 = {"b1", "b2"};
        String[] node08 = {"a1", "a2", "b1", "b2"};
        String[] node09 = {"a1", "a2", "b1", "b2", "c1"};
        String[] node10 = {"a1", "a2", "b1", "b2", "c1", "c2"};

        // correct tree rates were calculated by hand
        double rate00 = 0.9601001;
        double rate01 = 1.5991994;
        double rate02 = 1.8408454;
        double rate03 = 1.0190794;
        double rate04 = 1.7478792;
        double rate05 = 1.5396831;
        double rate06 = 1.7694276;
        double rate07 = 1.3241802;
        double rate08 = 1.1936728;
        double rate09 = 1.5037381;
        double rate10 = 1.5;

        assertEquals(rate00, getRate(testTree, node00), allowedError);
        assertEquals(rate01, getRate(testTree, node01), allowedError);
        assertEquals(rate02, getRate(testTree, node02), allowedError);
        assertEquals(rate03, getRate(testTree, node03), allowedError);
        assertEquals(rate04, getRate(testTree, node04), allowedError);
        assertEquals(rate05, getRate(testTree, node05), allowedError);
        assertEquals(rate06, getRate(testTree, node06), allowedError);
        assertEquals(rate07, getRate(testTree, node07), allowedError);
        assertEquals(rate08, getRate(testTree, node08), allowedError);
        assertEquals(rate09, getRate(testTree, node09), allowedError);
        assertEquals(rate10, getRate(testTree, node10), allowedError);
    }

    private boolean setRate(double rate, TreeParser tree, String[] target) {
        final Node targetNode = findNode(tree, target);
        if (targetNode == null) {
            return false;
        } else {
            branchRatesParameter.setValue(targetNode.getNr(), rate);
            return true;
        }
    }

    private double getRate(TreeParser tree, String[] target) {
        final Node targetNode = findNode(tree, target);
        return clockModel.getRateForBranch(targetNode);
    }

    private Node findNode(final TreeParser tree, final String[] targetArray) {
        final Node[] treeNodes = tree.getNodesAsArray();
        final Set<String> targetSet = new HashSet<>();
        for (int i = 0; i < targetArray.length; i++) {
            targetSet.add(targetArray[i]);
        }

        for (Node node: treeNodes) {
            Set<String> nodeSet = new HashSet<>();

            if (node.isLeaf()) {
                nodeSet.add(node.getID());
            } else {
                final List<Node> leafNodes = node.getAllLeafNodes();
                for (Node leaf: leafNodes) {
                    nodeSet.add(leaf.getID());
                }
            }

            if (targetSet.equals(nodeSet)) {
                return node;
            }
        }

        return null;
    }

    public void initializeTree() throws Exception {
        testTree = new TreeParser();
        testTree.initByName("newick", newickTree, "IsLabelledNewick", true);
    }
}
