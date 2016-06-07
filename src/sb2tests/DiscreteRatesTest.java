package sb2tests;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.math.distributions.LogNormalDistributionModel;
import beast.util.TreeParser;
import starbeast2.UncorrelatedRates;

public class DiscreteRatesTest {
    private String newickTree = "((((a1:0.3,a2:0.3):1.6,(b1:1.8,b2:1.8):0.1):0.5,c1:2.4):0.6,c2:3.0)";
    private TreeParser testTree;

    private final double meanRate = 1.5;
    private final int initialBranchRate = 1;

    final double allowedError = 10e-6;

    private RealParameter meanRateParameter;
    private IntegerParameter branchRatesParameter;

    private UncorrelatedRates clockModel;
    private LogNormalDistributionModel branchRateDistribution;

    @Test
    public void testRates() throws Exception {
        initializeTree();

        meanRateParameter = new RealParameter();
        branchRatesParameter = new IntegerParameter();

        meanRateParameter.initByName("value", String.valueOf(meanRate));
        branchRatesParameter.initByName("value", String.valueOf(initialBranchRate));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", meanRateParameter);
        state.initByName("stateNode", branchRatesParameter);
        state.initialise();

        branchRateDistribution = new LogNormalDistributionModel();
        branchRateDistribution.initByName("M", "1.0", "S", "1.0", "meanInRealSpace", true);

        clockModel = new UncorrelatedRates();
        clockModel.initByName("tree", testTree, "rates", branchRatesParameter, "distr", branchRateDistribution, "estimateRoot", false, "clock.rate", meanRateParameter);

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

        int rate00 = 0;
        int rate01 = 4;
        int rate02 = 8;
        int rate03 = 5;
        int rate04 = 3;
        int rate05 = 8;
        int rate06 = 9;
        int rate07 = 9;
        int rate08 = 3;
        int rate09 = 1;

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
        double rate00 = 0.1756278;
        double rate01 = 0.8023613;
        double rate02 = 2.5648462;
        double rate03 = 1.0316159;
        double rate04 = 0.6188729;
        double rate05 = 2.5648462;
        double rate06 = 4.7129721;
        double rate07 = 4.7129721;
        double rate08 = 0.6188729;
        double rate09 = 0.3227206;
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

    private boolean setRate(int rate, TreeParser tree, String[] target) {
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
