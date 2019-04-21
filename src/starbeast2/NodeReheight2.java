package starbeast2;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

@Description("Tree operator which randomly changes the height of a node, " +
        "then reconstructs the tree from node heights.")
public class NodeReheight2 extends Operator {
    public final Input<SpeciesTree> treeInput = new Input<>("tree", "the species tree", Validate.REQUIRED);
    public final Input<TaxonSet> taxonSetInput = new Input<>("taxonset", "taxon set describing species tree taxa and their gene trees", Validate.REQUIRED);
    public final Input<List<GeneTree>> geneTreesInput = new Input<>("geneTree", "list of gene trees that constrain species tree movement", new ArrayList<>());

    final private double window = 1.0;

    private int lastIndex;
    private int trueBifurcationCount;
    private Node[] canonicalOrder;
    private int[] trueBifurcations;
    private double[] nodeHeights;

    @Override
    public void initAndValidate() {
        final int nodeCount = treeInput.get().getNodeCount();
        canonicalOrder = new Node[nodeCount];
        trueBifurcations = new int[nodeCount];
        nodeHeights = new double[nodeCount];
    }

    @Override
    public double proposal() {
        final Node rootNode = treeInput.get().getRoot();

        // chooseCanonicalOrder also fills in nodeHeights and trueBifurcations
        // the lastIndex will be the last and right-most node index
        lastIndex = -1;
        trueBifurcationCount = 0;
        chooseCanonicalOrder(rootNode);

        // no nodes can be changed by this operator
        if (trueBifurcationCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // pick a bifurcation at random and change the height
        final int chosenNode = trueBifurcations[Randomizer.nextInt(trueBifurcationCount)];
        final double minHeight = Double.min(nodeHeights[chosenNode - 1], nodeHeights[chosenNode + 1]);

        double newHeight = window * (Randomizer.nextDouble() - 0.5);
        if (newHeight < minHeight) { // reflection
            newHeight = minHeight + minHeight - newHeight;
        }

        nodeHeights[chosenNode] = newHeight;
        canonicalOrder[chosenNode].setHeight(newHeight);

        rebuildTree(0, lastIndex);

        int[] visitedCounts = new int[lastIndex + 1];
        checkConsistency(rootNode, visitedCounts);
        for (int i = 0; i <= lastIndex; i++) {
            assert visitedCounts[i] == 1;
        }

        return 0.0;
    }

    /* Performs an in-order traversal of the species tree, randomly shuffling left and right nodes, to produce
       a canonical order in the sense of Mau et al 1999. Also identify which nodes are true bifurcations
       (not fake nodes used for sampled ancestors) */
    private double chooseCanonicalOrder(final Node node) {
        final double height = node.getHeight();

        if (node.isLeaf()) {
            lastIndex++;

            canonicalOrder[lastIndex] = node;
            nodeHeights[lastIndex] = height;

            return height;
        }

        Node canonicalLeft;
        Node canonicalRight;

        if (Randomizer.nextBoolean()) {
            canonicalLeft = node.getLeft();
            canonicalRight = node.getRight();
        } else {
            canonicalLeft = node.getRight();
            canonicalRight = node.getLeft();
        }

        final double leftChildHeight = chooseCanonicalOrder(canonicalLeft);

        lastIndex++;
        final int thisIndex = lastIndex;

        canonicalOrder[thisIndex] = node;

        final double rightChildHeight = chooseCanonicalOrder(canonicalRight);

        if (height != leftChildHeight && height != rightChildHeight) {
            trueBifurcations[trueBifurcationCount] = thisIndex;
            trueBifurcationCount++;
        }

        return height;
    }

    /* from and to are inclusive */
    private int rebuildTree(final int from, final int to) {
        double height = 0.0;
        int nodeIndex = -1;

        // only check internal nodes, which are odd numbered (leaves are even numbered)
        for (int i = from + 1; i < to; i = i + 2) {
            if (nodeHeights[i] > height) {
                height = nodeHeights[i];
                nodeIndex = i;
            }
        }

        int leftNodeIndex;
        if (from == nodeIndex - 1) {
            leftNodeIndex = from;
        } else {
            leftNodeIndex = rebuildTree(from, nodeIndex - 1);
        }

        canonicalOrder[leftNodeIndex].setParent(canonicalOrder[nodeIndex]);
        canonicalOrder[nodeIndex].setLeft(canonicalOrder[leftNodeIndex]);

        int rightNodeIndex;
        if (nodeIndex + 1 == to) {
            rightNodeIndex = to;
        } else {
            rightNodeIndex = rebuildTree(nodeIndex + 1, to);
        }

        canonicalOrder[rightNodeIndex].setParent(canonicalOrder[nodeIndex]);
        canonicalOrder[nodeIndex].setRight(canonicalOrder[rightNodeIndex]);

        return nodeIndex;
    }

    private void checkConsistency(Node node, int[] visitedCounts) {
        visitedCounts[node.getNr()]++;
        final double height = node.getHeight();
        final List<Node> children = node.getChildren();
        if (children.size() != 0) {
            assert children.size() == 2;
            Node leftChild = children.get(0);
            Node rightChild = children.get(1);
            assert leftChild.getParent() == node;
            assert rightChild.getParent() == node;
            double leftHeight = leftChild.getHeight();
            double rightHeight = rightChild.getHeight();
            assert leftHeight <= height;
            assert rightHeight <= height;
            if (leftHeight == rightHeight) {
                assert leftHeight != height;
            }
            checkConsistency(leftChild, visitedCounts);
            checkConsistency(rightChild, visitedCounts);
        }
    }
}
