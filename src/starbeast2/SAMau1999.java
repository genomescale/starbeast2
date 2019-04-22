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
public class SAMau1999 extends Operator {
    public final Input<SpeciesTree> treeInput = new Input<>("tree", "the species tree", Validate.REQUIRED);
    public final Input<TaxonSet> taxonSetInput = new Input<>("taxonset", "taxon set describing species tree taxa and their gene trees", Validate.REQUIRED);
    public final Input<List<GeneTree>> geneTreesInput = new Input<>("geneTree", "list of gene trees that constrain species tree movement", new ArrayList<>());

    private int lastIndex;
    private int trueBifurcationCount;
    private Node[] canonicalOrder;
    private int[] trueBifurcations;
    private double[] nodeHeights;
    private boolean superimposedAncestors;

    @Override
    public void initAndValidate() {
        final int nodeCount = treeInput.get().getNodeCount();
        canonicalOrder = new Node[nodeCount];
        trueBifurcations = new int[nodeCount];
        nodeHeights = new double[nodeCount];
    }

    @Override
    public double proposal() {
        final SpeciesTree tree = treeInput.get();
        final Node originalRoot = tree.getRoot();

        // chooseCanonicalOrder also fills in nodeHeights and trueBifurcations
        // the lastIndex will be the last and right-most node index
        lastIndex = -1;
        trueBifurcationCount = 0;
        chooseCanonicalOrder(originalRoot);

        // no nodes can be changed by this operator
        if (trueBifurcationCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // pick a bifurcation at random and change the height
        final int chosenNode = trueBifurcations[Randomizer.nextInt(trueBifurcationCount)];
        final double minHeight = Double.max(nodeHeights[chosenNode - 1], nodeHeights[chosenNode + 1]);

        double newHeight = nodeHeights[chosenNode] + Randomizer.nextDouble() - 0.5;
        if (newHeight < minHeight) { // reflection
            newHeight = minHeight + minHeight - newHeight;
        }

        nodeHeights[chosenNode] = newHeight;
        canonicalOrder[chosenNode].setHeight(newHeight);

        superimposedAncestors = false;
        final int rootIndex = rebuildTree(0, lastIndex);
        if (superimposedAncestors) {
            return Double.NEGATIVE_INFINITY;
        }

        final Node newRoot = canonicalOrder[rootIndex];
        if (newRoot != originalRoot) {
            newRoot.setParent(null);
            tree.setRoot(newRoot);
        }

        assert checkVisitedCounts(tree);

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
        nodeHeights[thisIndex] = height;

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

        /* Only check internal nodes, which are odd numbered (leaves are even numbered). Reject move if multiple
           internal nodes in the range have the same height, as they are likely fake bifurcations, and connecting
           them will result in multiple sampled ancestors at the same point in time along the same lineage. This
           matches the following behaviour of LeafToSampledAncestorJump (see lines 68-70):
           if (getOtherChild(parent, leaf).getHeight() >= leaf.getHeight()) return Double.NEGATIVE_INFINITY; */
        for (int i = from + 1; i < to; i = i + 2) {
            if (nodeHeights[i] > height) {
                height = nodeHeights[i];
                nodeIndex = i;
            } else if (nodeHeights[i] == height) {
                superimposedAncestors = true;
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

    private boolean checkVisitedCounts(SpeciesTree tree) {
        int[] visitedCounts = new int[lastIndex + 1];
        recurseVisitedCounts(tree.getRoot(), visitedCounts);
        for (int i = 0; i <= lastIndex; i++) {
            if (visitedCounts[i] != 1) {
                return false;
            }
        }
        return true;
    }

    private void recurseVisitedCounts(Node node, int[] visitedCounts) {
        visitedCounts[node.getNr()]++;
        final List<Node> children = node.getChildren();
        if (!node.isLeaf()) {
            assert children.size() == 2;
            final Node leftChild = children.get(0);
            final Node rightChild = children.get(1);
            assert leftChild.getParent() == node;
            assert rightChild.getParent() == node;
            assert leftChild.getHeight() <= node.getHeight();
            assert rightChild.getHeight() <= node.getHeight();
            recurseVisitedCounts(leftChild, visitedCounts);
            recurseVisitedCounts(rightChild, visitedCounts);
        }
    }
}
