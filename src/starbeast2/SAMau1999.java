package starbeast2;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * @author Huw Ogilvie
 */

@Description("Tree operator which randomly changes the height of a node, " +
        "then reconstructs the tree from node heights.")
public class SAMau1999 extends Operator {
    public final Input<Tree> treeInput = new Input<>("tree", "the tree", Validate.REQUIRED);
    public final Input<Double> windowInput = new Input<>("window", "size of the random walk window", 10.0);
    public final Input<RealParameter> originInput = new Input<RealParameter>("origin", "The time when the process started", (RealParameter) null);

    private int nextIndex;
    private int nodeCount;
    private int trueBifurcationCount;
    private Node[] canonicalOrder;
    private int[] trueBifurcations;
    private double[] nodeHeights;
    private boolean superimposedAncestors;
    private double window;
    private boolean originSpecified;

    @Override
    public void initAndValidate() {
        final Tree tree = treeInput.get();
        nodeCount = tree.getNodeCount();
        canonicalOrder = new Node[nodeCount];
        trueBifurcations = new int[nodeCount];
        nodeHeights = new double[nodeCount];
        window = windowInput.get();
        originSpecified = originInput.get() != null;
    }

    /* This improves the proposal suggested by Mau (1999). It has been made compatible with sampled ancestors by
    enforcing a minimum height and disallowing superimposed sampled ancestor nodes. */
    @Override
    public double proposal() {
        final Tree tree = treeInput.get();
        final Node originalRoot = tree.getRoot();

        double maxHeight;
        if (originSpecified) {
            maxHeight = originInput.get().getValue();
        } else {
            maxHeight = Double.POSITIVE_INFINITY;
        }

        // chooseCanonicalOrder also fills in nodeHeights and trueBifurcations
        // the lastIndex will be the last and right-most node index
        trueBifurcationCount = 0;
        nextIndex = 0;
        chooseCanonicalOrder(originalRoot);

        // no nodes can be changed by this operator
        if (trueBifurcationCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // pick a bifurcation at random and change the height
        final int chosenNode = trueBifurcations[Randomizer.nextInt(trueBifurcationCount)];
        // as long as the height is above the tips (or sampled ancestors) immediately either side in the canonical
        // order, the species tree seems buildable from the new times
        final double minHeight = Double.max(nodeHeights[chosenNode - 1], nodeHeights[chosenNode + 1]);

        final double heightDelta = window * (Randomizer.nextDouble() - 0.5);

        // use reflection to avoid invalid heights
        double newHeight = nodeHeights[chosenNode] + heightDelta;
        while (newHeight < minHeight || newHeight > maxHeight) {
            if (newHeight < minHeight) {
                newHeight = minHeight + minHeight - newHeight;
            }
            if (newHeight > maxHeight) {
                newHeight = maxHeight + maxHeight - newHeight;
            }
        }

        nodeHeights[chosenNode] = newHeight;
        canonicalOrder[chosenNode].setHeight(newHeight);

        superimposedAncestors = false;
        final int rootIndex = rebuildTree(0, nodeCount - 1);
        if (superimposedAncestors) {
            return Double.NEGATIVE_INFINITY;
        }

        final Node newRoot = canonicalOrder[rootIndex];
        if (newRoot != originalRoot) {
            newRoot.setParent(null);
            tree.setRoot(newRoot);
        }

        return 0.0;
    }

    /* Performs an in-order traversal of the species tree, randomly shuffling left and right nodes, to produce
       a canonical order in the sense of Mau et al 1999. Also identify which nodes are true bifurcations
       (not fake nodes used for sampled ancestors) */
    private double chooseCanonicalOrder(final Node node) {
        Node canonicalLeft;
        Node canonicalRight;

        if (Randomizer.nextBoolean()) {
            canonicalLeft = node.getLeft();
            canonicalRight = node.getRight();
        } else {
            canonicalLeft = node.getRight();
            canonicalRight = node.getLeft();
        }

        double leftChildHeight;
        if (canonicalLeft.isLeaf()) {
            final int leftChildIndex = nextIndex;
            nextIndex++;

            canonicalOrder[leftChildIndex] = canonicalLeft;

            leftChildHeight = canonicalLeft.getHeight();
            nodeHeights[leftChildIndex] = leftChildHeight;
        } else {
            leftChildHeight = chooseCanonicalOrder(canonicalLeft);
        }

        final int thisIndex = nextIndex;
        nextIndex++;

        canonicalOrder[thisIndex] = node;

        final double thisHeight = node.getHeight();
        nodeHeights[thisIndex] = thisHeight;

        double rightChildHeight;
        if (canonicalRight.isLeaf()) {
            final int rightChildIndex = nextIndex;
            nextIndex++;

            canonicalOrder[rightChildIndex] = canonicalRight;

            rightChildHeight = canonicalRight.getHeight();
            nodeHeights[rightChildIndex] = rightChildHeight;
        } else {
            rightChildHeight = chooseCanonicalOrder(canonicalRight);
        }

        if (thisHeight > leftChildHeight && thisHeight > rightChildHeight) {
            trueBifurcations[trueBifurcationCount] = thisIndex;
            trueBifurcationCount++;
        }

        return thisHeight;
    }

    /* from and to are inclusive */
    private int rebuildTree(final int from, final int to) {
        double thisHeight = 0.0;
        int nodeIndex = -1;

        /* Only check internal nodes, which are odd numbered (leaves are even numbered). Reject move if multiple
           internal nodes in the range have the same height, as they are likely fake bifurcations, and connecting
           them will result in multiple sampled ancestors at the same point in time along the same lineage. This
           matches the following behaviour of LeafToSampledAncestorJump (see lines 68-70):
           if (getOtherChild(parent, leaf).getHeight() >= leaf.getHeight()) return Double.NEGATIVE_INFINITY; */
        for (int i = from + 1; i < to; i = i + 2) {
            if (nodeHeights[i] > thisHeight) {
                thisHeight = nodeHeights[i];
                nodeIndex = i;
            } else if (nodeHeights[i] == thisHeight) {
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
}
