package starbeast2;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;


/**
 * @author Huw Ogilvie
 */

@Description("Tree operator which randomly changes the height of a node, " +
        "then reconstructs the tree from node heights.")
public class SAMau1999 extends Operator {
    public final Input<SpeciesTree> treeInput = new Input<>("tree", "the species tree", Validate.REQUIRED);
    public final Input<Double> windowInput = new Input<>("window", "size of the random walk window", 10.0);
    public final Input<RealParameter> originInput = new Input<RealParameter>("origin", "The time when the process started", (RealParameter) null);

    private int nextIndex;
    private int nodeCount;
    private int trueBifurcationCount;
    private Node[] canonicalOrder;
    private int[] trueBifurcations;
    private double[] nodeHeights;
    private Node[] leftChildren;
    private Node[] rightChildren;
    private Node[] parents;
    private boolean superimposedAncestors;
    private double window;
    private boolean originSpecified;

    @Override
    public void initAndValidate() {
        final SpeciesTree speciesTree = treeInput.get();
        nodeCount = speciesTree.getNodeCount();
        canonicalOrder = new Node[nodeCount];
        trueBifurcations = new int[nodeCount];
        nodeHeights = new double[nodeCount];
        leftChildren = new Node[nodeCount];
        rightChildren = new Node[nodeCount];
        parents = new Node[nodeCount];
        window = windowInput.get();
        originSpecified = originInput.get() != null;
    }

    /* This proposal improves TREE SLIDE, developed by Joseph Heled. See section 3.4.1 of Heled's 2011 PhD thesis
    "Bayesian Computational Inference of Species Trees and Population Sizes". TREE SLIDE was developed for ultrametric
    binary species trees, this proposal has been made compatible with sampled ancestors by enforcing a minimum height
    and disallowing superimposed sampled ancestor nodes. Also uses a random walk window with reflection in order to
    sample the heights of nodes without maximum height constraints. */
    @Override
    public double proposal() {
        final SpeciesTree tree = treeInput.get();
        final Node originalRoot = tree.getRoot();

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
        final double originalHeight = nodeHeights[chosenNode];

        // as long as the height is above the tips (or sampled ancestors) immediately either side in the canonical
        // order, the species tree seems buildable from the new times (except in the case of superimposed
        // sampled ancestors, as discussed below)
        final double minHeight = Double.max(nodeHeights[chosenNode - 1], nodeHeights[chosenNode + 1]);

        double maxHeight;
        if (originSpecified) {
            maxHeight = originInput.get().getValue();
        } else {
            maxHeight = Double.POSITIVE_INFINITY;
        }

        // Use reflection to avoid invalid heights. Height returns to original position every 2 * (max - min) units,
        // so modulus is used to avoid unnecessary looping if the difference between window size and the tree scale
        // is extreme.
        final double heightDelta = (window * (Randomizer.nextDouble() - 0.5)) % (2.0 * (maxHeight - minHeight));
        double newHeight = originalHeight + heightDelta;
        while (newHeight < minHeight || newHeight > maxHeight) {
            if (newHeight < minHeight) {
                newHeight = minHeight + minHeight - newHeight;
            }
            if (newHeight > maxHeight) {
                newHeight = maxHeight + maxHeight - newHeight;
            }
        }

        nodeHeights[chosenNode] = newHeight;

        superimposedAncestors = false;
        final int rootIndex = rebuildTree(0, nodeCount - 1);
        parents[rootIndex] = null;
        if (superimposedAncestors) {
            return Double.NEGATIVE_INFINITY;
        }

        // wait until after checking for superimposed ancestors before modifying tree
        canonicalOrder[chosenNode].setHeight(newHeight);

        for (int i = 0; i < nodeCount; i++) {
            canonicalOrder[i].setParent(parents[i]);

            if (i % 2 == 1) { // internal node
                canonicalOrder[i].setLeft(leftChildren[i]);
                canonicalOrder[i].setRight(rightChildren[i]);
            }
        }

        final Node newRoot = canonicalOrder[rootIndex];
        // for some reason if the root is not reset - even if the root node is the same node as before! - the
        // morphological likelihood will be radically wrong (idk why)
        tree.setRoot(newRoot);

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

        /* Only check internal nodes, which are odd numbered (leaves are even numbered). If there are multiple highest
           internal nodes in the range, they are likely fake bifurcations, and connecting
           them will result in multiple sampled ancestors at the same point in time along the same lineage.
           In this case we repeat changing the height of the chosen node until this no longer occurs.
           This is similar to the following behaviour of LeafToSampledAncestorJump (see lines 68-70):
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

        parents[leftNodeIndex] = canonicalOrder[nodeIndex];
        leftChildren[nodeIndex] = canonicalOrder[leftNodeIndex];

        int rightNodeIndex;
        if (nodeIndex + 1 == to) {
            rightNodeIndex = to;
        } else {
            rightNodeIndex = rebuildTree(nodeIndex + 1, to);
        }

        parents[rightNodeIndex] = canonicalOrder[nodeIndex];
        rightChildren[nodeIndex] = canonicalOrder[rightNodeIndex];

        return nodeIndex;
    }
}
