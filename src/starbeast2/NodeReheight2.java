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

/**
 * @author Huw Ogilvie
 */

@Description("Tree operator which randomly changes the height of a node, " +
        "then reconstructs the tree from node heights.")
public class NodeReheight2 extends Operator {
    public final Input<SpeciesTree> treeInput = new Input<>("tree", "the species tree", Validate.REQUIRED);
    public final Input<TaxonSet> taxonSetInput = new Input<>("taxonset", "taxon set describing species tree taxa and their gene trees", Validate.REQUIRED);
    public final Input<List<GeneTree>> geneTreesInput = new Input<>("geneTree", "list of gene trees that constrain species tree movement", new ArrayList<>());
    public final Input<Double> windowInput = new Input<>("window", "size of the random walk window", 10.0);

    private enum RelativePosition {LEFT, RIGHT, BOTH};

    private int nextIndex;
    private int nodeCount;
    private int geneTreeCount;
    private int[][] leafNodeMaps;
    private RelativePosition[][] leafPositionArrays;
    private int trueBifurcationCount;
    private Node[] canonicalOrder;
    private int[] canonicalMap;
    private int[] trueBifurcations;
    private double[] nodeHeights;
    private boolean superimposedAncestors;
    private double maxHeight;
    private double window;

    @Override
    public void initAndValidate() {
        final SpeciesTree speciesTree = treeInput.get();
        nodeCount = speciesTree.getNodeCount();
        canonicalOrder = new Node[nodeCount];
        canonicalMap = new int[nodeCount];
        trueBifurcations = new int[nodeCount];
        nodeHeights = new double[nodeCount];
        window = windowInput.get();

        final List<GeneTree> geneTrees = geneTreesInput.get();
        geneTreeCount = geneTrees.size();
        leafNodeMaps = new int[geneTreeCount][];
        leafPositionArrays = new RelativePosition[geneTreeCount][];
        for (int i = 0; i < geneTreeCount; i++) {
            leafNodeMaps[i] = geneTrees.get(i).getTipNumberMap();
            leafPositionArrays[i] = new RelativePosition[leafNodeMaps[i].length];
        }
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
        // as long as the height is above the tips (or sampled ancestors) immediately either side in the canonical
        // order, the species tree seems buildable from the new times
        final double minHeight = Double.max(nodeHeights[chosenNode - 1], nodeHeights[chosenNode + 1]);

        recalculateMaxHeight(chosenNode);

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

        assert checkVisitedCounts(tree);

        return 0.0;
    }

    private void recalculateMaxHeight(final int centerIndex) {
        maxHeight = Double.POSITIVE_INFINITY;

        final List<GeneTree> geneTrees = geneTreesInput.get();
        for (int i = 0; i < geneTreeCount; i++) {
            final int[] leafNodeMap = leafNodeMaps[i];
            final RelativePosition[] leafPositions = leafPositionArrays[i];
            for (int j = 0; j < leafNodeMap.length; j++) {
                final int speciesNodeNumber = leafNodeMap[j];
                final int speciesIndex = canonicalMap[speciesNodeNumber];
                if (speciesIndex < centerIndex) {
                    leafPositions[j] = RelativePosition.LEFT;
                } else {
                    leafPositions[j] = RelativePosition.RIGHT;
                }
            }

            final Node geneTreeRoot = geneTrees.get(i).getRoot();
            recurseMaxHeight(geneTreeRoot, leafPositions);
        }
    }

    private RelativePosition recurseMaxHeight(final Node node, final RelativePosition[] leafPositions) {
        final Node leftChild = node.getLeft();
        final Node rightChild = node.getRight();

        RelativePosition leftDescendantPosition;
        if (leftChild.isLeaf()) {
            leftDescendantPosition = leafPositions[leftChild.getNr()];
        } else {
            leftDescendantPosition = recurseMaxHeight(leftChild, leafPositions);
        }

        RelativePosition rightDescendantPosition;
        if (rightChild.isLeaf()) {
            rightDescendantPosition = leafPositions[rightChild.getNr()];
        } else {
            rightDescendantPosition = recurseMaxHeight(rightChild, leafPositions);
        }

        if (leftDescendantPosition == rightDescendantPosition) {
            return leftDescendantPosition;
        } else {
            // if all descendants of one child are on the left, and all descendants of the other child are on the right
            if (leftDescendantPosition != RelativePosition.BOTH && rightDescendantPosition != RelativePosition.BOTH) {
                maxHeight = Double.min(maxHeight, node.getHeight());
            }
            return RelativePosition.BOTH;
        }
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

            canonicalMap[canonicalLeft.getNr()] = leftChildIndex;
            canonicalOrder[leftChildIndex] = canonicalLeft;

            leftChildHeight = canonicalLeft.getHeight();
            nodeHeights[leftChildIndex] = leftChildHeight;
        } else {
            leftChildHeight = chooseCanonicalOrder(canonicalLeft);
        }

        final int thisIndex = nextIndex;
        nextIndex++;

        canonicalMap[node.getNr()] = thisIndex;
        canonicalOrder[thisIndex] = node;

        final double thisHeight = node.getHeight();
        nodeHeights[thisIndex] = thisHeight;

        double rightChildHeight;
        if (canonicalRight.isLeaf()) {
            final int rightChildIndex = nextIndex;
            nextIndex++;

            canonicalMap[canonicalRight.getNr()] = rightChildIndex;
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

    // for debugging, only called when assertions are enabled
    private boolean checkVisitedCounts(SpeciesTree tree) {
        int[] visitedCounts = new int[nodeCount];
        recurseVisitedCounts(tree.getRoot(), visitedCounts);
        for (int i = 0; i < nodeCount; i++) {
            if (visitedCounts[i] != 1) {
                return false;
            }
        }
        return true;
    }

    // for debugging, only called when assertions are enabled
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
