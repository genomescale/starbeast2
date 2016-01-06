// Based on Exchange.java in BEAST2, which in turn was based on ExchangeOperator.java in BEAST.

package starbeast2;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.SetMultimap;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

/**
* @author Remco Bouckaert
* @author Alexei Drummond
* @author Huw Ogilvie
* @author Andrew Rambaut
 */

@Description("Implements the co-ordinated species and gene tree operator described in Yang & Rannala 2015. "
        + "This performs a narrow exchange operation, then prunes-and-regrafts gene tree nodes made invalid "
        + "by the operation onto a valid contemporary branch (without changing node heights). "
        + "See http://doi.org/10.1093/molbev/msu279 for full details.")
public class CoordinatedExchange extends CoordinatedOperator {
    private Node brother;
    private Node parent;
    private Node uncle;

    private List<SortedMap<Node, Node>> forwardMovedNodes;
    private SetMultimap<Integer, Node> forwardGraftNodes;

    final static Comparator<Node> nhc = new NodeHeightComparator().reversed();

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        double fLogHastingsRatio = narrow();

        // only rearrange gene trees if the species tree has changed
        if (fLogHastingsRatio != Double.NEGATIVE_INFINITY) {
            fLogHastingsRatio += rearrangeGeneTrees();
        }

        return fLogHastingsRatio;
    }

    private int isg(final Node n) {
      return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    private int sisg(final Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }

    public double narrow() {
        final TreeInterface speciesTree = speciesTreeInput.get().getTree();
        final int nInternalNodes = speciesTree.getInternalNodeCount();
        if (nInternalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        Node grandparent = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        while (grandparent.getLeft().isLeaf() && grandparent.getRight().isLeaf()) {
            grandparent = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        }

        parent = grandparent.getLeft();
        uncle = grandparent.getRight();
        if (parent.getHeight() < uncle.getHeight()) {
            parent = grandparent.getRight();
            uncle = grandparent.getLeft();
        }

        if (parent.isLeaf()) {
            return Double.NEGATIVE_INFINITY;
        }

        brother = Randomizer.nextBoolean() ? parent.getLeft() : parent.getRight();

        // must be done before any changes are made to the gene or species trees
        forwardMovedNodes = getMovedPairs(brother);
        forwardGraftNodes = getGraftBranches(uncle);

        int validGP = 0;
        {
            for(int i = nInternalNodes + 1; i < 1 + 2*nInternalNodes; ++i) {
                validGP += isg(speciesTree.getNode(i));
            }
        }

        final int c2 = sisg(parent) + sisg(uncle);

        exchangeNodes(parent, grandparent, brother, uncle);

        final int validGPafter = validGP - c2 + sisg(parent) + sisg(uncle);

        return Math.log((float)validGP/validGPafter);
    }

    // for testing purposes
    public void manipulateSpeciesTree(Node argBrother) {
        brother = argBrother;
        parent = brother.getParent();

        final Node grandparent = parent.getParent();
        if (grandparent.getLeft().getNr() == parent.getNr()) {
            uncle = grandparent.getRight();
        } else {
            uncle = grandparent.getLeft();
        }

        forwardMovedNodes = getMovedPairs(brother);
        forwardGraftNodes = getGraftBranches(uncle);

        exchangeNodes(parent, grandparent, brother, uncle);
    }

    public double rearrangeGeneTrees() {
        final List<GeneTree> geneTrees = geneTreeInput.get();
        final int nGeneTrees = geneTrees.size();
        final List<Map<Node, Integer>> forwardGraftCounts = new ArrayList<>();

        for (int j = 0; j < nGeneTrees; j++) {
            final Map<Node, Integer> jForwardGraftCounts = new HashMap<>();
            final Set<Node> jForwardGraftNodes = forwardGraftNodes.get(j);
            final SortedMap<Node, Node> jForwardMovedNodes = forwardMovedNodes.get(j);
            for (final Entry<Node, Node> nodePair: jForwardMovedNodes.entrySet()) {
                final Node movedNode = nodePair.getKey();
                final Node disownedChild = nodePair.getValue();
                final double movedNodeHeight = movedNode.getHeight();

                final List<Node> validGraftBranches = new ArrayList<>();
                int forwardGraftCount = 0;
                for (Node potentialGraft: jForwardGraftNodes) {
                    final double potentialGraftBottom = potentialGraft.getHeight();
                    double potentialGraftTop;

                    if (potentialGraft.isRoot()) {
                        potentialGraftTop = Double.POSITIVE_INFINITY;
                    } else {
                        potentialGraftTop = potentialGraft.getParent().getHeight();
                    }

                    if (movedNodeHeight > potentialGraftBottom && movedNodeHeight < potentialGraftTop) {
                        forwardGraftCount++;
                        validGraftBranches.add(potentialGraft);
                    }
                }

                // no compatible branches to graft this node on to
                // this only occurs when there is missing data and the gene tree root is in the "parent" branch
                // or if two gene tree nodes which need moving are of equal height
                if (forwardGraftCount == 0) {
                    return Double.NEGATIVE_INFINITY;
                } else {
                    final Node chosenGraft = validGraftBranches.get(Randomizer.nextInt(forwardGraftCount));
                    pruneAndRegraft(movedNode, chosenGraft, disownedChild);
                    movedNode.makeDirty(Tree.IS_FILTHY);
                    jForwardGraftCounts.put(movedNode, forwardGraftCount);
                }
            }

            forwardGraftCounts.add(jForwardGraftCounts);
        }

        SetMultimap<Integer, Node> reverseGraftNodes = getGraftBranches(brother);

        double logHastingsRatio = 0.0;
        for (int j = 0; j < nGeneTrees; j++) {
            final Map<Node, Integer> jForwardGraftCounts = forwardGraftCounts.get(j);
            final Set<Node> jReverseGraftNodes = reverseGraftNodes.get(j);
            for (Node movedNode: jForwardGraftCounts.keySet()) {
                final double movedNodeHeight = movedNode.getHeight();
                final int forwardGraftCount = jForwardGraftCounts.get(movedNode);
                int reverseGraftCount = 0;

                for (Node potentialGraft: jReverseGraftNodes) {
                    final double potentialGraftBottom = potentialGraft.getHeight();
                    double potentialGraftTop;
                    if (potentialGraft.isRoot()) {
                        potentialGraftTop = Double.POSITIVE_INFINITY;
                    } else {
                        potentialGraftTop = potentialGraft.getParent().getHeight();
                    }

                    if (movedNodeHeight > potentialGraftBottom && movedNodeHeight < potentialGraftTop) {
                        reverseGraftCount++;
                    }
                }

                logHastingsRatio += Math.log(forwardGraftCount) - Math.log(reverseGraftCount);
            }
        }

        return logHastingsRatio;
    }

    protected static void pruneAndRegraft(final Node nodeToMove, final Node newChild, final Node disownedChild) {
        final Node sourceParent = nodeToMove.getParent();
        final Node destinationParent = newChild.getParent();

        // debug string
        // System.out.println(String.format("%d-%d-%d > %d-%d", sourceParent.getNr(), nodeToMove.getNr(), disownedChild.getNr(), destinationParent.getNr(), newChild.getNr()));

        nodeToMove.removeChild(disownedChild);
        sourceParent.removeChild(nodeToMove);
        destinationParent.removeChild(newChild);

        nodeToMove.addChild(newChild);
        sourceParent.addChild(disownedChild);
        destinationParent.addChild(nodeToMove);

        nodeToMove.makeDirty(Tree.IS_FILTHY);
        newChild.makeDirty(Tree.IS_FILTHY);
        disownedChild.makeDirty(Tree.IS_FILTHY);
        sourceParent.makeDirty(Tree.IS_FILTHY);
        destinationParent.makeDirty(Tree.IS_FILTHY);
    }

    protected static void exchangeNodes(final Node parent, final Node grandparent, final Node brother, final Node uncle) {
        grandparent.removeChild(uncle);
        parent.addChild(uncle);

        parent.removeChild(brother);
        grandparent.addChild(brother);

        parent.makeDirty(Tree.IS_FILTHY);
        grandparent.makeDirty(Tree.IS_FILTHY);        
        brother.makeDirty(Tree.IS_FILTHY);
        uncle.makeDirty(Tree.IS_FILTHY);
    }

    // identify nodes that can serve as graft branches as part of a coordinated exchange move
    protected SetMultimap<Integer, Node> getGraftBranches(Node uncleNode) {
        final int uncleNodeNumber = uncleNode.getNr();
        final Set<String> uncleDescendants = findDescendants(uncleNode, uncleNodeNumber);

        final SetMultimap<Integer, Node> allGraftBranches = HashMultimap.create();
        final List<GeneTree> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final Node geneTreeRootNode = geneTrees.get(j).getRoot();
            final Set<Node> jGraftBranches = new HashSet<Node>();
            findGraftBranches(geneTreeRootNode, jGraftBranches, uncleDescendants);
            allGraftBranches.putAll(j, jGraftBranches);
        }

        return allGraftBranches;
    }

    // identify nodes that can serve as graft branches as part of a coordinated exchange move
    private boolean findGraftBranches(Node geneTreeNode, Set<Node> graftNodes, Set<String> branchDescendants) {
        if (geneTreeNode.isLeaf()) {
            final String descendantName = geneTreeNode.getID();
            return branchDescendants.contains(descendantName);
        }

        final Node leftChild = geneTreeNode.getLeft();
        final Node rightChild = geneTreeNode.getRight();
        final boolean leftOverlaps = findGraftBranches(leftChild, graftNodes, branchDescendants);
        final boolean rightOverlaps = findGraftBranches(rightChild, graftNodes, branchDescendants);

        // subtree defined by a child node overlaps species subtree defined by branch
        if (leftOverlaps || rightOverlaps) {
            if (leftOverlaps) {
                graftNodes.add(leftChild);
            }

            if (rightOverlaps) {
                graftNodes.add(rightChild);
            }

            return true;
        }

        return false;
    }

    // identify nodes to be grafted in a narrow move, and children to be "disowned" (joined directly to their grandparent)
    private List<SortedMap<Node, Node>> getMovedPairs(Node brotherNode) {
        final int brotherNodeNumber = brotherNode.getNr();
        final Set<String> brotherDescendants = findDescendants(brotherNode, brotherNodeNumber);

        final double lowerHeight = brotherNode.getParent().getHeight(); // parent height (bottom of parent branch)
        final double upperHeight = brotherNode.getParent().getParent().getHeight(); // grandparent height (top of parent branch)

        final List<SortedMap<Node, Node>> allMovedNodes = new ArrayList<>();
        final List<GeneTree> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final Node geneTreeRootNode = geneTrees.get(j).getRoot();
            final SortedMap<Node, Node> jMovedNodes = new TreeMap<>(nhc);
            findMovedPairs(geneTreeRootNode, jMovedNodes, brotherDescendants, lowerHeight, upperHeight);
            allMovedNodes.add(jMovedNodes);
        }

        return allMovedNodes;
    }

    // identify nodes to be moved as part of a coordinated exchange move
    private boolean findMovedPairs(Node geneTreeNode, Map<Node, Node> movedNodes, Set<String> brotherDescendants, double lowerHeight, double upperHeight) {
        if (geneTreeNode.isLeaf()) {
            final String descendantName = geneTreeNode.getID();
            return brotherDescendants.contains(descendantName);
        }

        final Node leftChild = geneTreeNode.getLeft();
        final Node rightChild = geneTreeNode.getRight();

        final boolean leftOverlapsBrother = findMovedPairs(leftChild, movedNodes, brotherDescendants, lowerHeight, upperHeight);
        final boolean rightOverlapsBrother = findMovedPairs(rightChild, movedNodes, brotherDescendants, lowerHeight, upperHeight);

        final double nodeHeight = geneTreeNode.getHeight();
        if (nodeHeight >= lowerHeight && nodeHeight < upperHeight) {
            if (leftOverlapsBrother ^ rightOverlapsBrother) {
                if (leftOverlapsBrother) {
                    movedNodes.put(geneTreeNode, leftChild);
                } else {
                    movedNodes.put(geneTreeNode, rightChild);
                }
            }
        }

        return leftOverlapsBrother || rightOverlapsBrother;
    }
}
