package starbeast2;

import java.util.ArrayList;
import java.util.Arrays;
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
* @author Huw Ogilvie
 */

@Description("Implements the co-ordinated species and gene tree operator described in Rannala & Yang 2015. "
        + "This performs a wide exchange operation, then prunes-and-regrafts gene tree nodes made invalid "
        + "by the operation onto a valid contemporary branch (without changing node heights). "
        + "See http://arxiv.org/abs/1512.03843 for full details.")
public class WideExchange extends CoordinatedOperator {
    private Node aNode; // naming follows Rannala & Yang 2015
    private Node brotherNode;
    private Node parentNode;
    private Node cousinNode;
    private Node mrcaNode;

    private List<List<SortedMap<Node, Node>>> forwardMovedNodes;
    private List<SetMultimap<Integer, Node>> forwardGraftNodes;

    final static Comparator<Node> nhc = new NodeHeightComparator().reversed();

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        double fLogHastingsRatio = pickNodes();

        // only rearrange gene trees if the species tree has changed
        if (fLogHastingsRatio != Double.NEGATIVE_INFINITY) {
            findForwardNodes();
            pruneAndRegraft(parentNode, cousinNode, brotherNode);
            fLogHastingsRatio += rearrangeGeneTrees();
        }

        return fLogHastingsRatio;
    }

    public double pickNodes() {
        final TreeInterface speciesTree = speciesTreeInput.get().getTree();
        final Node[] speciesTreeNodes = speciesTree.getNodesAsArray();

        final int nLeafNodes = speciesTree.getLeafNodeCount();
        final int nInternalNodes = speciesTree.getInternalNodeCount() - 1; // excludes the root
        final int nNodesExceptRoot = nLeafNodes + nInternalNodes;
        final Node rootNode = speciesTreeNodes[nNodesExceptRoot];

        // pick an internal node at random (excluding the root)
        final int parentNodeNumber = nLeafNodes + Randomizer.nextInt(nInternalNodes);
        parentNode = speciesTreeNodes[parentNodeNumber];
        final double parentNodeHeight = parentNode.getHeight();
        if (Randomizer.nextBoolean()) {
            aNode = parentNode.getLeft();
            brotherNode = parentNode.getRight();
        } else {
            aNode = parentNode.getRight();
            brotherNode = parentNode.getLeft();
        }

        // for all internal nodes (excluding the root)
        final int[] mrcaNodeNumbers = new int[nNodesExceptRoot];
        // -1 means the branch defined by this node does not overlap the height of the parentNode 
        // that is it is not a valid graft branch
        Arrays.fill(mrcaNodeNumbers, -1);

        findCousinNodes(parentNode, rootNode, parentNodeHeight, Double.NEGATIVE_INFINITY, mrcaNodeNumbers);

        // pick a cousin from the available candidates
        int cousinNodeNumber = Randomizer.nextInt(nNodesExceptRoot);
        int mrcaNodeNumber = mrcaNodeNumbers[cousinNodeNumber];
        while (mrcaNodeNumber == -1) {
            cousinNodeNumber = Randomizer.nextInt(nNodesExceptRoot);
            mrcaNodeNumber = mrcaNodeNumbers[cousinNodeNumber];
        }

        cousinNode = speciesTreeNodes[cousinNodeNumber];
        mrcaNode = speciesTreeNodes[mrcaNodeNumber];
        
        return 0.0;
    }
    
    public void findForwardNodes() {
        // must be done before any changes are made to the gene or species trees
        Node currentNode = cousinNode;
        forwardMovedNodes = new ArrayList<>();
        forwardGraftNodes = new ArrayList<>();
        while (currentNode != mrcaNode) {
            final Node nextNode = currentNode.getParent();
            final double currentNodeHeight = currentNode.getHeight();
            final double nextNodeHeight = nextNode.getHeight();
            final List<SortedMap<Node, Node>> currentMovedNodes = getMovedPairs(aNode, currentNodeHeight, nextNodeHeight);
            final SetMultimap<Integer, Node> currentGraftNodes = getGraftBranches(currentNode);
            forwardMovedNodes.add(currentMovedNodes);
            forwardGraftNodes.add(currentGraftNodes);
            currentNode = nextNode;
        }
    }

    protected static void pruneAndRegraft(final Node nodeToMove, final Node newChild, final Node disownedChild) {
        final Node parent = nodeToMove.getParent();
        final Node destinationParent = newChild.getParent();

        // debug string
        // System.out.println(String.format("%d-%d-%d > %d-%d", parent.getNr(), nodeToMove.getNr(), disownedChild.getNr(), destinationParent.getNr(), newChild.getNr()));

        nodeToMove.removeChild(disownedChild);
        parent.removeChild(nodeToMove);
        destinationParent.removeChild(newChild);

        nodeToMove.addChild(newChild);
        parent.addChild(disownedChild);
        destinationParent.addChild(nodeToMove);

        nodeToMove.makeDirty(Tree.IS_FILTHY);
        parent.makeDirty(Tree.IS_FILTHY);
        destinationParent.makeDirty(Tree.IS_FILTHY);
        newChild.makeDirty(Tree.IS_FILTHY);
        disownedChild.makeDirty(Tree.IS_FILTHY);
        // markPathFilthy(disownedChild, newChild);
    }

    // mark nodes from two nodes to their MRCA as "filthy"
    /* private static void markPathFilthy(final Node childOne, final Node childTwo) {
        Node iup = childOne;
        Node jup = childTwo;
        iup.makeDirty(Tree.IS_FILTHY);
        jup.makeDirty(Tree.IS_FILTHY);
        while (iup != jup) {
            if( iup.getHeight() < jup.getHeight() ) {
                assert !iup.isRoot();
                iup = iup.getParent();
                iup.makeDirty(Tree.IS_FILTHY);
            } else {
                assert !jup.isRoot();
                jup = jup.getParent();
                jup.makeDirty(Tree.IS_FILTHY);
            }
        }
    } */

    private List<Integer> findCousinNodes(final Node parentNode, final Node currentNode, final double parentNodeHeight, final double branchTopHeight, final int[] mrcaNodeNumbers) {
        // valid graft nodes (nodes defining branches which include the height of the parent node)
        final List<Integer> candidateList = new ArrayList<>();
        final int currentNodeNumber = currentNode.getNr();
        final double currentNodeHeight = currentNode.getHeight();

        if (parentNode == currentNode) {
            return null;
        } else if (parentNodeHeight < branchTopHeight) {
            // this is a candidate node (would be a valid choice to graft parentNode to)
            candidateList.add(currentNode.getNr());
            return candidateList;
        } else {
            final List<Integer> leftCandidateNodeNumbers = findCousinNodes(parentNode, currentNode.getLeft(), parentNodeHeight, currentNodeHeight, mrcaNodeNumbers);
            final List<Integer> rightCandidateNodeNumbers = findCousinNodes(parentNode, currentNode.getRight(), parentNodeHeight, currentNodeHeight, mrcaNodeNumbers);

            if (leftCandidateNodeNumbers == null) { // parent is the left child or descendant of the left child
                // therefore the current node is the most recent common ancestor connecting the parent and right candidates
                for (final Integer candidateNodeNumber: rightCandidateNodeNumbers) {
                    mrcaNodeNumbers[candidateNodeNumber] = currentNodeNumber;
                }
                return null;
            } else if (rightCandidateNodeNumbers == null) { // parent is the right child or descendant of the right child
                // therefore the current node is the most recent common ancestor connecting the parent and left candidates
                for (final Integer candidateNodeNumber: leftCandidateNodeNumbers) {
                    mrcaNodeNumbers[candidateNodeNumber] = currentNodeNumber;
                }
                return null;
            } else {
                candidateList.addAll(leftCandidateNodeNumbers);
                candidateList.addAll(rightCandidateNodeNumbers);
                return candidateList;
            }
        }
    }

    // for testing purposes
    /* public void manipulateSpeciesTree(Node argBrother) {
        brother = argBrother;
        parent = brother.getParent();

        final Node grandparent = parent.getParent();
        if (grandparent.getLeft().getNr() == parent.getNr()) {
            cousin = grandparent.getRight();
        } else {
            cousin = grandparent.getLeft();
        }

        forwardMovedNodes = getMovedPairs(brother);
        forwardGraftNodes = getGraftBranches(cousin);

        exchangeNodes(parent, grandparent, brother, cousin);
    } */

    public double rearrangeGeneTrees() {
        double logHastingsRatio = 0.0;

        for (int i = 0; i < 5; i++) {
            final List<Map<Node, Integer>> forwardGraftCounts = new ArrayList<>();
            for (int j = 0; j < nGeneTrees; j++) {
                final Map<Node, Integer> jForwardGraftCounts = new HashMap<>();
                final Set<Node> jForwardGraftNodes = forwardGraftNodes.get(i).get(j);
                final SortedMap<Node, Node> jForwardMovedNodes = forwardMovedNodes.get(i).get(j);
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
    
            SetMultimap<Integer, Node> reverseGraftNodes = getGraftBranches(brotherNode);
    
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
        }

        return logHastingsRatio;
    }

    // identify nodes that can serve as graft branches as part of a coordinated exchange move
    protected SetMultimap<Integer, Node> getGraftBranches(Node cousinNode) {
        final int cousinNodeNumber = cousinNode.getNr();
        final Set<String> cousinDescendants = findDescendants(cousinNode, cousinNodeNumber);

        final SetMultimap<Integer, Node> allGraftBranches = HashMultimap.create();
        final List<TreeInterface> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final Node geneTreeRootNode = geneTrees.get(j).getRoot();
            final Set<Node> jGraftBranches = new HashSet<Node>();
            findGraftBranches(geneTreeRootNode, jGraftBranches, cousinDescendants);
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
    private List<SortedMap<Node, Node>> getMovedPairs(final Node chosenNode, final double lowerHeight, final double upperHeight) {
        final int chosenNodeNumber = chosenNode.getNr();
        final Set<String> chosenDescendants = findDescendants(chosenNode, chosenNodeNumber);

        final List<SortedMap<Node, Node>> allMovedNodes = new ArrayList<>();
        final List<TreeInterface> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final Node geneTreeRootNode = geneTrees.get(j).getRoot();
            final SortedMap<Node, Node> jMovedNodes = new TreeMap<>(nhc);
            findMovedPairs(geneTreeRootNode, jMovedNodes, chosenDescendants, lowerHeight, upperHeight);
            allMovedNodes.add(jMovedNodes);
        }

        return allMovedNodes;
    }

    // identify nodes to be moved as part of a coordinated exchange move
    private boolean findMovedPairs(Node geneTreeNode, Map<Node, Node> movedNodes, Set<String> chosenDescendants, double lowerHeight, double upperHeight) {
        if (geneTreeNode.isLeaf()) {
            final String descendantName = geneTreeNode.getID();
            // returns true if this leaf is a descendant of the chosen species
            return chosenDescendants.contains(descendantName);
        }

        final Node leftChild = geneTreeNode.getLeft();
        final Node rightChild = geneTreeNode.getRight();

        // left child descendants are exclusively descendants of the chosen species tree node
        final boolean leftExclusive = findMovedPairs(leftChild, movedNodes, chosenDescendants, lowerHeight, upperHeight);
        // right child descendants are exclusively descendants of the chosen species tree node
        final boolean rightExclusive = findMovedPairs(rightChild, movedNodes, chosenDescendants, lowerHeight, upperHeight);

        final double nodeHeight = geneTreeNode.getHeight();
        if (nodeHeight >= lowerHeight && nodeHeight < upperHeight) {
            if (leftExclusive ^ rightExclusive) {
                if (leftExclusive) { // leave right child attached to original parent
                    movedNodes.put(geneTreeNode, rightChild);
                } else { // leaf left child attached to original parent
                    movedNodes.put(geneTreeNode, leftChild);
                }
            }
        }

        // returns true if all descendants of this gene tree node are also descendants of the chosen species tree node
        return leftExclusive && rightExclusive;
    }
}
