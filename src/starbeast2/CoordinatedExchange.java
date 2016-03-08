package starbeast2;

import java.util.ArrayList;
import java.util.Comparator;
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
import beast.core.Input;
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

@Description("Implements the co-ordinated species and gene tree operator similar to those described in "
        + "Yang & Rannala 2015 and Rannala & Yang 2015. This performs an NNI or SPR exchange operation, "
        + "then prunes-and-regrafts gene tree nodes made invalid by the operation onto a valid "
        + "contemporary branch (without changing node heights). "
        + "See http://doi.org/10.1093/molbev/msu279 and http://arxiv.org/abs/1512.03843 for full details.")
public class CoordinatedExchange extends CoordinatedOperator {
    final public Input<Boolean> isNarrowInput = new Input<>("isNarrow", "if true (default) a narrow exchange is performed, otherwise a wide exchange", true);
    final public Input<Boolean> isTestInput = new Input<>("testing", "for performing unit tests, do not pick species tree nodes", false);

    private Node[] speciesTreeNodes;

    // public for unit testing
    public Node aNode; // naming follows Rannala & Yang 2015
    public Node bNode;
    public Node yNode;
    public Node cNode;
    public Node zNode;

    private int nLeafNodes;
    private int nInternalNodes; // excludes the root node
    private int nSpeciesNodes;
    private int czBranchCount;

    private List<List<SortedMap<Node, Node>>> movedNodes;
    private List<SetMultimap<Integer, Node>> graftNodes;
    
    private boolean testing = false;

    // reversed because nodes must be grafted oldest to youngest
    private final static Comparator<Node> nhc = new NodeHeightComparator().reversed();

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        testing = isTestInput.get();
        TreeInterface speciesTree = speciesTreeInput.get().getTree();
        speciesTreeNodes = speciesTree.getNodesAsArray();
        nLeafNodes = speciesTree.getLeafNodeCount();
        nInternalNodes = speciesTree.getInternalNodeCount();
        nSpeciesNodes = speciesTree.getNodeCount();

        boolean isNarrow = isNarrowInput.get();
        double logHastingsRatio = 0.0;
        if (isNarrow) {
            // only proceed to rearrange gene trees if the species tree can be changed
            // doesn't execute if testing
            if (!testing && !pickNarrow()) return Double.NEGATIVE_INFINITY;

            int validGP = 0;
            for(int i = nLeafNodes; i < nSpeciesNodes; ++i) {
                validGP += isg(speciesTree.getNode(i));
            }

            final int c2 = sisg(yNode) + sisg(cNode);

            fillNodes(); // fills in movedNodes and graftNodes
            pruneAndRegraft(yNode, cNode, bNode);

            final int validGPafter = validGP - c2 + sisg(yNode) + sisg(cNode);

            logHastingsRatio += Math.log(validGP) - Math.log(validGPafter);
        } else {
            // only proceed to rearrange gene trees if the species tree can be changed
            // doesn't execute if testing
            if (!testing && !pickWide()) return Double.NEGATIVE_INFINITY;

            fillNodes(); // fills in movedNodes and graftNodes
            pruneAndRegraft(yNode, cNode, bNode);
        }

        for (int i = 0; i < czBranchCount; i++) {
            final List<SortedMap<Node, Node>> perBranchMovedNodes = movedNodes.get(i);
            final SetMultimap<Integer, Node> perBranchGraftNodes = graftNodes.get(i);
            final double logForward = rearrangeGeneTrees(perBranchMovedNodes, perBranchGraftNodes, true);
            assert logForward != Double.NEGATIVE_INFINITY;
            if (logForward == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
            else logHastingsRatio += logForward;
        }

        // compute reverse move (Hastings ratio denominator)
        final Node bNodeTmp = bNode;
        final Node cNodeTmp = cNode;

        bNode = cNodeTmp;
        cNode = bNodeTmp;

        fillNodes(); // fills in movedNodes and graftNodes for reverse move

        for (int i = 0; i < czBranchCount; i++) {
            final List<SortedMap<Node, Node>> perBranchMovedNodes = movedNodes.get(i);
            final SetMultimap<Integer, Node> perBranchGraftNodes = graftNodes.get(i);
            final double logReverse = rearrangeGeneTrees(perBranchMovedNodes, perBranchGraftNodes, false);
            assert logReverse != Double.NEGATIVE_INFINITY;
            if (logReverse == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
            else logHastingsRatio -= logReverse;
        }

        return logHastingsRatio;
    }

    private boolean pickNarrow() {
        zNode = speciesTreeNodes[nLeafNodes + Randomizer.nextInt(nInternalNodes)];
        while (zNode.getLeft().isLeaf() && zNode.getRight().isLeaf()) {
            zNode = speciesTreeNodes[nLeafNodes + Randomizer.nextInt(nInternalNodes)];
        }

        yNode = zNode.getLeft();
        cNode = zNode.getRight();
        if (yNode.getHeight() < cNode.getHeight()) {
            yNode = zNode.getRight();
            cNode = zNode.getLeft();
        }

        if (yNode.isLeaf()) {
            return false;
        } else if (Randomizer.nextBoolean()) {
            aNode = yNode.getLeft();
            bNode = yNode.getRight();
        } else {
            aNode = yNode.getRight();
            bNode = yNode.getLeft();
        }

        return true;
    }

    private boolean pickWide() {
        final int nNodesExceptRoot = nSpeciesNodes - 1;
        final Node rootNode = speciesTreeNodes[nNodesExceptRoot];

        // pick an internal node at random (excluding the root)
        final int yNodeNumber = nLeafNodes + Randomizer.nextInt(nInternalNodes - 1);
        yNode = speciesTreeNodes[yNodeNumber];
        final double yNodeHeight = yNode.getHeight();
        if (Randomizer.nextBoolean()) {
            aNode = yNode.getLeft();
            bNode = yNode.getRight();
        } else {
            aNode = yNode.getRight();
            bNode = yNode.getLeft();
        }

        // for all internal nodes (excluding the root)
        final Node[] zNodes = new Node[nNodesExceptRoot];

        czNodeFinder(yNode, rootNode, yNodeHeight, zNodes);

        // pick a cousin from the available candidates
        int cousinNodeNumber = Randomizer.nextInt(nNodesExceptRoot);
        zNode = zNodes[cousinNodeNumber];
        while (zNode == null) {
            cousinNodeNumber = Randomizer.nextInt(nNodesExceptRoot);
            //System.out.println(String.format("%d/%d", cousinNodeNumber, nNodesExceptRoot));
            zNode = zNodes[cousinNodeNumber];
        }

        cNode = speciesTreeNodes[cousinNodeNumber];

        return true;
    }

    private List<Integer> czNodeFinder(final Node parentNode, final Node currentNode, final double parentNodeHeight, final Node[] zNodes) {
        // valid graft nodes (nodes defining branches which include the height of the parent node)
        final List<Integer> candidateList = new ArrayList<>();
        final double currentNodeHeight = currentNode.getHeight();

        if (parentNode == currentNode) {
            return null;
        } else if (parentNodeHeight >= currentNodeHeight) {
            // this is a candidate node (would be a valid choice to graft parentNode to)
            candidateList.add(currentNode.getNr());
            return candidateList;
        } else {
            final List<Integer> leftCandidateNodeNumbers = czNodeFinder(parentNode, currentNode.getLeft(), parentNodeHeight, zNodes);
            final List<Integer> rightCandidateNodeNumbers = czNodeFinder(parentNode, currentNode.getRight(), parentNodeHeight, zNodes);

            if (leftCandidateNodeNumbers == null) { // parent is the left child or descendant of the left child
                // therefore the current node is the most recent common ancestor connecting the parent and right candidates
                for (final Integer candidateNodeNumber: rightCandidateNodeNumbers) {
                    zNodes[candidateNodeNumber] = currentNode;
                }
                return null;
            } else if (rightCandidateNodeNumbers == null) { // parent is the right child or descendant of the right child
                // therefore the current node is the most recent common ancestor connecting the parent and left candidates
                for (final Integer candidateNodeNumber: leftCandidateNodeNumbers) {
                    zNodes[candidateNodeNumber] = currentNode;
                }
                return null;
            } else {
                candidateList.addAll(leftCandidateNodeNumbers);
                candidateList.addAll(rightCandidateNodeNumbers);
                return candidateList;
            }
        }
    }

    // fills forward nodes by destination branch (c through z) 
    private void fillNodes() {
        // must be done before any changes are made to the gene or species trees
        Node childNode = cNode;
        movedNodes = new ArrayList<>();
        graftNodes = new ArrayList<>();
        czBranchCount = 0;
        while (childNode != zNode) {
            czBranchCount++;
            final Node parentNode = childNode.getParent();
            final double childNodeHeight = childNode.getHeight();
            final double parentNodeHeight = parentNode.getHeight();
            final List<SortedMap<Node, Node>> perBranchMovedNodes = getMovedPairs(childNodeHeight, parentNodeHeight);
            final SetMultimap<Integer, Node> perBranchGraftNodes = getGraftBranches(childNode);
            movedNodes.add(0, perBranchMovedNodes); // needs to be added in reverse order
            graftNodes.add(0, perBranchGraftNodes); // because nodes must be grafted oldest to youngest
            childNode = parentNode;
        }
    }

    private static void pruneAndRegraft(final Node nodeToMove, final Node newChild, final Node disownedChild) {
        final Node sourceParent = nodeToMove.getParent();
        final Node destinationParent = newChild.getParent();

        // debug string
        // System.out.println(String.format("%d-%d-%d > %d-%d", parent.getNr(), nodeToMove.getNr(), disownedChild.getNr(), destinationParent.getNr(), newChild.getNr()));

        nodeToMove.removeChild(disownedChild);
        sourceParent.removeChild(nodeToMove);
        destinationParent.removeChild(newChild);

        nodeToMove.addChild(newChild);
        sourceParent.addChild(disownedChild);
        destinationParent.addChild(nodeToMove);

        nodeToMove.makeDirty(Tree.IS_FILTHY);
        sourceParent.makeDirty(Tree.IS_FILTHY);
        destinationParent.makeDirty(Tree.IS_FILTHY);
        newChild.makeDirty(Tree.IS_FILTHY);
        disownedChild.makeDirty(Tree.IS_FILTHY);
    }

    private double rearrangeGeneTrees(List<SortedMap<Node, Node>> branchMovedNodes, SetMultimap<Integer, Node> branchGraftNodes, boolean forwardMove) {
        double logHastingsRatio = 0.0;

        for (int j = 0; j < nGeneTrees; j++) {
            final Set<Node> jForwardGraftNodes = branchGraftNodes.get(j);
            final SortedMap<Node, Node> jForwardMovedNodes = branchMovedNodes.get(j);
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
                    logHastingsRatio += Math.log(forwardGraftCount);
                    if (forwardMove) { // only actually change gene trees if this is the forward move
                        final Node chosenGraft = validGraftBranches.get(Randomizer.nextInt(forwardGraftCount));
                        pruneAndRegraft(movedNode, chosenGraft, disownedChild);
                    }
                }
            }
        }

        return logHastingsRatio;
    }

    // identify nodes that can serve as graft branches as part of a coordinated exchange move
    protected SetMultimap<Integer, Node> getGraftBranches(Node yNode) {
        final int yNumber = yNode.getNr();
        final Set<String> cousinDescendants = findDescendants(yNode, yNumber);

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
    private List<SortedMap<Node, Node>> getMovedPairs(final double lowerHeight, final double upperHeight) {
        final int aNodeNumber = aNode.getNr();
        final Set<String> chosenDescendants = findDescendants(aNode, aNodeNumber);

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

    private int isg(final Node n) {
      return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    private int sisg(final Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }
}
