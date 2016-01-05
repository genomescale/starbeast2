package starbeast2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.SetMultimap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;

/**
* @author Remco Bouckaert
* @author Joseph Heled
* @author Huw Ogilvie
 */

@Description("Calculates probability of gene trees conditioned on a species tree (the multi-species coalescent).")
public class MultispeciesCoalescent extends TreeDistribution {
    public Input<List<GeneTreeWithinSpeciesTree>> geneTreeInput = new Input<>("geneTree", "Gene tree within the species tree.", new ArrayList<>());
    public Input<TaxonSet> taxonSuperSetInput = new Input<>("taxonSuperSet", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);
    public Input<MultispeciesPopulationModel> populationFunctionInput = new Input<>("populationModel", "The species tree population model.", Validate.REQUIRED);

    private boolean needsUpdate;
    private boolean allTreesCompatible;
    private TreeInterface speciesTree;
    private TaxonSet taxonSuperSet;
    private MultispeciesPopulationModel populationModel;
    private int nGeneTrees;
    private int speciesTreeNodeCount;
    private double[] perGenePloidy;

    private double[] speciesStartTimes;
    private double[] speciesEndTimes;

    final private Map<String, Integer> tipNumberMap = new HashMap<>();
    final private Map<Node, double[]> speciesOccupancy = new HashMap<>();
    final private Multimap<Integer, String> numberTipMap = HashMultimap.create();
    final private List<int[]> allLineageCounts = new ArrayList<>();
    final private List<int[]> allEventCounts = new ArrayList<>();
    final private List<List<Double[]>> allCoalescentTimes = new ArrayList<>();

    final static Comparator<Node> nhc = new NodeHeightComparator().reversed();

    private enum descendsThrough {
       LEFT_ONLY, RIGHT_ONLY, BOTH, NEITHER
    }

    @Override
    public void initAndValidate() throws Exception {
        final HashMap<String, Integer> speciesNumberMap = new HashMap<>();

        TreeInterface speciesTree = treeInput.get();
        List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();

        nGeneTrees = geneTrees.size();

        taxonSuperSet = taxonSuperSetInput.get();
        populationModel = populationFunctionInput.get();

        speciesTreeNodeCount = speciesTree.getNodeCount();

        speciesStartTimes = new double[speciesTreeNodeCount]; // the earlier date (rootward end)
        speciesEndTimes = new double[speciesTreeNodeCount]; // the later date (tipward end)

        // generate map of species tree tip node names to node numbers
        Node speciesTreeRoot = speciesTree.getRoot();
        for (Node leafNode: speciesTreeRoot.getAllLeafNodes()) {
            final String speciesName = leafNode.getID();
            final int speciesNumber = leafNode.getNr();

            speciesNumberMap.put(speciesName, speciesNumber);
        }

        // generate map of gene tree tip node names to species tree tip node numbers
        /** replacement line for compatibility with BEAST 2.3.0 **/
        // final Set<Taxon> speciesSet = taxonSuperSet.getTaxonSet();
        final Set<Taxon> speciesSet = new HashSet<>(taxonSuperSet.taxonsetInput.get());

        for (Taxon species: speciesSet) {
            final String speciesName = species.getID();
            final int speciesNumber = speciesNumberMap.get(speciesName);
            final TaxonSet speciesTaxonSet = (TaxonSet) species;
            /** replacement line for compatibility with BEAST 2.3.0 **/
            // final Set<Taxon> tipSet = speciesTaxonSet.getTaxonSet();
            final Set<Taxon> tipSet = new HashSet<>(speciesTaxonSet.taxonsetInput.get());

            for (Taxon tip: tipSet) {
                final String tipName = tip.getID();
                tipNumberMap.put(tipName, speciesNumber);
                numberTipMap.put(speciesNumber, tipName);
            }
        }

        // initialize gene trees and store ploidy
        perGenePloidy = new double[nGeneTrees];
        for (int i = 0; i < nGeneTrees; i++) {
            final GeneTreeWithinSpeciesTree geneTreeI = geneTrees.get(i);
            perGenePloidy[i] = geneTreeI.ploidy;
        }

        populationModel.initPopSizes(speciesTreeNodeCount);
        
        needsUpdate = true;
    }

    void update() {
        speciesTree = treeInput.get();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        logP = 0.0;

        assert checkTreeSanity(speciesTree.getRoot()); // species tree should not be insane

        for (int i = 0; i < speciesTreeNodeCount; i++) {
            final Node speciesNode = speciesTree.getNode(i);
            final Node parentNode = speciesNode.getParent();

            speciesEndTimes[i] = speciesNode.getHeight();
            
            if (parentNode == null) {
                speciesStartTimes[i] = Double.POSITIVE_INFINITY;
            } else {
                speciesStartTimes[i] = parentNode.getHeight();
            }
        }

        allLineageCounts.clear();
        allEventCounts.clear();
        allCoalescentTimes.clear();
        speciesOccupancy.clear();

        for (int i = 0; i < speciesTreeNodeCount; i++) {
            allLineageCounts.add(new int[nGeneTrees]);
            allEventCounts.add(new int[nGeneTrees]);
            allCoalescentTimes.add(new ArrayList<>());
        }

        // transpose gene-branch list of lists to branch-gene list of lists
        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            assert checkTreeSanity(geneTree.getRoot()); // gene trees should not be insane either
            if (geneTree.computeCoalescentTimes(speciesTree, tipNumberMap)) {
                geneTree.addSpeciesOccupancy(speciesOccupancy);
                for (int i = 0; i < speciesTreeNodeCount; i++) { // for each species tree node/branch "i"
                    final List<Double> timesView = geneTree.coalescentTimes.get(i);
                    final int geneBranchEventCount = timesView.size();
                    final Double[] geneBranchCoalescentTimes = new Double[geneBranchEventCount];
                    timesView.toArray(geneBranchCoalescentTimes);
                    Arrays.sort(geneBranchCoalescentTimes);

                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts.count(i);

                    final Double[] coalescentTimesIJ = new Double[geneBranchEventCount + 2];
                    coalescentTimesIJ[0] = speciesEndTimes[i];
                    for (int k = 0; k < geneBranchEventCount; k++) {
                        coalescentTimesIJ[k + 1] = geneBranchCoalescentTimes[k];
                    }
                    coalescentTimesIJ[geneBranchEventCount + 1] = speciesStartTimes[i];

                    allLineageCounts.get(i)[j] = geneBranchLineageCount;
                    allEventCounts.get(i)[j] = geneBranchEventCount;
                    allCoalescentTimes.get(i).add(coalescentTimesIJ);
                }
            } else { // this gene tree IS NOT compatible with the species tree
                allTreesCompatible = false;
                needsUpdate = false;
                return;
            }
        }

        allTreesCompatible = true;
        needsUpdate = false;
    }

    public double calculateLogP() {
        if (needsUpdate) {
            update();
        }

        if (!allTreesCompatible) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        logP = 0.0;
        for (int i = 0; i < speciesTreeNodeCount; i++) {
            final Node speciesTreeNode = speciesTree.getNode(i); 
            final List<Double[]> branchCoalescentTimes = allCoalescentTimes.get(i);
            final int[] branchLineageCounts = allLineageCounts.get(i);
            final int[] branchEventCounts = allEventCounts.get(i);
            logP += populationModel.branchLogP(i, speciesTreeNode, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
        }

        return logP;
    }

    public double getRootHeight() {
        double tallestGeneTreeHeight = 0.0;

        for (final GeneTreeWithinSpeciesTree geneTree: geneTreeInput.get()) {
            tallestGeneTreeHeight = Math.max(tallestGeneTreeHeight, geneTree.getTreeHeight());
        }

        return tallestGeneTreeHeight;
    }

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = true;
        return true;
    }

    public void restore() {
        needsUpdate = true;
        super.restore();
    }

    @Override
    public boolean canHandleTipDates() {
        return false;
    }

    public Tree getSpeciesTree() {
        return (Tree) treeInput.get();
    }

    public List<GeneTreeWithinSpeciesTree> getGeneTrees() {
        return geneTreeInput.get();
    }

    // for testing purposes (called by assert statement)
    public boolean computeCoalescentTimes() {
        final TreeInterface speciesTree = treeInput.get();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        for (GeneTreeWithinSpeciesTree geneTree: geneTrees) {
            if (!geneTree.computeCoalescentTimes(speciesTree, tipNumberMap)) {
                // this gene tree IS NOT compatible with the species tree
                return false;
            }
        }

        return true;
    }

    // identify nodes to be grafted in a narrow move, and children to be "disowned" (joined directly to their grandparent)
    protected List<SortedMap<Node, Node>> getMovedPairs(Node brotherNode) {
        final int brotherNodeNumber = brotherNode.getNr();
        final Set<String> brotherDescendants = findDescendants(brotherNode, brotherNodeNumber);

        final double lowerHeight = brotherNode.getParent().getHeight(); // parent height (bottom of parent branch)
        final double upperHeight = brotherNode.getParent().getParent().getHeight(); // grandparent height (top of parent branch)

        final List<SortedMap<Node, Node>> allMovedNodes = new ArrayList<>();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final Node geneTreeRootNode = geneTrees.get(j).getRoot();
            final SortedMap<Node, Node> jMovedNodes = new TreeMap<>(nhc);
            findMovedPairs(geneTreeRootNode, jMovedNodes, brotherDescendants, lowerHeight, upperHeight);
            allMovedNodes.add(jMovedNodes);
        }

        return allMovedNodes;
    }

    // identify nodes that can serve as graft branches as part of a coordinated exchange move
    protected SetMultimap<Integer, Node> getGraftBranches(Node uncleNode) {
        final int uncleNodeNumber = uncleNode.getNr();
        final Set<String> uncleDescendants = findDescendants(uncleNode, uncleNodeNumber);

        final SetMultimap<Integer, Node> allGraftBranches = HashMultimap.create();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final Node geneTreeRootNode = geneTrees.get(j).getRoot();
            final Set<Node> jGraftBranches = new HashSet<Node>();
            findGraftBranches(geneTreeRootNode, jGraftBranches, uncleDescendants);
            allGraftBranches.putAll(j, jGraftBranches);
        }

        return allGraftBranches;
    }

    // identify gene tree nodes which descend through both (and also descend exclusively through)
    // the left and right children of the species tree node of interest
    protected SetMultimap<Integer, Node> getConnectingNodes(Node speciesTreeNode, MinimumDouble tipwardFreedom, MinimumDouble rootwardFreedom) {
        final Node leftChildNode = speciesTreeNode.getLeft();
        final Node rightChildNode = speciesTreeNode.getRight();
        final int leftChildNodeNumber = leftChildNode.getNr();
        final int rightChildNodeNumber = rightChildNode.getNr();
        final Set<String> leftChildDescendants = findDescendants(leftChildNode, leftChildNodeNumber);
        final Set<String> rightChildDescendants = findDescendants(rightChildNode, rightChildNodeNumber);

        final SetMultimap<Integer, Node> allConnectingNodes = HashMultimap.create();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final Node geneTreeRootNode = geneTrees.get(j).getRoot();
            final Set<Node> jConnectingNodes = new HashSet<Node>();
            findConnectingNodes(geneTreeRootNode, jConnectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom, rootwardFreedom);
            allConnectingNodes.putAll(j, jConnectingNodes);
        }

        return allConnectingNodes;
    }

    private Set<String> findDescendants(Node speciesTreeNode, int speciesTreeNodeNumber) {
        final Set<String> descendantNames = new HashSet<>();

        if (speciesTreeNode.isLeaf()) {
            descendantNames.addAll(numberTipMap.get(speciesTreeNodeNumber));
        } else {
            final Node leftChild = speciesTreeNode.getLeft();
            final Node rightChild = speciesTreeNode.getRight();
            final int leftChildNumber = leftChild.getNr();
            final int rightChildNumber = rightChild.getNr();

            descendantNames.addAll(findDescendants(leftChild, leftChildNumber));
            descendantNames.addAll(findDescendants(rightChild, rightChildNumber));
        }

        return descendantNames;
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

    private descendsThrough findConnectingNodes(Node geneTreeNode, Set<Node> connectingNodes, Set<String> leftChildDescendants, Set<String> rightChildDescendants, MinimumDouble tipwardFreedom, MinimumDouble rootwardFreedom) {
        if (geneTreeNode.isLeaf()) {
            final String descendantName = geneTreeNode.getID();
            if (leftChildDescendants.contains(descendantName)) {
                return descendsThrough.LEFT_ONLY;
            } else if (rightChildDescendants.contains(descendantName)) {
                return descendsThrough.RIGHT_ONLY;
            } else {
                return descendsThrough.NEITHER;
            }
        }

        final Node leftChild = geneTreeNode.getLeft();
        final Node rightChild = geneTreeNode.getRight();
        final descendsThrough leftDescent = findConnectingNodes(leftChild, connectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom, rootwardFreedom);
        final descendsThrough rightDescent = findConnectingNodes(rightChild, connectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom, rootwardFreedom);

        if (leftDescent == rightDescent) {
            if (leftDescent == descendsThrough.BOTH) {
                connectingNodes.add(geneTreeNode);
            }

            return leftDescent;
        }

        // this code only executes when the left and right gene tree child nodes descend through different species tree node of interest children
        final double geneTreeNodeHeight = geneTreeNode.getHeight();
        if (leftDescent == descendsThrough.BOTH) { // the gene tree node left child is a member of a connected component
            if (rightDescent == descendsThrough.NEITHER) { // the gene tree node left child is the root node of a connected component
                final double connectedComponentRootFreedom = geneTreeNodeHeight - leftChild.getHeight();
                rootwardFreedom.set(connectedComponentRootFreedom);
                return descendsThrough.NEITHER;
            } else { // the gene tree node right child descends exclusively through the left XOR right child of the species tree node of interest
                // so the current gene tree node is part of a connected component but the right child is not
                final double connectedComponentDescendantBranchLength = geneTreeNodeHeight - rightChild.getHeight();
                tipwardFreedom.set(connectedComponentDescendantBranchLength);
                connectingNodes.add(geneTreeNode);
                return descendsThrough.BOTH;
            }
        } else if (rightDescent == descendsThrough.BOTH) { // the gene tree node right child is a member of a connected component
            if (leftDescent == descendsThrough.NEITHER) { // the gene tree node right child is the root node of a connected component
                final double connectedComponentRootFreedom = geneTreeNodeHeight - rightChild.getHeight();
                rootwardFreedom.set(connectedComponentRootFreedom);
                return descendsThrough.NEITHER;
            } else { // the gene tree node left child descends exclusively through the left XOR right child of the species tree node of interest
             // so the current gene tree node is part of a connected component but the left child is not
                final double connectedComponentTipFreedom = geneTreeNodeHeight - leftChild.getHeight();
                tipwardFreedom.set(connectedComponentTipFreedom);
                connectingNodes.add(geneTreeNode);
                return descendsThrough.BOTH;
            }
        } else if (leftDescent == descendsThrough.NEITHER || rightDescent == descendsThrough.NEITHER) {
            return descendsThrough.NEITHER; // the current gene tree node does not descend exclusively through the species tree node of interest
        } else { // this is a tip node of a connected component
            final double leftChildBranchLength = geneTreeNodeHeight - leftChild.getHeight();
            final double rightChildBranchLength = geneTreeNodeHeight - rightChild.getHeight();
            tipwardFreedom.set(leftChildBranchLength);
            tipwardFreedom.set(rightChildBranchLength);
            connectingNodes.add(geneTreeNode);
            return descendsThrough.BOTH;
        }
    }

    private boolean checkTreeSanity(Node node) {
        final List<Node> children = node.getChildren();
        final int nChildren = children.size();

        for (Node childNode: children) {
            assert childNode.getParent() == node;
            assert childNode.getHeight() <= node.getHeight();
            if (!node.isLeaf()) {
                checkTreeSanity(childNode);
            }
        }

        if (node.isLeaf()) {
            assert nChildren == 0;
        } else {
            assert nChildren == 2;
        }

        return true;
    }

    public double[] getOccupancy(Node node) {
        if (needsUpdate) {
            update();
        }

        if (speciesOccupancy.get(node) == null) {
            System.out.println(String.format("%s: %d %d", node.getTree().getID(), node.getNr(), node.hashCode()));
            for (Node n: speciesOccupancy.keySet()) {
                if (n.getTree() == node.getTree()) {
                    System.out.println(String.format("%s: %d %d", n.getTree().getID(), n.getNr(), n.hashCode()));
                }
            }
        }

        return speciesOccupancy.get(node);
    }
}
