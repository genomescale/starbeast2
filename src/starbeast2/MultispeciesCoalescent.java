package starbeast2;

import java.util.ArrayList;
import java.util.Arrays;
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

    private TaxonSet taxonSuperSet;
    private MultispeciesPopulationModel populationModel;
    private int nGeneTrees;
    private int speciesTreeNodeCount;
    private double[] perGenePloidy;

    private Double[] speciesStartTimes;
    private Double[] speciesEndTimes;

    final private Map<String, Integer> tipNumberMap = new HashMap<>();
    final private Multimap<Integer, String> numberTipMap = HashMultimap.create();

    final static NodeHeightComparator nhc = new NodeHeightComparator();

    @Override
    public void initAndValidate() throws Exception {
        final HashMap<String, Integer> speciesNumberMap = new HashMap<>();

        TreeInterface speciesTree = treeInput.get();
        List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();

        nGeneTrees = geneTrees.size();

        taxonSuperSet = taxonSuperSetInput.get();
        populationModel = populationFunctionInput.get();

        speciesTreeNodeCount = speciesTree.getNodeCount();

        speciesStartTimes = new Double[speciesTreeNodeCount]; // the earlier date (rootward end)
        speciesEndTimes = new Double[speciesTreeNodeCount]; // the later date (tipward end)

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
    }

    public double calculateLogP() {
        final TreeInterface speciesTree = treeInput.get();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        logP = 0.0;

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

        final List<int[]> allLineageCounts = new ArrayList<>();
        final List<int[]> allEventCounts = new ArrayList<>();
        final List<List<Double[]>> allCoalescentTimes = new ArrayList<>();
        for (int i = 0; i < speciesTreeNodeCount; i++) {
            allLineageCounts.add(new int[nGeneTrees]);
            allEventCounts.add(new int[nGeneTrees]);
            allCoalescentTimes.add(new ArrayList<>());
        }

        // transpose gene-branch list of lists to branch-gene list of lists
        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            if (geneTree.computeCoalescentTimes(speciesTree, tipNumberMap)) {
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
                logP = Double.NEGATIVE_INFINITY;
                return logP;
            }
        }

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

    public String serializePopulation(Node node) {
        return populationModel.serialize(node);
    }

    @Override
    public boolean requiresRecalculation() {
        return true;
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

    // identify children to be moved with parent nodes as part of a coordinated exchange move
    protected List<SortedMap<Node, Node>> getMovedChildren(Node brotherNode) {
        final int brotherNodeNumber = brotherNode.getNr();
        final Set<String> brotherDescendants = findDescendants(brotherNode, brotherNodeNumber);

        final double lowerHeight = brotherNode.getParent().getHeight(); // parent height (bottom of parent branch)
        final double upperHeight = brotherNode.getParent().getParent().getHeight(); // grandparent height (top of parent branch)

        final List<SortedMap<Node, Node>> allBranchNodes = new ArrayList<>();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            final Node geneTreeRootNode = geneTree.getRoot();
            final SortedMap<Node, Node> movedNodes = new TreeMap<>(nhc);
            geneTree.findMovedChildren(geneTreeRootNode, movedNodes, brotherDescendants, lowerHeight, upperHeight);
            allBranchNodes.add(movedNodes);
        }

        return allBranchNodes;
    }

    // identify nodes that can serve as graft branches as part of a coordinated exchange move
    protected SetMultimap<Integer, Node> getGraftBranches(Node uncleNode) {
        final int uncleNodeNumber = uncleNode.getNr();
        final Set<String> uncleDescendants = findDescendants(uncleNode, uncleNodeNumber);

        final SetMultimap<Integer, Node> allGraftBranches = HashMultimap.create();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        for (int j = 0; j < nGeneTrees; j++) {
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            final Node geneTreeRootNode = geneTree.getRoot();
            final Set<Node> graftBranches = new HashSet<Node>();
            geneTree.findGraftBranches(geneTreeRootNode, graftBranches, uncleDescendants);
            allGraftBranches.putAll(j, graftBranches);
        }

        return allGraftBranches;
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
}

