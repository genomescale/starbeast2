package starbeast2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
    public Input<TaxonSet> taxonSuperSetInput = new Input<TaxonSet>("taxonSuperSet", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);
    public Input<MultispeciesPopulationModel> populationFunctionInput = new Input<MultispeciesPopulationModel>("populationModel", "The species tree population model.", Validate.REQUIRED);

    private TaxonSet taxonSuperSet;
    private MultispeciesPopulationModel populationModel;
    private int nGeneTrees;
    private int nSpeciesBranches;
    private double[] perGenePloidy;

    private Double[] speciesStartTimes;
    private Double[] speciesEndTimes;

    final private HashMap<String, Integer> tipNumberMap = new HashMap<String, Integer>();

    @Override
    public void initAndValidate() throws Exception {
        final HashMap<String, Integer> speciesNumberMap = new HashMap<String, Integer>();

        TreeInterface speciesTree = treeInput.get();
        List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();

        nGeneTrees = geneTrees.size();

        taxonSuperSet = taxonSuperSetInput.get();
        populationModel = populationFunctionInput.get();

        nSpeciesBranches = speciesTree.getNodeCount();

        speciesStartTimes = new Double[nSpeciesBranches]; // the earlier date (rootward end)
        speciesEndTimes = new Double[nSpeciesBranches]; // the later date (tipward end)

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
            }
        }

        // initialize gene trees and store ploidy
        perGenePloidy = new double[nGeneTrees];
        for (int i = 0; i < nGeneTrees; i++) {
            final GeneTreeWithinSpeciesTree geneTreeI = geneTrees.get(i);
            perGenePloidy[i] = geneTreeI.ploidy;
        }

        populationModel.initPopSizes(nSpeciesBranches);
    }

    public double calculateLogP() {
        final TreeInterface speciesTree = treeInput.get();
        final List<GeneTreeWithinSpeciesTree> geneTrees = geneTreeInput.get();
        logP = 0.0;

        for (int i = 0; i < nSpeciesBranches; i++) {
            final Node speciesNode = speciesTree.getNode(i);
            final Node parentNode = speciesNode.getParent();

            speciesEndTimes[i] = speciesNode.getHeight();
            
            if (parentNode == null) {
                speciesStartTimes[i] = Double.POSITIVE_INFINITY;
            } else {
                speciesStartTimes[i] = parentNode.getHeight();
            }
        }

        final List<int[]> allLineageCounts = new ArrayList<int[]>();
        final List<int[]> allEventCounts = new ArrayList<int[]>();
        final List<List<Double[]>> allCoalescentTimes = new ArrayList<List<Double[]>>();
        for (int i = 0; i < nSpeciesBranches; i++) {
            allLineageCounts.add(new int[nGeneTrees]);
            allEventCounts.add(new int[nGeneTrees]);
            allCoalescentTimes.add(new ArrayList<Double[]>());
        }

        // transpose gene-branch list of lists to branch-gene list of lists
        for (int j = 0; j < nGeneTrees; j++) { // for each gene "j"
            final GeneTreeWithinSpeciesTree geneTree = geneTrees.get(j);
            if (geneTree.computeCoalescentTimes(speciesTree, tipNumberMap)) {
                for (int i = 0; i < nSpeciesBranches; i++) { // for each species tree node/branch "i"
                    final List<Double> branchCoalescentTimes = geneTree.coalescentTimes.get(i);
                    final int geneBranchEventCount = branchCoalescentTimes.size();;
                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts[i];

                    final List<Double> geneBranchTimes = new ArrayList<>();
                    geneBranchTimes.add(speciesStartTimes[i]); // add the start time for each branch
                    geneBranchTimes.addAll(branchCoalescentTimes); // add the coalescent event times for each branch
                    geneBranchTimes.add(speciesEndTimes[i]); // add the end time for each branch

                    final Double[] coalescentTimesIJ = new Double[geneBranchEventCount + 2];
                    geneBranchTimes.toArray(coalescentTimesIJ);
                    Arrays.sort(coalescentTimesIJ);

                    allLineageCounts.get(i)[j] = geneBranchLineageCount;
                    allEventCounts.get(i)[j] = geneBranchEventCount;
                    allCoalescentTimes.get(i).add(coalescentTimesIJ);
                }
            } else { // this gene tree IS NOT compatible with the species tree
                logP = Double.NEGATIVE_INFINITY;
                return logP;
            }
        }

        for (int i = 0; i < nSpeciesBranches; i++) {
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
        boolean allCompatible = true;
        for (GeneTreeWithinSpeciesTree geneTree: geneTrees) {
            allCompatible = allCompatible && geneTree.computeCoalescentTimes(speciesTree, tipNumberMap);
        }
        
        return allCompatible;
    }
}
