package starbeast2;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

/**
* @author Remco Bouckaert
* @author Huw Ogilvie
 */

public class GeneTree extends Distribution {
    public Input<Tree> treeInput = new Input<>("tree", "The gene tree.", Validate.REQUIRED);
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "Species tree for embedding the gene tree.", Validate.REQUIRED);
    public Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for this gene, typically a whole number or half (default is 2).", 2.0);
    public Input<PopulationModel> popModelInput = new Input<>("populationModel", "Population model used to infer the multispecies coalescent probability for this gene");

    private double ploidy;
    private int geneTreeLeafNodeCount;
    private int geneTreeNodeCount;
    private int speciesLeafNodeCount;
    private int speciesNodeCount;
    private boolean needsUpdate;
    private boolean needsLeafMapping;

    // the following are matrices associated with each branch of the species tree
    // they are flattened to arrays for optimal java performance
    protected double[] coalescentTimes; // the coalescent event times for this gene tree for all species tree branches
    protected double[] storedCoalescentTimes; // the coalescent event times for this gene tree for all species tree branches
    int coalescentTimesLength; // length of coalescentTimes array
    protected int[] coalescentCounts; // the number of coalescent events in each branch
    protected int[] storedCoalescentCounts; // stored version of coalescentCounts
    final static int DELTA_BLOCK_SIZE = 4;
    private int blocksize = DELTA_BLOCK_SIZE; // size of blocks for storing coalescentTimes, may grow (and shrink) throughout the MCMC
    int maxCoalescentCounts, storedMaxCoalescentCounts; // maximum number of coalescent events in a branch -- blocksize must always be at least as large

    protected int[] coalescentLineageCounts; // the number of lineages at the tipward end of each branch
    protected int[] storedCoalescentLineageCounts; // the number of lineages at the tipward end of each branch

    protected int[] geneNodeSpeciesAssignment;
    protected int[] storedGeneNodeSpeciesAssignment;
    protected double[] speciesOccupancy;
    protected double[] storedSpeciesOccupancy;
    protected boolean geneTreeCompatible;
    protected boolean storedGeneTreeCompatible;

    // pre-calculated lineage counts and node assignments for gene tree leaf nodes
    int[] leafCoalescentLineageCounts;
    int[] leafGeneNodeSpeciesAssignment;

    int updateCount = 0;
    boolean stopPopping = false;

    private boolean[] speciesBranchIsDirty;

    private SpeciesTreeInterface spTree;
    private Tree geneTree;
    private PopulationModel popModel;

    private double[] perBranchLogP;
    private double[] storedPerBranchLogP;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = true;
        return needsUpdate;
    }

    @Override
    public void store() {
        super.store();

        System.arraycopy(coalescentCounts, 0, storedCoalescentCounts, 0, coalescentCounts.length);
        System.arraycopy(coalescentTimes, 0, storedCoalescentTimes, 0, coalescentTimesLength);
        System.arraycopy(coalescentLineageCounts, 0, storedCoalescentLineageCounts, 0, coalescentLineageCounts.length);

        System.arraycopy(geneNodeSpeciesAssignment, 0, storedGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment.length);
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);
        System.arraycopy(perBranchLogP, 0, storedPerBranchLogP, 0, perBranchLogP.length);

        storedGeneTreeCompatible = geneTreeCompatible;
        storedMaxCoalescentCounts = maxCoalescentCounts;
    }

    @Override
    public void restore() {
        super.restore();

        double[] tmpCoalescentTimes = coalescentTimes;
        int[] tmpCoalescentCounts = coalescentCounts;
        int[] tmpCoalescentLineageCounts = coalescentLineageCounts;
        int[] tmpGeneNodeSpeciesAssignment = geneNodeSpeciesAssignment;
        double[] tmpSpeciesOccupancy = speciesOccupancy;
        double[] tmpPerBranchLogP = perBranchLogP;
        boolean tmpGeneTreeCompatible = geneTreeCompatible;

        coalescentTimes = storedCoalescentTimes;
        coalescentCounts = storedCoalescentCounts;
        coalescentLineageCounts = storedCoalescentLineageCounts;
        speciesOccupancy = storedSpeciesOccupancy;
        geneNodeSpeciesAssignment = storedGeneNodeSpeciesAssignment;
        perBranchLogP = storedPerBranchLogP;
        geneTreeCompatible = storedGeneTreeCompatible;

        storedCoalescentTimes = tmpCoalescentTimes;
        storedCoalescentCounts = tmpCoalescentCounts;
        storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
        storedSpeciesOccupancy = tmpSpeciesOccupancy;
        storedGeneNodeSpeciesAssignment = tmpGeneNodeSpeciesAssignment;
        storedPerBranchLogP = tmpPerBranchLogP;
        storedGeneTreeCompatible = tmpGeneTreeCompatible;

        maxCoalescentCounts = storedMaxCoalescentCounts;
    }

    public void initAndValidate() {
        if (popModelInput.get() != null) popModel = popModelInput.get().getBaseModel();

        ploidy = ploidyInput.get();
        geneTree = treeInput.get();
        spTree = speciesTreeInput.get();

        geneTreeLeafNodeCount = treeInput.get().getLeafNodeCount();
        geneTreeNodeCount = geneTree.getNodeCount();
        speciesLeafNodeCount = spTree.getLeafNodeCount();
        speciesNodeCount = spTree.getNodeCount();

        geneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        storedGeneNodeSpeciesAssignment = new int[geneTreeNodeCount];

        geneTreeCompatible = false;
        storedGeneTreeCompatible = false;

        coalescentLineageCounts = new int[speciesNodeCount];
        storedCoalescentLineageCounts = new int[speciesNodeCount];
        
        coalescentCounts = new int[speciesNodeCount];
        storedCoalescentCounts = new int[speciesNodeCount];
        coalescentTimesLength = speciesNodeCount * blocksize;
        coalescentTimes = new double[coalescentTimesLength + geneTreeNodeCount];
        storedCoalescentTimes = new double[coalescentTimesLength + geneTreeNodeCount];

        speciesOccupancy = new double[geneTreeNodeCount * speciesNodeCount];
        storedSpeciesOccupancy = new double[geneTreeNodeCount * speciesNodeCount];

        perBranchLogP = new double[speciesNodeCount];
        storedPerBranchLogP = new double[speciesNodeCount];

        speciesBranchIsDirty = new boolean[speciesNodeCount];
        Arrays.fill(speciesBranchIsDirty, true);

        leafCoalescentLineageCounts = new int[speciesLeafNodeCount];
        leafGeneNodeSpeciesAssignment = new int[geneTreeLeafNodeCount];

        needsUpdate = true;
        needsLeafMapping = true;

        logP = 0.0;
    }

    @Override
    public double calculateLogP() {
        assert SanityChecks.checkTreeSanity(speciesTreeInput.get().getRoot());
        if (needsUpdate) update();

        if (!geneTreeCompatible) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        logP = 0.0;
        // if using analytical integration no need to specify a population model
        if (popModel == null || popModel instanceof DummyModel) return logP;

        final Node[] speciesTreeNodes = spTree.getNodesAsArray();
        for (int nodeI = 0; nodeI < speciesNodeCount; nodeI++) {
            Node speciesNode = speciesTreeNodes[nodeI];
            if (isDirtyBranch(nodeI) || popModel.isDirtyBranch(speciesNode)) {
                final int branchLineageCount = coalescentLineageCounts[nodeI];
                final int branchEventCount = coalescentCounts[nodeI];
                final double[] branchCoalescentTimes = getCoalescentTimes(nodeI);
                perBranchLogP[nodeI] = popModel.branchLogP(nodeI, speciesNode, ploidy, branchCoalescentTimes, branchLineageCount, branchEventCount);
            }

            // System.out.println(String.format("%s-%d: %f", getID(), nodeI, logP));
            logP += perBranchLogP[nodeI];
        }

        // System.out.println(String.format("%s-%d: %f", getID(), speciesNodeCount, logP));
        return logP;
    }

    void update() {
    	synchronized (this) {
			if (needsUpdate) {
				updateCount++;

                final TreeInterface geneTree = treeInput.get();

                // generate map of species tree tip node names to node numbers
                // and count up the number of gene copies for each species
                if (needsLeafMapping) {
                    final Map<String, Integer> tipNumberMap = spTree.getTipNumberMap();

                    Arrays.fill(leafCoalescentLineageCounts, 0);

                    for (int i = 0; i < geneTreeLeafNodeCount; i++) {
                        final Node geneTreeLeafNode = geneTree.getNode(i);
                        final String geneTreeLeafName = geneTreeLeafNode.getID();
                        final int speciesTreeLeafNumber = tipNumberMap.get(geneTreeLeafName);

                        leafGeneNodeSpeciesAssignment[i] = speciesTreeLeafNumber;
                        leafCoalescentLineageCounts[speciesTreeLeafNumber]++;
                    }

                    needsLeafMapping = false;
                }

				// shrink memory reservation for coalescent times?
				if (! stopPopping &&  (updateCount & 0x7fff) == 0 && maxCoalescentCounts < blocksize - 4) {
					// ensure stored coalescent times are valid, so that a restore gives proper times
	            	double [] stmp = new double[speciesNodeCount * (blocksize - 4) + geneTreeNodeCount];
	            	for (int i = 0; i < speciesNodeCount; i++) {
	            		System.arraycopy(storedCoalescentTimes, i * blocksize, stmp, i * (blocksize - 4), blocksize - 4);
	            	}
            		System.arraycopy(stmp, 0, storedCoalescentTimes, 0, speciesNodeCount * (blocksize - 4));
	            	
					blocksize -= 4;
	            	coalescentTimesLength = speciesNodeCount * blocksize;
	            	// System.err.print("pop");
				}

                // reset arrays as these values need to be recomputed after any changes to the species or gene tree
                System.arraycopy(leafGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment, 0, geneTreeLeafNodeCount);
                System.arraycopy(leafCoalescentLineageCounts, 0, coalescentLineageCounts, 0, speciesLeafNodeCount);

                // -1 means no species assignment for that gene tree node has been made yet
                Arrays.fill(geneNodeSpeciesAssignment, geneTreeLeafNodeCount, geneTreeNodeCount, -1);
                Arrays.fill(coalescentLineageCounts, speciesLeafNodeCount, speciesNodeCount, 0);

                Arrays.fill(speciesOccupancy, 0);
                Arrays.fill(coalescentCounts, 0);
                Arrays.fill(speciesBranchIsDirty, false);

                for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
                    final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
                    final int speciesTreeLeafNumber = leafGeneNodeSpeciesAssignment[geneTreeLeafNumber];
                    final Node speciesTreeLeafNode = spTree.getNode(speciesTreeLeafNumber);
                    final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
                    final int firstCoalescenceNumber = firstCoalescenceNode.getNr();
                    final double lastHeight = 0.0;

                    if (!collateCoalescenceEvents(geneTreeLeafNumber, lastHeight,
                            firstCoalescenceNode, firstCoalescenceNumber,
                            speciesTreeLeafNode, speciesTreeLeafNumber)) {
                        // this gene tree IS NOT compatible with the species tree
                        geneTreeCompatible = false;
                        needsUpdate = false;
                        return;
                    }
                }

                maxCoalescentCounts = 0;
                for (int j : coalescentCounts) {
                    if (j > maxCoalescentCounts) {maxCoalescentCounts = j;}
                }
                if (maxCoalescentCounts > blocksize) {
                    // grow memory reservation for coalescent times
                    int DELTA_BLOCK_SIZE = 4*((maxCoalescentCounts+3)/4) - blocksize;
                    coalescentTimesLength = speciesNodeCount * (blocksize + DELTA_BLOCK_SIZE);
                    double [] tmp = new double[coalescentTimesLength + geneTreeNodeCount];
                    double [] stmp = new double[coalescentTimesLength + geneTreeNodeCount];
                    for (int i = 0; i < speciesNodeCount; i++) {
                        //System.arraycopy(coalescentTimes, i * blocksize, tmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
                        System.arraycopy(storedCoalescentTimes, i * blocksize, stmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
                    }
                    coalescentTimes = tmp;
                    storedCoalescentTimes = stmp;
                    blocksize += DELTA_BLOCK_SIZE;
                    // System.err.print("blocksize = " + blocksize + " ");

                    // do calculation again, this time with properly sized array
                    // we only get here very occasionally (only when blocksize is updated)
                    update();

                    if (updateCount > 0x7fff) {
                        stopPopping = true;
                    }
                    return;
                }

                // determine which species tree branch is dirty for this gene tree
                for (int i = 0; i < speciesNodeCount; i++) {
                    if (coalescentLineageCounts[i] != storedCoalescentLineageCounts[i] ||
                        coalescentCounts[i] != storedCoalescentCounts[i]) {
                        speciesBranchIsDirty[i] = true;
                    } else {
                        Node node = spTree.getNode(i);
                        if (node.isDirty() != Tree.IS_CLEAN ||
                            (!node.isRoot() && node.getParent().isDirty() != Tree.IS_CLEAN) ||
                            coalescentTimesChanged(i)) speciesBranchIsDirty[i] = true;
                    }
                }

                geneTreeCompatible = true;
                needsUpdate = false;
			}
    	}
    }

    private boolean coalescentTimesChanged(int i) {
    	int k = i * blocksize;
    	for (int j = 0; j < coalescentLineageCounts[i]; j++) {
    		if (coalescentTimes[k] != storedCoalescentTimes[k]) {
    			return true;
    		}
    		k++;
    	}
		return false;
	}

    // non-recursive version of recurseCoalescenceEvents
    private boolean collateCoalescenceEvents(int lastGeneTreeNodeNumber, double lastHeight, Node geneTreeNode, int geneTreeNodeNumber, Node speciesTreeNode, int speciesTreeNodeNumber) {
        while (true) {
            final double geneTreeNodeHeight = geneTreeNode.getHeight();

            // check if the next coalescence event occurs in an ancestral branch
            while (!speciesTreeNode.isRoot() && geneTreeNodeHeight >= speciesTreeNode.getParent().getHeight()) {
                /* if (geneTreeNode.isDirty() != Tree.IS_CLEAN )
                    speciesBranchIsDirty[speciesTreeNodeNumber] = true;
                */
                final Node speciesTreeParentNode = speciesTreeNode.getParent();
                final double speciesTreeParentHeight = speciesTreeParentNode.getHeight();
                final int speciesTreeParentNodeNumber = speciesTreeParentNode.getNr();

                speciesOccupancy[lastGeneTreeNodeNumber * speciesNodeCount + speciesTreeNodeNumber] = speciesTreeParentHeight - lastHeight;
                coalescentLineageCounts[speciesTreeParentNodeNumber]++;

                speciesTreeNode = speciesTreeParentNode;
                speciesTreeNodeNumber = speciesTreeParentNodeNumber;
                lastHeight = speciesTreeParentHeight;
            }

            // this code executes if the next coalescence event occurs within the current branch
            speciesOccupancy[lastGeneTreeNodeNumber * speciesNodeCount + speciesTreeNodeNumber] = geneTreeNodeHeight - lastHeight;
            final int existingSpeciesAssignment = geneNodeSpeciesAssignment[geneTreeNodeNumber];
            if (existingSpeciesAssignment == -1) {
                geneNodeSpeciesAssignment[geneTreeNodeNumber] = speciesTreeNodeNumber;

                coalescentTimes[speciesTreeNodeNumber * blocksize + coalescentCounts[speciesTreeNodeNumber]++] = geneTreeNodeHeight;

                final Node nextGeneTreeNode = geneTreeNode.getParent();
                if (nextGeneTreeNode == null) {
                    // this is the root of the gene tree and no incompatibilities were detected
                    return true;
                } else {
                    // if this is not the root of the gene tree, check the subsequent (back in time) coalescence event
                    lastGeneTreeNodeNumber = geneTreeNodeNumber;
                    lastHeight = geneTreeNodeHeight;
                    geneTreeNode = nextGeneTreeNode;
                    geneTreeNodeNumber = nextGeneTreeNode.getNr();
                    //return collateCoalescenceEvents(geneTreeNodeNumber, geneTreeNodeHeight, nextGeneTreeNode, nextGeneTreeNodeNumber, speciesTreeNode, speciesTreeNodeNumber);
                }
            } else if (existingSpeciesAssignment == speciesTreeNodeNumber) {
                return true; // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
            } else {
                return false; // this gene tree IS NOT compatible with the species tree
            }
        }
    }

    public double[] getSpeciesOccupancy() {
        if (needsUpdate) update();

        return speciesOccupancy;
    }

	public double[] getCoalescentTimes(int nodeI) {
        if (needsUpdate) update();

        final Node speciesNode = spTree.getNode(nodeI);
        final Node parentNode = speciesNode.getParent();

        final double speciesEndTime = speciesNode.getHeight();
        final double speciesStartTime = (parentNode == null) ? Double.POSITIVE_INFINITY : parentNode.getHeight();
        final int branchEventCount = coalescentCounts[nodeI];

		final double[] branchCoalescentTimes = new double[branchEventCount + 2];
		branchCoalescentTimes[0] = speciesEndTime;
        branchCoalescentTimes[branchEventCount + 1] = speciesStartTime;

		System.arraycopy(coalescentTimes, nodeI * blocksize, branchCoalescentTimes, 1, branchEventCount);
		Arrays.sort(branchCoalescentTimes);

		return branchCoalescentTimes;
	}

	protected boolean isDirtyBranch(int nodeNr) {
		return speciesBranchIsDirty[nodeNr];
	}

	int[] getTipNumberMap() {
		return leafGeneNodeSpeciesAssignment;
	}

    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO Auto-generated method stub
        
    }

    public int getNodeCount() {
        return geneTreeNodeCount;
    }

    public Node getRoot() {
        return treeInput.get().getRoot();
    }

    public double getPloidy() {
        return ploidy;
    }
}
