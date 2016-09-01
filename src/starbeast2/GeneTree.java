package starbeast2;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
 */

public class GeneTree extends TreeWrapper {
    public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "Species tree for embedding the gene tree.");
    public Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for this gene, typically a whole number or half (default is 2).", 2.0);
    protected double ploidy;

    private int geneTreeLeafNodeCount;
    private int geneTreeNodeCount;
    private int speciesTreeNodeCount;
    private boolean needsUpdate;

    // the following are matrices associated with each branch of the species tree
    // they are flattened to arrays for optimal java performance
    protected double[] coalescentTimes; // the coalescent event times for this gene tree for all species tree branches
    protected double[] storedCoalescentTimes; // the coalescent event times for this gene tree for all species tree branches
    int coalescentTimesLength; // lenght of coalescentTimes array
    protected int [] coalescentCounts; // the number of coalescent events in each branch
    protected int [] storedCoalescentCounts; // stored version of coalescentCounts
    final static int DELTA_BLOCK_SIZE = 4;
    private int blocksize = DELTA_BLOCK_SIZE; // size of blocks for storing coalescentTimes, may grow (and shrink) throughout the MCMC
    int maxCoalescentCounts, storedMaxCoalescentCounts; // maximum number of coalescent events in a branch -- blocksize must always be at least as large
    
    
    protected int [] coalescentLineageCounts; // the number of lineages at the tipward end of each branch
    protected int [] storedCoalescentLineageCounts; // the number of lineages at the tipward end of each branch

    protected int[] geneNodeSpeciesAssignment;
    protected int[] storedGeneNodeSpeciesAssignment;
    protected double[] speciesOccupancy;
    protected double[] storedSpeciesOccupancy;
    protected boolean geneTreeCompatible;
    protected boolean storedGeneTreeCompatible;

    // pre-calculated lineage counts and node assignments for gene tree leaf nodes
    int [] leafCoalescentLineageCounts;
    int [] leafGeneNodeSpeciesAssignment;

    
    int updateCount = 0;
    boolean stopPopping = false;
   
    // maps gene tree tip numbers to species tree tip number
    private int [] localTipNumberMap;

    private boolean [] speciesBranchIsDirty;
    //private boolean [] storedSpeciesBranchIsDirty;

    Tree spTree, geneTree;
    //List<Node> dirtyList = new ArrayList<>();
    
    @Override
    public boolean requiresRecalculation() {
        needsUpdate = spTree.somethingIsDirty() || geneTree.somethingIsDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
    }

    private void doStore() {
        System.arraycopy(coalescentCounts, 0, storedCoalescentCounts, 0, coalescentCounts.length);
        System.arraycopy(coalescentTimes, 0, storedCoalescentTimes, 0, coalescentTimesLength);
        System.arraycopy(coalescentLineageCounts, 0, storedCoalescentLineageCounts, 0, coalescentLineageCounts.length);

        System.arraycopy(geneNodeSpeciesAssignment, 0, storedGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment.length);
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);
        
        storedGeneTreeCompatible = geneTreeCompatible;
        storedMaxCoalescentCounts = maxCoalescentCounts;
        
        super.store();
    }

    @Override
    public void restore() {
    	if (needsUpdate) {
    		return;
    	}
        double [] tmpCoalescentTimes = coalescentTimes;
        int [] tmpCoalescentCounts = coalescentCounts;
        int[] tmpCoalescentLineageCounts = coalescentLineageCounts;
        int[] tmpGeneNodeSpeciesAssignment = geneNodeSpeciesAssignment;
        double[] tmpSpeciesOccupancy = speciesOccupancy;
        boolean tmpGeneTreeCompatible = geneTreeCompatible;

        coalescentTimes = storedCoalescentTimes;
        coalescentCounts = storedCoalescentCounts;
        coalescentLineageCounts = storedCoalescentLineageCounts;
        speciesOccupancy = storedSpeciesOccupancy;
        geneNodeSpeciesAssignment = storedGeneNodeSpeciesAssignment;
        geneTreeCompatible = storedGeneTreeCompatible;

        storedCoalescentTimes = tmpCoalescentTimes;
        storedCoalescentCounts = tmpCoalescentCounts;
        storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
        storedSpeciesOccupancy = tmpSpeciesOccupancy;
        storedGeneNodeSpeciesAssignment = tmpGeneNodeSpeciesAssignment;
        storedGeneTreeCompatible = tmpGeneTreeCompatible;

        maxCoalescentCounts = storedMaxCoalescentCounts;
        
        super.restore();
    }

    
    public void initAndValidate() {
        ploidy = ploidyInput.get();

        geneTree = treeInput.get();
        geneTreeNodeCount = geneTree.getNodeCount();
        geneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        storedGeneNodeSpeciesAssignment = new int[geneTreeNodeCount];

        geneTreeLeafNodeCount = treeInput.get().getLeafNodeCount();

        // generate map of species tree tip node names to node numbers
        final SpeciesTree speciesTree = speciesTreeInput.get();
        spTree = speciesTree.treeInput.get();

        
        final SpeciesTree speciesTreeWrapper = speciesTreeInput.get();
        final Map<String, Integer> tipNumberMap = speciesTreeWrapper.getTipNumberMap();
        TreeInterface geneTree = treeInput.get();
        localTipNumberMap = new int[geneTree.getLeafNodeCount()];
        for (int i = 0; i < geneTree.getLeafNodeCount(); i++) {
        	Node geneTreeLeafNode = geneTree.getNode(i);
        	localTipNumberMap[geneTreeLeafNode.getNr()] = tipNumberMap.get(geneTreeLeafNode.getID());
        }

        geneTreeCompatible = false;
        storedGeneTreeCompatible = false;

        needsUpdate = true;
        speciesTreeNodeCount = speciesTree.treeInput.get().getNodeCount();
        coalescentLineageCounts = new int[speciesTreeNodeCount];
        storedCoalescentLineageCounts = new int[speciesTreeNodeCount];
        
        coalescentCounts = new int[speciesTreeNodeCount];
        storedCoalescentCounts = new int[speciesTreeNodeCount];
        coalescentTimesLength = speciesTreeNodeCount * blocksize;
        coalescentTimes = new double[coalescentTimesLength + geneTreeNodeCount];
        storedCoalescentTimes = new double[coalescentTimesLength + geneTreeNodeCount];

        speciesOccupancy = new double[geneTreeNodeCount * speciesTreeNodeCount];
        storedSpeciesOccupancy = new double[geneTreeNodeCount * speciesTreeNodeCount];
        
        speciesBranchIsDirty = new boolean[speciesTreeNodeCount];
        
        
        leafCoalescentLineageCounts = new int[speciesTreeNodeCount];
        leafGeneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        Arrays.fill(leafGeneNodeSpeciesAssignment, -1);
        
        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
            final int speciesTreeLeafNumber = localTipNumberMap[geneTreeLeafNode.getNr()];
            leafCoalescentLineageCounts[speciesTreeLeafNumber]++;            
            leafGeneNodeSpeciesAssignment[geneTreeLeafNumber] = speciesTreeLeafNumber;
        }
    }

    protected boolean computeCoalescentTimes() {
        if (needsUpdate) {
            update();
        }

        // the number of coalescent times should equal the number of internal gene tree nodes (each a coalescent event)
        if (geneTreeCompatible) {
            //assert coalescentTimes.size() == geneTreeNodeCount - geneTreeLeafNodeCount;
            // this gene tree IS compatible with the species tree
            return true;
        } else {
            return false;
        }
    }

    void update() {
    	synchronized (this) {
			if (needsUpdate) {
				updateCount++;

				// shrink memory reservation for coalescent times?
				if (! stopPopping &&  (updateCount & 0x7fff) == 0 && maxCoalescentCounts < blocksize - 4) {
					// ensure stored coalescent times are valid, so that a restore gives proper times
	            	double [] stmp = new double[speciesTreeNodeCount * (blocksize - 4) + geneTreeNodeCount];
	            	for (int i = 0; i < speciesTreeNodeCount; i++) {
	            		System.arraycopy(storedCoalescentTimes, i * blocksize, stmp, i * (blocksize - 4), blocksize - 4);
	            	}
            		System.arraycopy(stmp, 0, storedCoalescentTimes, 0, speciesTreeNodeCount * (blocksize - 4));
	            	
					blocksize -= 4;
	            	coalescentTimesLength = speciesTreeNodeCount * blocksize;
	            	System.err.print("pop");
				}
	    	doStore();

	        Arrays.fill(speciesOccupancy, 0);
	        
	        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
	        //Arrays.fill(geneNodeSpeciesAssignment, -1); // -1 means no species assignment for that gene tree node has been made yet
	        System.arraycopy(leafGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment, 0, geneTreeNodeCount);
	
	        
	        // Arrays.fill(coalescentLineageCounts, 0);
	        System.arraycopy(leafCoalescentLineageCounts, 0, coalescentLineageCounts, 0, speciesTreeNodeCount);
	        Arrays.fill(coalescentCounts, 0);
	        
        	Arrays.fill(speciesBranchIsDirty, false);
	
	        final TreeInterface geneTree = treeInput.get();
	        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
	            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
	            final int speciesTreeLeafNumber = localTipNumberMap[geneTreeLeafNode.getNr()];
	            final Node speciesTreeLeafNode = spTree.getNode(speciesTreeLeafNumber);	
	            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
	            final int firstCoalescenceNumber = firstCoalescenceNode.getNr();
	            final double lastHeight = 0.0;
	
	            if (!recurseCoalescenceEvents(geneTreeLeafNumber, lastHeight, 
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
            	coalescentTimesLength = speciesTreeNodeCount * (blocksize + DELTA_BLOCK_SIZE);
            	double [] tmp = new double[coalescentTimesLength + geneTreeNodeCount];
            	double [] stmp = new double[coalescentTimesLength + geneTreeNodeCount];
            	for (int i = 0; i < speciesTreeNodeCount; i++) {
            		//System.arraycopy(coalescentTimes, i * blocksize, tmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
            		System.arraycopy(storedCoalescentTimes, i * blocksize, stmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
            	}
            	coalescentTimes = tmp;
            	storedCoalescentTimes = stmp;
            	blocksize += DELTA_BLOCK_SIZE;
            	System.err.print("blocksize = " + blocksize + " ");
            	
            	// do calculation again, this time with properly sized array
            	// we only get here very occasionally (only when blocksize is updated)
            	update();
            	
            	if (updateCount > 0x7fff) {
            		stopPopping = true;
            	}
            	return;
            }

	        
            // determine which species tree branch is dirty for this gene tree
	        for (int i = 0; i < speciesTreeNodeCount; i++) {
	        	if (coalescentLineageCounts[i] != storedCoalescentCounts[i] ||
	        		coalescentCounts[i] != storedCoalescentCounts[i]) {
	        		speciesBranchIsDirty[i] = true;
	        	} else {
	        		Node node = spTree.getNode(i);
	        		if (node.isDirty() != Tree.IS_CLEAN ||
	        			(node.isRoot() && node.getParent().isDirty() != Tree.IS_CLEAN)) {
		        		speciesBranchIsDirty[i] = true;
	        		} else if (coalescentTimesChanged(i)) {
		        		speciesBranchIsDirty[i] = true;
	        		}

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

    // non-recursive version of recurseCoalescenceEventsOrg
    private boolean recurseCoalescenceEvents(int lastGeneTreeNodeNumber, double lastHeight, Node geneTreeNode, int geneTreeNodeNumber, Node speciesTreeNode, int speciesTreeNodeNumber) {
        while (true) {
	    	final double geneTreeNodeHeight = geneTreeNode.getHeight();
	
	        // check if the next coalescence event occurs in an ancestral branch
	    	while (!speciesTreeNode.isRoot() && geneTreeNodeHeight >= speciesTreeNode.getParent().getHeight()) {
//	            if (geneTreeNode.isDirty() != Tree.IS_CLEAN ) {
//	            	speciesBranchIsDirty[speciesTreeNodeNumber] = true;
//	            }
    			final Node speciesTreeParentNode = speciesTreeNode.getParent();
                speciesOccupancy[lastGeneTreeNodeNumber * speciesTreeNodeCount + speciesTreeNodeNumber] = speciesTreeParentNode.getHeight() - lastHeight;
                final int speciesTreeParentNodeNumber = speciesTreeParentNode.getNr();
                coalescentLineageCounts[speciesTreeParentNodeNumber]++;
                speciesTreeNode = speciesTreeParentNode;
                speciesTreeNodeNumber = speciesTreeParentNodeNumber;
           }
	
	        // this code executes if the next coalescence event occurs within the current branch
	        speciesOccupancy[lastGeneTreeNodeNumber * speciesTreeNodeCount + speciesTreeNodeNumber] = geneTreeNodeHeight - lastHeight;
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
	                //return recurseCoalescenceEvents(geneTreeNodeNumber, geneTreeNodeHeight, nextGeneTreeNode, nextGeneTreeNodeNumber, speciesTreeNode, speciesTreeNodeNumber);
	            }
	        } else if (existingSpeciesAssignment == speciesTreeNodeNumber) {
	            return true; // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
	        } else {
	            return false; // this gene tree IS NOT compatible with the species tree
	        }
	    }
    }
    

    public double[] getSpeciesOccupancy() {
        if (needsUpdate) {
            update();
        }

        return speciesOccupancy;
    }

	public double[] getCoalescentTimes(int i) {
		double [] geneBranchCoalescentTimes = new double[coalescentCounts[i]];
		System.arraycopy(coalescentTimes, i * blocksize, geneBranchCoalescentTimes, 0, geneBranchCoalescentTimes.length);
		Arrays.sort(geneBranchCoalescentTimes);
		return geneBranchCoalescentTimes;
	}

	protected boolean isDirtyBranch(int nodeNr) {
		return speciesBranchIsDirty[nodeNr];
	}

	int [] getTipNumberMap() {
		return localTipNumberMap;
	}
	
    /* public double[] getOccupancy(Node node) {
        if (needsUpdate) {
            update();
        }

        final int geneTreeNodeNumber = node.getNr();
        return speciesOccupancy[geneTreeNodeNumber];
    } */
}
