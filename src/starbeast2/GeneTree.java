package starbeast2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
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

    protected double[] coalescentTimes; // the coalescent event times for this gene tree for all species tree branches
    protected double[] storedCoalescentTimes; // the coalescent event times for this gene tree for all species tree branches
    protected int [] coalescentCounts; // the number of coalescent events in each branch
    protected int [] storedCoalescentCounts; // stored version of coalescentCounts
    final static int DELTA_BLOCK_SIZE = 4;
    private int blocksize = DELTA_BLOCK_SIZE; // size of blocks for storing coalescentTimes, may grow throughout the MCMC
    
    protected int [] coalescentLineageCounts; // the number of lineages at the tipward end of each branch
    protected int [] storedCoalescentLineageCounts; // the number of lineages at the tipward end of each branch

    protected int[] geneNodeSpeciesAssignment;
    protected int[] storedGeneNodeSpeciesAssignment;
    protected double[] speciesOccupancy;
    protected double[] storedSpeciesOccupancy;
    protected boolean geneTreeCompatible;
    protected boolean storedGeneTreeCompatible;
    
    // maps gene tree tip numbers to species tree tip number
    private int [] localTipNumberMap;

    private boolean [] speciesBranchIsDirty;
    //private boolean [] storedSpeciesBranchIsDirty;

    Tree spTree;
    //List<Node> dirtyList = new ArrayList<>();
    
    @Override
    public boolean requiresRecalculation() {
        needsUpdate = treeInput.isDirty() || speciesTreeInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
    }

    private void doStore() {
        System.arraycopy(coalescentCounts, 0, storedCoalescentCounts, 0, coalescentCounts.length);
        System.arraycopy(coalescentTimes, 0, storedCoalescentTimes, 0, coalescentTimes.length);
        System.arraycopy(coalescentLineageCounts, 0, storedCoalescentLineageCounts, 0, coalescentLineageCounts.length);

        //storedSpeciesOccupancy = new double[speciesOccupancy.length][speciesOccupancy[0].length];
        System.arraycopy(geneNodeSpeciesAssignment, 0, storedGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment.length);
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);
        
        //System.arraycopy(speciesBranchIsDirty, 0, storedSpeciesBranchIsDirty, 0, speciesBranchIsDirty.length);

        storedGeneTreeCompatible = geneTreeCompatible;

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

        //boolean [] tmp = speciesBranchIsDirty;
        //speciesBranchIsDirty = storedSpeciesBranchIsDirty;
        //storedSpeciesBranchIsDirty = tmp;
        
        super.restore();
    }

    public void initAndValidate() {
        ploidy = ploidyInput.get();

        geneTreeNodeCount = treeInput.get().getNodeCount();
        geneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        storedGeneNodeSpeciesAssignment = new int[geneTreeNodeCount];

        geneTreeLeafNodeCount = treeInput.get().getLeafNodeCount();

        // generate map of species tree tip node names to node numbers
        final SpeciesTree speciesTree = speciesTreeInput.get();
        final HashMap<String, Integer> speciesNumberMap = new HashMap<>();

        Node speciesTreeRoot = speciesTree.getRoot();
        spTree = speciesTree.treeInput.get();
        for (Node leafNode: speciesTreeRoot.getAllLeafNodes()) {
            final String speciesName = leafNode.getID();
            final int speciesNumber = leafNode.getNr();

            speciesNumberMap.put(speciesName, speciesNumber);
        }

        
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
        coalescentTimes = new double[speciesTreeNodeCount * blocksize];
        storedCoalescentTimes = new double[speciesTreeNodeCount * blocksize];

        speciesOccupancy = new double[geneTreeNodeCount * speciesTreeNodeCount];
        storedSpeciesOccupancy = new double[geneTreeNodeCount * speciesTreeNodeCount];
        
        speciesBranchIsDirty = new boolean[speciesTreeNodeCount];
        //storedSpeciesBranchIsDirty = new boolean[speciesTreeNodeCount];
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
	    	doStore();
	        final SpeciesTree speciesTreeWrapper = speciesTreeInput.get();
	        final TreeInterface speciesTree = speciesTreeWrapper.getTree();
	        //final Map<String, Integer> tipNumberMap = speciesTreeWrapper.getTipNumberMap();
	
	//        final int speciesTreeNodeCount = speciesTree.getNodeCount();
	        Arrays.fill(speciesOccupancy, 0);
	        
	        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
	        Arrays.fill(geneNodeSpeciesAssignment, -1); // -1 means no species assignment for that gene tree node has been made yet
	
	        Arrays.fill(coalescentLineageCounts, 0);
	        Arrays.fill(coalescentCounts, 0);
	        
        	Arrays.fill(speciesBranchIsDirty, false);
	        if (speciesTreeInput.isDirty()) {
		        for (int i = 0; i < speciesTreeNodeCount; i++) {
		        	final Node node = spTree.getNode(i);
		        	if (node.isDirty() != Tree.IS_CLEAN) {
		        		speciesBranchIsDirty[i] = true;
		        		for (Node child : node.getChildren()) {
			        		speciesBranchIsDirty[child.getNr()] = true;
		        		}
		        	}
		        }
	        }
	        //dirtyList.clear();
	
	        final TreeInterface geneTree = treeInput.get();
	        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
	            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
	            // final int speciesTreeLeafNumber = tipNumberMap.get(geneTreeLeafNode.getID());
	            final int speciesTreeLeafNumber = localTipNumberMap[geneTreeLeafNode.getNr()];
	            final Node speciesTreeLeafNode = speciesTree.getNode(speciesTreeLeafNumber);
	            coalescentLineageCounts[speciesTreeLeafNumber]++;
	            
	            geneNodeSpeciesAssignment[geneTreeLeafNumber] = speciesTreeLeafNumber;
	
	            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
	            final int firstCoalescenceNumber = firstCoalescenceNode.getNr();
	            final double lastHeight = 0.0;
	
	            if (!recurseCoalescenceEvents(geneTreeLeafNumber, lastHeight, firstCoalescenceNode, firstCoalescenceNumber, speciesTreeLeafNode, speciesTreeLeafNumber)) {
	                // this gene tree IS NOT compatible with the species tree
	                geneTreeCompatible = false;
	                needsUpdate = false;
	                return;
	            }
	        }
	
	        for (int i = 0; i < speciesTreeNodeCount; i++) {
	        	if (coalescentLineageCounts[i] != storedCoalescentCounts[i]) {
	        		speciesBranchIsDirty[i] = true;
	        	}
	        }
//        	for (Node geneTreeNode : dirtyList) {
//            	// deal with parent node
//            	Node node = spTree.getNode(geneNodeSpeciesAssignment[geneTreeNode.getNr()]);
//            	if (!geneTreeNode.isRoot()) {
//	            	double h = geneTreeNode.getParent().getHeight();
//	            	while (node != null && h > node.getHeight()) {
//		            	speciesBranchIsDirty[node.getNr()] = true;
//		            	node = node.getParent();
//	            	}
//            	}
//
//            	// deal with child nodes
////            	for (Node child : geneTreeNode.getChildren()) {
////            		double h = child.getHeight();
////            		node = spTree.getNode(geneNodeSpeciesAssignment[child.getNr()]);
////	            	while (node != null && h > node.getHeight()) {
////		            	speciesBranchIsDirty[node.getNr()] = true;
////		            	node = node.getParent();
////	            	}
////	            }
//            }

	        geneTreeCompatible = true;
	        needsUpdate = false;
			}
    	}
    }

    private boolean recurseCoalescenceEventsOrg(final int lastGeneTreeNodeNumber, final double lastHeight, final Node geneTreeNode, final int geneTreeNodeNumber, final Node speciesTreeNode, final int speciesTreeNodeNumber) {
    	final double geneTreeNodeHeight = geneTreeNode.getHeight();

        // check if the next coalescence event occurs in an ancestral branch
        if (!speciesTreeNode.isRoot()) {
            final Node speciesTreeParentNode = speciesTreeNode.getParent();
            final double speciesTreeParentHeight = speciesTreeParentNode.getHeight();
            if (geneTreeNodeHeight >= speciesTreeParentHeight) {
                speciesOccupancy[lastGeneTreeNodeNumber * speciesTreeNodeCount + speciesTreeNodeNumber] = speciesTreeParentHeight - lastHeight;
                final int speciesTreeParentNodeNumber = speciesTreeParentNode.getNr();
                coalescentLineageCounts[speciesTreeParentNodeNumber]++;
                return recurseCoalescenceEvents(lastGeneTreeNodeNumber, speciesTreeParentHeight, geneTreeNode, geneTreeNodeNumber, speciesTreeParentNode, speciesTreeParentNodeNumber);
            }
        }

        // this code executes if the next coalescence event occurs within the current branch
        speciesOccupancy[lastGeneTreeNodeNumber * speciesTreeNodeCount + speciesTreeNodeNumber] = geneTreeNodeHeight - lastHeight;
        final int existingSpeciesAssignment = geneNodeSpeciesAssignment[geneTreeNodeNumber];
        if (existingSpeciesAssignment == -1) {
            geneNodeSpeciesAssignment[geneTreeNodeNumber] = speciesTreeNodeNumber;
            
            if (coalescentCounts[speciesTreeNodeNumber] == blocksize) {
            	// grow memory reservation for coalescent times
            	int speciesTreeNodeCount = coalescentTimes.length / blocksize;
            	double [] tmp = new double[speciesTreeNodeCount * (blocksize + DELTA_BLOCK_SIZE)];
            	double [] stmp = new double[speciesTreeNodeCount * (blocksize + DELTA_BLOCK_SIZE)];
            	for (int i = 0; i < speciesTreeNodeCount; i++) {
            		System.arraycopy(coalescentTimes, i * blocksize, tmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
            		System.arraycopy(storedCoalescentTimes, i * blocksize, stmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
            	}
            	coalescentTimes = tmp;
            	storedCoalescentTimes = stmp;
            	blocksize += DELTA_BLOCK_SIZE;
            	System.err.print("blocksize = " + blocksize + " ");
            }
            coalescentTimes[speciesTreeNodeNumber * blocksize + coalescentCounts[speciesTreeNodeNumber]++] = geneTreeNodeHeight;
            
            final Node nextGeneTreeNode = geneTreeNode.getParent();
            if (nextGeneTreeNode == null) {
                // this is the root of the gene tree and no incompatibilities were detected
                return true;
            } else {
                // if this is not the root of the gene tree, check the subsequent (back in time) coalescence event
                final int nextGeneTreeNodeNumber = nextGeneTreeNode.getNr();               
                return recurseCoalescenceEvents(geneTreeNodeNumber, geneTreeNodeHeight, nextGeneTreeNode, nextGeneTreeNodeNumber, speciesTreeNode, speciesTreeNodeNumber);
            }
        } else if (existingSpeciesAssignment == speciesTreeNodeNumber) {
            return true; // gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
        } else {
            return false; // this gene tree IS NOT compatible with the species tree
        }
    }
    
    // non-recursive version of recurseCoalescenceEventsOrg
    private boolean recurseCoalescenceEvents(int lastGeneTreeNodeNumber, double lastHeight, Node geneTreeNode, int geneTreeNodeNumber, Node speciesTreeNode, int speciesTreeNodeNumber) {
        while (true) {
	    	final double geneTreeNodeHeight = geneTreeNode.getHeight();
	
	        // check if the next coalescence event occurs in an ancestral branch
	    	while(!speciesTreeNode.isRoot() && geneTreeNodeHeight >= speciesTreeNode.getParent().getHeight()) {
	            if (geneTreeNode.isDirty() != Tree.IS_CLEAN ) {
	            	speciesBranchIsDirty[speciesTreeNodeNumber] = true;
	            }
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
	            
            	// grow memory reservation for coalescent times?
	            if (coalescentCounts[speciesTreeNodeNumber] == blocksize) {
	            	int speciesTreeNodeCount = coalescentTimes.length / blocksize;
	            	double [] tmp = new double[speciesTreeNodeCount * (blocksize + DELTA_BLOCK_SIZE)];
	            	double [] stmp = new double[speciesTreeNodeCount * (blocksize + DELTA_BLOCK_SIZE)];
	            	for (int i = 0; i < speciesTreeNodeCount; i++) {
	            		System.arraycopy(coalescentTimes, i * blocksize, tmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
	            		System.arraycopy(storedCoalescentTimes, i * blocksize, stmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
	            	}
	            	coalescentTimes = tmp;
	            	storedCoalescentTimes = stmp;
	            	blocksize += DELTA_BLOCK_SIZE;
	            	System.err.print("blocksize = " + blocksize + " ");
	            }
	            coalescentTimes[speciesTreeNodeNumber * blocksize + coalescentCounts[speciesTreeNodeNumber]++] = geneTreeNodeHeight;
	            
	            if (geneTreeNode.isDirty() != Tree.IS_CLEAN) {
	            	speciesBranchIsDirty[speciesTreeNodeNumber] = true;
//	            	dirtyList.add(geneTreeNode);
	            	
	            	// deal with parent node
	            	Node node = speciesTreeNode;
	            	if (!geneTreeNode.isRoot()) {
		            	double h = geneTreeNode.getParent().getHeight();
		            	while (node != null && h > node.getHeight()) {
			            	speciesBranchIsDirty[node.getNr()] = true;
			            	node = node.getParent();
		            	}
	            	}

	            }
	            
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

    /* public double[] getOccupancy(Node node) {
        if (needsUpdate) {
            update();
        }

        final int geneTreeNodeNumber = node.getNr();
        return speciesOccupancy[geneTreeNodeNumber];
    } */
}
