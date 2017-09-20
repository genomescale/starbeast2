/*
 * Copyright (C) 2016 Nicola Felix Mueller (nicola.felix.mueller@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package starbeast2;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

import java.text.DecimalFormat;
import java.util.*;


/**
 * @author Nicola Felix Mueller
 */

@Description("Species tree that contains the information about things such as the speciation times and migration rates")
public class ConstantWithGeneFlow extends CalculationNode implements PopulationModel {
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.", Validate.REQUIRED);
    public Input<RealParameter> NeInput = new Input<RealParameter>("Ne","contains the Ne of each branch", Validate.REQUIRED);
    public Input<RealParameter> mInput  = new Input<>("m","relative migration rates between branches", Validate.REQUIRED);
    public Input<RealParameter> mClockInput  = new Input<>("mClock","absolute migration rates", Validate.REQUIRED);
    
    private SpeciesTreeInterface speciesTree;
    
    private boolean needsUpdate;
    private int leafNodeCount;
    
	public int nrSamples;
	
	// save all the species (sample) names
	private HashMap<String, Integer> speciesNames;
	
	// maps states to node numbers
	public ArrayList<ArrayList<Integer>> stateToNodeMap;
	private ArrayList<ArrayList<Integer>> storedStateToNodeMap;
	
	// maps branches migration rates
	private ArrayList<Integer[]> migrationMap;
	private ArrayList<Integer[]> storedMigrationMap;
	
	//  get the indices of the daughter lineages at the last coalescent event
	private int daughterIndex1, daughterIndex2;	

	
    protected double[] intervals;
    protected double[] storedIntervals;
    protected boolean[] isCoalescent;
    protected boolean[] intervalIsDirty;
    private boolean lastIntervalDirty;

    final static NodeHeightComparator nhc = new NodeHeightComparator();
    
    /**
     * parent and daugheter lineages at each coalescent event
     */
    protected List<Node>[] lineagesAdded;
    protected List<Node>[] lineagesRemoved;
    protected List<Node>[] storedLineagesAdded;
    protected List<Node>[] storedLineagesRemoved;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = true;
        return needsUpdate;
    }
    
    
    protected boolean checkDirty(){
    	// check if tree is dirty
		calculateIntervals();    			
    	for (int i = 0; i < intervals.length; i++){
    		if (speciesTreeInput.get().getNode(i).isDirty()>0){
    			needsUpdate = true;
    			stateToNodeMap();
    			return true;
    		}
    	}
   	
    	boolean isDirty = false;
    	
    	
    	// set all intervals clean. This point is only reached if only migration or
    	// Ne is changed
    	setIntervalsClean();  	

		// check what is dirty
		for (int i = 0; i < NeInput.get().getDimension(); i++){
			if (NeInput.get().isDirty(i)){
				isDirty = true;
				for (int a = 0; a < stateToNodeMap.size(); a++){
					for (int b = 0; b < stateToNodeMap.get(a).size(); b++ )
						if (i==stateToNodeMap.get(a).get(b))
							setIntervalIsDirty(nrSamples+a);
				}
			}
		}
		
		// TODO adjust for symmetric and all the same rates
		// TODO make sure the correct interval is used
		int diffMigCount = nrSamples*(nrSamples-1);
		int intCount = 0;
		int counti = 0;
		for (int i = 0; i < mInput.get().getDimension(); i++){
			if (counti == diffMigCount){
				counti = 0;
				intCount++;
				diffMigCount = (nrSamples-intCount)*(nrSamples-intCount-1) - (nrSamples-intCount-1)*(nrSamples-intCount-2);
			}
			if (mInput.get().isDirty(i)){
				isDirty = true;
				setIntervalIsDirty(nrSamples + intCount);
			}
			counti++;
		}		
		
		return isDirty;
    }

    

    @Override
    public void initAndValidate() {
        speciesTree = speciesTreeInput.get();
        final int speciesNodeCount = speciesTree.getNodeCount();
        leafNodeCount = speciesTree.getLeafNodeCount(); // also the number of "tip" population sizes
        NeInput.get().setDimension(2*leafNodeCount-1);
        mInput.get().setDimension((leafNodeCount-1)*(leafNodeCount-1)*2);
        
        needsUpdate = true;
        
    	// Calculate the tree intervals (time between events, which nodes participate at a event etc.)
        calculateIntervals();        
        // initialize the state to node mapping
        stateToNodeMap();     
    }
    
    
    @SuppressWarnings({ "unchecked", "deprecation" })
    protected void calculateIntervals() {
    	Node[] speciesNodesTmp = speciesTree.getNodesAsArray();
    	Node[] speciesNodes = new Node[speciesNodesTmp.length];
    	System.arraycopy(speciesNodesTmp, 0, speciesNodes, 0, speciesNodesTmp.length);

    	
    	Arrays.sort(speciesNodes, nhc);
    	
    	
 	
        intervals = new double[speciesNodes.length];
        isCoalescent = new boolean[speciesNodes.length];
        lineagesAdded = new List[speciesNodes.length];
        lineagesRemoved = new List[speciesNodes.length];
        
        ArrayList<Node> initAL = new ArrayList<>(); 
        intervals[0] = speciesNodes[0].getHeight();
        lineagesAdded[0] = new ArrayList<>(); 
		lineagesAdded[0].add(speciesNodes[0]);
		isCoalescent[0] = false;
    	for (int i = 1; i < speciesNodes.length; i++){
    		intervals[i] = speciesNodes[i].getHeight()-speciesNodes[i-1].getHeight();
            lineagesAdded[i] = new ArrayList<>(); 
    		lineagesAdded[i].add(speciesNodes[i]);
    		if (!speciesNodes[i].isLeaf()){
    			isCoalescent[i] = true;
    			lineagesRemoved[i] = new ArrayList<>(); 
	    		lineagesRemoved[i].add(speciesNodes[i].getLeft());
	    		lineagesRemoved[i].add(speciesNodes[i].getRight());  
    		}else{
    			isCoalescent[i] = false;
    			lineagesRemoved[i] = new ArrayList<>(); 
	    		lineagesRemoved[i].add(null);
	    		lineagesRemoved[i].add(null);    		  			
    		}
    	}   
    	stateToNodeMap();
    }
    

	// builds the map from state to node number
    public void stateToNodeMap(){
    	stateToNodeMap = new ArrayList<>();
    	ArrayList<Integer> activeStates = new ArrayList<>();
	
    	boolean firstCoalescent = true;
    	for (int i = 0; i < intervals.length; i++){
    		if (isCoalescent[i]){
    			if(firstCoalescent){
    				// sort sampled lineages
    				Collections.sort(activeStates);
    				ArrayList<Integer> addList = new ArrayList<>(activeStates);
    				stateToNodeMap.add(addList);
    				firstCoalescent = false;
    			}
			
    			// get the added and removed lineages
    			List<Node> lins = lineagesRemoved[i];
    			activeStates.remove(activeStates.indexOf(lins.get(0).getNr()));
    			activeStates.remove(activeStates.indexOf(lins.get(1).getNr()));
    			
    			if (lins.get(1).getParent().getNr() != lins.get(0).getParent().getNr()){
    				System.err.println("lineages don't have the same parent");
    				System.exit(0);
    			}
    			if (lins.get(1).getParent().getNr() != lineagesAdded[i].get(0).getNr()){
    				System.err.println("wrong lineage added");  				
    			}
    			
    			activeStates.add(lins.get(1).getParent().getNr());
				ArrayList<Integer> addList = new ArrayList<>(activeStates);
				stateToNodeMap.add(addList);
    		}else{
    			activeStates.add(lineagesAdded[i].get(0).getNr());
    		}
    	}
    	migrationMap = new ArrayList<>();
    	
    	// get the migration map for each interval
    	activeStates = new ArrayList<>();
    	int a = 0;
    	while (!isCoalescent[a]){
			activeStates.add(lineagesAdded[a].get(0).getNr());
    		a++;
    	}    		
    	
    	for (int i = 0; i < activeStates.size(); i++){
    		for (int j = 0; j < activeStates.size(); j++){
    			if (i!=j){
    				Integer[] migroute = new Integer[]{activeStates.get(i), activeStates.get(j)};
    				migrationMap.add(migroute);
    			}
    		}
    	}
    	
    	
    	while (a < intervals.length){
    		int parent = lineagesAdded[a].get(0).getNr();
    		activeStates.add(parent);
    		int d1 = activeStates.indexOf(lineagesRemoved[a].get(0).getNr());
    		activeStates.remove(d1);
    		int d2 = activeStates.indexOf(lineagesRemoved[a].get(1).getNr());
    		activeStates.remove(d2);	
        	for (int j = 0; j < activeStates.size()-1; j++){
				Integer[] migroute = new Integer[]{activeStates.get(j), parent};
				migrationMap.add(migroute);
        	}
        	for (int j = 0; j < activeStates.size()-1; j++){
				Integer[] migroute = new Integer[]{parent, activeStates.get(j)};
				migrationMap.add(migroute);
        	}
        	a++;
    	}
    }

    
    // get the time of the next speciation event
    protected double getNextSpeciationTime(int currTreeInterval){
    	if (currTreeInterval >= intervals.length)
    		return Double.POSITIVE_INFINITY;
    	else
    		return intervals[currTreeInterval];
    }
	
	
	private double[] savedNes;
	private double[] storedNes;	
    
    
    // return the number of internal nodes
    protected int getIntNodes(){
    	return nrSamples-1;
    }
    

    // get the effective population size of a state in the current interval
    public double getPopulationSize(int currentInterval, int state){
    	return NeInput.get().getArrayValue(
    			stateToNodeMap.get(currentInterval-getNumberOfSpecies()).get(state));
    }
    
    //slow to return migration rates
	public double getMigrationRates(int currentInterval, int state1, int state2) {
		for (int i = 0; i < migrationMap.size(); i++){
			if (migrationMap.get(i)[0]==stateToNodeMap.get(currentInterval-getNumberOfSpecies()).get(state1) 
					&& migrationMap.get(i)[1]==stateToNodeMap.get(currentInterval-getNumberOfSpecies()).get(state2))
				return mInput.get().getArrayValue(i);
		}
		return 0.0;
	}
	
	// return all migration rates from a state
	public ArrayList<Double> getAllMigrationRates(int state){
		ArrayList<Double> migRates = new ArrayList<>();

		for (int i = 0; i < migrationMap.size(); i++){
			if (migrationMap.get(i)[0]==state){
				migRates.add(mInput.get().getArrayValue(i));
				migRates.add((double) migrationMap.get(i)[1]);
			}
		}

		return migRates;
	}
	
    protected void coalesce(int currTreeInterval) {
    	
    	// Get the daughter lineages
		List<Node> coalLines = lineagesRemoved[currTreeInterval];
    	if (coalLines.size() > 2) {
			System.err.println("Unsupported coalescent at non-binary node");
			System.exit(0);
		}
    	if (coalLines.size() < 2) {
    		System.err.println();
    		System.err.println("WARNING: Less than two lineages found at coalescent event!");
    		System.err.println();
			System.exit(0);
		}
		
    	//  get the indices of the daughter lineages
    	daughterIndex1 = coalLines.get(0).getNr();
		daughterIndex2 = coalLines.get(1).getNr();	
		
		
		if (daughterIndex1 < 0)
			System.out.println("daughter not found");
		else if(daughterIndex2 < 0)
			System.out.println("daughter not found");
    }
    
    protected int getDaughter1(int currTreeInterval){
    	int index = stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()).indexOf(daughterIndex1);
    	if (index < 0){
    		System.out.println(daughterIndex1);
    		System.out.println(stateToNodeMap);
    		System.out.println(stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()));
    		System.out.println(daughterIndex1 + " " + daughterIndex2);
    		System.out.println(lineagesRemoved[currTreeInterval]);
    	}
    	return index;
    }
    
    protected int getDaughter2(int currTreeInterval){
    	int index = stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()).indexOf(daughterIndex2);

//    	System.out.println(daughterIndex2);
//    	System.out.println(currTreeInterval-getNumberOfSpecies());
//    	System.out.println(stateToNodeMap);
//    	System.out.println("daughter2 " + index);
    	return stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()).indexOf(daughterIndex2);
    }

	public int getCurrentNumberOfStates(int speciesInterval) {	
		if (speciesInterval < (intervals.length+1)/2)
			return (intervals.length+1)/2;
		else
			return intervals.length-speciesInterval;
	}
	
	// gets all the species names from the species Tree. Gets it ones such that the 
	// order is fixed for the MCMC analysis
    private void getSpeciesNames() {
		int currInt = 0; 
		speciesNames = new HashMap<>();
		do{
			List<Node> incomingLines = lineagesAdded[currInt];
			for (Node l:incomingLines){
				speciesNames.put(l.getID(), l.getNr());
			}
			currInt++;
		}while (!isCoalescent[currInt]);
	}

    // get the state of a sample, i.e. the int corresponding to a sampled species
    public Integer getSampleState(String species){
    	for (Node leaf: speciesTreeInput.get().getNodesAsArray()) {
    		if (leaf.getID().equals(species)) return leaf.getNr();
    	}
    	return -1;
    }
    
    public int getNumberOfSpecies(){
    	return speciesTreeInput.get().getLeafNodeCount();	
    }
    
	public double getNodeNe(int nr) {
    	return NeInput.get().getArrayValue(nr);
	}
	
     
	protected int getSpeciesState(int currentInterval, int state){
		return stateToNodeMap.get(currentInterval - speciesTreeInput.get().getLeafNodeCount()).get(state);
	}
	
	
	@Override
	protected void store(){		
		storedStateToNodeMap = new ArrayList<>();
		for (int i= 0; i < stateToNodeMap.size(); i++){
			ArrayList<Integer> add = new ArrayList<>();
			for (int j = 0; j < stateToNodeMap.get(i).size(); j++)
				add.add(stateToNodeMap.get(i).get(j));
			storedStateToNodeMap.add(add);			
		}
		
		storedMigrationMap = new ArrayList<>();
		for (Integer[] m : migrationMap)
			storedMigrationMap.add(m);
		
//		storedNes = Arrays.copyOf(savedNes, savedNes.length);
		super.store();
	}
	
	@Override
	protected void restore(){
		stateToNodeMap = new ArrayList<>();
		for (int i= 0; i < storedStateToNodeMap.size(); i++){
			ArrayList<Integer> add = new ArrayList<>();
			for (int j = 0; j < storedStateToNodeMap.get(i).size(); j++)
				add.add(storedStateToNodeMap.get(i).get(j));
			stateToNodeMap.add(add);			
		}
		
		migrationMap = new ArrayList<>();
		for (Integer[] m : storedMigrationMap)
			migrationMap.add(m);
		
//		savedNes = Arrays.copyOf(storedNes, storedNes.length);
		super.restore();
	}

	@Override
	public double branchLogP(int speciesTreeNodeNumber, Node speciesTreeNode, double ploidy,
			double[] branchCoalescentTimes, int branchLineageCount, int branchEventCount) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void initPopSizes(double initialPopSizes) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean isDirtyBranch(Node speciesTreeNode) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public PopulationModel getBaseModel() {
		// TODO Auto-generated method stub
		return null;
	}
	
    protected boolean intervalIsDirty(int i) {
        if (needsUpdate) {
            calculateIntervals();
        }
        //TODO change to correct
        return true;
//        if (i < intervals.length)
//        	return intervalIsDirty[i];
//        else
//        	return lastIntervalDirty;
        	
    }    
    
    protected void setIntervalIsDirty(int i) {
        if (needsUpdate) {
            calculateIntervals();
        }
        if (i < intervals.length){
        	intervalIsDirty[i] = true;
        }else{
        	lastIntervalDirty = true;
        }
    }    
  
    protected void setIntervalsClean(){
        intervalIsDirty = new boolean[intervals.length];
        lastIntervalDirty = false;
    }


}
