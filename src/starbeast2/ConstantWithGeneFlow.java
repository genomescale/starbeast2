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
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

import java.text.DecimalFormat;
import java.util.*;

/**
 * @author Nicola Felix Mueller
 */

@Description("Species tree that contains the information about things such as the speciation times and migration rates")
@Citation("Nicola F. Müller, Huw A. Ogilvie, Chi Zhang, Alexei J. Drummond and Tanja Stadler(2018)\n  Inference of species histories in the presence of gene flow.\n  bioRxiv doi: 10.1101/348391")
public class ConstantWithGeneFlow extends CalculationNode implements PopulationModel {
    public Input<RealParameter> NeInput = new Input<RealParameter>("Ne","contains the Ne of each branch",Input.Validate.REQUIRED);   
    public Input<RealParameter> NeMeanInput = new Input<RealParameter>("NeMean","contains the Ne of each branch");   
    public Input<RealParameter> mInput  = new Input<>("m","relative migration rates between branches",Input.Validate.REQUIRED);
    public Input<BooleanParameter> indicatorInput  = new Input<>("indicator","indicator if rate is not 0");
    public Input<MigrationModel> migrationModelInput  = new Input<>("migrationModel","input of model of migration",Input.Validate.REQUIRED);
    public Input<Boolean> rateIsForwardInput  = new Input<>("rateIsForward","input of direction of migration",false);
   
    SpeciesTreeInterface speciesTree;
    private MigrationModel migModel;
    
    private boolean needsUpdate = true;
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
	
	private ArrayList<Boolean[]> connected;
	private ArrayList<Boolean[]> storedConnected;
	
	private ArrayList<DoubleArrayForList> migrationRates;
	private ArrayList<IntegerArrayForList> indicators;
	private ArrayList<DoubleArrayForList> storedMigrationRates;
	private ArrayList<IntegerArrayForList> storedIndicators;
	
	

	
    protected double[] intervals;
    protected double[] storedIntervals;
    protected boolean[] isCoalescent;
    protected boolean[] storedIsCoalescent;

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
    
    

    @Override
    public void initAndValidate() {
    	migModel = migrationModelInput.get();
        speciesTree = migModel.speciesTreeInput.get();
        final int speciesNodeCount = speciesTree.getNodeCount();
        leafNodeCount = speciesTree.getLeafNodeCount(); // also the number of "tip" population sizes
        NeInput.get().setDimension(2*leafNodeCount-1);
        mInput.get().setDimension((leafNodeCount-1)*(leafNodeCount-1)*2);
        if (indicatorInput.get()!=null)
        	indicatorInput.get().setDimension((leafNodeCount-1)*(leafNodeCount-1)*2);
        
        needsUpdate = true;
        
    	// Calculate the tree intervals (time between events, which nodes participate at a event etc.)
        calculateIntervals();        
    }
    
    public void checkConnections(){
//    	System.out.println(stateToNodeMap);
//    	System.out.println(speciesTree);
    	connected = new ArrayList<>();
        // check the connectivity for 
        ArrayList<Boolean[][]> isConnectedList = new ArrayList<>();
        
    	ArrayList<Integer> activeStates = new ArrayList<>();
    	
    	int currentInterval = 0;
    	boolean firstCoalescent = true;
    	Boolean[][] isConnected = new Boolean[0][0];
    	for (int i = 0; i < intervals.length; i++){
    		if (isCoalescent[i]){
    			if(firstCoalescent){
    		    	isConnected = new Boolean[activeStates.size()][activeStates.size()];
    		    	for (int a = 0; a < activeStates.size(); a++)
    		    		for (int b = 0; b < activeStates.size(); b++)
    		    			isConnected[a][b] = false;
    		    	getInitialConnectedStates(activeStates, currentInterval, isConnected);
    		    	getConnectedStates(isConnected, activeStates);
    				firstCoalescent = false;
    			}
    			// get the added and removed lineages
    			List<Node> lins = lineagesRemoved[i];
    			int ind1 = activeStates.indexOf(lins.get(0).getNr());
    			int ind2 = activeStates.indexOf(lins.get(1).getNr());
    			// get the correct order of indices
    			int min_ind = Math.min(ind1,  ind2);
    			int max_ind = Math.max(ind1,  ind2);
		    	Boolean[][] new_isConnected = new Boolean[activeStates.size()-1][activeStates.size()-1];
		    	for (int a = 0; a < activeStates.size()-1; a++)
		    		for (int b = 0; b < activeStates.size()-1; b++)
		    			new_isConnected[a][b] = false;

		    	
		    	// add all connections of the preexisting and the daughter species
		    	int a_val,b_val;
		    	for (int a = 0; a < activeStates.size(); a++){
    				if (a<min_ind){
    					a_val = a;
    				}else if(a==min_ind){
    					a_val = activeStates.size()-2;
    				}else if(a>min_ind && a<max_ind){
    					a_val = a-1;
    				}else if(a==max_ind){
    					a_val = activeStates.size()-2;   					
    				}else{
    					a_val = a-2;
    				}
		    		for (int b = 0; b < activeStates.size(); b++){
	    				if (b<min_ind){
	    					b_val = b;
	    				}else if(b==min_ind){
	    					b_val = activeStates.size()-2;
	    				}else if(b>min_ind && b<max_ind){
	    					b_val = b-1;
	    				}else if(b==max_ind){
	    					b_val = activeStates.size()-2;   					
	    				}else{
	    					b_val = b-2;
	    				}
	    				if (isConnected[a][b])
	    					new_isConnected[a_val][b_val] = true;
		    		}
		    	}
		    	// remove daughter lineages
    			activeStates.remove(max_ind);
    			activeStates.remove(min_ind);
		    	
		    	// get the parent lineage
		    	int par_lin = lins.get(1).getParent().getNr();
		    	// find potential connections of the parent lineage
	    		for (int a = 0; a < activeStates.size(); a++){
    				if (a!=par_lin){					
    					for (int c = 0; c < migrationMap.size(); c++){
    						if (migrationMap.get(c)[0]==activeStates.get(a)
    								&& migrationMap.get(c)[1]==par_lin){
    							if(indicatorInput.get().getArrayValue(c) > 0.5)
    								new_isConnected[a][activeStates.size()] = true;
    						}
    					}
    					for (int c = 0; c < migrationMap.size(); c++){
    						if (migrationMap.get(c)[0]==par_lin
    								&& migrationMap.get(c)[1]==activeStates.get(a)){
    							if(indicatorInput.get().getArrayValue(c) > 0.5)
    								new_isConnected[activeStates.size()][a] = true;
    						}
    					}

    				}
	    			
	    		}
				isConnectedList.add(new_isConnected);
				isConnected = new_isConnected;    			    			
    			activeStates.add(par_lin);
    			
//		    	System.out.println(Arrays.toString(isConnected[0]));
//		    	System.out.println(Arrays.toString(isConnected[1]));

    			
		    	getConnectedStates(isConnected, activeStates);
    			
    			currentInterval++;

    			
    		}else{
    			activeStates.add(lineagesAdded[i].get(0).getNr());
    		}
    	}
    	
//    	System.out.println("");
//    	for (int i = 0; i < connected.size(); i++)
//    		System.out.println(Arrays.toString(connected.get(i)));

    	
    }
    
    private void getInitialConnectedStates(ArrayList activeStates, int currentInterval, Boolean[][] isConnected){
    	
		// check migration between any two states
		for (int a = 0; a < activeStates.size(); a++){
			for (int b = 0; b < activeStates.size(); b++){
				if (a!=b){					
					for (int c = 0; c < migrationMap.size(); c++){
						if (migrationMap.get(c)[0]==activeStates.get(a)
								&& migrationMap.get(c)[1]==activeStates.get(b)){
							if(indicatorInput.get().getArrayValue(c) > 0.5){
								isConnected[a][b] = true;
							}
						}
					}
				}
			}
		}
    }
    
    
    private void getConnectedStates(Boolean[][] isConnected, ArrayList<Integer> activeStates){
    	Boolean[] con = new Boolean[speciesTree.getNodeCount()];
    	for (int a = 0; a < con.length; a++)
    		con[a] = false;
    	for (int a = 0; a < isConnected.length; a++){
    		for (int b = 0; b < isConnected.length; b++){
    			if (a!=b){
	    			if (isConnected[a][b]){
	    				con[activeStates.get(a)] = true;
	    				con[activeStates.get(b)] = true;
	    			}
    			}
    		}
    	}
    	connected.add(con);
    	
    }
    
    @SuppressWarnings({ "unchecked", "deprecation" })
    private void calculateIntervals() {
    	Node[] speciesNodesTmp = speciesTree.getNodesAsArray();
    	Node[] speciesNodes = new Node[speciesNodesTmp.length];
    	// make a deep copy of every node
    	for (int i = 0; i < speciesNodesTmp.length; i++)
    		speciesNodes[i] = speciesNodesTmp[i].copy();

    	
    	Arrays.sort(speciesNodes, nhc);
    	
        intervals = new double[speciesNodes.length];
        isCoalescent = new boolean[speciesNodes.length];
        lineagesAdded = new List[speciesNodes.length];
        lineagesRemoved = new List[speciesNodes.length];
        
        ArrayList<Node> initAL = new ArrayList<>(); 
        intervals[0] = speciesNodes[0].getHeight();
        lineagesAdded[0] = new ArrayList<>(); 
		lineagesAdded[0].add(speciesNodes[0].copy());
		isCoalescent[0] = false;
		lineagesRemoved[0] = new ArrayList<>(); 
		lineagesRemoved[0].add(null);
		lineagesRemoved[0].add(null);    		  			
		for (int i = 1; i < speciesNodes.length; i++){
    		intervals[i] = speciesNodes[i].getHeight()-speciesNodes[i-1].getHeight();
            lineagesAdded[i] = new ArrayList<>(); 
    		lineagesAdded[i].add(speciesNodes[i].copy());
    		if (!speciesNodes[i].isLeaf()){
    			isCoalescent[i] = true;
    			lineagesRemoved[i] = new ArrayList<>(); 
	    		lineagesRemoved[i].add(speciesNodes[i].getLeft().copy());
	    		lineagesRemoved[i].add(speciesNodes[i].getRight().copy()); 
	    		lineagesRemoved[i].get(0).setParent(lineagesAdded[i].get(0));
	    		lineagesRemoved[i].get(1).setParent(lineagesAdded[i].get(0));
    		}else{
    			isCoalescent[i] = false;
    			lineagesRemoved[i] = new ArrayList<>(); 
	    		lineagesRemoved[i].add(null);
	    		lineagesRemoved[i].add(null);    		  			
    		}
    	}   
    	
    	stateToNodeMap();
    	preComputeMigrationRates();
        if (indicatorInput.get()!=null)
        	checkConnections();
        
        needsUpdate = false;
        

    }
    
	// builds the map from state to node number
    private void stateToNodeMap(){
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
    				System.out.println(lins.get(1).getParent().getNr() + " " + lineagesAdded[i].get(0).getNr());
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
    	migrationRates = new ArrayList<>();
    	
    	int migrationRatePoint = 0;
    	
    	// get the migration map for each interval
    	activeStates = new ArrayList<>();
    	int a = 0;
    	while (!isCoalescent[a]){
			activeStates.add(lineagesAdded[a].get(0).getNr());
    		a++;
    	}    	
    	// sorting ensures the correct order of migration rate elements
    	Collections.sort(activeStates);
    	
    	
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
    
    private void preComputeMigrationRates(){
    	indicators = new ArrayList<>();
    	migrationRates = new ArrayList<>();
    	// get the correct mapping of migration rates
    	//TODO speed up that thing
    	for (int i = 0; i < stateToNodeMap.size(); i++){
    		int[][] mRateMapping = new int[stateToNodeMap.get(i).size()][stateToNodeMap.get(i).size()];
    		for (int a = 0; a < mRateMapping.length; a++){
    			for (int b = 0; b < mRateMapping.length; b++){
    				for (int k = 0; k< migrationMap.size(); k++){
    					if (migrationMap.get(k)[0]==stateToNodeMap.get(i).get(a) 
    							&& migrationMap.get(k)[1]==stateToNodeMap.get(i).get(b)){
    						mRateMapping[a][b] = k;
    					}
    				}
    			}
    		}
    		
    		// build the indicators list
    		if (indicatorInput.get()!=null){
	        	ArrayList<Integer[]> indicatorList = new ArrayList<>();
	        	for (int a = 0; a < mRateMapping.length; a++){
	    			for (int b = 0; b < mRateMapping.length; b++){
	    				if (a!=b){
		    				if(indicatorInput.get().getArrayValue(mRateMapping[a][b]) > 0.5){
		    					Integer[] add={a,b};
		    					indicatorList.add(add);
		    				}
	    				}
	    			}
	    		}
	        	IntegerArrayForList indicators_array = new IntegerArrayForList(indicatorList.size());
	    		for (int a = 0; a < indicatorList.size(); a++){
	    			indicators_array.setValue(a, indicatorList.get(a)[0], indicatorList.get(a)[1]);
	    		}  
	    		indicators.add(indicators_array);
	    		
	    		// also add migration rates
	    		DoubleArrayForList migration_array = new DoubleArrayForList(mRateMapping.length);
	    		for (int a = 0; a < indicatorList.size(); a++){
	    			int x_c = indicators_array.getValue(a, 0);
	    			int y_c = indicators_array.getValue(a, 1);
	    			migration_array.setValue(x_c, y_c, migModel.getMigration(stateToNodeMap.get(i).get(x_c) , stateToNodeMap.get(i).get(y_c))
	    					* mInput.get().getArrayValue(mRateMapping[x_c][y_c]));
	    		}
	    		migrationRates.add(migration_array);
	    				
    		}else{
	    		DoubleArrayForList migration_array = new DoubleArrayForList(mRateMapping.length);
	        	for (int a = 0; a < mRateMapping.length; a++){
	    			for (int b = 0; b < mRateMapping.length; b++){
	    				if (a!=b){
	    					migration_array.setValue(a, b, migModel.getMigration(stateToNodeMap.get(i).get(a) , stateToNodeMap.get(i).get(b))
	    							* mInput.get().getArrayValue(mRateMapping[a][b]));
	    				}
	    			}
	    		}    
	    		migrationRates.add(migration_array);
    		}
    	}
    	
    }
    
    // get the time of the next speciation event
    protected double getNextSpeciationTime(int currTreeInterval){
    	if (needsUpdate)
    		calculateIntervals();
    	
    	if (currTreeInterval >= intervals.length)
    		return Double.POSITIVE_INFINITY;
    	else
    		return intervals[currTreeInterval];
    }

   
    // return the number of internal nodes
    protected int getIntNodes(){
    	if (needsUpdate)
    		calculateIntervals();

    	return nrSamples-1;
    }
    

    // get the effective population size of a state in the current interval
    public double getPopulationSize(int currentInterval, int state){
    	if (needsUpdate)
    		calculateIntervals();
    	
    	if (NeMeanInput.get()!=null)
    		return NeMeanInput.get().getValue()*NeInput.get().getArrayValue(
        			stateToNodeMap.get(currentInterval-getNumberOfSpecies()).get(state));
    	else
    		return NeInput.get().getArrayValue(
    			stateToNodeMap.get(currentInterval-getNumberOfSpecies()).get(state));
    }
    
    public boolean getIsConnected(int currentInterval, int state){
    	if (needsUpdate)
    		calculateIntervals();

        if (indicatorInput.get()!=null)      
        	return connected.get(currentInterval-getNumberOfSpecies())[
    			stateToNodeMap.get(currentInterval-getNumberOfSpecies()).get(state)];
        else
        	return true;
    }
    
    public int getIsNodeConnected(int nodeNr){
    	if (needsUpdate)
    		calculateIntervals();

    	for (int i =0; i < connected.size();i++)
    		if (connected.get(i)[nodeNr])
    			return 1;
    			
    	return 0;
    }
    
    
//    migrationRates
//    
//    private void computeAllMigrationRates(){
//    	
//    }
//    
//    private double[][] getAllMigrationRates(){
//    	
//    }
    
    // faster way to return migration rates
    
    public double[][] getMigrationRates(int currentInterval){
    	if (needsUpdate)
    		calculateIntervals();
    	
    	
    	return migrationRates.get(currentInterval-getNumberOfSpecies()).getArray();
    }
    
    public int[][] getIndicatorsRates(int currentInterval){
    	if (needsUpdate)
    		calculateIntervals();

    	return indicators.get(currentInterval-getNumberOfSpecies()).getArray();
    }

    
    //slow to return migration rates
	public double getMigrationRates(int currentInterval, int state1, int state2) {
    	if (needsUpdate)
    		calculateIntervals();

		int interval = currentInterval-getNumberOfSpecies();
//		double Nesink = NeInput.get().getArrayValue(
//    			stateToNodeMap.get(currentInterval-getNumberOfSpecies()).get(state1));
//		double Nesource = NeInput.get().getArrayValue(
//    			stateToNodeMap.get(currentInterval-getNumberOfSpecies()).get(state2));
		
		double migration = migModel.getMigration(stateToNodeMap.get(interval).get(state1) , stateToNodeMap.get(interval).get(state2));
//		migration *= Nesource;
//		migration /= Nesink;
		for (int i = 0; i < migrationMap.size(); i++){
			if (migrationMap.get(i)[0]==stateToNodeMap.get(interval).get(state1) 
					&& migrationMap.get(i)[1]==stateToNodeMap.get(interval).get(state2)){
					if (indicatorInput.get()!=null){
						if(indicatorInput.get().getArrayValue(i) > 0.5)
							return migration*mInput.get().getArrayValue(i);
						else
							return 0.0;
					}else{
						return migration*mInput.get().getArrayValue(i);						
					}
			}
		}
		return 0.0;
	}
	
    //get migration rates between nodes
	
	public double getMigrationRates(int node1, int node2) {	
    	if (needsUpdate)
    		calculateIntervals();

		double migration = migModel.getMigration(node1 , node2);
		for (int i = 0; i < migrationMap.size(); i++){
			if (migrationMap.get(i)[0]==node1 
					&& migrationMap.get(i)[1]==node2){
				if (!rateIsForwardInput.get()){
					double NeRatio = NeInput.get().getArrayValue(node2)/NeInput.get().getArrayValue(node1);
					return migration*mInput.get().getArrayValue(i) * NeRatio;
				}else{
					return migration*mInput.get().getArrayValue(i);
				}
			}
		}
		return 0.0;
	}
	
	
	public ArrayList<ArrayList<Integer>> getStateToNodeMap(){
    	if (needsUpdate)
    		calculateIntervals();

		ArrayList<ArrayList<Integer>> returnList = new ArrayList<>();
		for (int i= 0; i < stateToNodeMap.size(); i++){
			ArrayList<Integer> add = new ArrayList<>();
			for (int j = 0; j < stateToNodeMap.get(i).size(); j++)
				add.add(stateToNodeMap.get(i).get(j));
			returnList.add(add);		
		}		
		return returnList;
	}

	
	// return all migration rates from a node
	public String getAllMigrationRates(int nodeNr){
    	if (needsUpdate)
    		calculateIntervals();

		boolean isEmpty = true;
		String migRates = new String();	
		String migTo = new String();
		migRates ="abcdef";
		migTo ="abcdef";
		ArrayList<Integer> visited = new ArrayList<>();
		for (int i = 0; i < (stateToNodeMap.size()-1); i++){
			for (int j = 0; j < stateToNodeMap.get(i).size(); j++){
				if (stateToNodeMap.get(i).get(j) == nodeNr){
					for (int k = 0; k < stateToNodeMap.get(i).size(); k++){
						if (k!=j && visited.indexOf(stateToNodeMap.get(i).get(k))==-1){
							visited.add(stateToNodeMap.get(i).get(k));
							double migration = migModel.getMigration(stateToNodeMap.get(i).get(j) , stateToNodeMap.get(i).get(k));
							for (int l = 0; l < migrationMap.size(); l++){
								if (migrationMap.get(l)[0]==stateToNodeMap.get(i).get(j) 
										&& migrationMap.get(l)[1]==stateToNodeMap.get(i).get(k)){
										if (indicatorInput.get()!=null){
											if (indicatorInput.get().getArrayValue(l) > 0.5){
												migration*=mInput.get().getArrayValue(l);											
											}else{
												migration*=0.0;											
											}											
										}else{
											migration*=mInput.get().getArrayValue(l);
										}
								}
							}
							migRates = migRates + "," + migration;
							migTo = migTo + "," + stateToNodeMap.get(i).get(k);
							isEmpty = false;
						}
					}
				}
			}
		}		
		migRates = migRates.replace("abcdef,", ",rates={");
		migTo = migTo.replace("abcdef,", ",to={");
		
		migRates = migRates + "}";
		migTo = migTo + "}";
		String returnString = new String();
		returnString =  migTo + migRates;

		if (isEmpty)
			return "";
		else
			return returnString;
	}
	
    
	public String getAllMigrationRatesLong(int nodeNr){
    	if (needsUpdate)
    		calculateIntervals();

		boolean isEmpty = true;
		String migRates = new String();	
		String migInd = new String();
		migRates ="abcdef";
		migInd ="abcdef";
		ArrayList<Integer> visited = new ArrayList<>();
		for (int i = 0; i < (stateToNodeMap.size()-1); i++){
			for (int j = 0; j < stateToNodeMap.get(i).size(); j++){
				if (stateToNodeMap.get(i).get(j) == nodeNr){
					for (int k = 0; k < stateToNodeMap.get(i).size(); k++){
						if (k!=j && visited.indexOf(stateToNodeMap.get(i).get(k))==-1){
							visited.add(stateToNodeMap.get(i).get(k));
							double migration = migModel.getMigration(stateToNodeMap.get(i).get(j) , stateToNodeMap.get(i).get(k));
							for (int l = 0; l < migrationMap.size(); l++){
								if (migrationMap.get(l)[0]==stateToNodeMap.get(i).get(j) 
										&& migrationMap.get(l)[1]==stateToNodeMap.get(i).get(k)){
										if (indicatorInput.get()!=null){
											if (indicatorInput.get().getArrayValue(l) > 0.5){
												migration=mInput.get().getArrayValue(l)*migModel.getEM();											
											}else{
												migration*=0.0;											
											}											
										}else{
											migration=mInput.get().getArrayValue(l)*migModel.getEM();
										}
								}
							}
							// get the node numbers of all children
							ArrayList<Integer> node_nr = new ArrayList<>();
							for (Node n : speciesTree.getNode(stateToNodeMap.get(i).get(k)).getAllLeafNodes())
								node_nr.add(n.getNr());
							// sort the node numbers
							Collections.sort(node_nr);	
							
							migRates = migRates + ",r";
							migInd = migInd + ",i";
							
							for (Integer leaf_nr : node_nr){
								migRates = migRates + "." + (leaf_nr+1);
								migInd = migInd + "." + (leaf_nr+1);
							}
							
							if (node_nr.size()==0){
								migRates = migRates + "." + (stateToNodeMap.get(i).get(k)+1);
								migInd = migInd + "." + (stateToNodeMap.get(i).get(k)+1);
							}
								
							migRates = migRates + "=" + migration;
							if (migration==0)
								migInd = migInd + "=0";
							else
								migInd = migInd + "=1";
							
							isEmpty = false;
						}
					}
				}
			}
		}		
		migRates = migRates.replace("abcdef,", ",");
		migInd = migInd.replace("abcdef,", ",");
		
		String returnString = new String();
		returnString =  migInd + "" + migRates;

		if (isEmpty)
			return "";
		else
			return returnString;
	}
	

	
	protected void coalesce(int currTreeInterval) {
    	if (needsUpdate)
    		calculateIntervals();
    	
    	// Get the daughter lineages
		List<Node> coalLines = lineagesRemoved[currTreeInterval];
    	if (coalLines.size() > 2) {
			System.err.println("Unsupported coalescent at non-binary node");
			System.exit(0);
		}
    	if (coalLines.size() < 2) {
    		System.out.println(coalLines);
    		System.out.println(lineagesAdded[currTreeInterval]);
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
    	if (needsUpdate)
    		calculateIntervals();

    	int index = stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()).indexOf(daughterIndex1);
    	if (index < 0){
    		System.out.println(daughterIndex1);
    		System.out.println(stateToNodeMap);
    		System.out.println(stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()));
    		System.out.println(daughterIndex1 + " " + daughterIndex2);
    		System.out.println(lineagesRemoved[currTreeInterval]);
    		System.out.println("daughter lineage 1 not found");
    		for (int i = 0; i < lineagesAdded.length; i++)
    			System.out.println(lineagesAdded[i]);
    		System.out.println(lineagesAdded.length);
    		System.out.println(speciesTree);
    		calculateIntervals();
    		System.out.println("");
    		for (int i = 0; i < lineagesAdded.length; i++)
    			System.out.println(lineagesAdded[i]);
    		System.out.println(lineagesAdded.length);
    		System.out.println(lineagesAdded.length);
    		System.exit(0);
    	}
    	return index;
    }
    
    protected int getDaughter2(int currTreeInterval){
    	if (needsUpdate)
    		calculateIntervals();

    	int index = stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()).indexOf(daughterIndex2);
    	if (index < 0){
    		System.out.println(daughterIndex1);
    		System.out.println(stateToNodeMap);
    		System.out.println(stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()));
    		System.out.println(daughterIndex1 + " " + daughterIndex2);
    		System.out.println(lineagesRemoved[currTreeInterval]);
    		System.out.println("daughter lineage 2 not found");
    	}

    	return stateToNodeMap.get(currTreeInterval-getNumberOfSpecies()).indexOf(daughterIndex2);
    }

	public int getCurrentNumberOfStates(int speciesInterval) {	
    	if (needsUpdate)
    		calculateIntervals();

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
    	for (Node leaf: speciesTree.getNodesAsArray()) {
    		if (leaf.getID().equals(species)) return leaf.getNr();
    	}
    	return -1;
    }
    
    public int getNumberOfSpecies(){
    	return speciesTree.getLeafNodeCount();	
    }
    
	public double getNodeNe(int nr) {
		if (NeMeanInput.get()!=null)
			return NeMeanInput.get().getValue()*NeInput.get().getArrayValue(nr);
		else
			return NeInput.get().getArrayValue(nr);
	}
	
     
	protected int getSpeciesState(int currentInterval, int state){
    	if (needsUpdate)
    		calculateIntervals();

		return stateToNodeMap.get(currentInterval - speciesTree.getLeafNodeCount()).get(state);
	}
	
	
	@SuppressWarnings("unchecked")
	@Override
	protected void store(){	
//		System.out.println("store..");

	    storedIntervals = new double[intervals.length];
		System.arraycopy(intervals,0,storedIntervals,0,intervals.length);
	    		
		
		
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
		
		// store the lineages added and removed
		storedLineagesAdded = new List[lineagesAdded.length];
		storedLineagesRemoved = new List[lineagesRemoved.length];
		
//		System.out.println("");
//		for (int i = 0; i <lineagesAdded.length; i++)
//			System.out.println(lineagesAdded[i] + " " + storedLineagesAdded[i]);

		
		for (int i = 0; i < lineagesAdded.length; i++){
			ArrayList<Node> add = new ArrayList<>();
			for (Node n : lineagesAdded[i]){
				add.add(n.copy());
			}
			storedLineagesAdded[i] = new ArrayList<>(add);
		}
			
		for (int i = 0; i < lineagesRemoved.length; i++){
			ArrayList<Node> remove = new ArrayList<Node>();
			if (lineagesRemoved[i] != null){
				for (Node n : lineagesRemoved[i]){
					if (n!=null){
						remove.add(n.copy());
					}else{
						remove.add(null);
					}						
				}
				storedLineagesRemoved[i] = new ArrayList<>(remove);
			}else{
				storedLineagesRemoved[i] = null;
			}

		}	
		
		storedIsCoalescent = new boolean[isCoalescent.length];
		System.arraycopy(isCoalescent,0,storedIsCoalescent,0,isCoalescent.length);
		
		storedMigrationRates = new ArrayList<>();
		for (DoubleArrayForList dafl : migrationRates){
			DoubleArrayForList add  = new DoubleArrayForList(dafl.getArray());
			storedMigrationRates.add(add);
		}
		
		storedIndicators = new ArrayList<>();
		for (IntegerArrayForList iafl : indicators){
			IntegerArrayForList add  = new IntegerArrayForList(iafl.getArray());
			storedIndicators.add(add);
		}
//		System.exit(0);

		storedConnected = new ArrayList<>();
		for (Boolean[] bool : connected){
			Boolean[] con = new Boolean[bool.length];
			for (int i = 0; i < con.length; i++)
				con[i] = bool[i];
			storedConnected.add(con);
		}
		
//		System.out.println("");
//		for (int i = 0; i <lineagesAdded.length; i++)
//			System.out.println(lineagesAdded[i]);
//		System.out.println("");
//		for (int i = 0; i <lineagesAdded.length; i++)
//			System.out.println(storedLineagesAdded[i]);

		
//		storedNes = Arrays.copyOf(savedNes, savedNes.length);
		super.store();
	}
	
	@SuppressWarnings("unchecked")
	@Override
	protected void restore(){
//		System.out.println("restore..");
		intervals = new double[storedIntervals.length];
		System.arraycopy(storedIntervals,0,intervals,0,storedIntervals.length);

		
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
		
		// restore the lineages added and removed
		lineagesAdded = new List[storedLineagesAdded.length];
		lineagesRemoved = new List[storedLineagesRemoved.length];
		for (int i = 0; i < storedLineagesAdded.length; i++){
			ArrayList<Node> add = new ArrayList<Node>();
			for (Node n : storedLineagesAdded[i])
				add.add(n.copy());
			lineagesAdded[i] = new ArrayList<>(add);			
		}
		
		for (int i = 0; i < storedLineagesRemoved.length; i++){
			ArrayList<Node> remove = new ArrayList<Node>();
			if (storedLineagesRemoved[i]!=null){
				for (Node n : storedLineagesRemoved[i])
					if (n!=null){
						remove.add(n.copy());
					}else{
						remove.add(null);
					}

				lineagesRemoved[i] = new ArrayList<>(remove);
			}else{
				lineagesRemoved[i] = null;
			}

		}				

		isCoalescent = new boolean[storedIsCoalescent.length];
		System.arraycopy(storedIsCoalescent,0,isCoalescent,0,storedIsCoalescent.length);

		
		migrationRates = new ArrayList<>();
		for (DoubleArrayForList dafl : storedMigrationRates){
			DoubleArrayForList add  = new DoubleArrayForList(dafl.getArray());
			migrationRates.add(add);
		}
		
		indicators = new ArrayList<>();
		for (IntegerArrayForList iafl : storedIndicators){
			IntegerArrayForList add  = new IntegerArrayForList(iafl.getArray());
			indicators.add(add);
		}

		connected = new ArrayList<>();
		for (Boolean[] bool : storedConnected){
			Boolean[] con = new Boolean[bool.length];
			for (int i = 0; i < con.length; i++)
				con[i] = bool[i];
			connected.add(con);
		}
		
//		System.out.println("");
//		for (int i = 0; i < lineagesAdded.length; i++)
//			System.out.println(lineagesAdded[i]);


		
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

//	
//    protected boolean intervalIsDirty(int i) {
//        if (needsUpdate) {
//            calculateIntervals();
//        }
//        //TODO change to correct
//        return true;
////        if (i < intervals.length)
////        	return intervalIsDirty[i];
////        else
////        	return lastIntervalDirty;
//        	
//    }    
//    
//    protected void setIntervalIsDirty(int i) {
//        if (needsUpdate) {
//            calculateIntervals();
//        }
//        if (i < intervals.length){
//        	intervalIsDirty[i] = true;
//        }else{
//        	lastIntervalDirty = true;
//        }
//    }    
//  
//    protected void setIntervalsClean(){
//        intervalIsDirty = new boolean[intervals.length];
//        lastIntervalDirty = false;
//    }

    // helper class to have primitiv classes
    private class DoubleArrayForList{
    	double[][] array;
    	
    	DoubleArrayForList(int states){
    		this.array = new double[states][states];
    	}
    	
    	DoubleArrayForList(double[][] array){
    		this.array = new double[array.length][array.length];
    		for (int i = 0; i < this.array.length; i++)
    			System.arraycopy(array[i], 0, this.array[i], 0, array[i].length);
    	}
    	
    	void setValue(int i, int j, double val){
    		this.array[i][j] = val;
    	}
    	
    	double[][] getArray(){
    		return array;
    	}
    }
    
    // helper class to have primitiv classes
    private class IntegerArrayForList{
    	int[][] array;
    	
    	IntegerArrayForList(int states){
    		this.array = new int[states][2];
    	}
    	
    	IntegerArrayForList(int[][] array){
    		this.array = new int[array.length][2];
    		for (int i = 0; i < this.array.length; i++)
    			System.arraycopy(array[i], 0, this.array[i], 0, array[i].length);
//    		for (int i = 0; i < this.array.length; i++){
//    			System.out.println(Arrays.toString(this.array[i]));
//    			System.out.println(Arrays.toString(array[i]));
//    		}
    	}

    	
    	void setValue(int i, int val0, int val1){
    		this.array[i][0] = val0;
    		this.array[i][1] = val1;
    	}
    	
    	int getValue(int i, int j){
    		return this.array[i][j];
    	}
    	
    	int[][] getArray(){
    		return array;
    	}
    }


}
