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

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.*;


/**
 * @author Nicola Felix Mueller
 */
// Gene tree type trait nested in species tree as taxon superset
@Description("Calculate the probability of a tree under the structured coalescent assuming lineage independence with constant rates" +
		" as described in Mueller et. al.,2016. The Input rates are backwards in time migration rates" +
		" and pairwise coalescent rates translating to 1/Ne for m=1 in the Wright Fisher model")
public class GeneTreeWithMigration extends Distribution {
    public Input<Tree> treeInput = new Input<>("tree", "The gene tree.", Validate.REQUIRED);
    public Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for this gene, typically a whole number or half (default is 2).", 2.0);
    public Input<ConstantWithGeneFlow> popModelInput = new Input<>("populationModel", "Population model used to infer the multispecies coalescent probability for this gene");

//    public Input<TraitSet> typeTraitInput = new Input<>(
//            "species", "Type trait set.", Input.Validate.REQUIRED);     
//            		
//    public Input<SpeciesTree> speciesTreeInput = new Input<>(
//    		"speciesTree",
//    		"species tree input", Input.Validate.REQUIRED);    
//        
//    public Input<Double> inheritanceInput = new Input<>(
//    		"inheritance",
//    		"inheritance scaler", 4.0);
    
   
	public int samples;
	public int nrSamples;
	
	
	private double[] coalescentRates;
	private double[][] migrationRates;
	private double[] linProbs;
	private int[] multiplicator;
	private int[] mostLikelyState;
	
    public int states;    
    private boolean recalculateLogP;    
    private double maxTolerance;            
    private int nr_lineages;  
    
    // store the linProbs, multiplicators and logP's at coalescent points in jagged arrays from last time
    private double[][] coalLinProbs;
    private int[][] coalMultiplicator;
    private double[] coalLogP;
    private int[] coalSpeciesInterval;
    private ArrayList<ArrayList<Integer>> coalActiveLineages;
    private ArrayList<ArrayList<Integer>> coalSampleState;
    
    // deep store the things above for MCMC
    private double[][] storeLinProbs;
    private int[][] storeMultiplicator;
    private double[] storeLogP;
    private int[] storeSpeciesInterval;
    private ArrayList<ArrayList<Integer>> storeActiveLineages;
    private ArrayList<ArrayList<Integer>> storeSampleState;
    private boolean storedFirst;
    
        
    // Set up for lineage state probabilities
    private ArrayList<Integer> activeLineages;
    private ArrayList<Integer> sampleState;
    
    // node count
    private int nodeCount;
    
    // check if this is the first calculation
    private boolean first = true;
    
    @Override
    public boolean requiresRecalculation(){
    	return true;   		
    }
    
    @Override
    public void initAndValidate(){    
    	popModelInput.get().calculateIntervals();
    	// Calculate the tree intervals (time between events, which nodes participate at a event etc.)
    	calculateIntervals(); 
                
    	// initialize the maximum relative tolerance 
    	maxTolerance = 1e-3;    	
    	    	
//    	// initialize storing arrays and ArrayLists
//    	coalLinProbs = new double[nodeCount+1][];
//    	coalMultiplicator = new int[nodeCount+1][];
//    	coalLogP = new double[nodeCount+1];
//    	coalSpeciesInterval = new int[nodeCount+1];
//    	coalActiveLineages = new ArrayList<>();
//    	coalSampleState = new ArrayList<>();
//    	
//    	
//    	ArrayList<Integer> emptyList = new ArrayList<>();
//    	for (int i = 0; i <= nodeCount; i++){
//    		coalActiveLineages.add(emptyList);
//    		coalSampleState.add(emptyList);
//    	}
    	
        calculateLogP();
    }
    
    public double calculateLogP() {  
//    	System.out.println("new tree");
    	popModelInput.get().calculateIntervals();
    	recalculateLogP = false;
		calculateIntervals();
    	
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 
        sampleState = new ArrayList<Integer>(); 
        
        // set back logP
        logP = 0;
        
        // initialize the counter for the current interval for the gene and the species trees
        int geneInterval = 0;
        int speciesInterval = popModelInput.get().getNumberOfSpecies();
        
        // initialize the number of state as the number of leaves in the species tree
        states = speciesInterval;   
        // initialize the number of leafs for the gene tree
        nr_lineages = treeInput.get().getLeafNodeCount() + 1;        
   	
    	double nextGeneTime = intervals[geneInterval];
    	double nextSpeciesTime = popModelInput.get().getNextSpeciationTime(speciesInterval);  
//    	if (first) {	        
	        // initialize the array that saves the most likely state of each node
	        mostLikelyState = new int[2*nr_lineages+1];	        
	        
	        // initialize the coalescent and migration rates
	        updateRatesList(speciesInterval);
	        
	        int linsAdded = 0;
	        while (!isCoalescent[geneInterval]){
	        	addLineages(geneInterval, linsAdded);
	        	geneInterval++;
	        	linsAdded++;
	        }
	        // initialize the lineage state probs array
	        initializeP();
	        
			// store the node
			storeNode(geneInterval - nodeCount - 1, linProbs,
						multiplicator, logP, speciesInterval,
						activeLineages, sampleState);
			nextGeneTime = intervals[geneInterval];

//    	} else {	
//    		// Check if the first interval is already dirty
//			if (popModelInput.get().intervalIsDirty(speciesInterval) || intervalIsDirty(geneInterval)){ 
//				restoreNode(0);				
//        		updateRatesList(speciesInterval);
//			}else{   	
//				double lastNextSpeciesTime = nextSpeciesTime;
//				int lastSpeciesInterval = speciesInterval;
//				
//	    		while(geneInterval < Integer.MAX_VALUE){
//	    			// proceed to the next interval
//		    		if (nextSpeciesTime < nextGeneTime){  
//		    			// check if next interval is dirty, if yes stay in that gene interval
//		    			if (popModelInput.get().intervalIsDirty(speciesInterval+1)){
//		    				// use the next species time from the point of the last gene coalescent event
//		    				nextSpeciesTime = lastNextSpeciesTime;
//			    			nextGeneTime = intervals[geneInterval];
//			    			speciesInterval = lastSpeciesInterval;
//		    				// update the number of states
//		    				states = 2*popModelInput.get().getNumberOfSpecies() - speciesInterval;
//
//		    				// restore from last (gene) coalescent event
//		    				restoreNode(geneInterval-treeInput.get().getLeafNodeCount()-1);
//		    				updateRatesList(speciesInterval);
//		    				break;
//		    			}
//		    			states--;
//		        		speciesInterval++; 
//		        		nextGeneTime -= nextSpeciesTime;
//		        		nextSpeciesTime = popModelInput.get().getNextSpeciationTime(speciesInterval);
//		    		}else{
//		        		geneInterval++;
//		        		// Check if the last interval was reached
//		        		if (geneInterval == intervals.length){
//		        			logP = coalLogP[coalLogP.length-1];
//		        			return logP;
//		        		}			 
//		        		
//		        		
//		        		nr_lineages--;
//		        		nextSpeciesTime -= nextGeneTime;
//		        		lastNextSpeciesTime = nextSpeciesTime;
//		        		lastSpeciesInterval = speciesInterval;
//		        		
//		    			nextGeneTime = intervals[geneInterval];
//		    			// check if this gene interval is dirty
//		    			if (intervalIsDirty(geneInterval)){
//		    				// if the interval is dirty, restore from last coalescent event
//		    				speciesInterval = restoreNode(geneInterval-treeInput.get().getLeafNodeCount()-1);
//		    				updateRatesList(speciesInterval);
//		    				break;
//		    			}
//		    		}
//	    		}
//			}
//		}	   	
//			if(geneInterval!=26){
//				System.exit(0);
//			}
				
//		System.out.println(logP);
//		System.out.println(treeInput.get().getID() + " " + geneInterval + " " + intervals[geneInterval] + " next speciation " + nextSpeciesTime + " species interval " + speciesInterval);
        // integrate until there are no more intervals
        do {			
        	// Length of the current interval
        	final double duration = Math.min(nextGeneTime, nextSpeciesTime);
        	
        	// if the current interval has a length greater than 0, integrate
        	if (duration > 0) {        		
                if(recalculateLogP){
                	System.out.println(Arrays.toString(linProbs));
    				System.err.println("ode calculation stuck, reducing tolerance");
    				System.err.println("new tolerance: " + maxTolerance);
    				maxTolerance *=0.9;
                	return calculateLogP();
                }              
                if (states>1){	        		
		        	double[] linProbs_tmp = new double[linProbs.length+1]; 
		        	double[] linProbs_tmpdt = new double[linProbs.length+1]; 
		        	double[] linProbs_tmpddt = new double[linProbs.length+1]; 
		        	double[] linProbs_tmpdddt = new double[linProbs.length+1]; 
	        		
	                Euler2ndOrder euler = new Euler2ndOrder(multiplicator, migrationRates, coalescentRates, multiplicator.length , states, 0.001, 0.2);
	                
	                
		        	for (int i = 0; i < linProbs.length; i++) linProbs_tmp[i] = linProbs[i];		        	
		        	
		        	linProbs[linProbs.length-1] = 0;
		        	euler.calculateValues(duration, linProbs_tmp, linProbs_tmpdt, linProbs_tmpddt, linProbs_tmpdddt);		        	
	        		
		            for (int i = 0; i < linProbs.length; i++) linProbs[i] = linProbs_tmp[i]; 
		            logP += linProbs_tmp[linProbs_tmp.length-1];
		            
                }else{
                	// use analytical solution instead
                	double nonCoalProb = Math.exp(-duration * coalescentRates[0] * activeLineages.size()*(activeLineages.size()-1));
	                for (int i = 0; i<linProbs.length; i++)
	                	linProbs[i] *= nonCoalProb;  
                }
                
        	}
        	
        	if (nextSpeciesTime < nextGeneTime){
                // update the number of states
        		states--;
                // update the coalescent and migration rates
        		updateRatesList(speciesInterval+1);
        		// combine the state probs of the newly combined states
        		combineStates(speciesInterval, logP, geneInterval);        		
        		speciesInterval++;  
        		nextGeneTime -= nextSpeciesTime;
        		nextSpeciesTime = popModelInput.get().getNextSpeciationTime(speciesInterval);
        	}else if (nextSpeciesTime >= nextGeneTime){
                logP += coalesce(geneInterval, speciesInterval); 
                
                if (logP == Double.NEGATIVE_INFINITY){
                	return logP;                
                }
        		geneInterval++;
        		nr_lineages--;
        		nextSpeciesTime -= nextGeneTime;
        		if (geneInterval<intervals.length)
        			nextGeneTime = intervals[geneInterval];
        		else
        			break;          		
        	}else{
        		System.err.println("speciation and coalescence coincide, return negative infinity");
        	}
        }while(activeLineages.size()>1);
        first = false;  
        return logP;  	
    }    
   
    // works: adds all sampled genes to the active lineages and saves their initial state
    private void addLineages(int currTreeInterval, int alreadyAdded) {
		List<Node> incomingLines = lineagesAdded[currTreeInterval];
		if (incomingLines.size()>1){
			System.err.println("too many lineages in sampling interval");
		}
		final int geneLeafNr = incomingLines.get(0).getNr();
		activeLineages.add(incomingLines.get(0).getNr());
		final String geneLeafName = incomingLines.get(0).getID();
		//String sampledSpecies = typeTraitInput.get().getStringValue(geneLeafName);
		final Map<String, Integer> tipNumberMap = popModelInput.get().speciesTreeInput.get().getTipNumberMap();
		final int sampledSpeciesNr = tipNumberMap.get(geneLeafName);
		final String sampledSpeciesName = popModelInput.get().speciesTreeInput.get().getNode(sampledSpeciesNr).getID();
		sampleState.add(popModelInput.get().getSampleState(sampledSpeciesName));
		mostLikelyState[incomingLines.get(0).getNr()] = popModelInput.get().getSampleState(sampledSpeciesName);
    }
 
    private void initializeP(){
    	final int nrSpecies = popModelInput.get().speciesTreeInput.get().getLeafNodeCount();
    	linProbs = new double[nrSpecies*nrSpecies];
    	multiplicator = new int[nrSpecies];
    	for (int i = 0; i < nrSpecies; i++){
    		linProbs[nrSpecies*i+i]=1; 
    	}
    	for (int i = 0; i < sampleState.size(); i++){
			multiplicator[sampleState.get(i)] += 1;
    	}    	
    }
    
    // works: normalizes the lineage state probabilities 
    private double normalizeLineages(){
    	double interval = 0.0;
    	for (int i = 0; i < linProbs.length/states; i++){
    		double lineProbs = 0.0;
    		for (int j = 0; j < states; j++)
    			if (linProbs[i*states+j]>=0.0){
    				lineProbs += linProbs[i*states+j];
    			}else{
    				// try recalculation after lowering the tolerance
    				recalculateLogP = true;
    				return Math.log(1.0);
    			}
    		for (int j = 0; j < states; j++){
    			linProbs[i*states+j] = linProbs[i*states+j]/lineProbs;
    		}    		
    		interval +=lineProbs;
    	}	
    	
		if (nr_lineages>0){
			// return mean P_t(T)
			return Math.log(interval/(multiplicator.length));
		}else{
			return 0;
		}

    }
    
    // seems to work:
    private double coalesce(int currTreeInterval, int speciesInterval) {
    	// the probability that no coalescent event was observed in the last interval
    	double intervalProbability = normalizeLineages();    	
    	
    	// Get the daughter lineages
		List<Node> coalLines = lineagesRemoved[currTreeInterval];
		List<Node> parentLineage = lineagesAdded[currTreeInterval];
		
    	if (coalLines.get(0) == null) {
    		// TODO throw exception
			System.err.println("daughter nodes not found");
			System.exit(0);
		}
    	
    	//  get the indices of the daughter lineages
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {			
			System.out.println("active lineages: "+ activeLineages + "\n lin1 " 
					+ coalLines.get(0).getNr() + " lin2 " + coalLines.get(1).getNr() + " parent " + parentLineage.get(0).getNr() 
					+ "\n " + coalLines.get(1).getParent().getNr());
			System.err.println("daughter lineages at coalescent event not found");
			first = true;
//			System.exit(0);
			return Double.NEGATIVE_INFINITY;
		}
		double[] coalProb = new double[states];
		
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */		
        for (int k = 0; k < states; k++) { 
        	// calculate the coalscent rate. The 2 originates from lineage i coalescing with lineage j
       		// AND vice versa.
			final double pairCoalRate = coalescentRates[k] * 2 *
					(linProbs[sampleState.get(daughterIndex1)*states + k] * linProbs[sampleState.get(daughterIndex2)*states + k]);	
			
			if (!Double.isNaN(pairCoalRate)){
				coalProb[k] = pairCoalRate;
			}
			else{
				return Double.NEGATIVE_INFINITY;
			}
        }
        
        int maxstate=0;
        double maxval = 0.0;
        for (int k = 0; k < states; k++){
        	if(coalProb[k]>maxval){
        		maxstate=k;
        		maxval = coalProb[k];
        	}
        }
        
        // save the most likely state of coalescence
        mostLikelyState[parentLineage.get(0).getNr()] = popModelInput.get().getSpeciesState(speciesInterval, maxstate);
        
        activeLineages.add(parentLineage.get(0).getNr());
                      
        boolean sameOri;
        if (multiplicator[sampleState.get(daughterIndex1)]
        		==multiplicator[sampleState.get(daughterIndex2)]) sameOri=true;
        else sameOri = false;
        boolean remove1,remove2;
        
        // check how the size of the linProbs vector changes
        if (multiplicator[sampleState.get(daughterIndex1)]==1){
        	remove1 = true;
        }else{
        	remove1 = false;
        	multiplicator[sampleState.get(daughterIndex1)] -= 1;
        }
        if (multiplicator[sampleState.get(daughterIndex2)]<1){
        	System.err.println("seems like a lineage coalesces with itself");
        }
        if (multiplicator[sampleState.get(daughterIndex2)]==1){
        	remove2 = true;
        }else{
        	remove2 = false;
        	multiplicator[sampleState.get(daughterIndex2)] -= 1;
        }
          
        
        // init the arrays that will hold the new lin probs and multiplicators
        double[] newLinProbs;
        int[] newMultiplicator;
        
        if (remove1 && remove2){
    		newLinProbs = new double[linProbs.length-states];
    		newMultiplicator = new int[multiplicator.length-1];
        }else if (remove1 || remove2){
			newLinProbs = new double[linProbs.length];
    		newMultiplicator = new int[multiplicator.length];       	
        }else{
			newLinProbs = new double[linProbs.length+states];
    		newMultiplicator = new int[multiplicator.length+1];       	
        }     
        
        
        // ad the old information about the 
		int newIndex = 0;
		for (int i = 0; i < linProbs.length/states; i++){
			// do not add this to the old states
			if ((i==sampleState.get(daughterIndex1)||i==sampleState.get(daughterIndex2)) 
					&& sameOri && remove2){
			}else if(i==sampleState.get(daughterIndex1) && remove1){
				
			}else if(i==sampleState.get(daughterIndex2) && remove2){
				
			}else {
				for (int j = 0; j < states; j++){
					newLinProbs[states*newIndex+j] = linProbs[states*i+j];
					newMultiplicator[newIndex] = multiplicator[i];
				}
				newIndex++;				
			}
		}
		
		// add the new lineage
		newMultiplicator[newMultiplicator.length-1] = 1;
		double sumProb = 0.0;
		for (int j = 0; j < states; j++)
			sumProb += coalProb[j];
		for (int j = 0; j < states; j++)
			newLinProbs[newLinProbs.length-states+j] += coalProb[j]/sumProb;

		
        // update the sampleStates
        sampleState.add(multiplicator.length);            
        
        // make state mapping consistent		
        for (int i = 0; i < sampleState.size(); i++){
        	if(i!=daughterIndex1 && i!=daughterIndex2){
        		if ((sampleState.get(i) > sampleState.get(daughterIndex1) && remove1)
        				&& (sampleState.get(i) > sampleState.get(daughterIndex2) && remove2)){
        			sampleState.set(i, sampleState.get(i)-2);
        		}else if ((sampleState.get(i) > sampleState.get(daughterIndex1) && remove1)
        				|| (sampleState.get(i) > sampleState.get(daughterIndex2) && remove2)){
        			sampleState.set(i, sampleState.get(i)-1);
        		}
        	}
        }

		
		//Remove daughter lineages
		if (daughterIndex1>daughterIndex2){
			activeLineages.remove(daughterIndex1);
			sampleState.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);
			sampleState.remove(daughterIndex2);
		}else{
			activeLineages.remove(daughterIndex2);
			sampleState.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);
			sampleState.remove(daughterIndex1);
		}
     
		if (sumProb<0.0)
			System.err.println("Coalescent probability smaller than 0");
					
		// set new to old
		linProbs = newLinProbs;
		multiplicator = newMultiplicator;		
		
		// store the node
		storeNode(currTreeInterval - nodeCount, linProbs,
					multiplicator, logP + Math.log(sumProb) + intervalProbability, speciesInterval,
					activeLineages, sampleState);
		
		
    	return Math.log(sumProb) + intervalProbability;
    }   
     
    private void updateRatesList(int speciesInterval){
    	// initialize coalescent and migration rates
    	coalescentRates = new double[states];
    	migrationRates = new double[states][states];
    	//TODO check ploidity
    	for (int i = 0; i < states; i++){
    		coalescentRates[i] =  
    				1/(ploidyInput.get()*popModelInput.get().getPopulationSize(speciesInterval, i));
    	}
    	
    	for (int i = 0; i < states; i++){
    		for (int j = 0; j < states; j++){
    			if (i != j)
    				migrationRates[i][j] = popModelInput.get().getMigrationRates(speciesInterval, i, j);
    		}
    	}    	
    }
    
    // Could work
    private void combineStates(int speciesInterval, double probability, int geneInterval){
    	popModelInput.get().coalesce(speciesInterval);
    	// Get the two states that will be combined after the species tree coalescent event
    	int state1 = popModelInput.get().getDaughter1(speciesInterval);
    	int state2 = popModelInput.get().getDaughter2(speciesInterval);
    	
    	// make sure that state2 is larger than state1
    	if (state2<state1){
    		int tmp = state1;
    		state1 = state2;
    		state2 = tmp;
    	}

    	// skip the states involved in speciation    	
    	double[] newLinProbs = new double[multiplicator.length*states];
    	for (int i = 0; i < multiplicator.length; i++){
        	int counter = 0;
    		for (int j = 0; j < (states+1); j++){
    			if (j!=state1 && j!=state2){
    				newLinProbs[states*i+counter] = linProbs[(states+1)*i+j];
    				counter++;
    			}  			
    		}
    	}   	
   	
    	// add the combined states at the end
    	for (int i = 0; i < multiplicator.length; i++){
    		if (((states+1)*i+state1)<0){
    			System.exit(0);
    			System.out.println("daughter lineages not found at speciation event");
    			System.out.println(states + " " + state1 + " " + state2);
    			
    		}

			newLinProbs[states*i+(states-1)] = 
					linProbs[(states+1)*i+state1]
							+linProbs[(states+1)*i+state2];
    	}    	
    	
    	
    	// make new to old
    	linProbs = newLinProbs;   	
//    	System.out.println("linProbs after speciation " + Arrays.toString(linProbs));
    	
   }
        
    private void storeNode(int storingInterval, double[] storeLinProbs,
    		int[] storeMultiplicator, double probability, int storeSpeciesInterval,
    		ArrayList<Integer> storeActiveLineages, ArrayList<Integer> storeSampleState){    	
    	
//    	coalLinProbs[storingInterval] = Arrays.copyOf(storeLinProbs, storeLinProbs.length);
//    	coalMultiplicator[storingInterval] = Arrays.copyOf(storeMultiplicator, storeMultiplicator.length);
//    	coalLogP[storingInterval] = probability;		
//    	coalSpeciesInterval[storingInterval] = storeSpeciesInterval;   	
//   	
//    	// copy the array lists
//    	ArrayList<Integer> tmp_activeLins = new ArrayList<>(storeActiveLineages);
//    	ArrayList<Integer> tmp_sampleState = new ArrayList<>(storeSampleState);
//    	coalActiveLineages.set(storingInterval, tmp_activeLins);
//    	coalSampleState.set(storingInterval, tmp_sampleState);
    }
        
    private int restoreNode(int restoringInterval){   
//    	System.out.println(restoringInterval);
//    	linProbs = Arrays.copyOf(coalLinProbs[restoringInterval], coalLinProbs[restoringInterval].length);
//    	multiplicator = Arrays.copyOf(coalMultiplicator[restoringInterval], coalMultiplicator[restoringInterval].length);
//    	logP = coalLogP[restoringInterval];    	
//    	activeLineages = new ArrayList<>(coalActiveLineages.get(restoringInterval));
//    	sampleState = new ArrayList<>(coalSampleState.get(restoringInterval));
//    	
    	return 0;
//    	return coalSpeciesInterval[restoringInterval];
    }
    
    @Override
	public void store(){
    	// store the intermediate results
//    	storeLinProbs = new double[coalLinProbs.length][];
//    	storeMultiplicator  = new int[coalMultiplicator.length][];
//    	for (int i = 0; i < coalLinProbs.length; i++)
//    		storeLinProbs[i] = Arrays.copyOf(coalLinProbs[i], coalLinProbs[i].length);
//    	for (int i = 0; i < coalMultiplicator.length; i++)
//    		storeMultiplicator[i] = Arrays.copyOf(coalMultiplicator[i], coalMultiplicator[i].length);
//    	
//        storeLogP = Arrays.copyOf(coalLogP, coalLogP.length);
//        storeSpeciesInterval = Arrays.copyOf(coalSpeciesInterval, coalSpeciesInterval.length);
//        
//        
//        storeActiveLineages = new ArrayList<>();
//		for (int i= 0; i < coalActiveLineages.size(); i++){
//			ArrayList<Integer> add = new ArrayList<>();
//			for (int j = 0; j < coalActiveLineages.get(i).size(); j++)
//				add.add(coalActiveLineages.get(i).get(j));
//			storeActiveLineages.add(add);			
//		}
//		
//		storeSampleState = new ArrayList<>();
//		for (int i= 0; i < coalSampleState.size(); i++){
//			ArrayList<Integer> add = new ArrayList<>();
//			for (int j = 0; j < coalSampleState.get(i).size(); j++)
//				add.add(coalSampleState.get(i).get(j));
//			storeSampleState.add(add);			
//		}        
//        
//        
//        storedFirst = first;
    	// store the 
    	super.store();
    }
        
    @Override
	public void restore(){
    	// store the intermediate results
//    	coalLinProbs = Arrays.copyOf(storeLinProbs, storeLinProbs.length);
//    	coalMultiplicator = Arrays.copyOf(storeMultiplicator, storeMultiplicator.length);
//    	coalLogP = Arrays.copyOf(storeLogP, storeLogP.length);
//    	coalSpeciesInterval = Arrays.copyOf(storeSpeciesInterval, storeSpeciesInterval.length);
//    	coalStates = Arrays.copyOf(storeStates, storeStates.length);
//    	coalActiveLineages = storeActiveLineages;
//    	coalSampleState = storeSampleState;
    	
    	
//    	coalActiveLineages = new ArrayList<>();
//		for (int i= 0; i < storeActiveLineages.size(); i++){
//			ArrayList<Integer> add = new ArrayList<>();
//			for (int j = 0; j < storeActiveLineages.get(i).size(); j++)
//				add.add(storeActiveLineages.get(i).get(j));
//			coalActiveLineages.add(add);			
//		}
//		
//		coalSampleState = new ArrayList<>();
//		for (int i= 0; i < storeSampleState.size(); i++){
//			ArrayList<Integer> add = new ArrayList<>();
//			for (int j = 0; j < storeSampleState.get(i).size(); j++)
//				add.add(storeSampleState.get(i).get(j));
//			coalSampleState.add(add);			
//		}        
//
//    	
//    	
//    	first = storedFirst;
//    	calculateIntervals();
    	// store the 
    	super.restore();
    }
        
	public String getType(){
   		return "state";    	
    }
            
	public int getMaxState(int nodeNr){
		return mostLikelyState[nodeNr];
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

	
	/*
	 * tree intervals
	 */
			
	
    protected double[] intervals;
    protected double[] storedIntervals;
    protected boolean[] isCoalescent;
    protected boolean[] intervalIsDirty;


    
    /**
     * parent and daugheter lineages at each coalescent event
     */
    protected List<Node>[] lineagesAdded;
    protected List<Node>[] lineagesRemoved;
    protected List<Node>[] storedLineagesAdded;
    protected List<Node>[] storedLineagesRemoved;

    final static NodeHeightComparator nhc = new NodeHeightComparator();

    
    @SuppressWarnings({ "unchecked", "deprecation" })
    private void calculateIntervals() {
//    	final Tree tree = treeInput.get()
    	Node[] treeNodesTmp = treeInput.get().getNodesAsArray();
    	Node[] treeNodes = new Node[treeNodesTmp.length];
    	System.arraycopy(treeNodesTmp, 0, treeNodes, 0, treeNodesTmp.length);
    	
    	
    	Arrays.sort(treeNodes, nhc);
    	
 	
        intervals = new double[treeNodes.length];
        isCoalescent = new boolean[treeNodes.length];
        lineagesAdded = new List[treeNodes.length];
        lineagesRemoved = new List[treeNodes.length];
        
        //TODO init
        
        intervals[0] = treeNodes[0].getHeight();
        lineagesAdded[0] = new ArrayList<>(); 
		lineagesAdded[0].add(treeNodes[0]);
		isCoalescent[0] = false;
    	for (int i = 1; i < treeNodes.length; i++){
    		intervals[i] = treeNodes[i].getHeight()-treeNodes[i-1].getHeight();
            lineagesAdded[i] = new ArrayList<>(); 
            lineagesAdded[i].add(treeNodes[i]);
    		if (!treeNodes[i].isLeaf()){
    			isCoalescent[i] = true;
    			lineagesRemoved[i] = new ArrayList<>(); 
	    		lineagesRemoved[i].add(treeNodes[i].getLeft());
	    		lineagesRemoved[i].add(treeNodes[i].getRight());    	
	    		if(treeNodes[i].getLeft().getParent().getNr() != treeNodes[i].getRight().getParent().getNr()){
	    			System.err.println("parent offspring error " + treeInput.get().getID());
	    			System.err.println(treeNodes[i].getLeft().getParent().getNr() + " " + treeNodes[i].getRight().getParent().getNr());
	    			System.err.println(treeInput.get().getNode(treeNodes[i].getNr()).getLeft().getParent().getNr());
	    			System.exit(0);
	    		}
    		}else{
    			isCoalescent[i] = false;
    			lineagesRemoved[i] = new ArrayList<>(); 
	    		lineagesRemoved[i].add(null);
	    		lineagesRemoved[i].add(null);    		
    			
    		}
    	}    	
    }
    
    protected boolean intervalIsDirty(int i) {
    	return true;
    	
    	//TODO change again
//    	Node node = treeInput.get().getNode(i);
//    	if (node.isRoot()){
//    		return false;
//    	}else{
//    		if(node.getParent().isDirty()>0)
//    			return true;
//    		else
//    			return false;
//    		
//    	}
    }    

    
}
