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

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;

import java.util.*;


/**
 * @author Nicola Felix Mueller
 */
// Gene tree type trait nested in species tree as taxon superset
@Description("Calculate the probability of a tree under the structured coalescent assuming lineage independence with constant rates" +
		" as described in Mueller et. al.,2017. The Input rates are backwards in time migration rates" +
		" and pairwise coalescent rates translating to 1/Ne for m=1 in the Wright Fisher model")
@Citation("Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2017)\n  The Structured Coalescent and its Approximations.\n  Mol Biol Evol 2017 msx186. doi: 10.1093/molbev/msx186")
@Citation("Nicola Felix Mueller, Huw Ogilvie, Chi Zhang, Alexei Drummond, Tanja Stadler (2018)\n  Inference of species histories in the presence of gene flow.\n  BioRxiv 348391; doi: https://doi.org/10.1101/348391")
public class GeneTreeWithMigration extends Distribution {
    public Input<Tree> treeInput = new Input<>("tree", "The gene tree.", Validate.REQUIRED);
    public Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for this gene, typically a whole number or half (default is 4).", 4.0);
    public Input<ConstantWithGeneFlow> popModelInput = new Input<>("populationModel", "Population model used to infer the multispecies coalescent probability for this gene");
   
	public int samples;
	public int nrSamples;
	
    SpeciesTreeInterface speciesTree;
	
//    private boolean hasIndicators;
    
	private double[] coalescentRates;
	private double[][] migrationRates;
	private int[][] indicators;
	private boolean[] isConnected;
	private double[] linProbs;
	private int[] multiplicator;
	private int[] mostLikelyState;
	
    public int states;    
    private boolean recalculateLogP;    
    private double maxTolerance;            
    private int nr_lineages;  
    
    // store the linProbs, multiplicators and logP's at coalescent points in jagged arrays from last time
//    private double[][] coalLinProbs;
//    private int[][] coalMultiplicator;
//    private double[] coalLogP;
//    private int[] coalSpeciesInterval;
//    private ArrayList<ArrayList<Integer>> coalActiveLineages;
//    private ArrayList<ArrayList<Integer>> coalSampleState;
//    
//    // deep store the things above for MCMC
//    private double[][] storeLinProbs;
//    private int[][] storeMultiplicator;
//    private double[] storeLogP;
//    private int[] storeSpeciesInterval;
//    private ArrayList<ArrayList<Integer>> storeActiveLineages;
//    private ArrayList<ArrayList<Integer>> storeSampleState;
//    private boolean storedFirst;
    
        
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
        speciesTree = popModelInput.get().migrationModelInput.get().speciesTreeInput.get();
//    	popModelInput.get().calculateIntervals();
    	// Calculate the tree intervals (time between events, which nodes participate at a event etc.)
    	calculateIntervals(); 
                
    	// initialize the maximum relative tolerance 
    	maxTolerance = 1e-3;
    	    	
    	
    	int MAX_SIZE_TMP = popModelInput.get().migrationModelInput.get().speciesTreeInput.get().getLeafNodeCount() *
    			treeInput.get().getLeafNodeCount();
    	int MAX_SIZE = Math.max(MAX_SIZE_TMP, 
    			popModelInput.get().migrationModelInput.get().speciesTreeInput.get().getLeafNodeCount()*
    			popModelInput.get().migrationModelInput.get().speciesTreeInput.get().getLeafNodeCount())
    			+1;

    	linProbs_tmp = new double[MAX_SIZE];
    	linProbs_tmpdt = new double[MAX_SIZE];
    	linProbs_tmpddt = new double[MAX_SIZE];
    	linProbs_tmpdddt = new double[MAX_SIZE];

        // make sure we are not in BEAUti
        final int speciesNodeCount = speciesTree.getNodeCount();
        if (speciesNodeCount != 1) {
            calculateLogP();
        }
    }

    double [] linProbs_tmp;
    double [] linProbs_tmpdt;
    double [] linProbs_tmpddt;
    double [] linProbs_tmpdddt;
    
    public double calculateLogP() { 
    	if (!popModelInput.get().checkMaxRates())
    		return Double.NEGATIVE_INFINITY;
    	
//    	if(treeInput.get().getID().contentEquals("Tree.t:chr3L-10165")){
//    		System.out.println(treeInput.get());
//        	System.exit(0);   		
//    	}
    	
    	
//    	System.out.println(popModelInput.get().speciesTree);
//    	System.out.println("new tree");
//    	popModelInput.get().calculateIntervals();
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
        // initialize the array that saves the most likely state of each node
        mostLikelyState = new int[2*nr_lineages+1];	        
        
        linProbs = new double[0];
        
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
		nextGeneTime = intervals[geneInterval];
		
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
	        		logP += doEuler(duration);	
                }else{
                	// use analytical solution instead
                	logP -= duration * coalescentRates[0] * activeLineages.size()*(activeLineages.size()-1)/2;	                
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
//        		System.out.println(nextSpeciesTime);
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
        	if (logP == Double.NEGATIVE_INFINITY)
        		return logP;

        }while(activeLineages.size()>1);
        first = false;  
//        System.exit(0);
        return logP;  	
    }   
    
	private double doEuler(double nextEventTime) {
		Euler2ndOrderAIM euler;
		
		double logVal = 0.0;
		int inactiveStates = 0;
		
		// add the prob of the unstructured parts
		int[] nrLins = new int[states];
		if (popModelInput.get().indicatorInput.get()!=null){
			for (int j = 0; j < states; j++){
				if (!isConnected[j]){
					inactiveStates++;
					for (int i = 0; i < multiplicator.length; i++){
						if (linProbs[states*i+j]>0.999999999999999){
							nrLins[j] += multiplicator[i];
						}else if (linProbs[states*i+j]>0.000000000000001){
							System.out.println("error in the lin probs of unconnected state, return negative infinity");
							logVal = Double.NEGATIVE_INFINITY;
						}	
					}
				}
			}
			for (int j = 0; j < states; j++){
				if (!isConnected[j]){
					logVal -= nextEventTime*coalescentRates[j]* nrLins[j]*(nrLins[j]-1)/2;
				}
			}
		}			

		double[] coaltmp = new double[coalescentRates.length];
		System.arraycopy(coalescentRates,0,coaltmp,0,coalescentRates.length);
		
		System.arraycopy(linProbs,0,linProbs_tmp,0,linProbs.length);
		linProbs_tmp[linProbs.length] = 0;

		if (popModelInput.get().indicatorInput.get()!=null){
			if (inactiveStates<states){
				euler = new Euler2ndOrderAIM(multiplicator, migrationRates, indicators, isConnected, coalescentRates, multiplicator.length , states, 0.001, 0.2);
				euler.calculateConnectedValues(nextEventTime, linProbs_tmp, linProbs_tmpdt, linProbs_tmpddt, linProbs_tmpdddt, linProbs.length + 1);
			}
		}else{
    		euler = new Euler2ndOrderAIM(multiplicator, migrationRates, coalescentRates, multiplicator.length , states, 0.001, 0.2);
			euler.calculateValues(nextEventTime, linProbs_tmp, linProbs_tmpdt, linProbs_tmpddt, linProbs_tmpdddt, linProbs.length + 1);
		}
		
		if (Double.isNaN(linProbs_tmp[linProbs.length])){
			return Double.NEGATIVE_INFINITY;
		}
		
		System.arraycopy(linProbs_tmp,0,linProbs,0,linProbs.length);

		
		return linProbs_tmp[linProbs.length] + logVal;
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
		final Map<String, Integer> tipNumberMap = speciesTree.getTipNumberMap();
		final int sampledSpeciesNr = tipNumberMap.get(geneLeafName);
		final String sampledSpeciesName =speciesTree.getNode(sampledSpeciesNr).getID();
		sampleState.add(popModelInput.get().getSampleState(sampledSpeciesName));
		mostLikelyState[incomingLines.get(0).getNr()] = popModelInput.get().getSampleState(sampledSpeciesName);
    }
 
    private void initializeP(){
    	final int nrSpecies = speciesTree.getLeafNodeCount();
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
    	// Get the daughter lineages
		List<Node> coalLines = lineagesRemoved[currTreeInterval];
		List<Node> parentLineage = lineagesAdded[currTreeInterval];
		
    	if (coalLines.get(0) == null) {
    		// TODO throw exception
    		System.err.println(treeInput.get());
			System.err.println("daughter nodes not found, possible multifurcations in the starting tree");
			throw new IllegalArgumentException();
		}	
    	
    	
    	//  get the indices of the daughter lineages
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {			
			System.out.println("active lineages: "+ activeLineages + "\n lin1 " 
					+ coalLines.get(0).getNr() + " lin2 " + coalLines.get(1).getNr() + " parent " + parentLineage.get(0).getNr() 
					+ "\n " + coalLines.get(1).getParent().getNr());
			System.err.println("daughter lineages at coalescent event not found");
			System.out.println(treeInput.get());
			first = true;
			return Double.NEGATIVE_INFINITY;
		}
		double[] coalProb = new double[states];
		
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */		
        for (int k = 0; k < states; k++) { 
        	// calculate the coalscent rate. 
        	Double pairCoalRate = coalescentRates[k] *
					(linProbs[sampleState.get(daughterIndex1)*states + k] * linProbs[sampleState.get(daughterIndex2)*states + k]);
        				
			if (!Double.isNaN(pairCoalRate)){
				coalProb[k] = pairCoalRate;
			}
			else{
				System.out.println("allalalala");
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
		
//		// store the node
//		storeNode(currTreeInterval - nodeCount, linProbs,
//					multiplicator, logP + Math.log(sumProb) + intervalProbability, speciesInterval,
//					activeLineages, sampleState);
		
		if (sumProb<=0){
//			System.out.println("sum prob");
			return Double.NEGATIVE_INFINITY;
		}else{
			return Math.log(sumProb);
		}
    }   
     
    private void updateRatesList(int speciesInterval){
    	// initialize coalescent and migration rates
    	coalescentRates = new double[states];
    	migrationRates = popModelInput.get().getMigrationRates(speciesInterval);
    	isConnected = new boolean[states];
    	
    	//TODO check ploidity
    	for (int i = 0; i < states; i++){
    		coalescentRates[i] =  
    				1/(ploidyInput.get()*popModelInput.get().getPopulationSize(speciesInterval, i));
    	}    	
    	
    	if (popModelInput.get().indicatorInput.get()!=null){
    		indicators = popModelInput.get().getIndicatorsRates(speciesInterval);
    		for (int i = 0; i < states; i++){
    			isConnected[i] = popModelInput.get().getIsConnected(speciesInterval, i);
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
    			System.out.println("daughter lineages not found at speciation event");
    			System.out.println(states + " " + state1 + " " + state2);
    			System.exit(0);    			
    		}

			newLinProbs[states*i+(states-1)] = 
					linProbs[(states+1)*i+state1]
							+linProbs[(states+1)*i+state2];
    	}    	
    	
    	
    	// make new to old
    	linProbs = newLinProbs;   	
//    	System.out.println("linProbs after speciation " + Arrays.toString(linProbs));
    	
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
     * parent and daughter lineages at each coalescent event
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
    }    

    
}
