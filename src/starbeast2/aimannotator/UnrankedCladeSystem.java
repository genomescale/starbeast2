package starbeast2.aimannotator;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.DiscreteStatistics;

import java.io.PrintStream;
import java.util.*;

/**
 * extracted from TreeAnnotator
 */
//TODO merge with CladeSet?
public class UnrankedCladeSystem {
	private boolean processSA = true;

    protected List<UnrankedTree> unrankedTrees = new ArrayList<>();
    protected List<UnrankedTree> orderedUnrankedTrees = new ArrayList<>();
    protected UnrankedTree newRankedTree;

    public UnrankedCladeSystem() { }

    public UnrankedCladeSystem(Tree targetTree, Set<String> attributeNames) {
        add(targetTree, true, attributeNames);
    }
    
    /**
     * adds all the clades in the tree
     */
    public int add(Tree tree, boolean includeTips, Set<String> attributeNames) {
    	getCurrentRankedTree(tree, includeTips, attributeNames);
    	// get all node heights for ranking the nodes
        // compare the new ranked tree to the current ranked trees
        int rtIndex = unrankedTrees.indexOf(newRankedTree);
        if (rtIndex==-1) {
        	unrankedTrees.add(newRankedTree);
        	unrankedTrees.get(unrankedTrees.size()-1).setCount(1);
        	unrankedTrees.get(unrankedTrees.size()-1).tree = tree.copy();
        	return unrankedTrees.size()-1;
        }else {
        	unrankedTrees.get(rtIndex).setCount(unrankedTrees.get(rtIndex).getCount()+1);
        	// add the attribute values
        	addAttributes(rtIndex);
        	return rtIndex;
        }       
    }
    
    private void getCurrentRankedTree(Tree tree, boolean includeTips, Set<String> attributeNames) {
    	Node[] nodes = tree.getNodesAsArray();
    	List<Double> nodeHeights = new ArrayList<>();
    	for (Node n : nodes) {
    		if (!n.isLeaf()) {
    			nodeHeights.add(n.getHeight());
    		}
    	}
    	
    	// puts the node heights in the correct order, otherwise orderred according to nr
    	Collections.sort(nodeHeights);
    	    	
        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
    	newRankedTree = new UnrankedTree();
        addClades(tree.getRoot(), includeTips, nodeHeights, attributeNames);

    }

    private BitSet addClades(Node node, boolean includeTips, List<Double> nodeHeights, Set<String> attributeNames) {

        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            bits.set(2*index);

            if (includeTips) {
                addClade(bits, -1, node, attributeNames);
            }

        } else {

            for (int i = 0; i < node.getChildCount(); i++) {

                Node node1 = node.getChild(i);

                bits.or(addClades(node1, includeTips, nodeHeights, attributeNames));
            }

            for (int i=1; i<bits.length(); i=i+2) {
                bits.set(i, false);
            }
            if (node.isFake() && processSA) {
                int index = getTaxonIndex(node.getDirectAncestorChild());
                bits.set(2 * index + 1);
            }
            int rank = nodeHeights.indexOf(node.getHeight());
            addClade(bits, rank, node, attributeNames);
        }

        return bits;
    }

    private void addClade(BitSet bits, int rank, Node node, Set<String> attributeNames) {
    	UnrankedBitSet rbs = new UnrankedBitSet(bits, node.getNr());
    	collectAttributesForClade(rbs, node, attributeNames);
    	newRankedTree.addRBS(rbs);
	}   

    public void calculateCladeCredibilities(int totalTreesUsed) {
        for (UnrankedTree rankedTree : unrankedTrees) {
            rankedTree.setCredibility(((double) rankedTree.getCount()) / (double) totalTreesUsed);
        }
    }
    
    public int[] orderRankedTrees() {
    	int[] order = new int[unrankedTrees.size()];
    	// get the order of trees by rank
    	for (int i = 0; i < unrankedTrees.size(); i++) {
    		int highestCount = -1;
    		int highestIndex = -1;
        	for (int j = 0; j < unrankedTrees.size(); j++) {
	    		if (unrankedTrees.get(j).count>highestCount) {
	    			highestCount = unrankedTrees.get(j).count;
	    			highestIndex = j;
	    		}
	    	}
        	order[i] = highestIndex;
        	unrankedTrees.get(highestIndex).count=-1;
    	}
    	return order;
    }
    
    private void collectAttributesForClade(UnrankedBitSet rbs, Node node, Set<String> attributeNames) {
        if (rbs.attributeValues == null) {
        	rbs.attributeValues = new ArrayList<>();
        }

        int i = 0;
        Object[] values = new Object[attributeNames.size()];
        for (String attributeName : attributeNames) {
            Object value;
            switch (attributeName) {
                case "height":
                    value = node.getHeight();
                    break;
                case "length":
                    value = getBranchLength(node);
                    break;
                default:
                    value = node.getMetaData(attributeName);
                    if (value instanceof String && ((String) value).startsWith("\"")) {
                        value = ((String) value).replaceAll("\"", "");
                    }
                    break;
            }

            values[i] = value;

            i++;
        }
        rbs.attributeValues.add(values); 
    }

    private void addAttributes(int rtIndex) {
    	for (UnrankedBitSet newRBS : newRankedTree.rankedSet) {
    		int oldRBSindex = unrankedTrees.get(rtIndex).rankedSet.indexOf(newRBS);
    		UnrankedBitSet rbs = unrankedTrees.get(rtIndex).rankedSet.get(oldRBSindex);
    		
	        if (rbs.attributeValues == null) {
	        	throw new IllegalArgumentException("attributes should not be empty");
	        }
	        rbs.attributeValues.addAll(newRBS.attributeValues);   
    	}
    }

    private Object getBranchLength(Node node) {
        if (node.isRoot()) {
            return 0;
        }
        return node.getParent().getHeight() - node.getHeight();
    }
        
    public int getTaxonIndex(Node node) {
        return node.getNr();
    }
    
    public Tree compactMetaData(int index, Set<String> attributeNames, boolean useMean) {
    	Tree tree = new Tree(unrankedTrees.get(index).tree.getRoot());
    	
    	// build the species number to clade mapping
    	int[][] oriSpeciesNr = new int[unrankedTrees.get(index).rankedSet.size()][];
    	for (int i = 0; i < oriSpeciesNr.length;i++) {
            List<Object[]> attributeValues = unrankedTrees.get(index).rankedSet.get(i).getAttributeValues();
            int j=0;
            for (String attributeName : attributeNames) {
         		if (attributeName.contentEquals("species")) {
         			oriSpeciesNr[i] = new int[attributeValues.size()];
         			for (int k = 0; k< attributeValues.size(); k++) {
         				Double val = (Double) attributeValues.get(k)[j];
         				oriSpeciesNr[i][k] = (int) Math.round(val);
         			}
         		}
         		j++;
            }
    	}
    	
    	
    	for (int i = 0; i < unrankedTrees.get(index).rankedSet.size();i++) {
    		int nodeNr = unrankedTrees.get(index).rankedSet.get(i).nodeNr;
    		Node n = tree.getNode(nodeNr);
    		n.metaDataString = "species=" + oriSpeciesNr[i][0];
    		
            int j=0;
            
            Integer[][] to = new Integer[0][0];
            Double[][] allRates = new Double[0][0];
            
            List<Object[]> attributeValues = unrankedTrees.get(index).rankedSet.get(i).getAttributeValues();
            for (String attributeName : attributeNames) {
            	// get the height
         		if (attributeName.contentEquals("height")) {
         			double[] vals = get95interval(j, attributeValues); 
         			if (useMean)
         				n.setHeight(vals[0]);
         			else
         				n.setHeight(vals[1]);
         			
         			n.metaDataString = n.metaDataString + ",height_mean=" + vals[0];
         			n.metaDataString = n.metaDataString + ",height_median=" + vals[1];
         			n.metaDataString = n.metaDataString + ",height_95%_HPD={" + vals[2] + "," + vals[3] + "}";
         		}
         		if (attributeName.contentEquals("Ne")) {
         			double[] vals = get95interval(j, attributeValues); 
         			n.metaDataString = n.metaDataString + ",Ne_mean=" + vals[0];
         			n.metaDataString = n.metaDataString + ",Ne_median=" + vals[1];
         			n.metaDataString = n.metaDataString + ",Ne_95%_HPD={" + vals[2] + "," + vals[3] + "}";
         		}    
         		if (attributeName.contentEquals("rates")) {
         			if (!n.isRoot()) {
         				allRates = getRates(j, attributeValues);
         			}
         		}  
         		if (attributeName.contentEquals("to")) {
         			if (!n.isRoot()) {
	         			to = getTo(j, attributeValues, oriSpeciesNr);
         			}
         		}
            	j++;	
            }
            
            // compute the mean rates etc. based on the "to" and "allRates" values
 			if (!n.isRoot()) {	
 				List<Integer> tovalues = getIntersectionToVals(to);
 				Collections.sort(tovalues);

     			// compute posterior support for rates
     			n.metaDataString = n.metaDataString + ",rates_mean={" + getMean(allRates, tovalues.get(0), to);
     			for (int k = 1; k < tovalues.size();k++)
     				n.metaDataString = n.metaDataString + "," + getMean(allRates, tovalues.get(k), to);
     			n.metaDataString = n.metaDataString + "}";

     			// compute posterior support for rates
     			n.metaDataString = n.metaDataString + ",rates_posterior={" + getPost(allRates, tovalues.get(0), to);
     			for (int k = 1; k < tovalues.size();k++)
     				n.metaDataString = n.metaDataString + "," + getPost(allRates, tovalues.get(k), to);
     			n.metaDataString = n.metaDataString + "}";
 				
     			
     			// compute posterior support for rates
     			n.metaDataString = n.metaDataString + ",rates_to={" + tovalues.get(0);
     			for (int k = 1; k < tovalues.size();k++)
     				n.metaDataString = n.metaDataString + "," + tovalues.get(k);
     			n.metaDataString = n.metaDataString + "}";
 			}
    	}
    	return tree;
    }
    
    private double[] get95interval(int j, List<Object[]> attributeValues) {
		double[] array = new double[attributeValues.size()]; 
		for (int i = 0; i < attributeValues.size(); i++) {
			try {
				array[i] = (double) attributeValues.get(i)[j];
			}
			catch(Exception e) {
				System.out.println(attributeValues.get(i)[j]);
			}
		}
		
		double[] vals = new double[4];
		
		vals[0] = DiscreteStatistics.mean(array);
		vals[1] = DiscreteStatistics.median(array);
        Arrays.sort(array);
        vals[2] = array[(int)(0.025 * array.length)];
        vals[3] = array[(int)(0.975 * array.length)];

		
		return vals;
    }
    
    private Double[][] getRates(int j, List<Object[]> attributeValues) {
    	// get the maximal size of the rates vector
    	int max_val = 0;
    	Double[] rates;
		for (int i = 0; i < attributeValues.size(); i++) {
			rates = (Double[]) attributeValues.get(i)[j];
			max_val = Math.max(max_val, rates.length);
		} 

		Double[][] allRates = new Double[max_val][attributeValues.size()];
    	
		for (int i = 0; i < attributeValues.size(); i++) {
	    	rates = (Double[]) attributeValues.get(i)[j];
	    	for (int k = 0; k < rates.length; k++) {
	    		allRates[k][i] = (Double) rates[k];
	    	}
		}
		return allRates;
    }
  
    private Integer[][] getTo(int j, List<Object[]> attributeValues, int[][] oriSpeciesNr) {
    	// get the maximal size of the rates vector
    	int max_val = 0;
    	Double[] to;
		for (int i = 0; i < attributeValues.size(); i++) {
			to = (Double[]) attributeValues.get(i)[j];
			max_val = Math.max(max_val, to.length);
		} 

		Integer[][] allTo = new Integer[max_val][attributeValues.size()];

		for (int i = 0; i < attributeValues.size(); i++) {
	    	to = (Double[]) attributeValues.get(i)[j];
	    	for (int k = 0; k < to.length; k++) {
	    		allTo[k][i] = getOriSpecies(oriSpeciesNr,i,(int) Math.round(to[k]));
	    	}
		}
		

		return allTo;
    }
    
    private int getOriSpecies(int[][] oriSpeciesNr, int i, int toval) {
    	for (int k = 0; k < oriSpeciesNr.length; k++) {
    		if (oriSpeciesNr[k][i]==toval)
    			return oriSpeciesNr[k][0];
    	}
    	return -1;    	
    }
  
    private List<Integer> getIntersectionToVals(Integer[][] to) {
    	List<Integer> intersectionVals = new ArrayList<>();
		for (int i = 0; i < to.length; i++) {
			if (to[i][0]!=null)
				intersectionVals.add(to[i][0]);
		}
		for (int j = 1; j < to[0].length; j++) {
	    	List<Integer> newVals = new ArrayList<>();

			for (int i = 0; i < to.length; i++) {
				if (to[i][j]!=null)
					newVals.add(to[i][j]);
			}
			// take the intersection
			for (int i= intersectionVals.size()-1; i>=0; i--) {
				if (!newVals.contains(intersectionVals.get(i))) {
					intersectionVals.remove(i);
				}
			}

		}    			
    	
    	return intersectionVals;
    }

    private double getMean(Double[][] allRates, Integer toval, Integer[][] to) {
        double m = 0;
        int l = 0;
    	for (int i=0; i < to[0].length; i++)
    		for (int j=0; j < to.length;j++)
    			if (to[j][i]==toval) {
    				l++;
    				m+=allRates[j][i];
    			}		
        

        return m / l;
    	
    }
    
    private double getPost(Double[][] allRates, Integer toval, Integer[][] to) {
        int nonzero = 0;
        int l = 0;
    	for (int i=0; i < to[0].length; i++)
    		for (int j=0; j < to.length;j++)
    			if (to[j][i]==toval) {
    				l++;
    				if (allRates[j][i]>0)
    					nonzero++;
    			}		
			   
	   return nonzero/((double) l);	   
   }
   
    public void printMetaData(PrintStream ps, int index,  Set<String> attributeNames) {
    	Tree tree = new Tree(unrankedTrees.get(index).tree.getRoot());
    	
    	
    	// build the species number to clade mapping
    	int[][] oriSpeciesNr = new int[unrankedTrees.get(index).rankedSet.size()][];
    	for (int i = 0; i < oriSpeciesNr.length;i++) {
            List<Object[]> attributeValues = unrankedTrees.get(index).rankedSet.get(i).getAttributeValues();
            int j=0;
            for (String attributeName : attributeNames) {
         		if (attributeName.contentEquals("species")) {
         			oriSpeciesNr[i] = new int[attributeValues.size()];
         			for (int k = 0; k< attributeValues.size(); k++) {
         				Double val = (Double) attributeValues.get(k)[j];
         				oriSpeciesNr[i][k] = (int) Math.round(val);
         			}
         		}
         		j++;
            }
    	}

    	// get the species names for tips and ancestral nodes
    	String[] names = new String[unrankedTrees.get(index).rankedSet.size()];
    	List<Double> species = new ArrayList<>();
    	for (int i = 0; i < names.length;i++) {
    		int nodeNr = unrankedTrees.get(index).rankedSet.get(i).nodeNr;
    		Node n = tree.getNode(nodeNr);
    		List<String> childnames = getChildNames(n);
    		String newName = childnames.get(0);
    		for (int j = 1; j < childnames.size();j++) {
    			newName = newName + ":" +  childnames.get(j);
    		}
    		names[i] = newName;
    		species.add((Double) n.getMetaData("species"));
    	}
    	
    	Integer[][][] allTo = new Integer[names.length][][];
    	
    	// print header for the log file
    	ps.print("sample");
    	for (String attributeName : attributeNames) {
    		for (int i = 0; i < unrankedTrees.get(index).rankedSet.size();i++) {             
    			List<Object[]> attributeValues = unrankedTrees.get(index).rankedSet.get(i).getAttributeValues();	        
	     		if (attributeName.contentEquals("height")) {
	     			ps.print("\theight_" +names[i]);
	     		}
	     		if (attributeName.contentEquals("Ne")) {
	     			ps.print("\tNe_"+names[i]);
	     		}
	     		if (attributeName.contentEquals("rates")) {
	     			// find the to vector
	     			int tovec = 0;
	     			for (String an : attributeNames) {
	     				if (an.contentEquals("to")){
	     					break;
	     				}
	     				tovec++;
	     			}	     			
	     			Double[] to = (Double[]) attributeValues.get(0)[tovec];
	     			if (to!=null) {
	     				allTo[i] = getTo(tovec, attributeValues, oriSpeciesNr);
	     				List<Integer> tovalues = getIntersectionToVals(allTo[i]);
	     				Collections.sort(tovalues);
		     			for (int k = 0; k < tovalues.size();k++) {	     
		     				int toInd = species.indexOf((double) tovalues.get(k));
		     				ps.print("\tbmig_"+ names[i] + "_to_" + names[toInd]);
		     			}	     		
	     			}
	     		}
	        } 
        }
    	
    	// print the Data
    	for (int sample=0;sample<unrankedTrees.get(index).rankedSet.get(0).getAttributeValues().size();sample++ ) {
	    	ps.print("\n"+sample);
	    	int j=0;
	    	
	    	for (String attributeName : attributeNames) {
	    		for (int i = 0; i < unrankedTrees.get(index).rankedSet.size();i++) {             
	    			List<Object[]> attributeValues = unrankedTrees.get(index).rankedSet.get(i).getAttributeValues();	        
		     		if (attributeName.contentEquals("height")) {
		     			ps.print("\t" + attributeValues.get(sample)[j]);
		     		}
		     		if (attributeName.contentEquals("Ne")) {
		     			ps.print("\t" + attributeValues.get(sample)[j]);
		     		}
		     		if (attributeName.contentEquals("rates")) {
		     			// find the to vector
		     			Double[] rates = (Double[]) attributeValues.get(sample)[j];
		     			if (rates!=null) {
		     				List<Integer> tovalues = getIntersectionToVals(allTo[i]);
		     				Collections.sort(tovalues);
			     			for (int k = 0; k < tovalues.size();k++) {
			     				for (int l = 0; l <allTo[i].length; l++) {
			     					if (allTo[i][l][sample]==tovalues.get(k))
			     						ps.print("\t" + rates[l]);
			     				}
			     			}
		     			}//if
		     		}//if
		        }// i
	    		j++;
	        }// attributes
    	}// sample   	
    }

    public List<String> getChildNames(Node n) {
    	List<String> names = new ArrayList<>();;
    	if (n.isLeaf()) {
    		names.add(n.getID());
    	}else {
    		for (Node child : n.getChildren()) {
    			names.addAll(getChildNames(child));
    		}
    	}    	
    	return names;
    }

    public class UnrankedTree {    	
    	
        public UnrankedTree() {
        	rankedSet = new ArrayList<>();
        	count=1;
        }
        
        public void addRBS(UnrankedBitSet rbs) {
        	rankedSet.add(rbs);        	
        }
        
        public int getCount() {
            return count;
        }

        public void setCount(int count) {
            this.count = count;
        }

        public double getCredibility() {
            return credibility;
        }

        public void setCredibility(double credibility) {
            this.credibility = credibility;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final UnrankedTree rt = (UnrankedTree) o;
            
            for (UnrankedBitSet rbs :  rt.rankedSet) {
            	// compare all vs. all
            	boolean equal=false;
            	for (UnrankedBitSet rbs2 : rankedSet) {
            		if (rbs.equals(rbs2))
            			equal=true;
            	}
            	if (!equal)
            		return false;
            }

            return true;

        }
        
        @Override
        public String toString() {
            return "treeSet " +  count;
        }

        List<UnrankedBitSet> rankedSet;
        int count;
        double credibility;
        Tree tree;
    }
 
    public class UnrankedBitSet {
    	
        public UnrankedBitSet(BitSet bits, int nodeNr) {
            this.bits = bits;
            this.nodeNr = nodeNr;
        }
        
        public List<Object[]> getAttributeValues() {
            return attributeValues;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final UnrankedBitSet rbs = (UnrankedBitSet) o;
            
           
            if (!bits.equals(rbs.bits))
            	return false;
            

            return !(bits != null ? !bits.equals(rbs.bits) : rbs.bits != null);

        }

        @Override
        public int hashCode() {
            return (bits != null ? bits.hashCode() : 0);
        }

        @Override
        public String toString() {
            return bits.toString();
        }

        
        int species;
        int nodeNr;
        BitSet bits;
        List<Object[]> attributeValues = null;
    }

}
