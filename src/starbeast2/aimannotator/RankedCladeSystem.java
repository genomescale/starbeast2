package starbeast2.aimannotator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.app.treeannotator.CladeSystem;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.statistic.DiscreteStatistics;

/**
 * extracted from TreeAnnotator
 */
//TODO merge with CladeSet?
public class RankedCladeSystem {
	private boolean processSA = true;

    protected List<RankedTree> rankedTrees = new ArrayList<>();
    protected List<RankedTree> orderedRankedTrees = new ArrayList<>();
    protected RankedTree newRankedTree;

    public RankedCladeSystem() { }

    public RankedCladeSystem(Tree targetTree, Set<String> attributeNames) {
        add(targetTree, true, attributeNames);
    }
    
    /**
     * adds all the clades in the tree
     */
    public void add(Tree tree, boolean includeTips, Set<String> attributeNames) {
    	getCurrentRankedTree(tree, includeTips, attributeNames);
    	// get all node heights for ranking the nodes
        // compare the new ranked tree to the current ranked trees
        int rtIndex = rankedTrees.indexOf(newRankedTree);
        if (rtIndex==-1) {
        	rankedTrees.add(newRankedTree);
        	rankedTrees.get(rankedTrees.size()-1).setCount(1);
        	rankedTrees.get(rankedTrees.size()-1).tree = tree.copy();
        }else {
        	rankedTrees.get(rtIndex).setCount(rankedTrees.get(rtIndex).getCount()+1);
        	// add the attribute values
        	addAttributes(rtIndex);
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
    	
        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
    	newRankedTree = new RankedTree();
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
    	RankedBitSet rbs = new RankedBitSet(bits, rank, node.getNr());
    	collectAttributesForClade(rbs, node, attributeNames);
    	newRankedTree.addRBS(rbs);
	}   

    public void calculateCladeCredibilities(int totalTreesUsed) {
        for (RankedTree rankedTree : rankedTrees) {

//            if (rankedTree.getCount() > totalTreesUsed) {
//
//                throw new AssertionError("clade.getCount=(" + rankedTree.getCount() +
//                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
//            }

            rankedTree.setCredibility(((double) rankedTree.getCount()) / (double) totalTreesUsed);
        }
    }
    
    public int[] orderRankedTrees() {
    	int[] order = new int[rankedTrees.size()];
    	// get the order of trees by rank
    	for (int i = 0; i < rankedTrees.size(); i++) {
    		int highestCount = -1;
    		int highestIndex = -1;
        	for (int j = 0; j < rankedTrees.size(); j++) {
	    		if (rankedTrees.get(j).count>highestCount) {
	    			highestCount = rankedTrees.get(j).count;
	    			highestIndex = j;
	    		}
	    	}
        	order[i] = highestIndex;
        	rankedTrees.get(highestIndex).count=-1;
    	}
    	return order;
    }
    
    private void collectAttributesForClade(RankedBitSet rbs, Node node, Set<String> attributeNames) {
    	
    	
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
    	for (RankedBitSet newRBS : newRankedTree.rankedSet) {
    		int oldRBSindex = rankedTrees.get(rtIndex).rankedSet.indexOf(newRBS);
    		RankedBitSet rbs = rankedTrees.get(rtIndex).rankedSet.get(oldRBSindex);
    		
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
    	Tree tree = new Tree(rankedTrees.get(index).tree.getRoot());
    	
    	// build the species number to clade mapping
    	int[][] oriSpeciesNr = new int[rankedTrees.get(index).rankedSet.size()][];
    	for (int i = 0; i < oriSpeciesNr.length;i++) {
            List<Object[]> attributeValues = rankedTrees.get(index).rankedSet.get(i).getAttributeValues();
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
    	
    	
    	for (int i = 0; i < rankedTrees.get(index).rankedSet.size();i++) {
    		int nodeNr = rankedTrees.get(index).rankedSet.get(i).nodeNr;
    		Node n = tree.getNode(nodeNr);
    		n.metaDataString = "species=" + oriSpeciesNr[i][0];
    		
            int j=0;
            
            List<Object[]> attributeValues = rankedTrees.get(index).rankedSet.get(i).getAttributeValues();
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
         				double[][] allRates = getRates(j, attributeValues);
	         			// compute posterior support for rates
	         			n.metaDataString = n.metaDataString + ",rates_mean={" + DiscreteStatistics.mean(allRates[0]);
	         			for (int k = 1; k < allRates.length;k++)
	         				n.metaDataString = n.metaDataString + "," +  DiscreteStatistics.mean(allRates[k]);
	         			n.metaDataString = n.metaDataString + "}";
	         			
	         			// compute posterior support for rates
	         			n.metaDataString = n.metaDataString + ",rates_posterior={" + getPost(allRates[0]);
	         			for (int k = 1; k < allRates.length;k++)
	         				n.metaDataString = n.metaDataString + "," + getPost(allRates[k]);
	         			n.metaDataString = n.metaDataString + "}";
         			}
         		}  
         		if (attributeName.contentEquals("to")) {
         			if (!n.isRoot()) {
	         			Double[] to = (Double[]) attributeValues.get(0)[j];
	         			// compute posterior support for rates
	         			n.metaDataString = n.metaDataString + ",rates_to={" + Math.round(to[0]);
	         			for (int k = 1; k < to.length;k++)
	         				n.metaDataString = n.metaDataString + "," + Math.round(to[k]);
	         			n.metaDataString = n.metaDataString + "}";
         			}
         		}
            	j++;	
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
    
    private double[][] getRates(int j, List<Object[]> attributeValues) {
    	Double[] rates = (Double[]) attributeValues.get(0)[j];

    	double[][] allRates = new double[rates.length][attributeValues.size()];
    	
		for (int i = 0; i < attributeValues.size(); i++) {
	    	rates = (Double[]) attributeValues.get(i)[j];
	    	for (int k = 0; k < rates.length; k++) {
	    		allRates[k][i] = (double) rates[k];
	    	}
		}
		return allRates;
    }
  
    private double getPost(double[] rates) {
	   double c = 0.0;
	   for (int i = 0; i < rates.length; i++)
		   if (rates[i]>0)
			   c++;
			   
	   return c/rates.length;	   
   }
   
    public void printMetaData(PrintStream ps, int index,  Set<String> attributeNames) {
    	Tree tree = new Tree(rankedTrees.get(index).tree.getRoot());

    	// get the species names for tips and ancestral nodes
    	String[] names = new String[rankedTrees.get(index).rankedSet.size()];
    	List<Double> species = new ArrayList<>();
    	for (int i = 0; i < names.length;i++) {
    		int nodeNr = rankedTrees.get(index).rankedSet.get(i).nodeNr;
    		Node n = tree.getNode(nodeNr);
    		List<String> childnames = getChildNames(n);
    		String newName = childnames.get(0);
    		for (int j = 1; j < childnames.size();j++) {
    			newName = newName + ":" +  childnames.get(j);
    		}
    		names[i] = newName;
    		species.add((Double) n.getMetaData("species"));
    	}
    	// print header for the log file
    	ps.print("sample");
    	for (String attributeName : attributeNames) {
    		for (int i = 0; i < rankedTrees.get(index).rankedSet.size();i++) {             
    			List<Object[]> attributeValues = rankedTrees.get(index).rankedSet.get(i).getAttributeValues();	        
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
		     			for (int k = 0; k < to.length;k++) {	     
		     				int toInd = species.indexOf(to[k]);
		     				ps.print("\tbmig_"+ names[i] + "_to_" + names[toInd]);
		     			}	     		
	     			}
	     		}
	        } 
        }
    	
    	// print the Data
    	for (int sample=0;sample<rankedTrees.get(index).rankedSet.get(0).getAttributeValues().size();sample++ ) {
	    	ps.print("\n"+sample);
	    	int j=0;
	    	for (String attributeName : attributeNames) {
	    		for (int i = 0; i < rankedTrees.get(index).rankedSet.size();i++) {             
	    			List<Object[]> attributeValues = rankedTrees.get(index).rankedSet.get(i).getAttributeValues();	        
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
			     			for (int k = 0; k < rates.length;k++) {	     
			     				ps.print("\t" + rates[k]);
			     			}//k	     		
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

    public class RankedTree {    	
    	
        public RankedTree() {
        	rankedSet = new ArrayList<>();
        	count=1;
        }
        
        public void addRBS(RankedBitSet rbs) {
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

            final RankedTree rt = (RankedTree) o;
            
            for (RankedBitSet rbs :  rt.rankedSet) {
            	if (rankedSet.indexOf(rbs)==-1)
            		return false;
            }

            return true;

        }
        
        @Override
        public String toString() {
            return "treeSet " +  count;
        }

        List<RankedBitSet> rankedSet;
        int count;
        double credibility;
        Tree tree;
    }
    

    public class RankedBitSet {
    	
        public RankedBitSet(BitSet bits, int rank, int nodeNr) {
            this.bits = bits;
            this.rank = rank;
            this.nodeNr = nodeNr;
        }
        
        public List<Object[]> getAttributeValues() {
            return attributeValues;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final RankedBitSet rbs = (RankedBitSet) o;
            
            // compare ranks
            if (rbs.rank!=rank)
            	return false;

            return !(bits != null ? !bits.equals(rbs.bits) : rbs.bits != null);

        }

        @Override
        public int hashCode() {
            return (bits != null ? bits.hashCode() : 0);
        }

        @Override
        public String toString() {
            return bits.toString() + "r" + rank;
        }

        
        int rank;
        int species;
        int nodeNr;
        BitSet bits;
        List<Object[]> attributeValues = null;
    }

}
