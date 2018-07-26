package starbeast2;

import java.util.ArrayList;
import java.util.Collections;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

public class MinimalBranchLength extends MigrationModel {
	
    public Input<RealParameter> effectiveMigrantsInput  = new Input<>("effectiveMigrants","absolute migration rates",Input.Validate.REQUIRED);
    public Input<Double> minimalBranchLengthInput  = new Input<>("minimalBranchLength","minimum bound to stop arbitrarily large prior on migration rate", 0.0001);

    
	@Override
	public void initAndValidate() {
	}

	@Override
	public double getMigration(int sourceNode, int sinkNode) {
		Node n1 = speciesTreeInput.get().getNode(sourceNode);
		Node n2 = speciesTreeInput.get().getNode(sinkNode);
		
		ArrayList<Integer> parents1 = new ArrayList<>();
		Node p1 = n1.getParent();
		while (p1 != null){
			parents1.add(p1.getNr());
			p1 = p1.getParent();
		}
		ArrayList<Integer> parents2 = new ArrayList<>();
		Node p2 = n2.getParent();
		while (p2 != null){
			parents2.add(p2.getNr());
			p2 = p2.getParent();
		}
		parents2.retainAll(parents1);

		// get the node in parents2 with the lowest height
		ArrayList<Double> heights = new ArrayList<>();
		for (int i = 0; i < parents2.size(); i++)
			heights.add(speciesTreeInput.get().getNode(parents2.get(i)).getHeight());
		
		Collections.sort(heights);
		if (heights.size()==0){
			System.out.println(n1 + "\t" + n2);
			System.out.println(n1.getTree());
		}
		double b = Math.max(Math.min((heights.get(0)-n1.getHeight()), (heights.get(0)-n2.getHeight())), minimalBranchLengthInput.get());
		
		return effectiveMigrantsInput.get().getValue()/(b) ; 		
	}	

		
	
	
}
