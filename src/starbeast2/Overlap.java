package starbeast2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

public class Overlap extends MigrationModel {
	
    public Input<RealParameter> effectiveMigrantsInput  = new Input<>("effectiveMigrants","absolute migration rates",Input.Validate.REQUIRED);
    public Input<Double> minimalBranchLengthInput  = new Input<>("minimalBranchLength","absolute migration rates", 0.0);

    public Input<String> excludeInput = new Input<>("exclude", "nodes with no migration");
    
    private boolean hasExcluded = false;
    private ArrayList<Integer> exclNode;
    
	@Override
	public void initAndValidate() {
		if (excludeInput.get()!=null){
			exclNode = new ArrayList<Integer>();
			String[] splitStr = excludeInput.get().split("\\s+");
			for (int i = 0; i < splitStr.length; i++){
//				System.out.println(splitStr[i]);
				for (int j=0; j < speciesTreeInput.get().getLeafNodeCount(); j++){
					if(splitStr[i].equals(speciesTreeInput.get().getNode(j).getID() )){
						exclNode.add(j);
						break;
					}
				}
			}
				
			
			hasExcluded = true;
		}
	}

	@Override
	public double getMigration(int sourceNode, int sinkNode) {
		Node n1 = speciesTreeInput.get().getNode(sourceNode);
		Node n2 = speciesTreeInput.get().getNode(sinkNode);
		
		// check if some nodes are exluded
		if(hasExcluded){
			
			
			List<Node> al1 = n1.getAllLeafNodes();
			List<Node> al2 = n2.getAllLeafNodes();
						
			if (al1.size()>0){
				for (int j = 0; j < al1.size(); j++){
					if (exclNode.indexOf(al1.get(j).getNr())!=-1){
						return 0.0;							
					}
				}
			}else{
				if (exclNode.indexOf(n1.getNr())!=-1){
					return 0.0;
				}
			}
			
			if (al2.size()>0)
				for (int j = 0; j < al2.size(); j++)
					if (exclNode.indexOf(al2.get(j).getNr())!=-1)
						return 0.0;
			if (exclNode.indexOf(n2.getNr())!=-1)
				return 0.0;
				
		}
		
//		System.out.println(n1.getAllLeafNodes());
			
				
		double lower = Math.max(n1.getHeight(), n2.getHeight());
		double upper = Math.min(n1.getParent().getHeight(), n2.getParent().getHeight());

		double b = Math.max(upper-lower, minimalBranchLengthInput.get());
		
		return effectiveMigrantsInput.get().getValue()/(b); 		
	}

	@Override
	public double getEM() {
		return effectiveMigrantsInput.get().getValue();
	}	

		
	
	
}
