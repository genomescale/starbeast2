package starbeast2;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;

import java.util.ArrayList;
import java.util.List;

public class AllEqual extends MigrationModel {	

    public Input<RealParameter> effectiveMigrantsInput  = new Input<>("effectiveMigrants","absolute migration rates",Input.Validate.OPTIONAL);

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
		if (effectiveMigrantsInput.get()==null)
			return 1.0;
		else
			return effectiveMigrantsInput.get().getValue(); 		
	}

	@Override
	public double getEM() {
		if (effectiveMigrantsInput.get()==null)
			return 1.0;
		else
			return effectiveMigrantsInput.get().getValue();
	}	
		
}
