package starbeast2;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class AllEqual extends MigrationModel {	

    public Input<RealParameter> effectiveMigrantsInput  = new Input<>("effectiveMigrants","absolute migration rates",Input.Validate.REQUIRED);

	@Override
	public void initAndValidate() {
	}
    
	@Override
	public double getMigration(int sourceNode, int sinkNode) {		
		return effectiveMigrantsInput.get().getValue(); 		
	}	
		
}
