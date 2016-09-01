package starbeast2;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.CompoundDistribution;

public class MultispeciesCoalescent2 extends CompoundDistribution {
    public Input<MultispeciesPopulationModel> populationModelInput = new Input<>("populationModel", "The species tree population model.", Validate.REQUIRED);

    MultispeciesPopulationModel popModel;
    
    @Override
    public void initAndValidate() {
    	popModel = populationModelInput.get();
    	super.initAndValidate();
    }
    
    @Override
    public double calculateLogP() {
    	logP = super.calculateLogP();
    	logP += popModel.commonContributionLogP();
    	return logP;
    }
    
    @Override
    public void store() {
    	popModel.doStore();
    	super.store();
    }
    
    @Override
    public void restore() {
    	popModel.doRetore();
    	super.restore();
    }
}
