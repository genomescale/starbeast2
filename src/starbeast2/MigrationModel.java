package starbeast2;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public abstract class MigrationModel extends CalculationNode {
	
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.", Validate.REQUIRED);

	public abstract double getMigration(int sourceNode, int sinkNode);
	
}
