package starbeast2;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;

public abstract class MigrationModel extends CalculationNode {
	
    public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree", "The species tree this model applies to.", Validate.REQUIRED);

	public abstract double getMigration(int sourceNode, int sinkNode);

	public abstract double getEM();

}
