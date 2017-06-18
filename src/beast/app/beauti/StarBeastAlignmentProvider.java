package beast.app.beauti;

import java.io.File;
import java.util.List;

import beast.core.BEASTInterface;
import beast.evolution.operators.DeltaExchangeOperator;

public class StarBeastAlignmentProvider extends BeautiAlignmentProvider {

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		doc.autoSetClockRate = false;
		doc.beauti.autoSetClockRate.setSelected(false);

		System.out.println(String.format("N_ALIGNMENTS = %d", doc.alignments.size()));
		// initialize delta exchange operator in order to increase weight to something more sensible
		DeltaExchangeOperator operator = (DeltaExchangeOperator) doc.pluginmap.get("FixMeanMutationRatesOperator");
    	if (operator == null) {
    		operator = new DeltaExchangeOperator();
    		try {
    			operator.setID("FixMeanMutationRatesOperator");
				operator.initByName("weight", (double) files.length, "delta", 0.75);
			} catch (Throwable e1) {
				// ignore initAndValidate exception
			}
    		doc.addPlugin(operator);
    	} else {
    		final double updatedWeight = doc.alignments.size() + files.length;
    		operator.setInputValue("weight", updatedWeight);
    		System.out.println("HERE WE ARE");
    	}

		return super.getAlignments(doc, files);
	}
	
}
