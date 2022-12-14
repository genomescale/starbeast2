package starbeast2.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.inference.operator.DeltaExchangeOperator;
import beastfx.app.inputeditor.BeautiAlignmentProvider;
import beastfx.app.inputeditor.BeautiDoc;

import java.io.File;
import java.util.List;

public class StarBeastAlignmentProvider extends BeautiAlignmentProvider {

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		final List<BEASTInterface> newAlignments = super.getAlignments(doc, files);
		final int alignmentCount = newAlignments.size();

		doc.autoSetClockRate = false;
		doc.beauti.autoSetClockRate(false);

		System.out.println(String.format("N_ALIGNMENTS = %d", doc.alignments.size()));
		// initialize delta exchange operator in order to increase weight to something more sensible
		DeltaExchangeOperator operator = (DeltaExchangeOperator) doc.pluginmap.get("FixMeanMutationRatesOperator");
		if (operator == null) {
			operator = new DeltaExchangeOperator();
			try {
				operator.setID("FixMeanMutationRatesOperator");
				operator.initByName("weight", (double) alignmentCount, "delta", 0.75);
			} catch (Throwable e1) {
				// ignore initAndValidate exception
			}
			doc.addPlugin(operator);
		} else {
			final double updatedWeight = doc.alignments.size() + alignmentCount;
			operator.setInputValue("weight", updatedWeight);
		}

		return newAlignments;
	}
}
