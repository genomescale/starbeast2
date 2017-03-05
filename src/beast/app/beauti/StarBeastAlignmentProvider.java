package beast.app.beauti;

import java.io.File;
import java.util.List;

import beast.core.BEASTInterface;

public class StarBeastAlignmentProvider extends BeautiAlignmentProvider {

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		doc.autoSetClockRate = false;
		doc.beauti.autoSetClockRate.setSelected(false);
		return super.getAlignments(doc, files);
	}
	
}
