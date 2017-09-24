package beast.app.beauti;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.evolution.alignment.Alignment;

import java.util.*;
import starbeast2.SpeciesTree;
import starbeast2.StarBeastTaxonSet;

@Description("Class for creating new partitions for morphological data to be edited by AlignmentListInputEditor")
public class StarBeastMorphModelAlignmentProvider extends BeautiMorphModelAlignmentProvider {
    @Override
    public void processAlignment(Alignment alignment, List<BEASTInterface> filteredAlignments, boolean ascertained, BeautiDoc doc) throws Exception {
        StarBeastTaxonSet ts = (StarBeastTaxonSet) doc.pluginmap.get("taxonsuperset");
        ts.alignmentInput.set(alignment);
        ts.initAndValidate();

        SpeciesTree st = (SpeciesTree) doc.pluginmap.get("Tree.t:Species");
        st.m_taxonset.set(ts);
        st.makeCaterpillar(0, 1, false);

        super.processAlignment(alignment, filteredAlignments, ascertained, doc);
    }
}
