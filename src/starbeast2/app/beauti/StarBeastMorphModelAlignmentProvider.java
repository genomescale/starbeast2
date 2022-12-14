package starbeast2.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.evolution.alignment.Alignment;
import beastfx.app.inputeditor.BeautiDoc;
import morphmodels.app.beauti.BeautiMorphModelAlignmentProvider;
import starbeast2.SpeciesTree;
import starbeast2.StarBeastTaxonSet;

import java.util.List;

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
