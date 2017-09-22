package starbeast2;


import beast.core.Description;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;

import java.util.*;


@Description("A TaxonSet is an ordered set of taxa. The order on the taxa is provided at the time of construction" +
        " either from a list of taxon objects or an alignment.")
public class StarBeastTaxonSet extends TaxonSet {

    @Override
    public void initAndValidate() {
        taxaNames = new ArrayList<>();

        if (taxonsetInput.get() == null && alignmentInput.get() == null) {
            throw new IllegalArgumentException("Need a taxonset and/or an alignment as input");
        }

        if (taxonsetInput.get() != null) {
            for (Taxon t: taxonsetInput.get()) {
                final String taxonName = t.getID();
                taxaNames.add(taxonName);
            }
        }

        // Add taxon names in morphology but not in molecular data
        if (alignmentInput.get() != null) {
            for (String taxonName: alignmentInput.get().getTaxaNames()) {
                if (!taxaNames.contains(taxonName)) taxaNames.add(taxonName);
            }
        }
    }
}
