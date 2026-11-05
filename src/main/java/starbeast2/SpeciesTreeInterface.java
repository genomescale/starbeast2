package starbeast2;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import com.google.common.collect.Multimap;

import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

public interface SpeciesTreeInterface extends TreeInterface {
    Map<String, Integer> getTipNumberMap();
    Multimap<Integer, String> getNumberTipMap();

    void init(PrintStream out);
    void close(PrintStream out);

    StateNode getCurrent();

    default void makeMaps() {
        // generate map of species tree tip node names to node numbers
        final Map<String, Integer> speciesNumberMap = new LinkedHashMap<>();

        Node speciesTreeRoot = getRoot();
        for (Node leafNode: speciesTreeRoot.getAllLeafNodes()) {
            final String speciesName = leafNode.getID();
            final int speciesNumber = leafNode.getNr();
    
            speciesNumberMap.put(speciesName, speciesNumber);
        }

        // generate map of gene tree tip node names to species tree tip node numbers
        final Map<String, Integer> tipNumberMap = getTipNumberMap();
        final Multimap<Integer, String> numberTipMap = getNumberTipMap();
        final TaxonSet taxonSuperSet = getTaxonset();
        final Set<Taxon> speciesSet = new LinkedHashSet<>(taxonSuperSet.taxonsetInput.get());

        for (Taxon species: speciesSet) {
            final String speciesName = species.getID();
            int speciesNumber = 0;
            if (speciesNumberMap.containsKey(speciesName)) { // skipped for BEAUTi
                speciesNumber = speciesNumberMap.get(speciesName);
            }
            final TaxonSet speciesTaxonSet = (TaxonSet) species;
            final Set<Taxon> tipSet = new LinkedHashSet<>(speciesTaxonSet.taxonsetInput.get());

            for (Taxon tip: tipSet) {
                final String tipName = tip.getID();
                tipNumberMap.put(tipName, speciesNumber);
                numberTipMap.put(speciesNumber, tipName);
            }
        }
    }
}
