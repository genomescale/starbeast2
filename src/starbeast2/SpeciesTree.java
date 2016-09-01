package starbeast2;

import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
 */

public class SpeciesTree extends TreeWrapper {
    final private Map<String, Integer> tipNumberMap = new LinkedHashMap<>();
    final private Multimap<Integer, String> numberTipMap = HashMultimap.create();

    public void initAndValidate() {
        // generate map of species tree tip node names to node numbers
        final TreeInterface speciesTree = treeInput.get();
        final Map<String, Integer> speciesNumberMap = new LinkedHashMap<>();

        Node speciesTreeRoot = speciesTree.getRoot();
        for (Node leafNode: speciesTreeRoot.getAllLeafNodes()) {
            final String speciesName = leafNode.getID();
            final int speciesNumber = leafNode.getNr();

            speciesNumberMap.put(speciesName, speciesNumber);
        }

        // generate map of gene tree tip node names to species tree tip node numbers
        final TaxonSet taxonSuperSet = treeInput.get().getTaxonset();
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

    public Multimap<Integer, String> getNumberTipMap() {
        return numberTipMap;
    }

    public Map<String, Integer> getTipNumberMap() {
        return tipNumberMap;
    }
}

