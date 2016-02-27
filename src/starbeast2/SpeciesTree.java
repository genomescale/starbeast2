package starbeast2;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
 */

public class SpeciesTree extends TreeWrapper {
    public Input<TaxonSet> taxonSuperSetInput = new Input<>("taxonSuperSet", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);

    final private Map<String, Integer> tipNumberMap = new HashMap<>();
    final private Multimap<Integer, String> numberTipMap = HashMultimap.create();

    public void initAndValidate() {
        // generate map of species tree tip node names to node numbers
        final TreeInterface speciesTree = treeInput.get();
        final HashMap<String, Integer> speciesNumberMap = new HashMap<>();

        Node speciesTreeRoot = speciesTree.getRoot();
        for (Node leafNode: speciesTreeRoot.getAllLeafNodes()) {
            final String speciesName = leafNode.getID();
            final int speciesNumber = leafNode.getNr();

            speciesNumberMap.put(speciesName, speciesNumber);
        }

        // generate map of gene tree tip node names to species tree tip node numbers
        final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
        final Set<Taxon> speciesSet = new HashSet<>(taxonSuperSet.taxonsetInput.get());

        for (Taxon species: speciesSet) {
            final String speciesName = species.getID();
            final int speciesNumber = speciesNumberMap.get(speciesName);
            final TaxonSet speciesTaxonSet = (TaxonSet) species;
            final Set<Taxon> tipSet = new HashSet<>(speciesTaxonSet.taxonsetInput.get());

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

