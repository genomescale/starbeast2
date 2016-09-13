package starbeast2;

import java.util.LinkedHashMap;
import java.util.Map;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.util.TreeParser;

/**
* @author Huw Ogilvie
 */

// mainly for debugging... users who wish to use starting trees can more easily
// set the "newick" attribute of the StarBeastInitializer
public class SpeciesTreeParser extends TreeParser implements SpeciesTreeInterface {
    Map<String, Integer> tipNumberMap;
    Multimap<Integer, String> numberTipMap;

    public void initAndValidate() {
        super.initAndValidate();

        tipNumberMap = new LinkedHashMap<>();
        numberTipMap = HashMultimap.create();
        makeMaps(tipNumberMap, numberTipMap);
    }

    public Map<String, Integer> getTipNumberMap() {
        return tipNumberMap;
    }

    public Multimap<Integer, String> getNumberTipMap() {
        return numberTipMap;
    }
}
