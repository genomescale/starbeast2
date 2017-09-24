package starbeast2;

import java.util.LinkedHashMap;
import java.util.Map;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.evolution.tree.Tree;

/**
* @author Huw Ogilvie
 */

public class SpeciesTree extends Tree implements SpeciesTreeInterface {
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

	public void adjustTreeNodeHeights() {
		adjustTreeNodeHeights(root);
	}
}
