package starbeast2;

import beast.base.evolution.tree.Tree;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import java.util.LinkedHashMap;
import java.util.Map;

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

        makeMaps();
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
