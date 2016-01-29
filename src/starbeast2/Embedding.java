package starbeast2;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import beast.core.Input;
import beast.core.Loggable;
import beast.core.StateNode;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 *
 */

public class Embedding extends StateNode implements Loggable {

    /**
     * Map the reticulation network node number to the gene tree node map.
     * The gene tree node map maps the gene node number to a boolean indicating
     * left network node parent (true) or right (false).
     */
    private Map<Integer, Map> embeddedLineages = new HashMap<>();

    @Override
    public void initAndValidate() throws Exception {
    }

    protected void addReticulationNode(final Integer reticulationNodeNumber) {
        final Map<Integer, Boolean> orientationMap = new HashMap<>();
        embeddedLineages.put(reticulationNodeNumber, orientationMap);
    }

    protected void removeReticulationNode(final Integer reticulationNodeNumber) {
        embeddedLineages.remove(reticulationNodeNumber);
    }

    /**
     * @param reticulationNodeNumber reticulation node number
     * @return the set of gene node numbers for all lineages which traverse this reticulation node
     */
    protected Set<Integer> getGeneLineages(final Integer reticulationNodeNumber) {
        return embeddedLineages.get(reticulationNodeNumber).keySet();
    }

    protected Boolean traversesLeft(final Integer reticulationNodeNumber, final Integer geneNodeNumber) {
        Map<Integer, Boolean> geneLineageMap = embeddedLineages.get(reticulationNodeNumber);
        return geneLineageMap.get(geneNodeNumber);
    }

    protected void flipOrientation(final Integer reticulationNodeNumber, final Integer geneNodeNumber) {
        Map<Integer, Boolean> geneLineageMap = embeddedLineages.get(reticulationNodeNumber);
        Boolean flipped = !geneLineageMap.get(geneNodeNumber);
        geneLineageMap.put(geneNodeNumber, flipped);
    }
}
