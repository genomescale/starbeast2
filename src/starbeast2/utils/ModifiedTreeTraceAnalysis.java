package starbeast2.utils;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeTraceAnalysis;
import beast.util.FrequencySet;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * TreeTraceAnalysis class with some modifications to support sampled
 * ancestors and incremental sample addition.
 */
public class ModifiedTreeTraceAnalysis extends TreeTraceAnalysis {

    public ModifiedTreeTraceAnalysis() {
        super(new ArrayList<>(), 0.1);

        topologiesFrequencySet = new FrequencySet<>();
        topologiesFrequencySet.setCredSetProbability(1.0);
    }

    /**
     * Add topology corresponding to given tree to topolgoy frequency set.
     *
     * @param tree tree whose topology to add
     */
    public void addTree(Tree tree) {
        String topology = uniqueNewick(tree.getRoot());
        topologiesFrequencySet.add(topology, 1);
    }

    public void report(PrintStream oStream, int totalTrees, int burnin) {
        credibleSet = topologiesFrequencySet.getCredibleSet();
        this.totalTrees = totalTrees;
        this.burnin = burnin;

        super.report(oStream);
    }

    /**
     * Get tree topology in Newick that is sorted by taxon labels.
     *
     * @param node root of tree
     * @return newick string
     */
    String getSortedNewickWithSAs(Node node) {
        if (node.isLeaf()) {
            return String.valueOf(node.getID());
        } else {
            StringBuilder builder = new StringBuilder("(");

            List<String> subTrees = new ArrayList<>();
            for (Node child : node.getChildren()) {
                if (!child.isDirectAncestor())
                    subTrees.add(getSortedNewickWithSAs(child));
            }

            Collections.sort(subTrees);

            for (int i = 0; i < subTrees.size(); i++) {
                builder.append(subTrees.get(i));
                if (i < subTrees.size() - 1) {
                    builder.append(",");
                }
            }
            builder.append(")");

            if (node.isFake())
                builder.append(node.getDirectAncestorChild().getID());

            return builder.toString();
        }
    }

    @Override
    public String uniqueNewick(Node node) {
        return getSortedNewickWithSAs(node);
    }
}
