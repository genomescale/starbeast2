package starbeast2.utils;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Logger;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeTraceAnalysis;
import beast.util.FrequencySet;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 28/04/17.
 */
public class TreeTopologyDistLogger extends Logger {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree whose topology to consider",
            Input.Validate.REQUIRED);

    public Input<Integer> burninSamplesInput = new Input<>(
            "burninSamples",
           "Number of burn-in samples omitted from analysis.",
            Input.Validate.REQUIRED);

    public Input<Boolean> isSpeciesTreeInput = new Input<>(
            "isSpeciesTree",
            "Whether this is a species tree or not. (Default false.)", false);


    Tree tree;
    ModifiedTreeTraceAnalysis analysis;
    int burninSamples, nTreesTotal, nBurninTrees;

    public TreeTopologyDistLogger() {
        loggersInput.setRule(Input.Validate.OPTIONAL);
        loggersInput.setValue(new DummyLoggable(), this);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        tree = treeInput.get();

        burninSamples = burninSamplesInput.get();

        analysis = new ModifiedTreeTraceAnalysis(isSpeciesTreeInput.get());

    }

    @Override
    public void init() throws IOException {
        super.init();

        nTreesTotal = 0;
        nBurninTrees = 0;
    }


    @Override
    public void log(int sample) {
        nTreesTotal += 1;

        if (sample> burninSamples)
            analysis.addTree(tree);
        else
            nBurninTrees += 1;
    }

    @Override
    public void close() {
        analysis.report(getM_out(), nTreesTotal, nBurninTrees);
    }

    /**
     * TreeTraceAnalysis class with some modifications to support sampled
     * ancestors and incremental sample addition.
     */
    class ModifiedTreeTraceAnalysis extends TreeTraceAnalysis {

        boolean isSpeciesTree;

        public ModifiedTreeTraceAnalysis(boolean isSpeciesTree) {
            super(new ArrayList<>(), 0.1);

            this.isSpeciesTree = isSpeciesTree;
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

    class DummyLoggable extends BEASTObject implements Loggable {

        public DummyLoggable() { }

        @Override
        public void init(PrintStream out) { }

        @Override
        public void log(int sample, PrintStream out) { }

        @Override
        public void close(PrintStream out) { }

        @Override
        public void initAndValidate() { }
    }
}
