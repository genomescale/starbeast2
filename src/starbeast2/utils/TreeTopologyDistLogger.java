package starbeast2.utils;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeTraceAnalysis;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 28/04/17.
 */
public class TreeTopologyDistLogger extends BEASTObject implements Loggable {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree whose topology to consider",
            Input.Validate.REQUIRED);

    public Input<Double> burninInput = new Input<>(
            "burnInFrac",
           "Burn-in fraction omitted in analysis.",
            Input.Validate.REQUIRED);

    public Input<Boolean> isSpeciesTreeInput = new Input<>(
            "isSpeciesTree",
            "Whether this is a species tree or not. (Default false.)", false);


    Tree tree;
    List<Tree> treeList;
    boolean isSpeciesTree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        treeList = new ArrayList<>();
        isSpeciesTree = isSpeciesTreeInput.get();
    }

    void removeFakeNodes(Node node) {
        if (node.isFake()) {
            Node daChild = node.getDirectAncestorChild();
            node.setNr(daChild.getNr());
            node.setID(daChild.getID());
            node.removeChild(daChild);
        }

        for (Node child : node.getChildren())
            removeFakeNodes(child);
    }

    @Override
    public void init(PrintStream out) {
    }

    @Override
    public void log(int sample, PrintStream out) {
        treeList.add(tree.copy());
    }

    @Override
    public void close(PrintStream out) {
        TreeTraceAnalysis analysis;

        if (isSpeciesTree)
            analysis = new ModifiedTreeTraceAnalysis(treeList, burninInput.get());
        else
            analysis = new TreeTraceAnalysis(treeList, burninInput.get());

        analysis.analyze(1.0);
        analysis.report(out);
    }

    class ModifiedTreeTraceAnalysis extends TreeTraceAnalysis {

        public ModifiedTreeTraceAnalysis(List<Tree> posteriorTreeList, double burninFraction) {
            super(posteriorTreeList, burninFraction);
        }

        /**
         * get tree topology in Newick that is sorted by taxon labels.
         * @param node
         * @return
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
}
