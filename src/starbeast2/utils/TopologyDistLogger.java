package starbeast2.utils;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeTraceAnalysis;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 28/04/17.
 */
public class TopologyDistLogger extends BEASTObject implements Loggable {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree whose topology to consider",
            Input.Validate.REQUIRED);

    public Input<Double> burninInput = new Input<>(
            "burnInFrac",
           "Burn-in fraction omitted in analysis.",
            Input.Validate.REQUIRED);


    Tree tree;
    List<Tree> treeList;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        treeList = new ArrayList<>();
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
        TreeTraceAnalysis analysis = new TreeTraceAnalysis(treeList, burninInput.get());
        analysis.analyze(1.0);
        analysis.report(out);
    }
}
