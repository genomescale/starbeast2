package starbeast2.utils;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Logger;

import java.io.IOException;
import java.io.PrintStream;

/**
 * Computes the frequency distribution over topologies.
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

    Tree tree;
    TreeTraceAnalysisWithError analysis;
    int burninSamples, nTreesTotal, nBurninTrees;

    int sampleNr;

    public TreeTopologyDistLogger() {
        loggersInput.setRule(Input.Validate.OPTIONAL);
        loggersInput.setValue(new DummyLoggable(), this);
        modeInput.setRule(Input.Validate.OPTIONAL);
        modeInput.setValue("tree", this);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        tree = treeInput.get();

        burninSamples = burninSamplesInput.get();

        analysis = new TreeTraceAnalysisWithError();

    }

    @Override
    public void init() throws IOException {
        super.init();

        nTreesTotal = 0;
        nBurninTrees = 0;

        sampleNr = 0;
    }


    @Override
    public void log(long sample) {
        sampleNr += 1;

        if (sampleNr % everyInput.get() != 0)
            return;

        nTreesTotal += 1;

        if (sample> burninSamples)
            analysis.addTree(tree);
        else
            nBurninTrees += 1;
    }

    @Override
    public void close() {
        analysis.computeCredibleSet(1.0);
        analysis.report(getM_out());
    }

    class DummyLoggable extends BEASTObject implements Loggable {

        public DummyLoggable() { }

        @Override
        public void init(PrintStream out) { }

        @Override
        public void log(long sample, PrintStream out) { }

        @Override
        public void close(PrintStream out) { }

        @Override
        public void initAndValidate() { }
    }
}
