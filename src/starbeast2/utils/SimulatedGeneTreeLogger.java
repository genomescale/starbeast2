package starbeast2.utils;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

/**
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 28/04/17.
 */
public class SimulatedGeneTreeLogger extends BEASTObject implements Loggable {

    public Input<SimulatedGeneTree> geneTreeInput = new Input<>("simulatedGeneTree",
            "Simulated gene tree.", Input.Validate.REQUIRED);

    SimulatedGeneTree geneTree;

    @Override
    public void initAndValidate() {
        geneTree = geneTreeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        geneTree.init(out);
    }

    @Override
    public void log(long sample, PrintStream out) {
        geneTree.initAndValidate();
        geneTree.log(sample, out);
    }

    @Override
    public void close(PrintStream out) {
        geneTree.init(out);
    }
}
