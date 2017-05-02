package starbeast2.utils;

import beast.core.Input;
import beast.core.Runnable;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeTraceAnalysis;
import beast.util.NexusParser;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

/**
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 2/05/17.
 */
public class TreeTopologyDistGenerator extends Runnable {

    public Input<String> logFileNameInput = new Input<>(
            "logFileName",
            "Name of log file from which trees will be read.",
            Input.Validate.REQUIRED);

    public Input<Double> burninFracInput = new Input<>(
            "burninFrac",
            "Burn-in fraction. (Default 0.1.)", 0.1);

    public Input<String> reportFileNameInput = new Input<>(
            "reportFileName",
            "Name of file to which topology distribution report will be written.");

    String logFileName;
    public Tree tree;

    public TreeTopologyDistGenerator() { }

    public void initAndValidate() {
        logFileName = logFileNameInput.get();
    }

    @Override
    public void run() throws Exception {

        NexusParser nexusParser = new NexusParser();
        nexusParser.parseFile(new File(logFileName));

        int totalTrees = nexusParser.trees.size();
        int burnin = (int)Math.round(burninFracInput.get()*totalTrees);

        if (reportFileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(reportFileNameInput.get())) {
                ModifiedTreeTraceAnalysis analysis = new ModifiedTreeTraceAnalysis();

                for (int i=(burnin+1); i<nexusParser.trees.size(); i++)
                    analysis.addTree(nexusParser.trees.get(i));

//                analysis.analyze(1.0);
                analysis.report(ps, totalTrees, burnin);
            }
        }
    }
}
