package starbeast2.utils;

import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.TreeTraceAnalysis;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 27/04/17.
 */
public class GeneTreeSimulator extends Runnable {

    public Input<Tree> speciesTreeInput = new Input<>(
            "speciesTree",
            "Species tree on which to simulate gene trees",
            Input.Validate.REQUIRED);

    public Input<TraitSet> sampleCountsInput = new Input<>(
            "sampleCounts",
            "TraitSet defining number of  samples per node in species tree.",
            Input.Validate.REQUIRED);

    public Input<Integer> nSimsInput = new Input<>(
            "nSims",
            "Number of gene trees to simulate from given sample distribution.",
            Input.Validate.REQUIRED);

    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Name of file to which gene trees will be written.",
            Input.Validate.REQUIRED);

    public Input<String> reportFileNameInput = new Input<>(
            "reportFileName",
            "Name of file to which topology distribution report will be written.");

    public Input<Double> credibilityThresholdInput = new Input<>(
            "credibilityThreshold",
            "Maximum probability of topologies included in credible set written to report file.",
            0.95);

    public Tree speciesTree;
    public TraitSet sampleCounts;

    public GeneTreeSimulator() { }

    @Override
    public void initAndValidate() {
        speciesTree = speciesTreeInput.get();
        sampleCounts = sampleCountsInput.get();
    }

    @Override
    public void run() throws Exception {

        List<Tree> treeList = new ArrayList<>();

        try (PrintStream ps = new PrintStream(fileNameInput.get())) {
            for (int i = 0; i < nSimsInput.get(); i++) {
                Tree tree = new SimulatedGeneTree();
                tree.initByName(
                        "speciesTree", speciesTreeInput.get(),
                        "sampleCounts", sampleCountsInput.get());

                treeList.add(tree);
                ps.println(tree.toString() + ";");
            }
        }

        if (reportFileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(reportFileNameInput.get())) {
                TreeTraceAnalysis analysis = new TreeTraceAnalysis(treeList, 0.0);
                analysis.computeCredibleSet(credibilityThresholdInput.get());
                analysis.report(ps);
            }
        }
    }
}
