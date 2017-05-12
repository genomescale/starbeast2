package starbeast2.utils;

import beast.app.tools.SATreeTraceAnalysis;
import beast.core.util.ESS;
import beast.evolution.tree.Tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Simple extension to TTA which adds error estimates of topology probabilities to results
 */
public class TreeTraceAnalysisWithError extends SATreeTraceAnalysis {

    protected List<String> topologyList = new ArrayList<>();
    protected Map<String, Double> topologyESSs = new HashMap<>();

    @Override
    public void analyzeTree(Tree tree) {
        super.analyzeTree(tree);

        topologyList.add(uniqueNewick(tree));
    }

    public void computeErrors() {
        for (String topology : topologiesFrequencySet.getCredibleSet().credibleSetList) {
            System.out.println("Computing ESS for topology " + topology + "...");
            List<Double> presenceList = topologyList.stream().map(x -> x.equals(topology) ? 1.0 : 0.0).collect(Collectors.toList());
            topologyESSs.put(topology, ESS.calcESS(presenceList));
        }
    }

    @Override
    public void report(PrintStream ps, boolean verbose) {
        computeErrors();

        // prefix non-tabular lines with # so file can be read into R
        ps.println("# total number of trees used = " + String.valueOf(nTrees));

        // prefix non-tabular lines with # so file can be read into R
        ps.print("# \n# " + String.valueOf(topologiesFrequencySet.getCredSetProbability() * 100)
                + "% credible set");

        ps.println(" (" + String.valueOf(credibleSet.credibleSetList.size())
                + " unique tree topologies, "
                + String.valueOf(credibleSet.sumFrequency)
                + " trees in total)");

        if (verbose) {
            ps.println("Rank\tCount\tPercent\tRunning\tPercent_StdErr\tTree");
            double runningPercent = 0;
            for (int i = 0; i < credibleSet.credibleSetList.size(); i++) {
                double p = credibleSet.getFrequency(i, topologiesFrequencySet) / (double)nTrees;
                double percent = p * 100.0;
                String topology = credibleSet.credibleSetList.get(i);
                double percentSD = Math.sqrt(p*(1-p)/topologyESSs.get(topology)) * 100.0;
                runningPercent += percent;

                ps.print((i + 1) + "\t");
                ps.print(credibleSet.getFrequency(i, topologiesFrequencySet) + "\t");
                ps.format("%.5f%%\t", percent);
                ps.format("%.5f%%\t", runningPercent);
                ps.format("%.5f%%\t", percentSD);
                ps.format("%s\n", topology);
            }
        }
    }

}
