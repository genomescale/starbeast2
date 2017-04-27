package starbeast2.utils;

import beast.core.Input;
import beast.core.Runnable;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.*;

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

    public Tree speciesTree;
    public TraitSet sampleCounts;
    public int nSims;
    public String fileName;

    public GeneTreeSimulator() { }

    @Override
    public void initAndValidate() {
        speciesTree = speciesTreeInput.get();
        sampleCounts = sampleCountsInput.get();
        nSims = nSimsInput.get();
        fileName = fileNameInput.get();
    }

    int getTotalLineageCount(Map<Node, List<Node>> lineages) {
        int count = 0;

        for (List<Node> lineageList : lineages.values())
            count += lineageList.size();

        return count;
    }

    public Tree getSimulatedGeneTree() {

        List<Node> sortedSpeciesTreeNodes = new ArrayList<>(Arrays.asList(speciesTree.getNodesAsArray()));

        sortedSpeciesTreeNodes.sort((o1, o2) -> {
            if (o1.getHeight() < o2.getHeight())
                return -1;
            if (o1.getHeight() > o2.getHeight())
                return 1;

            return 0;
        });

        // Perform simulation

        Map<Node,List<Node>> activeLineages = new HashMap<>();
        int nextNodeNr = 0;
        double t = 0.0;

        while (getTotalLineageCount(activeLineages) > 1 || !sortedSpeciesTreeNodes.isEmpty()) {

            // Compute propensity

            double totalPropensity = 0;
            Map<Node, Double> propensities = new HashMap<>();
            for (Node speciesNode : activeLineages.keySet()) {
                int k=activeLineages.get(speciesNode).size();
                double thisProp = 0.5*k*(k-1);
                propensities.put(speciesNode, thisProp);
                totalPropensity += thisProp;
            }

            double dt = Randomizer.nextExponential(totalPropensity);

            if (!sortedSpeciesTreeNodes.isEmpty() && t + dt > sortedSpeciesTreeNodes.get(0).getHeight()) {
                Node speciesNode = sortedSpeciesTreeNodes.get(0);
                t = speciesNode.getHeight();

                activeLineages.put(speciesNode, new ArrayList<>());

                if (speciesNode.isLeaf()) {
                    int count = (int)Math.round(sampleCounts.getValue(speciesNode.getID()));

                    for (int i=0; i<count; i++) {
                        Node geneTreeSampleNode = new Node();
                        geneTreeSampleNode.setNr(nextNodeNr++);
                        geneTreeSampleNode.setHeight(speciesNode.getHeight());
                        activeLineages.get(speciesNode).add(geneTreeSampleNode);
                    }

                } else {
                    for (Node speciesChild : speciesNode.getChildren()) {
                        activeLineages.get(speciesNode).addAll(activeLineages.get(speciesChild));
                        activeLineages.get(speciesChild).clear();
                    }
                }

                sortedSpeciesTreeNodes.remove(0);

            } else {
                t += dt;

                // Coalesce a random pair of lineages

                double u = Randomizer.nextDouble()*totalPropensity;
                for (Node speciesNode : propensities.keySet()) {
                    u -= propensities.get(speciesNode);
                    if (u < 0) {
                        List<Node> lineageList = activeLineages.get(speciesNode);
                        int k = lineageList.size();

                        Node node1 = lineageList.get(Randomizer.nextInt(k));
                        Node node2;
                        do {
                            node2 = lineageList.get(Randomizer.nextInt(k));
                        } while (node2 == node1);

                        Node parent = new Node();
                        parent.setHeight(t);
                        parent.addChild(node1);
                        parent.addChild(node2);
                        lineageList.remove(node1);
                        lineageList.remove(node2);
                        lineageList.add(parent);

                        break;
                    }
                }
            }
        }

        // Return tree with remaining lineage as root
        return new Tree(activeLineages.get(speciesTree.getRoot()).get(0));
    }

    @Override
    public void run() throws Exception {
        try (PrintStream ps = new PrintStream(fileName)) {
            for (int i = 0; i < nSims; i++) {
                ps.println(getSimulatedGeneTree().toString() + ";");
            }
        }
    }
}
