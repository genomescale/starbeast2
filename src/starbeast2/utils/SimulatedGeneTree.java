package starbeast2.utils;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;

import java.util.*;

/**
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 27/04/17.
 */
public class SimulatedGeneTree extends Tree {

    public Input<Tree> speciesTreeInput = new Input<>(
            "speciesTree",
            "Species tree on which to simulate gene trees",
            Input.Validate.REQUIRED);

    public Input<TraitSet> sampleCountsInput = new Input<>(
            "sampleCounts",
            "TraitSet defining number of  samples per node in species tree.",
            Input.Validate.REQUIRED);

    public Tree speciesTree;
    public TraitSet sampleCounts;

    public SimulatedGeneTree() { }

    @Override
    public void initAndValidate() {
        speciesTree = speciesTreeInput.get();
        sampleCounts = sampleCountsInput.get();

        assignFromWithoutID(getSimulatedGeneTree());
    }

    int getTotalLineageCount(Map<Node, List<Node>> lineages) {
        int count = 0;

        for (List<Node> lineageList : lineages.values())
            count += lineageList.size();

        return count;
    }

    int getTotalSampleCount() {
        int count = 0;

        for (Node speciesNode : speciesTree.getExternalNodes())
            count += (int)Math.round(sampleCounts.getValue(speciesNode.getID()));

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
        int nextLeafNodeNr = 0;
        int nextIntNodeNr = getTotalSampleCount();
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
                        Node geneTreeSampleNode = new Node(speciesNode.getID() + String.valueOf(i + 1));
                        geneTreeSampleNode.setNr(nextLeafNodeNr++);
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

                        Node parent = new Node(String.valueOf(nextIntNodeNr));
                        parent.setNr(nextIntNodeNr++);
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

}
