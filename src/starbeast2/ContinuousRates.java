package starbeast2;

import java.util.HashMap;
import java.util.Map;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public class ContinuousRates extends SpeciesTreeRates {
    public Input<TreeInterface> speciesTreeInput = new Input<>("tree", "Species tree to apply per-branch rates to.");
    public Input<RealParameter> branchRatesInput = new Input<>("rates", "Per-branch rates.");

    final private Map<Node, Integer> nodeIndexMap = new HashMap<>();

    @Override
    public void initAndValidate() throws Exception {
        final RealParameter branchRates = branchRatesInput.get();
        final TreeInterface speciesTree = speciesTreeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();
        final int nNodes = speciesNodes.length;

        branchRates.setDimension(nNodes);

        for (int i = 0; i < nNodes; i++) {
            final Node nodeI = speciesNodes[i];
            final int stableNodeIndex = nodeI.getNr();
            nodeIndexMap.put(nodeI, stableNodeIndex);
        }
    }

    @Override
    double[] getRatesArray() {
        final Double[] branchRates = branchRatesInput.get().getValues();
        final TreeInterface speciesTree = speciesTreeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();
        final int nNodes = speciesNodes.length;

        final double[] newNodeRates = new double[nNodes];

        for (int i = 0; i < nNodes; i++) {
            final Node nodeI = speciesNodes[i];
            final int newNodeIndex = nodeI.getNr();
            final int stableNodeIndex = nodeIndexMap.get(nodeI);
            newNodeRates[newNodeIndex] = branchRates[stableNodeIndex];
        }

        return newNodeRates;
    }
}
