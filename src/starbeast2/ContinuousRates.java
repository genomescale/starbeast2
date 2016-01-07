package starbeast2;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public class ContinuousRates extends SpeciesTreeRates {
    public Input<RealParameter> branchRatesInput = new Input<>("rates", "Per-branch rates.", Input.Validate.REQUIRED);

    private Double[] ratesArray;
    private int nRates;
    private boolean needsUpdate;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = branchRatesInput.isDirty();
        return needsUpdate;
    }

    public void restore() {
        needsUpdate = true;
        super.restore();
    }

    @Override
    public void initAndValidate() throws Exception {
        final RealParameter branchRates = branchRatesInput.get();
        final TreeInterface speciesTree = speciesTreeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();

        nRates = speciesNodes.length;

        branchRates.setDimension(nRates);

        needsUpdate = true;
    }

    private void update() {
        ratesArray = branchRatesInput.get().getValues();

        needsUpdate = false;
    }

    @Override
    Double[] getRatesArray() {
        if (needsUpdate) {
            update();
        }

        return ratesArray;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate) {
            update();
        }

        return ratesArray[node.getNr()];
    }

    // for testing purposes
    public boolean setRate(String[] targetNames, double newRate) {
        final Set<String> s = new HashSet<>(Arrays.asList(targetNames));
        final RealParameter branchRates = branchRatesInput.get();
        final Node[] nodeArray = speciesTreeInput.get().getNodesAsArray();

        for (int i = 0; i < nodeArray.length; i++) {
            Node n = nodeArray[i];
            final HashSet<String> comparison = new HashSet<>();
            if (n.isLeaf()) {
                final String leafName = n.getID();
                comparison.add(leafName);
            } else {
                for (Node l: n.getAllLeafNodes()) {
                    final String leafName = l.getID();
                    comparison.add(leafName);
                }
            }

            if (s.equals(comparison)) {
                branchRates.setValue(i, newRate);;
                return true;
            }
        }
        
        return false;
    }
}
