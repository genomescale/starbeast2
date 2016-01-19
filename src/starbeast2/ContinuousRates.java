package starbeast2;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

@Description("Uses continuous rates stored in clock space, which should have a standard normal prior distribution." +
"These are scaled to have a spread of 'stdev' and a mean in real space of 1.")
public class ContinuousRates extends SpeciesTreeRates {
    public Input<RealParameter> stdevInput = new Input<>("stdev", "The standard deviation to apply to rates.", Input.Validate.REQUIRED);
    public Input<RealParameter> logRatesInput = new Input<>("logRates", "Per-branch rates.", Input.Validate.REQUIRED);

    private Double[] realRatesArray;
    private int nRates;
    private boolean needsUpdate;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = logRatesInput.isDirty() || stdevInput.isDirty();
        return needsUpdate;
    }

    public void restore() {
        needsUpdate = true;
        super.restore();
    }

    @Override
    public void initAndValidate() throws Exception {
        final RealParameter branchRates = logRatesInput.get();
        final TreeInterface speciesTree = speciesTreeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();

        nRates = speciesNodes.length;

        branchRates.setDimension(nRates);
        realRatesArray = new Double[nRates];

        needsUpdate = true;
    }

    private void update() {
        final Double realMean = meanRateInput.get().getValue();
        final Double stdev = stdevInput.get().getValue();
        final Double logMean = Math.log(realMean) - (0.5 * stdev * stdev);
        final Double[] logRatesArray = logRatesInput.get().getValues();
        for (int i = 0; i < logRatesArray.length; i++) {
            realRatesArray[i] = Math.exp((logRatesArray[i] * stdev) + logMean);
        }

        needsUpdate = false;
    }

    @Override
    Double[] getRatesArray() {
        if (needsUpdate) {
            update();
        }

        return realRatesArray;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate) {
            update();
        }

        return realRatesArray[node.getNr()];
    }

    // for testing purposes
    public boolean setRate(String[] targetNames, double newRate) {
        final Set<String> s = new HashSet<>(Arrays.asList(targetNames));
        final RealParameter branchRates = logRatesInput.get();
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
