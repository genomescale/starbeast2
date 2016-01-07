package starbeast2;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math.MathException;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.math.distributions.ParametricDistribution;

public class DiscreteRates extends SpeciesTreeRates {
    public Input<IntegerParameter> branchRatesInput = new Input<>("rates", "Discrete per-branch rates.", Input.Validate.REQUIRED);
    public Input<ParametricDistribution> rateDistributionInput = new Input<>("distr", "The distribution governing the rates among branches. Must have mean of 1.", Input.Validate.REQUIRED);

    private int nBins;
    private double[] binRates;
    private boolean binRatesNeedsUpdate;

    private Double[] ratesArray;
    private int nRates;
    private boolean needsUpdate;

    @Override
    public boolean requiresRecalculation() {
        binRatesNeedsUpdate = rateDistributionInput.isDirty();
        needsUpdate = binRatesNeedsUpdate || branchRatesInput.isDirty();
        return needsUpdate;
    }

    public void restore() {
        binRatesNeedsUpdate = true;
        needsUpdate = true;
        super.restore();
    }

    @Override
    public void initAndValidate() throws Exception {
        final IntegerParameter branchRates = branchRatesInput.get();
        final TreeInterface speciesTree = speciesTreeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();

        nBins = nRates = speciesNodes.length;
        binRates = new double[nBins];
        ratesArray = new Double[nRates];

        branchRates.setDimension(nRates);
        branchRates.setLower(0);
        branchRates.setUpper(nRates - 1);

        binRatesNeedsUpdate = true;
        needsUpdate = true;
    }

    private void update() {
        if (binRatesNeedsUpdate) {
            updateBinRates();
        }

        final Integer[] branchRatePointers = branchRatesInput.get().getValues();

        for (int i = 0; i < nRates; i++) {
            int b = branchRatePointers[i];
            //System.out.println(String.format("%d:%d/%d, %d:%d/%d", i, nRates, ratesArray.length, b, nBins, binRates.length));
            
            ratesArray[i] = binRates[b];
        }

        needsUpdate = false;
    }

    private void updateBinRates() {
        final ParametricDistribution rateDistribution = rateDistributionInput.get();

        try {
            for (int i = 0; i < nBins; i++) {
                binRates[i] = rateDistribution.inverseCumulativeProbability((i + 0.5) / nBins);
            }
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability!");
        }

        binRatesNeedsUpdate = false;
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
    public boolean setRate(String[] targetNames, int newBin) {
        final Set<String> s = new HashSet<>(Arrays.asList(targetNames));
        final IntegerParameter branchRates = branchRatesInput.get();
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
                branchRates.setValue(i, newBin);;
                return true;
            }
        }
        
        return false;
    }
}
