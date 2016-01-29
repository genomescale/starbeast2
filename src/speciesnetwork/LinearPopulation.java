package speciesnetwork;

import java.text.DecimalFormat;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

/**
* @author Huw Ogilvie
 */

public class LinearPopulation extends PopulationSizeModel {
    public Input<RealParameter> topPopSizesInput =
            new Input<>("topPopSizes", "Population sizes at the top (rootward) end of each branch.", Validate.REQUIRED);
    public Input<RealParameter> tipPopSizesInput =
            new Input<>("tipPopSizes", "Population sizes at the tips of leaf branches.", Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double branchLogP(int speciesNetworkPopNumber, NetworkNode speciesNetworkNode, double[] perGenePloidy,
                             List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final int nGenes = perGenePloidy.length;

        final double branchTopPopSize = topPopSizes.getValue(speciesNetworkPopNumber);

        final double branchTipPopSize;
        if (speciesNetworkNode.isLeaf()) {
            branchTipPopSize = tipPopSizes.getValue(speciesNetworkPopNumber);
        } else if (speciesNetworkNode.isReticulation()) {
            final int childNodeNumber;
            if (speciesNetworkNode.getLeftChild() != null)
                childNodeNumber = speciesNetworkNode.getLeftChild().getNr();
            else
                childNodeNumber = speciesNetworkNode.getRightChild().getNr();
            branchTipPopSize = topPopSizes.getValue(childNodeNumber*2);
        } else {
            final int leftChildNodeNumber = speciesNetworkNode.getLeftChild().getNr();
            final int rightChildNodeNumber = speciesNetworkNode.getRightChild().getNr();
            branchTipPopSize = topPopSizes.getValue(leftChildNodeNumber*2) +
                               topPopSizes.getValue(rightChildNodeNumber*2+1);
        }

        // set the root branch heights for each gene tree (to equal the highest gene tree root node)
        double tallestGeneTreeHeight = 0.0;
        if (speciesNetworkNode.isRoot()) {
            for (int j = 0; j < nGenes; j++) {
                tallestGeneTreeHeight = Math.max(tallestGeneTreeHeight, branchCoalescentTimes.get(j)[branchEventCounts[j]]);
            }
            for (int j = 0; j < nGenes; j++) {
                branchCoalescentTimes.get(j)[branchEventCounts[j] + 1] = tallestGeneTreeHeight;
            }
        }

        return linearLogP(branchTopPopSize, branchTipPopSize, perGenePloidy,
                          branchCoalescentTimes, branchLineageCounts, branchEventCounts);
    }

    // copied from *BEAST v2.3, with small modifications to work with starbeast2
    protected static double linearLogP(double topPopSize, double tipPopSize, double[] perGenePloidy, List<Double[]> branchCoalescentTimes,
                                       int[] branchLineageCounts, int[] branchEventCounts) {
        final int nGenes = perGenePloidy.length;

        double logP = 0.0;
        for (int j = 0; j < nGenes; j++) {
            final double fPopSizeTop = topPopSize * perGenePloidy[j];
            final double fPopSizeBottom = tipPopSize * perGenePloidy[j];
            final Double[] fTimes = branchCoalescentTimes.get(j);
            final int k = branchEventCounts[j];
            final int nLineagesBottom = branchLineageCounts[j];

            final double a = (fPopSizeTop - fPopSizeBottom) / (fTimes[k + 1] - fTimes[0]);
            // final double b = fPopSizeBottom;

            if (Math.abs(fPopSizeTop - fPopSizeBottom) < 1e-10) {
                // use approximation for small values to bypass numerical instability
                for (int i = 0; i <= k; i++) {
                    final double fTimeip1 = fTimes[i + 1];
                    final double fPopSize = a * (fTimeip1 - fTimes[0]) + fPopSizeBottom;
                    if( i < k ) {
                        logP += -Math.log(fPopSize);
                    }

                    // slope = 0, so population function is constant
                    final int i1 = nLineagesBottom - i;
                    logP -= (i1 * (i1 - 1.0) / 2.0) * (fTimeip1 - fTimes[i]) / fPopSize;
                }
            } else {
                final double vv = fPopSizeBottom - a * fTimes[0];
                for (int i = 0; i <= k; i++) {
                    final double fPopSize = a * fTimes[i + 1] + vv;
                    if( i < k ) {
                        logP += -Math.log(fPopSize);
                    }
                    final double f = fPopSize / (a * fTimes[i] + vv);

                    final int i1 = nLineagesBottom - i;
                    logP += -(i1 * (i1 - 1.0) / 2.0) * Math.log(f) / a;
                }
            }
        }

        return logP;
    }

    @Override
    public void initPopSizes(int nPopulation) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        //final int nSpecies = (nSpecies + 1) / 2;

        topPopSizes.setDimension(nPopulation);
        tipPopSizes.setDimension(nPopulation);
    }

    @Override
    public void initPopSizes(double popInitial) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();

        // TODO consider reticulation node ???
        for (int i = 0; i < topPopSizes.getDimension(); i++) {
            topPopSizes.setValue(i, popInitial / 2.0);
        }
        for (int i = 0; i < tipPopSizes.getDimension(); i++) {
            tipPopSizes.setValue(i, popInitial);
        }
    }

    @Override
    public void serialize(NetworkNode speciesNetworkNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter topPopSizes = topPopSizesInput.get();
        final RealParameter tipPopSizes = tipPopSizesInput.get();
        final int speciesNetworkNodeNumber = speciesNetworkNode.getNr();

        // TODO ???
        final double branchTopPopSize = topPopSizes.getValue(speciesNetworkNodeNumber);

        double branchTipPopSize;
        if (speciesNetworkNode.isLeaf()) {
            branchTipPopSize = tipPopSizes.getValue(speciesNetworkNodeNumber);
        } else {
            final int leftChildNodeNumber = speciesNetworkNode.getLeftChild().getNr();
            final int rightChildNodeNumber = speciesNetworkNode.getRightChild().getNr();
            branchTipPopSize = topPopSizes.getValue(leftChildNodeNumber) + topPopSizes.getValue(rightChildNodeNumber);
        }

        buf.append("dmv={");
        if (df == null) {
            buf.append(branchTopPopSize);
            buf.append(",");
            buf.append(branchTipPopSize);
        } else {
            buf.append(df.format(branchTopPopSize));
            buf.append(",");
            buf.append(df.format(branchTipPopSize));
        }
        buf.append("}");
    }
}
