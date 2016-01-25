package starbeast2;

import java.text.DecimalFormat;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulationIO extends PopulationSizeModel {
    public Input<RealParameter> invgammaShapeInput =
            new Input<>("alpha", "Shape of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);
    public Input<RealParameter> invgammaScaleInput =
            new Input<>("beta", "Scale of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double branchLogP(int speciesNetworkNodeNumber, NetworkNode speciesNetworkNode, double[] perGenePloidy,
                             List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter invgammaShape = invgammaShapeInput.get();
        final RealParameter invgammaScale = invgammaScaleInput.get();
        final double alpha = invgammaShape.getValue();
        final double beta = invgammaScale.getValue();
        
        final int nGenes = perGenePloidy.length;

        int branchQ = 0;
        double branchLogR = 0.0;
        double branchGamma = 0.0;

        for (int j = 0; j < nGenes; j++) {
            final int geneN = branchLineageCounts[j];
            final Double[] geneCoalescentTimes = branchCoalescentTimes.get(j);
            final int geneK = branchEventCounts[j];
            final double genePloidy = perGenePloidy[j]; 
            branchLogR -= geneK * Math.log(genePloidy);
            branchQ += geneK;

            double partialGamma = 0.0;
            for (int i = 0; i < geneK; i++) {
                partialGamma += (geneCoalescentTimes[i + 1] - geneCoalescentTimes[i])
                                * (geneN - i) * (geneN - i - 1.0) / 2.0;
            }
            
            if (geneN - geneK > 1) {
                partialGamma += (geneCoalescentTimes[geneK + 1] - geneCoalescentTimes[geneK])
                                * (geneN - geneK) * (geneN - geneK - 1.0) / 2.0;
            }

            branchGamma += partialGamma / genePloidy;
        }

        double logGammaRatio = 0.0;
        for (int i = 0; i < branchQ; i++) {
            logGammaRatio += Math.log(alpha + i);
        }

        return branchLogR + (alpha * Math.log(beta)) - ((alpha + branchQ) * Math.log(beta + branchGamma)) + logGammaRatio;
    }

    @Override
    public void initPopSizes(int nPopulation) {
        // do nothing
    }

    @Override
    public void initPopSizes(double popInitial) {
        // do nothing
    }

    @Override
    public void serialize(NetworkNode speciesNetworkNode, StringBuffer buf, DecimalFormat df) {
        // do nothing
    }
}
