package starbeast2;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulationIO extends MultispeciesPopulationModel {
    public Input<RealParameter> invgammaShapeInput = new Input<RealParameter>("alpha", "Shape of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);
    public Input<RealParameter> invgammaScaleInput = new Input<RealParameter>("beta", "Scale of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);
    public Input<Boolean> recordMathInput = new Input<Boolean>("recordMath", "Record the per-branch \"q\" and \"gamma\" intermediate calculations (default is false).", false);

    private boolean recordMath;
    private RealParameter invgammaShape;
    private RealParameter invgammaScale;
    private double[] perBranchQ;
    private double[] perBranchGamma;

    @Override
    public void initAndValidate() throws Exception {
        recordMath = recordMathInput.get();
        invgammaShape = invgammaShapeInput.get();
        invgammaScale = invgammaScaleInput.get();
    }

    @Override
    public double branchLogP(int speciesTreeNodeNumber, double[] perGenePloidy, List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
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
                partialGamma += (geneCoalescentTimes[i + 1] - geneCoalescentTimes[i]) * (geneN - i) * (geneN - (i + 1.0)) / 2.0;
            }
            
            if (geneN - geneK > 1) {
                partialGamma += (geneCoalescentTimes[geneK + 1] - geneCoalescentTimes[geneK]) * (geneN - geneK) * (geneN - (geneK + 1.0)) / 2.0;
            }

            branchGamma += partialGamma / genePloidy;
        }

        double logGammaRatio = 0.0;
        for (int i = 0; i < branchQ; i++) {
            logGammaRatio += Math.log(alpha + i);
        }

        final double logP = branchLogR + (alpha * Math.log(beta)) - ((alpha + branchQ) * Math.log(beta + branchGamma)) + logGammaRatio;

        return logP;
    }

    @Override
    public List<StateNode> initializePopSizes(TreeInterface speciesTree, double popInitial) {
        final int nBranches = speciesTree.getNodeCount();
        perBranchQ = new double[nBranches];
        perBranchGamma = new double[nBranches];

        Arrays.fill(perBranchQ, 0.0);
        Arrays.fill(perBranchGamma, 0.0);

        final List<StateNode> popSizeVectors = new ArrayList<StateNode>();
        return popSizeVectors;
    }

    @Override
    public String serialize(Node speciesTreeNode) {
        String dmv;

        if (recordMath) {
            final int speciesTreeNodeNumber = speciesTreeNode.getNr();
            dmv = String.format("q=%f,gamma=%f", perBranchQ[speciesTreeNodeNumber], perBranchGamma[speciesTreeNodeNumber]);
        } else {
            dmv = "";
        }

        return dmv;
    }
}
