package starbeast2;


import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;


@Description("calculates a prior from a distribution based on the ")
public class SpeciesNePrior extends Distribution {
    final public Input<Function> m_x = new Input<>("x", "point at which the density is calculated", Validate.REQUIRED);
    final public Input<ParametricDistribution> ancestralInput = new Input<>("ancestralDistribution", "distribution used to calculate prior, e.g. normal, beta, gamma.", Validate.REQUIRED);
    final public Input<ParametricDistribution> presentInput = new Input<>("presentDistribution", "distribution used to calculate prior, e.g. normal, beta, gamma.", Validate.REQUIRED);
    final public Input<Tree> speciesTreeInput = new Input<>("speciesTree", "input of species tree", Validate.REQUIRED);
    
    /**
     * shadows m_distInput *
     */
    protected ParametricDistribution ancestral;
    protected ParametricDistribution present;
    
    @Override
    public void initAndValidate() {
        ancestral = ancestralInput.get();
        present = presentInput.get();
        calculateLogP();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;
    	Tree speciesTree = speciesTreeInput.get();
    	// get all the leaf nodes and its daughters
    	ArrayList<Double> ratios = new ArrayList<>();
    	for (int i = 0; i < speciesTree.getNodeCount(); i++){
    		if (!speciesTree.getNode(i).isLeaf()){
    			double Npar = m_x.get().getArrayValue(speciesTree.getNode(i).getNr());
    			double N1 = m_x.get().getArrayValue(speciesTree.getNode(i).getLeft().getNr());
    			double N2 = m_x.get().getArrayValue(speciesTree.getNode(i).getRight().getNr());
    			double nom = (N1+N2)/2;
    			ratios.add(Npar/nom);
    		}else{
            	Double[] dP = new Double[1];
            	dP[0] = m_x.get().getArrayValue(speciesTree.getNode(i).getNr());        	
            	RealParameter dummyParam =  new RealParameter(dP);
                logP += present.calcLogP(dummyParam);   			
    		}
    			
    	}
    		
    	
        Function x = m_x.get();
        if (x instanceof RealParameter || x instanceof IntegerParameter) {
        	// test that parameter is inside its bounds
            double l = 0.0;
            double h = 0.0;
        	if (x instanceof RealParameter) {
                l = ((RealParameter) x).getLower();
                h = ((RealParameter) x).getUpper();
        	} else {
                l = ((IntegerParameter) x).getLower();
                h = ((IntegerParameter) x).getUpper();
        	}
            for (int i = 0; i < x.getDimension(); i++) {
            	double value = x.getArrayValue(i);
            	if (value < l || value > h) {
            		logP = Double.NEGATIVE_INFINITY;
            		return Double.NEGATIVE_INFINITY;
            	}
            }
        }
        for (int i = 0; i < ratios.size(); i++) {
        	Double[] dP = new Double[1];
        	dP[0] = ratios.get(i);        	
        	RealParameter dummyParam =  new RealParameter(dP);
            logP += ancestral.calcLogP(dummyParam);
        }
        return logP;
    }
    
    /** return name of the parameter this prior is applied to **/
    public String getParameterName() {
    	if (m_x.get() instanceof BEASTObject) {
    		return ((BEASTObject) m_x.get()).getID();
    	}
    	return m_x.get() + "";
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }
}
