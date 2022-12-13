package starbeast2.utils;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;



@Description("Dirichlet distribution.  p(x_1,...,x_n;alpha_1,...,alpha_n) = 1/B(alpha) prod_{i=1}^K x_i^{alpha_i - 1} " +
        "where B() is the beta function B(alpha) = prod_{i=1}^K Gamma(alpha_i)/ Gamma(sum_{i=1}^K alpha_i}. ")
public class BernoulliDistribution extends ParametricDistribution {
    final public Input<RealParameter> pInput = new Input<>("p", "success probability ", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public org.apache.commons.math.distribution.Distribution getDistribution() {
        return null;
    }

    @Override
    public double calcLogP(Function pX) {
        double p = pInput.get().getValue();
        double logP = 0;
        for (int i = 0; i < pX.getDimension(); i++) {
            double x = pX.getArrayValue(i);
            if (x==1){
            	logP += Math.log(p);
            }else if (x==0){
            	logP += Math.log(1-p);
            }else{
            	throw new IllegalArgumentException("value to calculate binomial is not 0 nor 1");
            }
        }
        return logP;
    }

	@Override
	public Double[][] sample(int size) {
		return null;
	}
}
