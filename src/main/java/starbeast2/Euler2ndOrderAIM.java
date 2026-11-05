package starbeast2;


import org.apache.commons.math3.util.FastMath;

public class Euler2ndOrderAIM {

	double epsilon;
	double max_step;
	
	double[][] migration_rates;
	int[] multiplicator;
	int[][] indicators;
	boolean[] isConnected;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
	double[] sumStates;
	boolean hasIndicators;
	boolean hasMultiplicator;
	
	int iterations;

	
	public Euler2ndOrderAIM(double[][] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = false;
    	hasMultiplicator = false;
    	
    	iterations=0;
	}
	
	public Euler2ndOrderAIM(double[][] migration_rates, int[][] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.indicators = indicators;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];
    	hasIndicators = true;
    	hasMultiplicator = false;
    	
    	iterations=0;
	}

	public Euler2ndOrderAIM(int[] multiplicator, double[][] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.multiplicator = multiplicator;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = false;
    	hasMultiplicator = true; 
    	
    	iterations=0;
	}
	
	public Euler2ndOrderAIM(int[] multiplicator, double[][] migration_rates, int[][] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.multiplicator = multiplicator;
        this.indicators = indicators;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = true;
    	hasMultiplicator = true;   
    	
    	iterations=0;
	}

	public Euler2ndOrderAIM(int[] multiplicator, double[][] migration_rates, int[][] indicators, boolean[] isConnected, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.multiplicator = multiplicator;
        this.indicators = indicators;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
        this.isConnected = isConnected;
    	sumStates = new double[states];    	
    	hasIndicators = true;
    	hasMultiplicator = true;   
    	
    	iterations=0;
	}

	
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
		clearArray(pDotDot, length);
		clearArray(pDotDotDot, length);

		while (duration > 0){
	    	iterations++;
			//pDot = new double[length];
			clearArray(pDot, length);
			computeDerivatives(p, pDot, pDotDot, pDotDotDot, length);
			computeSecondDerivate(p, pDot, pDotDot, length);
			approximateThirdDerivate(p, pDot, pDotDot, pDotDotDot, length);
			duration = updateP(duration, p,  pDot, pDotDot, pDotDotDot, length - 1);
						
			if (iterations>1000000000){
				System.err.println("too many iterations, return negative infinity");
				p[length-1] = Double.NEGATIVE_INFINITY;
			}
			
			if (p[length-1]==Double.NEGATIVE_INFINITY)
				break;
 
			
				
		}			
	}	
	
	public void calculateConnectedValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
		clearArray(pDotDot, length);
		clearArray(pDotDotDot, length);

		while (duration > 0){
	    	iterations++;
			//pDot = new double[length];
			clearArray(pDot, length);
			computeConnectedDerivatives(p, pDot, pDotDot, pDotDotDot, length);
			computeConnectedSecondDerivate(p, pDot, pDotDot, length);
			approximateConnectedThirdDerivate(p, pDot, pDotDot, pDotDotDot, length);
			duration = updateConnectedP(duration, p,  pDot, pDotDot, pDotDotDot, length - 1);
			
			if (iterations>1000000){
				System.err.println("too many iterations, return negative infinity");
				p[length-1] = Double.NEGATIVE_INFINITY;
			}
			
			if (p[length-1]==Double.NEGATIVE_INFINITY)
				break;		
				
		}			
	}	

	
	private void clearArray(double[] v, int n) {
		for (int i = 0; i < n; i++) {
			v[i] = 0.0;
		}		
	}

	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
		final double max_dotdotdot = maxAbs(pDotDotDot, length);	
		//double timeStep = FastMath.min(FastMath.pow(epsilon*6/max_dotdotdot, C), FastMath.min(duration, max_step));

		double timeStep = FastMath.min(FastMath.cbrt(epsilon*6/max_dotdotdot), FastMath.min(duration, max_step));

		double timeStepSquare = timeStep*timeStep*0.5;
		
		for (int i = 0; i < length; i++){
			double new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
			double diff = FastMath.abs(new_val - p[i]);
			int its = 0;
			while (new_val > 1 || new_val < 0 || diff>0.2){
				timeStep *= 0.9;
				timeStepSquare = timeStep*timeStep*0.5;
				new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
				diff = FastMath.abs(new_val - p[i]);
				its++;
				if (its>10000){
					System.err.println("cannot find proper time step, skip these parameter values");
					p[length] = Double.NEGATIVE_INFINITY;
					break;
				}

			}			
		}
		
		if (p[length]==Double.NEGATIVE_INFINITY)
			return 0.0;
		
		doUpdating(timeStep, timeStepSquare, p, pDot, pDotDot, length);
		p[length] += pDot[length]*timeStep;
		duration -= timeStep;
		return duration;
	}
	
	
	private double maxAbs(double[] pDotDotDot, int length) {
		double max_dotdotdot = 0.0;
		for (int i = 0; i < length; i++) {
			max_dotdotdot = FastMath.max(max_dotdotdot, FastMath.abs(pDotDotDot[i]));
		}
		return max_dotdotdot;
	}


	static final double C = 1.0/3.0;
	
	private void doUpdating(final double timeStep, final double timeStepSquare, double[] p, double[] pDot, double[] pDotDot, int length){
		updateP2(timeStep, timeStepSquare, p, length, pDot, pDotDot);
		
		// normalize to ensure stability
		for (int i = 0; i < lineages; i ++){
			normalise(i, p, length);
		}
	}
	    
	private void normalise(final int i, final double[] p, int length) {
		final int k = states*i;
		double linSum = 0;
		
		for (int j = 0; j < states; j++){
			if (p[k+j]>=0.0){
				linSum += p[k+j];
			}else{
//				System.err.println(Arrays.toString(p));
				p[length-1] = Double.NEGATIVE_INFINITY;
				return;
			}
		}
		for (int j = 0; j < states; j++){
			p[k+j] /= linSum;
		}
	}

	private void updateP2(final double timeStep, final double timeStepSquare, final double[] p, final int length, final double[] pDot,
			final double[] pDotDot) {
		for (int i = 0; i < length; i++)
			p[i] += pDot[i]*timeStep + pDotDot[i]*timeStepSquare;	
	}

	public void computeDerivatives (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		
    	double migrates;
    	// Compute the sum of line state probabilities for each state
     	sumStates = new double[states];
     	if (hasMultiplicator){
	    	for (int i = 0; i<lineages; i++) {
	    		int k = states * i;
	    		for (int j = 0; j<states; j++) {
					sumStates[j] += multiplicator[i]*p[k+j]; 
	    		}
	    	}
     	}else{
	    	for (int i = 0; i<lineages; i++) {
	    		int k = states * i;
	    		for (int j = 0; j<states; j++) {
					sumStates[j] += p[k+j];
	    		}
	    	}
     	}
    		
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		double[] tCR =  new double[states];
    		double sumCoal = 0;
    		int currlin = states*i;
    		for (int j = 0; j<states; j++){
    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
    			sumCoal += p[currlin+j]*tCR[j];
    		}
         	if (hasMultiplicator){
         		pDot[length-1] -= multiplicator[i]*sumCoal;
         	}else{
         		pDot[length-1] -= sumCoal;        		
         	}
    		for (int j = 0; j < states; j++){    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[currlin+j] = coal;
    			pDotDotDot[currlin+j] = coal;
    			pDot[currlin+j] +=	p[currlin+j] * coal;
    		}// j

    	}
    	
    	
    	// Calculate the probability of a lineage changing states
    	if (hasIndicators){
			for (int j = 0; j < indicators.length; j++){
				int source = indicators[j][0];
				int sink = indicators[j][1];
				double mrate = migration_rates[source][sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = p[states*i+source]*mrate;
					pDot[states*i+sink] += migrates;
					pDot[states*i+source] -= migrates;  			
		    	}
	    	}    
    	}else{
        	for (int i = 0; i<lineages; i++){
        		int currlin = states*i;
        		// Calculate the probability of a lineage changing states
        		for (int j = 0; j < states; j++){
        			double pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k][j] -
    							pj*migration_rates[j][k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j    			 
        		}// j
        	}// lineages       		
    	}
    	
		pDot[length-1]  /= 2;

    }
        
    public void computeSecondDerivate (double[] p, double[] pDot, double[] pDotDot, int length){  
    	double[] sumDotStates = new double[states];
    	if (hasMultiplicator){
	    	for (int i = 0; i<lineages; i++)
	    		for (int j = 0; j<states; j++)
	    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];   
    	}else{
	    	for (int i = 0; i<lineages; i++)
	    		for (int j = 0; j<states; j++)
	    			sumDotStates[j] += pDot[states*i+j];   
    	}
	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));     
    		
    		for (int j = 0; j<states; j++)
    			pDotDot[currlin+j] *= pDot[currlin+j];
    		
    		if (hasMultiplicator){
    			pDotDot[length-1] -= multiplicator[i]*pCoalRate;
    		}else{
    			pDotDot[length-1] -= pCoalRate;    			
    		}    		

    		for (int j = 0; j<states; j++){
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    		}// j    		
    	}// lineages 
    	
		double migrates;
		
		
		// Calculate the probability of a lineage changing states
		if (hasIndicators){
			for (int j = 0; j < indicators.length; j++){
				int source = indicators[j][0];
				int sink = indicators[j][1];
				double mrate = migration_rates[source][sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = pDot[states*i+source]*mrate;
					pDotDot[states*i+sink] += migrates;
					pDotDot[states*i+source] -= migrates;  			
		    	}
	    	}    
		}else{
        	for (int i = 0; i<lineages; i++){
        		int currlin = states*i;
        		// Calculate the probability of a lineage changing states
        		for (int j = 0; j < states; j++){
        			double pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k][j] -
    							pj*migration_rates[j][k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j    			 
        		}// j
        	}// lineages    
			
		}
		pDotDot[length-1] /= 2;

    }
    
    public void approximateThirdDerivate (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
    	
    	
    	double migrates;
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		// Calculate the probability of a lineage changing states
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pDotDotDot[currlin+j] *= pDotDot[currlin+j];
    	}
    	
		// Calculate the probability of a lineage changing states
    	if (hasIndicators){
			for (int j = 0; j < indicators.length; j++){
				int source = indicators[j][0];
				int sink = indicators[j][1];
				double mrate = migration_rates[source][sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = pDotDot[states*i+source]*mrate;
					pDotDotDot[states*i+sink] += migrates;
					pDotDotDot[states*i+source] -= migrates;  			
		    	}
	    	}    
    	}else{
			for (int j = 0; j < states; j++){
				for (int k = 0; k < states; k++){  
					double mrate = migration_rates[j][k];
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[states*i+j]*mrate;
						pDotDotDot[states*i+k] += migrates;
						pDotDotDot[states*i+j] -= migrates;  			
	    			} 			
				}
	    	}
    	}

    }


    
	private double updateConnectedP (double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
		final double max_dotdotdot = maxAbs(pDotDotDot, length);	
		//double timeStep = FastMath.min(FastMath.pow(epsilon*6/max_dotdotdot, C), FastMath.min(duration, max_step));

		double timeStep = FastMath.min(FastMath.cbrt(epsilon*6/max_dotdotdot), FastMath.min(duration, max_step));

		double timeStepSquare = timeStep*timeStep*0.5;
		
		for (int j = 0; j<states; j++) {
			if (isConnected[j]){
				for (int i = 0; i<lineages; i++) {
					int ind = states*i+j;
					double new_val = p[ind] + pDot[ind]*timeStep + pDotDot[ind]*timeStepSquare;
					double diff = FastMath.abs(new_val - p[ind]);
					int its = 0;
					while (new_val > 1 || new_val < 0 || diff>0.2){
						timeStep *= 0.9;
						timeStepSquare = timeStep*timeStep*0.5;
						new_val = p[ind] + pDot[ind]*timeStep + pDotDot[ind]*timeStepSquare;
						diff = FastMath.abs(new_val - p[ind]);
						its++;
						if (its>10000){
//							System.err.println("cannot find proper time step, skip these parameter values");
							p[length] = Double.NEGATIVE_INFINITY;
							break;
						}

					}			

				}
			}
		}
		
		if (p[length]==Double.NEGATIVE_INFINITY)
			return 0.0;
		
		doConnectedUpdating(timeStep, timeStepSquare, p, pDot, pDotDot, length);
		p[length] += pDot[length]*timeStep;
		duration -= timeStep;
		return duration;
	}
	
	private void doConnectedUpdating(final double timeStep, final double timeStepSquare, double[] p, double[] pDot, double[] pDotDot, int length){
		for (int j = 0; j<states; j++) {
			if (isConnected[j]){
				for (int i = 0; i<lineages; i++) {
					int ind = states*i+j;
					p[ind] += pDot[ind]*timeStep + pDotDot[ind]*timeStepSquare;
				}
			}
		}
		
		// normalize to ensure stability
		for (int i = 0; i < lineages; i ++){
			normalise(i, p, length);
		}
	}

	
	private void normaliseConnected(final int i, final double[] p, int length) {
		final int k = states*i;
		double linSum = 0;
		
		for (int j = 0; j < states; j++){
			if (isConnected[j]){
				if (p[k+j]>=0.0){
					linSum += p[k+j];
				}else{
	//				System.err.println(Arrays.toString(p));
					p[length-1] = Double.NEGATIVE_INFINITY;
					return;
				}
			}
		}
		for (int j = 0; j < states; j++){
			p[k+j] /= linSum;
		}
	}

	
	    

	   
	public void computeConnectedDerivatives (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		
    	double migrates;
    	// Compute the sum of line state probabilities for each state
     	sumStates = new double[states];
		for (int j = 0; j<states; j++) {
			if (isConnected[j]){
				for (int i = 0; i<lineages; i++) {
    				sumStates[j] += multiplicator[i]*p[states * i+j];
				}
    		}
    	}
    		
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i < lineages; i++){
    		double[] tCR =  new double[states];
    		double sumCoal = 0;
    		int currlin = states*i;
    		for (int j = 0; j<states; j++){
    			if (isConnected[j]){
	    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
	    			sumCoal += p[currlin+j]*tCR[j];
    			}
    		}
     		pDot[length-1] -= multiplicator[i]*sumCoal;
         
    		for (int j = 0; j < states; j++){    
    			if (isConnected[j]){
	    			// Calculate the Derivate of p:
	    			double coal = sumCoal - tCR[j];
	    			pDotDot[currlin+j] = coal;
	    			pDotDotDot[currlin+j] = coal;
	    			pDot[currlin+j] +=	p[currlin+j] * coal;
//    			}else{
//	    			pDotDot[currlin+j] = 0;
//	    			pDotDotDot[currlin+j] = 0;
//	    			pDot[currlin+j] = 0;   				
    			}
    		}// j
    	}
    	
    	
    	// Calculate the probability of a lineage changing states
    	if (indicators!=null){
			for (int j = 0; j < indicators.length; j++){
				int source = indicators[j][0];
				int sink = indicators[j][1];
				double mrate = migration_rates[source][sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = p[states*i+source]*mrate;
					pDot[states*i+sink] += migrates;
					pDot[states*i+source] -= migrates;  			
		    	}
	    	}    
    	}

    	
		pDot[length-1]  /= 2;
		

    }
        
    public void computeConnectedSecondDerivate (double[] p, double[] pDot, double[] pDotDot, int length){  
    	double[] sumDotStates = new double[states];
		for (int j = 0; j<states; j++){
			if (isConnected[j]){
				for (int i = 0; i<lineages; i++){
	    				sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
				}
			}
		}
	    			
	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			if (isConnected[j])
    				pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));     
    			
    		
    		for (int j = 0; j<states; j++)
    			if (isConnected[j]){
    				pDotDot[currlin+j] *= pDot[currlin+j];
    			}
    		
			pDotDot[length-1] -= multiplicator[i]*pCoalRate;

    		for (int j = 0; j<states; j++){
    			if (isConnected[j]){
    				pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    			}
    		}// j    		
    	}// lineages 
    	
		double migrates;
		
		
    	if (indicators!=null){
			// Calculate the probability of a lineage changing states
			for (int j = 0; j < indicators.length; j++){
				int source = indicators[j][0];
				int sink = indicators[j][1];
				double mrate = migration_rates[source][sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = pDot[states*i+source]*mrate;
					pDotDot[states*i+sink] += migrates;
					pDotDot[states*i+source] -= migrates;  			
		    	}
	    	}	
    	}
		pDotDot[length-1] /= 2;
		
    }
    
    public void approximateConnectedThirdDerivate (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
    	    	
    	double migrates;
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
		for (int j = 0; j<states; j++){
			if (isConnected[j]){
				for (int i = 0; i<lineages; i++){    		
	    		// Calculate the probability of a lineage changing states
	    		int currlin = states*i;
	    				pDotDotDot[currlin+j] *= pDotDot[currlin+j];
	    		}
			}
    	}
    	
    	if (indicators!=null){
			// Calculate the probability of a lineage changing states
			for (int j = 0; j < indicators.length; j++){
				int source = indicators[j][0];
				int sink = indicators[j][1];
				double mrate = migration_rates[source][sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = pDotDot[states*i+source]*mrate;
					pDotDotDot[states*i+sink] += migrates;
					pDotDotDot[states*i+source] -= migrates;  			
		    	}
	    	}    
    	}
    }

    
    
}

