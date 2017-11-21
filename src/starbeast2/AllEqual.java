package starbeast2;

public class AllEqual extends MigrationModel {	
    
	@Override
	public void initAndValidate() {
	}
    
	@Override
	public double getMigration(int sourceNode, int sinkNode) {		
		return 1;
	}	
		
}
