package starbeast2;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.StateNode;

import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;


/**
 * adapted by Nicola Felix Mueller from the tree logger
 */
@Description("log trees that also contain the node state probabilities")
public class SpeciesTreeLoggerWithGeneFlow extends Tree implements Loggable {
    public Input<ConstantWithGeneFlow> popModelInput = new Input<>("populationModel", "Population model used to infer the multispecies coalescent probability for this gene");
    public Input<BranchRateModel.Base> clockModelInput = new Input<BranchRateModel.Base>("branchratemodel", "rate to be logged with branches of the tree");
    public Input<List<Function>> parameterInput = new Input<List<Function>>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());
    public Input<Boolean> substitutionsInput = new Input<Boolean>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<Integer> decimalPlacesInput = new Input<Integer>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);
    
    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;
    boolean takeMax = true;
    boolean conditionals = true;

    private DecimalFormat df;	
	
    @Override
    public void initAndValidate() {
    	
        if (parameterInput.get().size() == 0 && clockModelInput.get() == null) {
        	someMetaDataNeedsLogging = false;
        	return;
            //throw new Exception("At least one of the metadata and branchratemodel inputs must be defined");
        }
    	someMetaDataNeedsLogging = true;
    	// without substitution model, reporting substitutions == reporting branch lengths 
        if (clockModelInput.get() != null) {
        	substitutions = substitutionsInput.get();
        }       

        int dp = decimalPlacesInput.get();

        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }       
        
    }

    @Override
    public void init(PrintStream out) {
    	popModelInput.get().migrationModelInput.get().speciesTreeInput.get().init(out);
    }

    @Override
    public void log(long nSample, PrintStream out) {
    	// build migration rate map	
    	popModelInput.get().calculateIntervals();
    	popModelInput.get().stateToNodeMap();
    	
        // make sure we get the current version of the inputs
        Tree tree = (Tree) popModelInput.get().migrationModelInput.get().speciesTreeInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
//        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), metadata, branchRateModel));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }

    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {

    	StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
//        if (!node.isLeaf()) {
	        buf.append("[&species=" + node.getNr() + ",Ne=" +  popModelInput.get().getNodeNe(node.getNr()) + "" + popModelInput.get().getAllMigrationRates(node.getNr()));		

	        
	        
	        buf.append("]");
//        }else{
//	        buf.append(String.format("[&Ne=%.3f", speciesTreeInput.get().getNodeNe(node.getNr())));		        
//	        buf.append(']');
//        }
        
        buf.append(":");
        if (substitutions) {
            appendDouble(buf, node.getLength() * branchRateModel.getRateForBranch(node));
        } else {
            appendDouble(buf, node.getLength());
        }
        return buf.toString();
    }


    @Override
    public void close(PrintStream out) {
    	popModelInput.get().migrationModelInput.get().speciesTreeInput.get().close(out);
    }

	
}
