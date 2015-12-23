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
import beast.core.BEASTObject;
import beast.core.Input.Validate;
import beast.core.parameter.Parameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Logs tree annotated with metadata and/or rates")
public class MultispeciesCoalescentLogger extends BEASTObject implements Loggable {
    public Input<MultispeciesCoalescent> multispeciesCoalescentInput = new Input<>("multispeciesCoalescent", "the multispecies coalescent (for the species tree and root branch length)", Validate.REQUIRED);
    public Input<List<Function>> parameterInput = new Input<>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());
    public Input<BranchRateModel.Base> clockModelInput = new Input<>("branchratemodel", "rate to be logged with branches of the tree");
    public Input<MultispeciesPopulationModel> populationModelInput = new Input<>("populationmodel", "population sizes to be logged with branches of the tree");
    public Input<Boolean> substitutionsInput = new Input<>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<Integer> decimalPlacesInput = new Input<>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);

    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;

    private DecimalFormat df;

    @Override
    public void initAndValidate() throws Exception {
        if (parameterInput.get().size() == 0 && clockModelInput.get() == null && populationModelInput.get() == null) {
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
    public void init(PrintStream out) throws Exception {
        final MultispeciesCoalescent msc = multispeciesCoalescentInput.get();
        msc.getSpeciesTree().init(out);
    }

    @Override
    public void log(int nSample, PrintStream out) {
        final MultispeciesCoalescent msc = multispeciesCoalescentInput.get();
        final MultispeciesPopulationModel populationModel = populationModelInput.get();
        // make sure we get the current version of the inputs
        Tree tree = (Tree) msc.getSpeciesTree().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), metadata, branchRateModel, populationModel, msc));
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

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel, MultispeciesPopulationModel populationModel, MultispeciesCoalescent msc) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel, populationModel, msc));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel, populationModel, msc));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
        if (someMetaDataNeedsLogging) {
	        buf.append("[&");
	        if (metadataList.size() > 0) {
	        	for (Function metadata : metadataList) {
		            buf.append(((BEASTObject)metadata).getID());
		            buf.append('=');
		            if (metadata instanceof Parameter<?>) {
		            	Parameter p = (Parameter) metadata;
		            	int dim = p.getMinorDimension1();
		            	if (dim > 1) {
			            	buf.append('{');
			            	for (int i = 0; i < dim; i++) {
				            	buf.append(p.getMatrixValue(node.getNr(), i));
				            	if (i < dim - 1) {
					            	buf.append(',');
				            	}
			            	}
			            	buf.append('}');
		            	} else {
			            	buf.append(metadata.getArrayValue(node.getNr()));
		            	}
		            } else {
		            	buf.append(metadata.getArrayValue(node.getNr()));
		            }
		            if (metadataList.indexOf(metadata) < metadataList.size() - 1) {
		            	buf.append(",");
		            }
	        	}
	            if (branchRateModel != null || populationModel != null) {
	                buf.append(",");
	            }
	        }

	        if (branchRateModel != null) {
	            buf.append("rate=");
                appendDouble(buf, branchRateModel.getRateForBranch(node));
                if (populationModel != null) {
                    buf.append(",");
                }
	        }

            if (populationModel != null) {
                populationModel.serialize(node, buf, df);
            }

	        buf.append(']');
        }
        buf.append(":");

        double nodeLength;
        if (node.isRoot()) {
            nodeLength = node.getHeight() - msc.getRootHeight();
        } else {
            nodeLength = node.getLength();
        }

        if (substitutions) {
            appendDouble(buf, nodeLength * branchRateModel.getRateForBranch(node));
        } else {
            appendDouble(buf, nodeLength);
        }

        return buf.toString();
    }


    @Override
    public void close(PrintStream out) {
        final MultispeciesCoalescent msc = multispeciesCoalescentInput.get();
        msc.getSpeciesTree().close(out);
    }

}

    