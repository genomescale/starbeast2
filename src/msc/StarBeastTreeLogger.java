package msc;


import java.io.PrintStream;
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
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
* @author Remco Bouckaert
* @author Joseph Heled
* @author Huw Ogilvie
 */

@Description("Logs tree annotated with metadata in StarBeast format")
public class StarBeastTreeLogger extends BEASTObject implements Loggable {
    public Input<Tree> treeInput = new Input<Tree>("tree", "tree to be logged", Validate.REQUIRED);
    public Input<MultispeciesPopulationModel> populationFunctionInput = new Input<MultispeciesPopulationModel>("populationModel", "The species tree population model.");
    public Input<List<Function>> metadataInput = new Input<List<Function>>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());

    private MultispeciesPopulationModel populationModel;

    String metaDataLabel;

    static final String dmv = "dmv";
    static final String dmt = "dmt";

    @Override
    public void initAndValidate() {
        metaDataLabel = "[&" + dmv + "=";
        populationModel = populationFunctionInput.get();
    }

    @Override
    public void init(final PrintStream out) throws Exception {
        treeInput.get().init(out);
    }

    @Override
    public void log(final int nSample, final PrintStream out) {
        // make sure we get the current version of the inputs
        final Tree tree = (Tree) treeInput.get().getCurrent();

        List<Function> metadataList = metadataInput.get();
        for (int i = 0; i < metadataList.size(); i++) {
        	if (metadataList.get(i) instanceof StateNode) {
        		metadataList.set(i, ((StateNode) metadataList.get(i)).getCurrent());
        	}
        }

        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), metadataList));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }


    String toNewick(final Node node, List<Function> metadataList) {
        final StringBuilder buf = new StringBuilder();
        final String dmv = populationModel.serialize(node);

        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr()+Tree.taxaTranslationOffset);
        }
        
        
        if (dmv.length() > 0 || metadataList.size() > 0) {
            buf.append("[&");
            buf.append(dmv);

            if (metadataList.size() > 0) {
            	for (Function metadata2 : metadataList) {
    	            if (metadataList.indexOf(metadata2) > 0 || buf.length() > 1) {
    	            	buf.append(",");
    	            }
    	            buf.append(((BEASTObject)metadata2).getID());
    	            buf.append('=');
    	            if (metadata2 instanceof Parameter<?>) {
    	            	Parameter p = (Parameter) metadata2;
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
    		            	buf.append(metadata2.getArrayValue(node.getNr()));
    	            	}
    	            } else {
    	            	buf.append(metadata2.getArrayValue(node.getNr()));
    	            }
            	}
            }
            buf.append(']');
        }
        if (!node.isRoot()) {
            buf.append(":").append(node.getLength());
        }
        return buf.toString();
    }

    double getMetaDataTopValue(final Node node, final Function metadataTop) {
        int nr = node.getNr();
        if (nr >= metadataTop.getDimension()) {
            nr = node.getTree().getRoot().getNr();
        }
        return metadataTop.getArrayValue(nr);
    }

    @Override
    public void close(final PrintStream out) {
        treeInput.get().close(out);
    }

}
