package starbeast2;


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
public class MultispeciesCoalescentLogger extends BEASTObject implements Loggable {
    public Input<MultispeciesCoalescent> multispeciesCoalescentInput = new Input<MultispeciesCoalescent>("multispeciesCoalescent", "The multispecies coalescent calculation node.", Validate.REQUIRED);
    public Input<List<Function>> metadataInput = new Input<List<Function>>("metadata", "Metadata to be logged with the tree nodes", new ArrayList<>());
    public Input<Boolean> logRootBranchLengthInput = new Input<Boolean>("logRootBranchLength", "Log the root branch length, based on the origin height (if available) or the implicit height.", false);

    private MultispeciesCoalescent msCoalescent;
    private boolean logRootBranchLength;

    @Override
    public void initAndValidate() {
        msCoalescent = multispeciesCoalescentInput.get();
        logRootBranchLength = logRootBranchLengthInput.get();
    }

    @Override
    public void init(final PrintStream out) throws Exception {
        msCoalescent.getSpeciesTree().init(out);
    }

    @Override
    public void log(final int nSample, final PrintStream out) {
        // make sure we get the current version of the inputs
        final Tree speciesTree = (Tree) msCoalescent.getSpeciesTree().getCurrent();

        List<Function> metadataList = metadataInput.get();
        for (int i = 0; i < metadataList.size(); i++) {
        	if (metadataList.get(i) instanceof StateNode) {
        		metadataList.set(i, ((StateNode) metadataList.get(i)).getCurrent());
        	}
        }

        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
        speciesTree.getRoot().sort();
        out.print(toNewick(speciesTree.getRoot(), metadataList));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }

    String toNewick(final Node node, List<Function> metadataList) {
        final StringBuilder buf = new StringBuilder();

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

        final String populationString = msCoalescent.serializePopulation(node);
        final boolean hasMetadata = (populationString.length() > 0 || metadataList.size() > 0);
        if (hasMetadata) {
            buf.append("[&");

            if (populationString.length() > 0) {
                buf.append(populationString);
                if (metadataList.size() > 0) {
                    buf.append(",");
                }
            }

            if (metadataList.size() > 0) {
                for (Function metadata2 : metadataList) {
                    if (metadataList.indexOf(metadata2) > 0) {
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

        if (node.isRoot()) {
            if (logRootBranchLength) {
                final double rootBranchLength = msCoalescent.getRootHeight() - node.getHeight();
                buf.append(":").append(rootBranchLength);
            }
        } else {
            final double branchLength = node.getLength();
            buf.append(":").append(branchLength);
        }

        return buf.toString();
    }

    @Override
    public void close(final PrintStream out) {
        msCoalescent.getSpeciesTree().close(out);
    }
}
