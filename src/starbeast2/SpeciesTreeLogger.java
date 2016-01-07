package starbeast2;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.StateNode;
import beast.core.parameter.Parameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Based on the TreeWithMetaDataLogger class, but with support for population sizes")
public class SpeciesTreeLogger extends BEASTObject implements Loggable {
    final public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree to be logged.", Validate.REQUIRED);
    final public Input<List<GeneTree>> geneTreeInput = new Input<>("geneTree", "Gene tree within the species tree.", new ArrayList<>());
    final public Input<MultispeciesPopulationModel> populationModelInput = new Input<>("populationmodel", "population sizes to be logged with branches of the tree");
    // TODO: make this input a list of valuables
    final public Input<List<Function>> parameterInput = new Input<>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());
    final public Input<BranchRateModel.Base> clockModelInput = new Input<>("branchratemodel", "rate to be logged with branches of the tree");
    final public Input<Boolean> substitutionsInput = new Input<>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    final public Input<Integer> decimalPlacesInput = new Input<>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);


    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;

    private DecimalFormat df;

    @Override
    public void initAndValidate() throws Exception {
        if (parameterInput.get().size() == 0 && clockModelInput.get() == null && populationModelInput.get() == null) {
            someMetaDataNeedsLogging = false;
            return;
            //throw new IllegalArgumentException("At least one of the metadata and branchratemodel inputs must be defined");
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
        Tree speciesTree = speciesTreeInput.get().getTree();
        speciesTree.init(out);
    }

    @Override
    public void log(int nSample, PrintStream out) {
        // make sure we get the current version of the inputs
        SpeciesTree speciesTree = speciesTreeInput.get();
        Tree tree = speciesTree.getCurrentTree();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
            if (metadata.get(i) instanceof StateNode) {
                metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
            }
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        MultispeciesPopulationModel populationModel = populationModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), metadata, branchRateModel, populationModel));
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

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel, MultispeciesPopulationModel populationModel) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel, populationModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel, populationModel));
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
                        Parameter<?> p = (Parameter<?>) metadata;
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
            nodeLength = node.getHeight() - getTreeHeight();
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

    // uses the height of the tallest species or gene tree
    private double getTreeHeight() {
        double speciesTreeHeight = speciesTreeInput.get().getTreeHeight();

        for (GeneTree gt: geneTreeInput.get()) {
            speciesTreeHeight = Double.max(speciesTreeHeight, gt.getTreeHeight());
        }

        return speciesTreeHeight;
    }

    @Override
    public void close(PrintStream out) {
        Tree speciesTree = speciesTreeInput.get().getTree();
        speciesTree.close(out);
    }

}

