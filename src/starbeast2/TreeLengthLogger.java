package starbeast2;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Logger to report total length of a tree")
public class TreeLengthLogger extends CalculationNode implements Loggable, Function {
    final public Input<Tree> treeInput = new Input<>("tree", "tree to report length for.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) {
        final Tree tree = treeInput.get();
        if (getID() == null || getID().matches("\\s*")) {
            out.print(tree.getID() + ".length\t");
        } else {
            out.print(getID() + "\t");
        }
    }

    @Override
    public void log(int sample, PrintStream out) {
        final Node treeRoot = treeInput.get().getRoot();
        final double treeHeight = treeRoot.getHeight();
        final double treeLength =  recurseLength(treeRoot, treeHeight);
        out.print(treeLength + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        final Node treeRoot = treeInput.get().getRoot();
        final double treeHeight = treeRoot.getHeight();
        return recurseLength(treeRoot, treeHeight);
    }

    @Override
    public double getArrayValue(int dim) {
        final Node treeRoot = treeInput.get().getRoot();
        final double treeHeight = treeRoot.getHeight();
        return recurseLength(treeRoot, treeHeight);
    }

    private double recurseLength(final Node treeNode, final double parentHeight) {
        if (treeNode.isLeaf()) {
            return parentHeight;
        } else {
            double subtreeLength = 0.0;

            final double nodeHeight = treeNode.getHeight();
            subtreeLength += parentHeight - nodeHeight;

            final double leftChildLength = recurseLength(treeNode.getLeft(), nodeHeight);
            final double rightChildLength = recurseLength(treeNode.getRight(), nodeHeight);
            subtreeLength += leftChildLength;
            subtreeLength += rightChildLength;

            return subtreeLength;
        }
    }
}

