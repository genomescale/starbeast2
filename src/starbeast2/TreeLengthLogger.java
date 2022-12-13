package starbeast2;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

@Description("Logger to report total length of a tree")
public class TreeLengthLogger extends CalculationNode implements Loggable, Function {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree to report length for.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) {
        final TreeInterface tree = treeInput.get();
        if (getID() == null || getID().matches("\\s*")) {
            out.print(tree.getID() + ".length\t");
        } else {
            out.print(getID() + "\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        final double treeLength = TreeStats.getLength(treeInput.get());
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
        return TreeStats.getLength(treeInput.get());
    }

    @Override
    public double getArrayValue(int dim) {
        return TreeStats.getLength(treeInput.get());
    }
}

