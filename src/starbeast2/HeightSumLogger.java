package starbeast2;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.core.Input.Validate;

@Description("Logger to report height of a tree")
public class HeightSumLogger extends CalculationNode implements Loggable, Function {
    public Input<Tree> treeInput = new Input<Tree>("tree", "Tree to report height for.", Validate.REQUIRED);
    public Input<Boolean> boolInput = new Input<Boolean>("logXform", "Log-transform node heights before summation.", false);

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) throws Exception {
        final Tree tree = treeInput.get();
        if (getID() == null || getID().matches("\\s*")) {
            out.print(tree.getID() + ".heightsum\t");
        } else {
            out.print(getID() + "\t");
        }
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(getArrayValue() + "\t");
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
        final Tree currentTree = (Tree) treeInput.get().getCurrent();
        final boolean logXform = boolInput.get();
        double heightSum = 0.0;
        for (Node treeNode: currentTree.getInternalNodes()) {
            if (logXform) {
                heightSum += Math.log(treeNode.getHeight());
            } else {
                heightSum += treeNode.getHeight();
            }
        }

        return heightSum;
    }

    @Override
    public double getArrayValue(int iDim) {
        return getArrayValue();
    }
}
