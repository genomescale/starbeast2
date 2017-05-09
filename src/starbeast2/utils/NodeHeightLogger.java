package starbeast2.utils;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.io.PrintStream;
import java.util.Comparator;
import java.util.List;

/**
 * Logs node heights in order of height.
 */
public class NodeHeightLogger extends BEASTObject implements Loggable {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree whose node ages to log.",
            Input.Validate.REQUIRED);

    Tree tree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        String prefix;

        if (getID() != null)
            prefix = getID() + ".";
        else
            prefix = tree.getID() + ".nodeHeight.";

        for (int i=0; i<tree.getInternalNodeCount(); i++)
            out.print("\t" + prefix + i);
    }

    @Override
    public void log(int sample, PrintStream out) {
        List<Node> sortedNodes = tree.getInternalNodes();
        sortedNodes.sort((o1, o2) -> {
            if (o1.getHeight() < o2.getHeight())
                return -1;

            if (o1.getHeight() > o2.getHeight())
                return 1;

            return 0;
        });

        for (Node node : sortedNodes)
            out.print("\t" + node.getHeight());
    }

    @Override
    public void close(PrintStream out) { }

}
