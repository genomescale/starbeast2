package starbeast2.utils;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.io.PrintStream;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Logs node heights in order of height.
 */
public class NodeHeightLogger extends BEASTObject implements Loggable {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree whose node ages to log.",
            Input.Validate.REQUIRED);

    public Input<Boolean>  excludeSAInput = new Input<>(
            "excludeSANodes",
            "If true, SA node heights will be recorded as NA.",
            false);

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
    public void log(long sample, PrintStream out) {
        List<Node> nonSANodes = tree.getInternalNodes().stream()
                .filter(x -> !x.isFake())
                .sorted((n1, n2) -> {
                    if (n1.getHeight() < n2.getHeight())
                        return -1;

                    if (n1.getHeight() > n2.getHeight())
                        return 1;

                    return 0;
                })
                .collect(Collectors.toList());

        for (Node node : nonSANodes)
            out.print("\t" + node.getHeight());

        for (int i=nonSANodes.size(); i<tree.getInternalNodeCount(); i++)
            out.print("\tNA");
    }

    @Override
    public void close(PrintStream out) { }

}
