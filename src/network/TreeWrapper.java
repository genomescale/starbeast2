package network;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public abstract class TreeWrapper extends CalculationNode {
    public Input<Tree> treeInput = new Input<>("tree", "Tree object for this wrapper.", Validate.REQUIRED);

    protected Node getRoot() {
        return treeInput.get().getRoot();
    }

    protected double getTreeHeight() {
        return treeInput.get().getRoot().getHeight();
    }

    protected Tree getTree() {
        return treeInput.get();
    }

    protected Tree getCurrentTree() {
        return (Tree) treeInput.get().getCurrent();
    }

    public String toString() {
        final String treeString = treeInput.get().getRoot().toNewick(); 
        return treeString;
    }
}
