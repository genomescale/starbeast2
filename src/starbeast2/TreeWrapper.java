package starbeast2;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public abstract class TreeWrapper extends CalculationNode {
    public Input<TreeInterface> treeInput = new Input<>("tree", "Tree object for this wrapper.");

    protected Node getRoot() {
        return treeInput.get().getRoot();
    }

    protected double getTreeHeight() {
        return treeInput.get().getRoot().getHeight();
    }

    protected TreeInterface getTree() {
        return treeInput.get();
    }

    public String toString() {
        final String treeString = treeInput.get().getRoot().toNewick(); 
        return treeString;
    }
}
