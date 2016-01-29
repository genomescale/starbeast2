package speciesnetwork;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Multimap;

import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public abstract class CoordinatedOperator extends Operator {
    public Input<TreeInterface> treeInput = new Input<>("tree", "The species tree state node.", Validate.REQUIRED);
    public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree calculation node.", Validate.REQUIRED);
    public Input<List<TreeInterface>> geneTreeInput = new Input<>("geneTree", "Gene tree within the species tree.", new ArrayList<>());

    protected int nGeneTrees;

    @Override
    public void initAndValidate() {
        nGeneTrees = geneTreeInput.get().size();
    }

    protected Set<String> findDescendants(Node speciesTreeNode, int speciesTreeNodeNumber) {
        final Multimap<Integer, String> numberTipMap = speciesTreeInput.get().getNumberTipMap();
        final Set<String> descendantNames = new HashSet<>();

        if (speciesTreeNode.isLeaf()) {
            descendantNames.addAll(numberTipMap.get(speciesTreeNodeNumber));
        } else {
            final Node leftChild = speciesTreeNode.getLeft();
            final Node rightChild = speciesTreeNode.getRight();
            final int leftChildNumber = leftChild.getNr();
            final int rightChildNumber = rightChild.getNr();

            descendantNames.addAll(findDescendants(leftChild, leftChildNumber));
            descendantNames.addAll(findDescendants(rightChild, rightChildNumber));
        }

        return descendantNames;
    }
}
