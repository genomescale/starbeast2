package starbeast2;

import com.google.common.collect.SetMultimap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
* @author Huw Ogilvie
 */

@Description("Implements a version of the co-ordinated species and gene tree operator described in Jones (2015)."
        + "Specifically, this operator moves a species tree node and a set of gene tree nodes related to the"
        + "species tree node by a uniform amount chosen from a range which preserves the topology of all trees."
        + "See http://dx.doi.org/10.1101/010199 for full details.")
public class CoordinatedShift extends Operator {
    public Input<MultispeciesCoalescent> mscInput = new Input<MultispeciesCoalescent>("multispeciesCoalescent", "Multispecies coalescent (gene trees contained within a species tree).");

    @Override
    public void initAndValidate() {
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final double fLogHastingsRatio = 0.0; // this move is uniform in both directions
        final MultispeciesCoalescent msc = mscInput.get();
        final Tree speciesTree = msc.getSpeciesTree();

        final int nInternalNodes = speciesTree.getInternalNodeCount();
        if (nInternalNodes == 1) { // if there are no internal nodes other than the root
            return Double.NEGATIVE_INFINITY;
        } // otherwise select a non-root internal node
        Node speciesTreeNode = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        while (speciesTreeNode.isRoot()) {
            speciesTreeNode = speciesTree.getNode(nInternalNodes + 1 + Randomizer.nextInt(nInternalNodes));
        }

        final double speciesTreeNodeHeight = speciesTreeNode.getHeight();

        final MinimumDouble tipwardFreedom = new MinimumDouble();
        final MinimumDouble rootwardFreedom = new MinimumDouble();
        final SetMultimap<Integer, Node> connectingNodes = msc.getConnectingNodes(speciesTreeNode, tipwardFreedom, rootwardFreedom);

        final double leftChildBranchLength = speciesTreeNodeHeight - speciesTreeNode.getLeft().getHeight();
        final double rightChildBranchLength = speciesTreeNodeHeight - speciesTreeNode.getRight().getHeight();
        final double speciesTreeNodeBranchLength = speciesTreeNode.getParent().getHeight() - speciesTreeNodeHeight;
        tipwardFreedom.set(leftChildBranchLength);
        tipwardFreedom.set(rightChildBranchLength);
        rootwardFreedom.set(speciesTreeNodeBranchLength);

        final double twf = tipwardFreedom.get();
        final double rwf = rootwardFreedom.get();
        final double uniformShift = (Randomizer.nextDouble() * (twf + rwf)) - twf;

        speciesTreeNode.setHeight(speciesTreeNode.getHeight() + uniformShift);
        for (Node geneTreeNode: connectingNodes.values()) {
            geneTreeNode.setHeight(geneTreeNode.getHeight() + uniformShift);
        }

        assert msc.computeCoalescentTimes(); // this move should always preserve gene-tree-within-species-tree compatibility

        return fLogHastingsRatio;
    }
}
