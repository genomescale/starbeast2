package starbeast2;

import java.util.List;

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
        + "specifically moves a species tree node, and equally all gene tree nodes within that species tree branch,"
        + "by a uniform amount chosen from a range which preserves the topology of all trees."
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
        final int nInternalNodes = speciesTree.getInternalNodeCount() - 1; // does not include root node
        final int speciesTreeNodeNumber = nInternalNodes + 2 + Randomizer.nextInt(nInternalNodes);
        final Node speciesTreeNode = speciesTree.getNode(speciesTreeNodeNumber);

        assert !speciesTreeNode.isRoot();
        assert !speciesTreeNode.isLeaf();

        final MinimumDouble tipwardFreedom = new MinimumDouble();
        final MinimumDouble rootwardFreedom = new MinimumDouble();
        final List<Node> parentalNodes = msc.getInternalBranchNodes(speciesTreeNodeNumber, tipwardFreedom, rootwardFreedom);

        final double twf = tipwardFreedom.get();
        final double rwf = rootwardFreedom.get();
        final double uniformShift = (Randomizer.nextDouble() * (twf + rwf)) - twf;

        speciesTreeNode.setHeight(speciesTreeNode.getHeight() + uniformShift);
        for (Node geneTreeNode: parentalNodes) {
            geneTreeNode.setHeight(geneTreeNode.getHeight() + uniformShift);
        }

        assert msc.computeCoalescentTimes(); // this move should always preserve gene-tree-within-species-tree compatibility

        return fLogHastingsRatio;
    }
}
