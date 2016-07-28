package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

/**
 * This proposal adds a reticulation event by connecting two existing branches (with length l1 and l2) with a new branch.
 * The cutting proportion of each picked branch by the connecting point, w1 and w2 ~ Uniform(0,1). The direction of the
 * new branch is determined by the two connecting points, the higher is speciation node, and the lower is reticulation node.
 * The gamma prob r = w3 ~ Uniform(0,1).
 * l11 = l1 * w1, l12 = l1 * (1-w1), l21 = l2 * w2, l22 = l2 * (1-w2), r = w3.
 * The Jacobian is l1 * l2.
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob. If there is no reticulation in the network,
 * the DeleteReticulation move is aborted.
 * Let k be the number of branches in the current network. The probability of adding this branch is (1/k)(1/(k-1))
 * Let m be the number of reticulation branches in the proposed network. The probability of selecting the same branch to
 * remove is (1/m).
 * The Hastings ratio is (1/m) / (1/k)(1/(k-1)) = k * (k-1) / m.
 *
 * @author Chi Zhang
 */

@Description("Relocate the source of an edge starting with speciation node, " +
        "or the destination of an edge ending with hybridization node.")
public class AddReticulation extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    public Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "list of gene trees embedded in species network", new ArrayList<>());
    public Input<List<IntegerParameter>> embeddingsInput =
            new Input<>("embedding", "The matrices to embed the gene trees in the species network.", new ArrayList<>());
    public Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonSuperset", "Super-set of taxon sets mapping lineages to species.", Input.Validate.REQUIRED);

    private enum Direction {LEFT, RIGHT}

    // empty constructor to facilitate construction by XML + initAndValidate
    public AddReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();

        // pick two branches randomly, including the root branch
        final int numBranches = speciesNetwork.getBranchCount();  // k
        int branchNr1 = Randomizer.nextInt(numBranches);
        int branchNr2 = Randomizer.nextInt(numBranches);
        while (branchNr2 == branchNr1)
            branchNr2 = Randomizer.nextInt(numBranches);  // get a different branch

        // get the node at the end of the branch, if it is a reticulation node, also get the direction
        NetworkNode netNode1, netNode2, parentNode1, parentNode2;
        Direction direction1, direction2;
        for (NetworkNode netNode: speciesNetwork.getAllNodesAsArray()) {
            if (netNode.isReticulation()) {
                if (netNode.getLeftBranchNumber() == branchNr1) {
                    netNode1 = netNode;
                    direction1 = Direction.LEFT;
                    parentNode1 = netNode.getLeftParent();
                } else if (netNode.getRightBranchNumber() == branchNr1) {
                    netNode1 = netNode;
                    direction1 = Direction.RIGHT;
                    parentNode1 = netNode.getRightParent();
                } else if (netNode.getLeftBranchNumber() == branchNr2) {
                    netNode2 = netNode;
                    direction2 = Direction.LEFT;
                    parentNode2 = netNode.getLeftParent();
                } else if (netNode.getRightBranchNumber() == branchNr2) {
                    netNode2 = netNode;
                    direction2 = Direction.RIGHT;
                    parentNode2 = netNode.getRightParent();
                }
            } else {
                if (netNode.getBranchNumber(0) == branchNr1 && netNode.getLeftParent() != null) {
                    netNode1 = netNode;
                    direction1 = Direction.LEFT;
                    parentNode1 = netNode.getLeftParent();
                } else if (netNode.getBranchNumber(0) == branchNr1 && netNode.getRightParent() != null) {
                    netNode1 = netNode;
                    direction1 = Direction.RIGHT;
                    parentNode1 = netNode.getRightParent();
                } else if (netNode.getBranchNumber(0) == branchNr2 && netNode.getLeftParent() != null) {
                    netNode2 = netNode;
                    direction2 = Direction.LEFT;
                    parentNode2 = netNode.getLeftParent();
                } else if (netNode.getBranchNumber(0) == branchNr2 && netNode.getRightParent() != null) {
                    netNode2 = netNode;
                    direction2 = Direction.RIGHT;
                    parentNode2 = netNode.getRightParent();
                }
            }
        }

        // add a branch joining the two picked branches
        // first create two new nodes


        return 0.0;
    }
}
