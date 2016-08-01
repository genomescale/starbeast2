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
import speciesnetwork.SanityChecks;

/**
 * This proposal adds a reticulation event by connecting two existing branches (with length l1 and l2) with a new branch.
 * The same branch can be picked twice (and forms a loop to that branch). The cutting proportion of each picked branch by
 * the connecting point, w1 and w2 ~ Uniform(0,1). Let l11 = l1 * w1, l12 = l1 * (1-w1), l21 = l2 * w2, l22 = l2 * (1-w2)
 * If the root branch is picked (with length l1 unknown), let l11 = w1 ~ exp(b), and l12 = unknown.
 * The direction of the new branch is determined by the two connecting points, the higher is speciation node, and the
 * lower is reticulation node. The gamma prob r = w3 ~ Uniform(0,1).
 * The Jacobian is l1 * l2. (li = 1 if the branch is a root branch. If the root branch is picked twice, the Jacobian = 1)
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob. If there is no reticulation in the network,
 * the DeleteReticulation move is aborted.
 * Let k be the number of branches in the current network. The probability of adding this branch is (1/k)(1/k)
 * Let m be the number of reticulation branches in the proposed network. The probability of selecting the same branch to
 * remove is (1/m).
 * The Hastings ratio is (1/m) / (1/k)(1/k)(f1)(f2)) = k^2 / (m * f1 * f2).
 * f1 = b * exp(-b * l11) if root, otherwise f1 = 1 (uniform density); f2 = b * exp(-b * l21) if root, otherwise f2 = 1.
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
    final static SanityChecks sc = new SanityChecks();

    // empty constructor to facilitate construction by XML + initAndValidate
    public AddReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();
        assert sc.checkNetworkSanity(speciesNetwork.getRoot()); // species network should not be insane

        // number of branches in the current network
        final int numBranches = speciesNetwork.getBranchCount();  // k

        // pick two branches randomly, including the root branch
        final int branchNr1 = Randomizer.nextInt(numBranches);
        final int branchNr2 = Randomizer.nextInt(numBranches);  // allow picking the same branch

        // get the node at the end of the branch, and the direction to the parent
        NetworkNode netNode1 = null, netNode2 = null, parentNode1 = null, parentNode2 = null;
        for (NetworkNode netNode: speciesNetwork.getAllNodesAsArray()) {
            if (netNode.isReticulation()) {
                if (netNode.getLeftBranchNumber() == branchNr1) {
                    netNode1 = netNode;
                    parentNode1 = netNode.getLeftParent();
                } else if (netNode.getRightBranchNumber() == branchNr1) {
                    netNode1 = netNode;
                    parentNode1 = netNode.getRightParent();
                }
                if (netNode.getLeftBranchNumber() == branchNr2) {
                    netNode2 = netNode;
                    parentNode2 = netNode.getLeftParent();
                } else if (netNode.getRightBranchNumber() == branchNr2) {
                    netNode2 = netNode;
                    parentNode2 = netNode.getRightParent();
                }
            } else {
                if (netNode.getBranchNumber(0) == branchNr1) {
                    netNode1 = netNode;
                    if (netNode.getLeftParent() != null) {
                        parentNode1 = netNode.getLeftParent();
                    } else {
                        parentNode1 = netNode.getRightParent();  // can be null if netNode is root
                    }
                }
                if (netNode.getBranchNumber(0) == branchNr2) {
                    netNode2 = netNode;
                    if (netNode.getLeftParent() != null) {
                        parentNode2 = netNode.getLeftParent();
                    } else {
                        parentNode2 = netNode.getRightParent();  // can be null if netNode is root
                    }
                }
            }
        }
        if (netNode1 == null || netNode2 == null )
            throw new RuntimeException("Developer ERROR: node not found!");

        // add a branch joining the two picked branches
        // propose the attaching position at each branch
        final double lambda = 1;
        final double l1, l2, l11, l21;
        double proposalRatio = 0.0;

        if (netNode1.isRoot()) {
            l1 = 1;
            l11 = Randomizer.nextExponential(lambda);
            proposalRatio += lambda * l11 - Math.log(lambda);
        } else {
            l1 = parentNode1.getHeight() - netNode1.getHeight();
            l11 = l1 * Randomizer.nextDouble();
        }
        if (netNode2.isRoot()) {
            l2 = 1;
            l21 = Randomizer.nextExponential(lambda);
            proposalRatio += lambda * l21 - Math.log(lambda);
        } else {
            l2 = parentNode2.getHeight() - netNode2.getHeight();
            l21 = l2 * Randomizer.nextDouble();
        }
        proposalRatio += Math.log(l1) + Math.log(l2);  // the Jacobian

        // create two new nodes
        NetworkNode middleNode1 = speciesNetwork.newNode();
        NetworkNode middleNode2 = speciesNetwork.newNode();
        middleNode1.setHeight(netNode1.getHeight() + l11);
        middleNode2.setHeight(netNode2.getHeight() + l21);

        // deal with the node numbers and relationships
        if (middleNode1.getHeight() > middleNode2.getHeight()) {
            // middleNode1 is bifurcation and middleNode2 is reticulation
            speciesNetwork.addSpeciationNode(middleNode1);
            speciesNetwork.addReticulationNode(middleNode2);

            if (netNode1 == netNode2 && parentNode1 == parentNode2) {
                // the two attaching points are on the same branch
                middleNode1.setLeftChild(middleNode2);
                middleNode1.setRightChild(middleNode2);
                middleNode2.setLeftParent(middleNode1);
                middleNode2.setRightParent(middleNode1);
                if (netNode1.getLeftParent() == parentNode1) {
                    middleNode1.setLeftParent(parentNode1);
                    middleNode2.setLeftChild(netNode1);
                } else {
                    middleNode1.setRightParent(parentNode1);
                    middleNode2.setRightChild(netNode1);
                }
            } else {
                // the two attaching points are on different branches
                if (netNode1.getLeftParent() == parentNode1) {
                    middleNode1.setLeftChild(netNode1);
                    middleNode1.setLeftParent(parentNode1);
                    middleNode1.setRightChild(middleNode2);
                } else {
                    middleNode1.setRightChild(netNode1);
                    middleNode1.setRightParent(parentNode1);
                    middleNode1.setLeftChild(middleNode2);
                }
                if (netNode2.getLeftParent() == parentNode2) {
                    middleNode2.setLeftChild(netNode2);
                    middleNode2.setLeftParent(parentNode2);
                    middleNode2.setRightParent(middleNode1);
                } else {
                    middleNode2.setRightChild(netNode2);
                    middleNode2.setRightParent(parentNode2);
                    middleNode2.setLeftParent(middleNode1);
                }

                // TODO: need to solve direction conflict!

            }
        } else {
            // middleNode1 is reticulation and middleNode2 is bifurcation
            speciesNetwork.addReticulationNode(middleNode1);
            speciesNetwork.addSpeciationNode(middleNode2);

            if (netNode1 == netNode2 && parentNode1 == parentNode2) {
                // the two attaching points are on the same branch
                middleNode1.setLeftParent(middleNode2);
                middleNode1.setRightParent(middleNode2);
                middleNode2.setLeftChild(middleNode1);
                middleNode2.setRightChild(middleNode1);
                if (netNode1.getLeftParent() == parentNode1) {
                    middleNode1.setLeftChild(netNode1);
                    middleNode2.setLeftParent(parentNode1);
                } else {
                    middleNode1.setRightChild(netNode1);
                    middleNode2.setRightParent(parentNode1);
                }
            } else {
                // the two attaching points are on different branches
                if (netNode1.getLeftParent() == parentNode1) {
                    middleNode1.setLeftChild(netNode1);
                    middleNode1.setLeftParent(parentNode1);
                    middleNode1.setRightParent(middleNode2);
                } else {
                    middleNode1.setRightChild(netNode1);
                    middleNode1.setRightParent(parentNode1);
                    middleNode1.setLeftParent(middleNode2);
                }
                if (netNode2.getLeftParent() == parentNode2) {
                    middleNode2.setLeftChild(netNode2);
                    middleNode2.setLeftParent(parentNode2);
                    middleNode2.setRightChild(middleNode1);
                } else {
                    middleNode2.setRightChild(netNode2);
                    middleNode2.setRightParent(parentNode2);
                    middleNode2.setLeftChild(middleNode1);
                }

                // TODO: need to solve direction conflict!

            }
        }

        // TODO: add the gamma prob.
        final double gamma = Randomizer.nextDouble();


        // number of reticulation branches in the proposed network
        final int numReticulationBranches = 2 * speciesNetwork.getReticulationNodeCount();  // m

        proposalRatio += 2 * Math.log(numBranches) - Math.log(numReticulationBranches);

        return proposalRatio;
    }
}
