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
 * This proposal delete a reticulation branch from the species network. If there is no reticulation, this is aborted.
 * The two branches at each connecting point are joined, resulting branches with length l1 and l2 respectively.
 * If the deleted node is root, its child becomes the new root. The gamma prob r is removed.
 * The Jacobian is 1/(l1*l2). (li = 1 if the deleted node is root. If the root branch is picked twice, the Jacobian = 1)
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob.
 * Let m' be the number of reticulation branches in the current network. The probability of selecting the this branch to
 * remove is (1/m').
 * Let k' be the number of branches in the proposed network. The probability of adding this branch is (1/k')(1/k')

 * The Hastings ratio is (1/k')(1/k')(f1)(f2) / (1/m') = m' * f1 * f2 / k^2.
 * f1 = b * exp(-b * l11) if root, otherwise f1 = 1 (uniform density); f2 = b * exp(-b * l21) if root, otherwise f2 = 1.
 * See also AddReticulation.
 *
 * @author Chi Zhang
 */

@Description("Relocate the source of an edge starting with speciation node, " +
        "or the destination of an edge ending with hybridization node.")
public class DeleteReticulation extends Operator {
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
    public DeleteReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();
        assert sc.checkNetworkSanity(speciesNetwork.getRoot()); // species network should not be insane

        List<NetworkNode> hybridNodes = speciesNetwork.getReticulationNodes();
        final int numHybridNodes = speciesNetwork.getReticulationNodeCount();

        //number of reticulation branches in the current network
        final int numReticulationBranches = 2 * numHybridNodes;  // m'

        // pick a reticulation branch randomly
        final NetworkNode netNode = hybridNodes.get(Randomizer.nextInt(numHybridNodes));
        final NetworkNode topNode;
        if (Randomizer.nextBoolean()) {
            topNode = netNode.getLeftParent();
        } else {
            topNode = netNode.getRightParent();
        }

        // get the child node and another parent node of netNode
        NetworkNode childNode1, parentNode1;
        if (netNode.getLeftChild() != null)
            childNode1 = netNode.getLeftChild();
        else
            childNode1 = netNode.getRightChild();
        if (netNode.getLeftParent() == topNode)
            parentNode1 = netNode.getRightParent();
        else
            parentNode1 = netNode.getLeftParent();

        // get the parent node and another child node of topNode
        NetworkNode childNode2, parentNode2;  // parentNode2 can be null if topNode is root
        if (topNode.getLeftParent() != null)
            parentNode2 = topNode.getLeftParent();
        else
            parentNode2 = topNode.getRightParent();
        if (topNode.getLeftChild() == netNode)
            childNode2 = netNode.getRightChild();
        else
            childNode2 = netNode.getLeftChild();

        final double lambda = 1;
        final double l1, l2, l11, l21;
        double proposalRatio = 0.0;

        // delete the reticulation branch, connect childNode1 and parentNode1, and connect childNode2 and parentNode2
        if (netNode == childNode2 && topNode == parentNode1) {
            // the two attaching points are on the same branch
            if (topNode.isRoot()) {
                l1 = l2 = 1;
                l11 = netNode.getHeight() - childNode1.getHeight();
                l21 = topNode.getHeight() - childNode1.getHeight();
                proposalRatio += 2 * Math.log(lambda) - lambda * (l11 + l21);
            } else {
                l1 = l2 = parentNode2.getHeight() - childNode1.getHeight();
            }

            if (childNode1.getLeftParent() == netNode) {
                childNode1.setLeftParent(parentNode2);
            } else {
                childNode1.setRightParent(parentNode2);
            }

            // TODO: need to solve direction conflict!
        } else {
            // the two attaching points are on different branches
            l1 = parentNode1.getHeight() - childNode1.getHeight();
            if (topNode.isRoot()) {
                l2 = 1;
                l21 = topNode.getHeight() - childNode2.getHeight();
                proposalRatio += Math.log(lambda) - lambda * l21;
            } else {
                l2 = parentNode2.getHeight() - childNode2.getHeight();
            }

            if (childNode1.getLeftParent() == netNode) {
                childNode1.setLeftParent(parentNode1);
            } else {
                childNode1.setRightParent(parentNode1);
            }

            if (childNode2.getLeftParent() == topNode) {
                childNode2.setLeftParent(parentNode2);
            } else {
                childNode2.setRightParent(parentNode2);
            }


            // TODO: need to solve direction conflict!
        }
        proposalRatio += - Math.log(l1) - Math.log(l2);  // the Jacobian

        // delete the two intermediate nodes
        speciesNetwork.deleteSpeciationNode(topNode);
        speciesNetwork.deleteReticulationNode(netNode);

        // TODO: add the gamma prob.


        // number of branches in the proposed network
        final int numBranches = speciesNetwork.getBranchCount();  // k'

        proposalRatio += Math.log(numReticulationBranches) - 2 * Math.log(numBranches);

        return proposalRatio;
    }
}
