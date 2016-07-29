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

    // empty constructor to facilitate construction by XML + initAndValidate
    public DeleteReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();

        // pick a reticulation branch randomly
        List<NetworkNode> hybridNodes = speciesNetwork.getReticulationNodes();
        final NetworkNode netNode = hybridNodes.get(Randomizer.nextInt(hybridNodes.size()));
        Direction direction;
        NetworkNode topNode;
        if (Randomizer.nextBoolean()) {
            direction = Direction.LEFT;
            topNode = netNode.getLeftParent();
        } else {
            direction = Direction.RIGHT;
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
        NetworkNode childNode2, parentNode2;
        if (topNode.getLeftParent() != null)
            parentNode2 = topNode.getLeftParent();
        else
            parentNode2 = topNode.getRightParent();  // can be null if topNode is root
        if (topNode.getLeftChild() == netNode)
            childNode2 = netNode.getRightChild();
        else
            childNode2 = netNode.getLeftChild();

        // delete the reticulation branch, connect childNode1 and parentNode1, and connect childNode2 and parentNode2




        return 0.0;
    }
}
