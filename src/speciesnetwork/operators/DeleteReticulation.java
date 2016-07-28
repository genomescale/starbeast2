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
 * The gamma prob r is removed. See also AddReticulation.
 * The Jacobian is 1/(l1*l2).
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob.
 * Let m' be the number of reticulation branches in the current network. The probability of selecting the this branch to
 * remove is (1/m').
 * Let k' be the number of branches in the proposed network. The probability of adding this branch is (1/k')(1/(k'-1))

 * The Hastings ratio is (1/k)(1/(k-1)) / (1/m) = m / (k*(k-1)).
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



        return 0.0;
    }
}
