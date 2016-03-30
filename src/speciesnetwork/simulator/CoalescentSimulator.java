package speciesnetwork.simulator;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.Operator;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import speciesnetwork.*;

/**
 * @author Chi Zhang
 */

@Description("Simulate gene trees given a species network (multispecies coalescent).")
public class CoalescentSimulator extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network for simulating the gene tree.", Validate.REQUIRED);
    public Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "Gene trees to initialize/simulate.", new ArrayList<>());
    public Input<PopulationSizeModel> populationInput =
            new Input<>("populationModel", "The species network population size model.", Validate.REQUIRED);

    private Network speciesNetwork;

    public void initStateNodes() {
        speciesNetwork = speciesNetworkInput.get();

        // population sizes
        final int nSpeciesBranches = speciesNetwork.getBranchCount();
        PopulationSizeModel populationModel = populationInput.get();
        populationModel.initPopSizes(nSpeciesBranches);

    }

    public void simulate() {
        NetworkNode networkRoot = speciesNetwork.getRoot();

        final List<Tree> geneTrees = geneTreesInput.get();
        for (Tree gtree : geneTrees) {
            // reset visited indicator
            networkRoot.recursiveResetVisited();
            simulate(networkRoot);
        }
    }

    // recursively simulate lineages coalescent in each population
    private void simulate(NetworkNode snNode) {
        if (snNode.isVisited())
            return;
        if (snNode.getLeftChild() != null) {
            simulate(snNode.getLeftChild());
        }
        if (snNode.getRightChild() != null) {
            simulate(snNode.getRightChild());
        }

        snNode.setVisited();
    }

}
