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
    public Input<PopulationSizeModel> populationInput =
            new Input<>("populationModel", "The species network population size model.", Validate.REQUIRED);
    public Input<Tree> geneTreeInput =
            new Input<>("geneTree", "The gene tree to be simulated.", Validate.REQUIRED);
    public Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonSuperset", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);
    public Input<IntegerParameter> embeddingInput =
            new Input<>("embedding", "The matrix to embed the gene tree within the species network.", Validate.REQUIRED);

    private Network speciesNetwork;

    @Override
    public void initAndValidate() {
        speciesNetwork = speciesNetworkInput.get();
        final int nSpeciesBranches = speciesNetwork.getBranchCount();

        // PopulationSizeModel populationModel = populationInput.get();
        // populationModel.initPopSizes(nSpeciesBranches);

    }

    @Override
    public double proposal() {
        NetworkNode networkRoot = speciesNetwork.getRoot();

        final Tree geneTree = geneTreeInput.get();

        // reset visited indicator
        networkRoot.recursiveResetVisited();
        simulate(networkRoot);



        return 0.0;
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
