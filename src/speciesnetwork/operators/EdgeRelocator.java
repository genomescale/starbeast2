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
 * @author Chi Zhang
 */

@Description("Relocate the source of an edge starting with speciation node, " +
        "or the destination of an edge ending with hybridization node.")
public class EdgeRelocator extends Operator {
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
    public EdgeRelocator() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();

        // pick an internal node randomly
        List<NetworkNode> intNodes = speciesNetwork.getInternalNodes();
        final NetworkNode snNode = intNodes.get(Randomizer.nextInt(intNodes.size()));

        final double newInterval, oldInterval;  // for calculating proposal ratio
        // start moving
        if (snNode.isReticulation()) {
            // move either the left parent or the right parent branch
            final Direction pickedDirection;
            final NetworkNode oldParent, anotherPa, oldChild;
            if (Randomizer.nextBoolean()) {
                pickedDirection = Direction.LEFT;
                oldParent = snNode.getLeftParent();
                anotherPa = snNode.getRightParent();
            } else {
                pickedDirection = Direction.RIGHT;
                oldParent = snNode.getRightParent();
                anotherPa = snNode.getLeftParent();
            }
            final double upperLimit = oldParent.getHeight();  // upper bound
            if (snNode.getLeftChild() != null)
                oldChild = snNode.getLeftChild();
            else
                oldChild = snNode.getRightChild();
            if (upperLimit > anotherPa.getHeight())
                oldInterval = anotherPa.getHeight() - oldChild.getHeight();
            else
                oldInterval = upperLimit - oldChild.getHeight();

            // look for all the candidate branches to attach to
            List<NetworkNode> candidates = new ArrayList<>();
            for (NetworkNode node : speciesNetwork.getNodes()) {
                if (node != snNode) {
                    NetworkNode lParent = node.getLeftParent();
                    NetworkNode rParent = node.getRightParent();
                    if (lParent != null && lParent != snNode && node.getHeight() < upperLimit)
                        candidates.add(node);
                    if (rParent != null && rParent != snNode && node.getHeight() < upperLimit)
                        candidates.add(node);
                }
            }
            // pick an candidate branch randomly
            final NetworkNode newNode = candidates.get(Randomizer.nextInt(candidates.size()));
            final Direction attachedDirection;
            final NetworkNode newParent;
            if (newNode.getRightParent() == null || newNode.getRightParent() == snNode) {
                attachedDirection = Direction.LEFT;
                newParent = newNode.getLeftParent();
            } else if (newNode.getLeftParent() == null || newNode.getLeftParent() == snNode) {
                attachedDirection = Direction.RIGHT;
                newParent = newNode.getRightParent();
            } else if (Randomizer.nextBoolean()) {
                attachedDirection = Direction.LEFT;
                newParent = newNode.getLeftParent();
            } else {
                attachedDirection = Direction.RIGHT;
                newParent = newNode.getRightParent();
            }

            // propose an attachment height
            final double upper;
            if (upperLimit > newParent.getHeight())
                upper = newParent.getHeight();
            else
                upper = upperLimit;
            final double lower = newNode.getHeight();
            final double newHeight = lower + (upper - lower) * Randomizer.nextDouble();
            newInterval = upper - lower;  // for calculating proposal ratio

            // deal with the left and right relationships


        } else {
            // move either the left child or the right child branch
            final Direction pickedDirection;
            final NetworkNode oldParent, anotherCh, oldChild;
            if (Randomizer.nextBoolean()) {
                pickedDirection = Direction.LEFT;
                oldChild = snNode.getLeftChild();
                anotherCh = snNode.getRightChild();
            }
            else {
                pickedDirection = Direction.RIGHT;
                oldChild = snNode.getRightChild();
                anotherCh = snNode.getLeftChild();
            }
            final double lowerLimit = oldChild.getHeight();  // lower bound
            if (snNode.getLeftParent() != null)
                oldParent = snNode.getLeftParent();
            else
                oldParent = snNode.getRightParent();
            final double oldParentHeight;
            if (oldParent != null)
                oldParentHeight = oldParent.getHeight();
            else  // oldParent is null when snNode is root, set as twice of the root height
                oldParentHeight = 2 * snNode.getHeight();
            if (lowerLimit < anotherCh.getHeight())
                oldInterval = oldParentHeight - anotherCh.getHeight();
            else
                oldInterval = oldParentHeight - lowerLimit;

            // look for all the candidate branches to attach to
            List<NetworkNode> candidates = new ArrayList<>();
            for (NetworkNode node : speciesNetwork.getNodes()) {
                if (node != snNode) {
                    NetworkNode lParent = node.getLeftParent();
                    NetworkNode rParent = node.getRightParent();
                    if (lParent != null && lParent != snNode && lParent.getHeight() > lowerLimit)
                        candidates.add(node);
                    if (rParent != null && rParent != snNode && rParent.getHeight() > lowerLimit)
                        candidates.add(node);
                    if (node.isRoot())
                        candidates.add(node);  // don't forget the root branch
                }
            }
            // pick an candidate branch randomly
            final NetworkNode newNode = candidates.get(Randomizer.nextInt(candidates.size()));
            final Direction attachedDirection;
            final NetworkNode newParent;  // newParent is null when newNode is root
            if (newNode.getRightParent() == null || newNode.getRightParent() == snNode) {
                attachedDirection = Direction.LEFT;
                newParent = newNode.getLeftParent();
            } else if (newNode.getLeftParent() == null || newNode.getLeftParent() == snNode) {
                attachedDirection = Direction.RIGHT;
                newParent = newNode.getRightParent();
            } else if (Randomizer.nextBoolean()) {
                attachedDirection = Direction.LEFT;
                newParent = newNode.getLeftParent();
            } else {
                attachedDirection = Direction.RIGHT;
                newParent = newNode.getRightParent();
            }

            // propose an attachment height



            // deal with the left and right relationships


        }


        return 0.0; //Math.log(newInterval/oldInterval);
    }
}
