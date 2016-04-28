package speciesnetwork.operators;

import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.StateNode;

/**
 * @author Chi Zhang
 */

@Description("Combine a gene tree operator with RebuildEmbedding.")
public class JointReembedding extends Operator {
    public Input<RebuildEmbedding> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene tree within species network.", Validate.REQUIRED);
    public Input<Operator> treeOperatorInput = new Input<>("operator",
            "Tree operator to combine into RebuildEmbedding.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        RebuildEmbedding reembedOp = rebuildEmbeddingInput.get();

        // count the number of alternative traversing choices for the current state (n0)
        final int oldChoices = reembedOp.getNumberOfChoices();
        if (oldChoices < 0)
            throw new RuntimeException("Developer ERROR: current embedding invalid!");

        // first make the operation
        Operator treeOp = treeOperatorInput.get();
        double logHR = treeOp.proposal();
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // Update calculation nodes as subsequent operators may depend on state nodes made dirty by this operation.
        if (!treeOp.listStateNodes().isEmpty())  // copied from JointOperator
            treeOp.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();

        // then rebuild the embedding
        if (!reembedOp.initializeEmbedding())
            return Double.NEGATIVE_INFINITY;

        // count the number of alternative traversing choices for the new state (n1)
        final int newChoices = reembedOp.getNumberOfChoices();
        if (newChoices < 0)
            throw new RuntimeException("Developer ERROR: new embedding invalid!");

        return logHR + (newChoices - oldChoices) * Math.log(2);
    }

    @Override
    public List<StateNode> listStateNodes() {
        List<StateNode> stateNodeList = new ArrayList<>();

        stateNodeList.addAll(treeOperatorInput.get().listStateNodes());
        stateNodeList.addAll(rebuildEmbeddingInput.get().listStateNodes());

        return stateNodeList;
    }
}
