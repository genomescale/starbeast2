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
        Operator treeOp = treeOperatorInput.get();

        /* IntegerParameter embedding = reembedOp.embeddingInput.get();
        Network speciesNetwork = reembedOp.speciesNetworkInput.get();
        Tree geneTree = reembedOp.geneTreeInput.get();
        // print matrix for debugging
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < embedding.getMinorDimension2(); i++) {
            for (int j = 0; j < embedding.getMinorDimension1(); j++) {
                sb.append(embedding.getMatrixValue(i, j));
                sb.append("\t");
            }
            sb.append("\n");
        }
        sb.append(geneTree.getRoot().toNewick());  sb.append("\n");
        sb.append(geneTree.getRoot().toString());  sb.append("\n");
        sb.append(speciesNetwork.toString());  //sb.append("\n");
        System.out.println(sb); */

        // count the number of alternative traversing choices for the current state (n0)
        final int oldChoices = reembedOp.getNumberOfChoices();
        if (oldChoices < 0)
            throw new RuntimeException("Developer ERROR: current embedding invalid!");

        // first make the operation
        double logHR = treeOp.proposal();
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // Update calculation nodes as subsequent operators may depend on state nodes made dirty by this operation.
        if (!treeOp.listStateNodes().isEmpty())  // copied from JointOperator
            treeOp.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();

        // then rebuild the embedding
        final int newChoices = reembedOp.initializeEmbedding();
        // Update calculation nodes as subsequent operators may depend on state nodes made dirty by this operation.
        if (!reembedOp.listStateNodes().isEmpty()) // copied from JointOperator
            reembedOp.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();

        if (newChoices < 0) return Double.NEGATIVE_INFINITY;

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
