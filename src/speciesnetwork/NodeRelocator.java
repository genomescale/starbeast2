package speciesnetwork;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;

import java.util.List;

/**
 * @author Chi Zhang
 */

@Description(".")
public class NodeRelocator extends Operator {

    // empty constructor to facilitate construction by XML + initAndValidate
    public NodeRelocator() {
    }

    @Override
    public void initAndValidate() {
    }

    /**
     *
     */
    @Override
    public double proposal() {

        return 0.0;
    }

}
