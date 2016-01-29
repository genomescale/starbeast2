package sb2tests;


import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import network.ConstantPopulation;
import network.PopulationSizeModel;

public class IncompatibleTreeTest extends PopulationTestHelper {
    private final int individualsPerSpecies = 2;
    private RealParameter popSizesParameter;

    public IncompatibleTreeTest() throws Exception {
        popSize = 0.3;
        ploidy = 2.0;
        nSpecies = 8;
        expectedLogP = Double.NEGATIVE_INFINITY;

        newickSpeciesTree = "((T4:5.982015323363934,((T1:1.9435075796666423,T2:1.9435075796666423):2.031076829347149,T3:3.9745844090137914):2.007430914350143):1.9848867018863912,(((T5:0.021716110420807247,T7:0.021716110420807247):0.005999952952004395,T6:0.027716063372811642):0.0043842395682471905,T8:0.03210030294105883):7.934801722309267)";
        newickGeneTrees.add("((((((T1_1:0.579713480872118,(T2_1:0.05611562075920323,T2_2:0.05611562075920323):0.5235978601129148):0.08787237243508017,T1_2:0.6675858533071982):0.09003817632305622,T3_2:0.7576240296302544):0.7947472706863555,T3_1:1.55237130031661):0.9175463345004531,((T4_1:0.07822345639893319,T4_2:0.07822345639893319):1.9276398321493908,T8_1:2.005863288548324):0.464054346268739):0.7187520201759066,((((T5_1:0.04820585227227065,T5_2:0.04820585227227065):0.44158170738649616,T8_2:0.4897875596587668):0.03103744358250654,((T6_1:0.03361079557684876,T7_1:0.03361079557684876):0.006140006250793695,T7_2:0.03975080182764246):0.4810742014136309):1.0088148036996385,T6_2:1.5296398069409118):1.6590298480520578)");
        newickGeneTrees.add("(((T1_1:1.1016256770243424,(T2_1:0.19828874072367977,T3_2:0.19828874072367977):0.9033369363006627):1.2800224449580238,((T2_2:1.1639549747191034,T3_1:1.1639549747191034):0.6916011764687042,((((T5_1:0.10938777466296375,T7_2:0.10938777466296375):0.4116904771092121,(T6_1:0.48594284383657343,T8_2:0.48594284383657343):0.03513540793560238):0.6520597260635511,(T5_2:0.14642483340635992,(T6_2:0.1325524781224538,T7_1:0.1325524781224538):0.013872355283906124):1.0267131444293671):0.12025146840841128,T8_1:1.2933894462441382):0.5621667049436694):0.5260919707945586):1.2502348026420713,(T1_2:2.4985948064682164,(T4_1:0.023970638819853857,T4_2:0.023970638819853857):2.4746241676483627):1.133288118156221)");
    }

    @Test
    public void testLogP() throws Exception {
        super.testLogP();
    }
    

    @Override
    public TaxonSet generateSuperset() throws Exception {
        List<Taxon> superSetList = new ArrayList<>();
        for (int i = 0; i < nSpecies; i++) {
            final String speciesName = String.format("T%d", i + 1);
            List<Taxon> taxonList = new ArrayList<>();
            for (int j = 0; j < individualsPerSpecies; j++) {
                final String taxonName = String.format("T%d_%d", i + 1, j + 1);
                taxonList.add(new Taxon(taxonName));
            }
            superSetList.add(new TaxonSet(speciesName, taxonList));
        }

        TaxonSet speciesSuperset = new TaxonSet(superSetList);
        
        return speciesSuperset;
    }

    @Override
    public PopulationSizeModel generatePopulationModel() throws Exception {
        popSizesParameter = new RealParameter();
        popSizesParameter.initByName("value", String.valueOf(popSize));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", popSizesParameter);
        state.initialise();

        PopulationSizeModel populationModel = new ConstantPopulation();
        populationModel.initByName("popSizes", popSizesParameter);
        
        return populationModel;
    }
}
