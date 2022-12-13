package sb2tests;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class MissingDataCoordinatedExchangeTest extends ExchangeTestHelper {
    public MissingDataCoordinatedExchangeTest() {
        bTipLabel = "B";
        bIsParent = false;

        newickSpeciesTree = "((A:1.0,B:1.0):2.0,C:3.0)";
        newickGeneTrees.add("((A1:1.5,B1:1.5):1.0,(A2:2.0,B2:2.0):0.5)");

        ploidy = 2.0;
        popSize = 0.3;
        expectedLogHR = Double.NEGATIVE_INFINITY;
    }

    public TaxonSet generateSuperset() throws Exception {
        List<Taxon> superSetList = new ArrayList<>();

        List<Taxon> taxonListA = new ArrayList<>();
        List<Taxon> taxonListB = new ArrayList<>();
        List<Taxon> taxonListC = new ArrayList<>();

        taxonListA.add(new Taxon("A1"));
        taxonListA.add(new Taxon("A2"));
        taxonListB.add(new Taxon("B1"));
        taxonListB.add(new Taxon("B2"));
        taxonListC.add(new Taxon("C1"));

        superSetList.add(new TaxonSet("A", taxonListA));
        superSetList.add(new TaxonSet("B", taxonListB));
        superSetList.add(new TaxonSet("C", taxonListC));


        TaxonSet speciesSuperSet = new TaxonSet(superSetList);

        return speciesSuperSet;
    }

    @Test
    public void testLogHR() throws Exception {
        super.testLogHR();
    }
}
