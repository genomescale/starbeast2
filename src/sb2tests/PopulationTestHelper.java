package sb2tests;

import java.util.ArrayList;
import java.util.List;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import starbeast2.GeneTreeWithinSpeciesTree;
import starbeast2.MultispeciesPopulationModel;

public class PopulationTestHelper {
    final int nSpecies = 4;
    final int nBranches = (nSpecies * 2) - 1;
    final double[] popSizes = new double[nBranches];

    final int individualsPerSpecies = 2;
    GeneTreeWithinSpeciesTree geneTreeTestA;
    GeneTreeWithinSpeciesTree geneTreeTestB;
    TaxonSet speciesSuperSet;
    final double allowedError = 10e-6;

    final TreeParser speciesTree = new TreeParser();
    final TreeParser geneTreeA = new TreeParser();
    final TreeParser geneTreeB = new TreeParser();
    
    final List<GeneTreeWithinSpeciesTree> geneTreeList = new ArrayList<GeneTreeWithinSpeciesTree>();
    
    MultispeciesPopulationModel populationModel;

    final double ploidy = 2.0;

    public void initializeTrees(final String newickSpeciesTree, final String newickGeneTreeA, final String newickGeneTreeB) throws Exception {
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
        geneTreeA.initByName("newick", newickGeneTreeA, "IsLabelledNewick", true);
        geneTreeB.initByName("newick", newickGeneTreeB, "IsLabelledNewick", true);

        List<Taxon> superSetList = new ArrayList<>();
        for (int i = 0; i < nSpecies; i++) {
            final String speciesName = String.format("s%d", i);
            List<Taxon> taxonList = new ArrayList<>();
            for (int j = 0; j < individualsPerSpecies; j++) {
                final String taxonName = String.format("s%d_tip%d", i, j);
                taxonList.add(new Taxon(taxonName));
            }
            superSetList.add(new TaxonSet(speciesName, taxonList));
        }
        speciesSuperSet = new TaxonSet(superSetList);

        geneTreeTestA = new GeneTreeWithinSpeciesTree();
        geneTreeTestB = new GeneTreeWithinSpeciesTree();

        geneTreeTestA.initByName("tree", geneTreeA, "ploidy", ploidy);
        geneTreeTestB.initByName("tree", geneTreeB, "ploidy", ploidy);
        
        geneTreeList.add(geneTreeTestA);
        geneTreeList.add(geneTreeTestB);
    }
}
