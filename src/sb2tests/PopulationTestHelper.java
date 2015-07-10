package sb2tests;

import java.util.ArrayList;
import java.util.List;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import starbeast2.GeneTreeWithinSpeciesTree;
import starbeast2.MultispeciesPopulationModel;

public class PopulationTestHelper {
    final String newickSpeciesTree = "((s0:0.32057156677143211,s3:0.32057156677143211):1.2653250035015629,(s1:0.56540722294658641,s2:0.56540722294658641):1.0204893473264085)";
    final String newickGeneTreeA = "((((s0_tip1:0.3416660303037105,s3_tip0:0.3416660303037105):0.024561190897159135,s0_tip0:0.36622722120086965):0.0643095990846464,s3_tip1:0.43053682028551604):1.4201019862262891,((s1_tip0:0.14473698225381706,s1_tip1:0.14473698225381706):0.5135479407233198,(s2_tip0:0.19897724687831703,s2_tip1:0.19897724687831703):0.4593076760988198):1.1923538835346683)";
    final String newickGeneTreeB = "(((s0_tip0:0.04173231934154758,s0_tip1:0.04173231934154758):0.7845256741090114,(s3_tip0:0.09482581277282173,s3_tip1:0.09482581277282173):0.7314321806777372):0.8703651781500925,((s1_tip0:0.33170960882423645,s1_tip1:0.33170960882423645):0.29497523293318856,(s2_tip0:0.2908611340994834,s2_tip1:0.2908611340994834):0.3358237076579416):1.0699383298432266)";

    final int n_species = 4;
    final int individuals_per_species = 2;
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

    public void initializeSpeciesTree() throws Exception {
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
    }

    public void initializeGeneTrees() throws Exception {
        speciesTree.initByName("newick", newickSpeciesTree, "IsLabelledNewick", true);
        geneTreeA.initByName("newick", newickGeneTreeA, "IsLabelledNewick", true);
        geneTreeB.initByName("newick", newickGeneTreeB, "IsLabelledNewick", true);

        List<Taxon> superSetList = new ArrayList<>();
        for (int i = 0; i < n_species; i++) {
            final String speciesName = String.format("s%d", i);
            List<Taxon> taxonList = new ArrayList<>();
            for (int j = 0; j < individuals_per_species; j++) {
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
