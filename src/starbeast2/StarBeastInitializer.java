package starbeast2;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.google.common.collect.Multimap;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.alignment.distance.JukesCantorDistance;
import beast.evolution.tree.Node;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.math.distributions.MRCAPrior;
import beast.util.ClusterTree;
import beast.util.TreeParser;

/**
* @author Joseph Heled
* @author Huw Ogilvie
 */

@Description("Set a starting point for a *BEAST analysis from gene alignment data.")
public class StarBeastInitializer extends Tree implements StateNodeInitialiser {

    static enum Method {
        POINT("point-estimate"),
        ALL_RANDOM("random");

        Method(final String name) {
            this.ename = name;
        }

        @Override
		public String toString() {
            return ename;
        }

        private final String ename;
    }
    final public Input<Method> initMethod = new Input<>("method", "Initialise either with a totally random " +
            "state or a point estimate based on alignments data (default point-estimate)",
            Method.POINT, Method.values());

    final public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree to initialize.", Input.Validate.REQUIRED);
    final public Input<String> newickInput = new Input<>("newick", "Newick string for a custom initial species tree.");

    final public Input<List<Tree>> genes = new Input<>("geneTree", "Gene trees to initialize", new ArrayList<>());

    final public Input<RealParameter> birthRate = new Input<>("birthRate",
            "Tree prior birth rate to initialize");

    final public Input<Function> muInput = new Input<>("baseRate",
            "Main clock rate used to scale trees (default 1).");

    final public Input<PopulationModel> populationFunctionInput = new Input<>("populationModel", "The species tree population model.");

    @Override
    public void initStateNodes() {
        final SpeciesTree speciesTree = speciesTreeInput.get();
        final Set<BEASTInterface> treeOutputs = speciesTreeInput.get().getOutputs();
        final Method method = initMethod.get();
        final String newick = newickInput.get();

        final List<MRCAPrior> calibrations = new ArrayList<>();
        for (final Object plugin : treeOutputs ) {
            if (plugin instanceof MRCAPrior) {
                calibrations.add((MRCAPrior) plugin);
            }
        }
        
        boolean geneTreesNeedInit = true;
        if (newick != null) {
            Log.info.println("StarBEAST2: using initFromNewick to initialize species tree.");
            initFromNewick(speciesTree, newick);
        } else if (method == Method.ALL_RANDOM) {
            Log.info.println("StarBEAST2: using randomInit to initialize species tree.");
            randomInit(speciesTree, calibrations);
        } else if (!checkSpeciesAlwaysRepresented())  {
            Log.info.println("StarBEAST2: using randomInit to initialize species tree (required by missing data)).");
            randomInit(speciesTree, calibrations);
        } else if (calibrations.size() > 0)  {
            Log.info.println("StarBEAST2: using randomInit to initialize species tree (required by calibrations)).");
            randomInit(speciesTree, calibrations);
        } else if (speciesTree.hasDateTrait()) {
            Log.info.println("StarBEAST2: using randomInit to initialize species tree (required by tip dates).");
            randomInit(speciesTree, calibrations);
        } else if (method == Method.POINT) {
            Log.info.println("StarBEAST2: using fullInit to initialize all trees.");
            fullInit(speciesTree);
            geneTreesNeedInit = false;
        }

        if (geneTreesNeedInit) {
            final double rootHeight = speciesTree.getRoot().getHeight();
            Log.info.println(String.format("StarBEAST2: using randomInitGeneTrees to initialize gene trees (%f).", rootHeight));

            final List<Tree> geneTrees = genes.get();
            for (final Tree gtree : geneTrees) {
                gtree.makeCaterpillar(rootHeight, rootHeight/gtree.getInternalNodeCount(), true);

            }
            


            // make sure the heights of all gene tree tips is equal to the height of corresponding species tree tips
            resetGeneTreeTipHeights();
        }

        // initialize population sizes to equal average branch length
        // this is equivalent to 2Ne = E[1/lambda]
        final double speciesTreeLength = TreeStats.getLength(speciesTree);
        final int nBranches = speciesTree.getNodeCount();
        final double averageBranchLength = speciesTreeLength / (nBranches - 1);

        final PopulationModel populationModel = populationFunctionInput.get();
        if (populationModel != null) populationModel.initPopSizes(averageBranchLength);
        
        final Set<String> tipNames = speciesTree.getTipNumberMap().keySet();
        for (Tree geneTree: genes.get()) {
        	for (Node geneNode: geneTree.getNodesAsArray()) {
        		if (geneNode.isLeaf()) {
        			final String tipName = geneNode.getID();
        			if (!tipNames.contains(tipName)) {
        	            throw new RuntimeException(String.format("ERROR: Gene tree tip name '%s' is missing from taxon map. "
        	            		+ "This typically occurs when a sequence or sample name is identical to a species name. "
        	            		+ "Make sure all species names are distinct from sequence or sample names.", tipName));
        			}
        		}
        	}
        }
    }

    // check if every gene tree has at least one specimen from every species
    private boolean checkSpeciesAlwaysRepresented() {
        final SpeciesTree speciesTree = speciesTreeInput.get();
        final int nSpecies = speciesTree.getLeafNodeCount();
        final Multimap<Integer, String> numberTipMap = speciesTree.getNumberTipMap();
        for (Tree geneTree: genes.get()) {
        	final String[] taxaNames = geneTree.getTaxaNames();
        	for (int speciesNr = 0; speciesNr < nSpecies; speciesNr++) {
        		final Collection<String> speciesTaxa = numberTipMap.get(speciesNr);
        		boolean speciesRepresented = false;
        		for (String geneTaxon: taxaNames)
        			speciesRepresented |= speciesTaxa.contains(geneTaxon);
        		if (!speciesRepresented) return false;
        	}
        }

        return true;
    }

    private void resetGeneTreeTipHeights() {
    	final SpeciesTreeInterface speciesTree = speciesTreeInput.get();
        final Map<String, Integer> tipNames = speciesTree.getTipNumberMap();
        final List<Tree> geneTrees = genes.get();

        for (final Tree gtree : geneTrees) {
        	final int leafNodeCount = gtree.getLeafNodeCount();
        	for (int geneLeafNr = 0; geneLeafNr < leafNodeCount; geneLeafNr++) {
        		final Node geneLeaf = gtree.getNode(geneLeafNr);
        		if (tipNames.containsKey(geneLeaf.getID())) {
                    final int speciesLeafNr = tipNames.get(geneLeaf.getID());
                    final Node speciesLeaf = speciesTree.getNode(speciesLeafNr);
                    geneLeaf.setHeight(speciesLeaf.getHeight());
                } else {
        		    throw new RuntimeException("The taxon " + geneLeaf.getID() + " is missing from the taxonsuperset!");
                }
        	}
        }
	}

	private double[] firstMeetings(final Tree gtree, final Map<String, Integer> tipName2Species, final int speciesCount) {
        final Node[] nodes = gtree.listNodesPostOrder(null, null);
        @SuppressWarnings("unchecked")
		final Set<Integer>[] tipsSpecies = new Set[nodes.length];
        for(int k = 0; k < tipsSpecies.length; ++k) {
            tipsSpecies[k] = new LinkedHashSet<>();
        }
        // d[i,j] = minimum height of node which has tips belonging to species i and j
        // d is is upper triangular
        final double[] dmin = new double[(speciesCount*(speciesCount-1))/2];
        Arrays.fill(dmin, Double.MAX_VALUE);

        for (final Node n : nodes) {
            if (n.isLeaf()) {
                tipsSpecies[n.getNr()].add(tipName2Species.get(n.getID()));
            } else {
                assert n.getChildCount() == 2;
                @SuppressWarnings("unchecked")
				final Set<Integer>[] sps = new Set[2];
                sps[0] = tipsSpecies[n.getChild(0).getNr()];
                sps[1] = tipsSpecies[n.getChild(1).getNr()];
                final Set<Integer> u = new LinkedHashSet<>(sps[0]);
                u.retainAll(sps[1]);
                sps[0].removeAll(u);
                sps[1].removeAll(u);

                for (final Integer s1 : sps[0]) {
                    for (final Integer s2 : sps[1]) {
                        final int i = getDMindex(speciesCount, s1, s2);
                        dmin[i] = min(dmin[i], n.getHeight());
                    }
                }
                u.addAll(sps[0]);
                u.addAll(sps[1]);
                tipsSpecies[n.getNr()] = u;
            }
        }
        return dmin;
    }

    private int getDMindex(final int speciesCount, final int s1, final int s2) {
        final int mij = min(s1,s2);
        return (mij*(2*speciesCount-1 - mij))/2 + (abs(s1-s2)-1);
    }


    private void fullInit(final SpeciesTree speciesTree) {
        // Build gene trees from  alignments
    	

        final Function muInput = this.muInput.get();
        final double mu =  (muInput != null )  ? muInput.getArrayValue() : 1;

        final TaxonSet species = speciesTree.m_taxonset.get();
        final List<String> speciesNames = species.asStringList();
        final int speciesCount = speciesNames.size();

        final List<Tree> geneTrees = genes.get();

        //final List<Alignment> alignments = genes.get();
        //final List<Tree> geneTrees = new ArrayList<>(alignments.size());
        double maxNsites = 0;
        //for( final Alignment alignment : alignments)  {
        for (final Tree gtree : geneTrees) {
            //final Tree gtree = new Tree();
            final Alignment alignment = gtree.m_taxonset.get().alignmentInput.get();

            final ClusterTree ctree = new ClusterTree();
            ctree.initByName("initial", gtree, "clusterType", "upgma", "taxa", alignment);
            gtree.scale(1 / mu);

            maxNsites = max(maxNsites, alignment.getSiteCount());
        }
        final Map<String, Integer> geneTips2Species = new LinkedHashMap<>();
        final List<Taxon> taxonSets = species.taxonsetInput.get();

        for(int k = 0; k < speciesNames.size(); ++k) {
            final Taxon nx = taxonSets.get(k);
            final List<Taxon> taxa = ((TaxonSet) nx).taxonsetInput.get();
            for( final Taxon n : taxa ) {
              geneTips2Species.put(n.getID(), k);
            }
        }
        final double[] dg = new double[(speciesCount*(speciesCount-1))/2];

        final double[][] genesDmins = new double[geneTrees.size()][];

        for( int ng = 0; ng < geneTrees.size(); ++ng ) {
            final Tree g = geneTrees.get(ng);
            final double[] dmin = firstMeetings(g, geneTips2Species, speciesCount);
            genesDmins[ng] = dmin;

            for(int i = 0; i < dmin.length; ++i) {
                dg[i] += dmin[i];
                if (dmin[i] == Double.MAX_VALUE) {
                	// this happens when a gene tree has no taxa for some species-tree taxon.
                	// TODO: ensure that if this happens, there will always be an "infinite"
                	// distance between species-taxon 0 and the species-taxon with missing lineages,
                	// so i < speciesCount - 1.
                	// What if lineages for species-taxon 0 are missing? Then all entries will be 'infinite'.
                	String id = (i < speciesCount - 1? speciesTree.getExternalNodes().get(i+1).getID() : "unknown taxon");
                	if (i == 0) {
                		// test that all entries are 'infinite', which implies taxon 0 has lineages missing 
                		boolean b = true;
                		for (int k = 1; b && k < speciesCount - 1; k++) {
                			b = (dmin[k] == Double.MAX_VALUE);
                		}
                		if (b) {
                			// if all entries have 'infinite' distances, it is probably the first taxon that is at fault
                			id = speciesTree.getExternalNodes().get(0).getID();
                		}
                	}
                	throw new RuntimeException("Gene tree " + g.getID() + " has no lineages for species taxon " + id + " ");
                }
            }
        }

        for(int i = 0; i < dg.length; ++i) {
            double d = dg[i] / geneTrees.size();
            if( d == 0 ) {
               d = (0.5/maxNsites) * (1/mu);
            } else {
                // heights to distances
                d *= 2;
            }
            dg[i] = d;
        }

        final ClusterTree ctree = new ClusterTree();
        final Distance distance = new Distance() {
            @Override
            public double pairwiseDistance(final int s1, final int s2) {
                final int i = getDMindex(speciesCount, s1,s2);
                return dg[i];
            }
        };
        ctree.initByName("initial", speciesTree, "taxonset", species,"clusterType", "upgma", "distance", distance);

        final Map<String, Integer> sptips2SpeciesIndex = new LinkedHashMap<>();
        for(int i = 0; i < speciesNames.size(); ++i) {
            sptips2SpeciesIndex.put(speciesNames.get(i), i);
        }
        final double[] spmin = firstMeetings(speciesTree, sptips2SpeciesIndex, speciesCount);

        for( int ng = 0; ng < geneTrees.size(); ++ng ) {
            final double[] dmin = genesDmins[ng];
            boolean compatible = true;
            for(int i = 0; i < spmin.length; ++i) {
                if( dmin[i] <= spmin[i] ) {
                    compatible = false;
                    break;
                }
            }
            if( ! compatible ) {
                final Tree gtree = geneTrees.get(ng);
                final TaxonSet gtreeTaxa = gtree.m_taxonset.get();
                final Alignment alignment = gtreeTaxa.alignmentInput.get();
                final List<String> taxaNames = alignment.getTaxaNames();
                final int taxonCount =  taxaNames.size();
                // speedup
                final Map<Integer,Integer> g2s = new LinkedHashMap<>();
                for(int i = 0; i < taxonCount; ++i) {
                    g2s.put(i, geneTips2Species.get(taxaNames.get(i)));
                }

                final JukesCantorDistance jc = new JukesCantorDistance();
                jc.setPatterns(alignment);
                final Distance gdistance = new Distance() {
                    @Override
                    public double pairwiseDistance(final int t1, final int t2) {
                        final int s1 = g2s.get(t1);
                        final int s2 = g2s.get(t2);
                        double d = jc.pairwiseDistance(t1,t2)/mu;
                        if( s1 != s2 ) {
                            final int i = getDMindex(speciesCount, s1,s2);
                            final double minDist = 2 * spmin[i];
                            if( d <= minDist ) {
                                d = minDist * 1.001;
                            }
                        }
                        return d;
                    }
                };
                final ClusterTree gtreec = new ClusterTree();
                gtreec.initByName("initial", gtree, "taxonset", gtreeTaxa,
                        "clusterType", "upgma", "distance", gdistance);
            }
        }

        final RealParameter lambda = birthRate.get();
        if (lambda != null && lambda instanceof StateNode) {
            final StateNode lambdaStateNode = (StateNode) lambda;

            // only change lambda if it is to be estimated
            if (lambdaStateNode.isEstimatedInput.get()) {
                final double rh = speciesTree.getRoot().getHeight();
                double l = 0;
                for(int i = 2; i < speciesCount+1; ++i) l += 1./i;
                lambda.setValue((1 / rh) * l);
            }
        }
    }

    private void randomInit(final SpeciesTree speciesTree, List<MRCAPrior> calibrations) {
    	final RealParameter birthRateParameter = birthRate.get();
    	final Double lambda = (birthRateParameter == null) ? 1.0 : birthRateParameter.getValue();
    	final Double initialPopSize = 1.0 / lambda; // scales coalescent tree height inverse to birth rate
    	final RealParameter popSize = new RealParameter(initialPopSize.toString());
        final ConstantPopulation pf = new ConstantPopulation();
        pf.setInputValue("popSize", popSize);

        final RandomTree rnd = new RandomTree();
        rnd.setInputValue("taxonset", speciesTree.getTaxonset());
        if (speciesTree.hasDateTrait()) rnd.setInputValue("trait", speciesTree.getDateTrait());

        for (final MRCAPrior cal: calibrations) rnd.setInputValue("constraint", cal);

        rnd.setInputValue("populationModel", pf);
        rnd.setInputValue("populationModel", pf);
        rnd.initAndValidate();

        speciesTree.assignFromWithoutID(rnd);
        // System.out.println("BEFORE = " + speciesTree.toString());

        /*final TraitSet speciesTipDates = speciesTree.getDateTrait();
        if (speciesTree.hasDateTrait()) {
        	for (Node node: speciesTree.getNodesAsArray()) {
        		if (node.isLeaf()) {
        			final String taxonName = node.getID();
        			final double taxonHeight = speciesTipDates.getValue(taxonName);
        			node.setHeight(taxonHeight);
        			
        		}
        	}
            System.out.println("AFTER = " + speciesTree.toString());
        }*/
    }

    private void initFromNewick(final Tree speciesTree, final String newick) {
        final TreeParser parser = new TreeParser(newick);
        speciesTree.assignFromWithoutID(parser);
        System.out.println(speciesTree.toString());
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
        stateNodes.add(speciesTreeInput.get());

        for (final Tree g : genes.get()) {
            stateNodes.add(g);
        }

        final RealParameter brate = birthRate.get();
        if (brate != null) {
            stateNodes.add(brate) ;
        }
    }
}
