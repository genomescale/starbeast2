package speciesnetwork;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.IntegerParameterList;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.alignment.distance.JukesCantorDistance;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.ClusterTree;

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
    final public Input<Method> initMethod = new Input<>("method", "Initialise either with a totally random state" +
            "or a point estimate based on alignments data (default point-estimate)", Method.POINT, Method.values());

    final public Input<Network> speciesNetworkInput = new Input<>("speciesNetwork", "Species network to initialize.");

    final public Input<List<Tree>> geneTreesInput = new Input<>("geneTrees", "Gene trees to initialize.", new ArrayList<>());

    final public Input<YuleHybridModel> hybridYuleInput = new Input<>("hybridYule",
            "The species network (with hybridization) to initialize.", Validate.XOR, speciesNetworkInput);

    final public Input<RealParameter> birthRateInput = new Input<>("birthRate", "Network prior birth rate to initialize.");

    final public Input<RealParameter> hybridRateInput = new Input<>("hybridRate", "Network hybridization rate to initialize.");

    final public Input<Function> muInput = new Input<>("baseRate", "Main clock rate used to scale trees (default 1).");

    final public Input<PopulationSizeModel> populationFunctionInput = new Input<>("populationModel",
            "The species network population size model.", Validate.REQUIRED);

    // private boolean hasCalibrations;  // don't deal with calibrations at the moment

    @Override
    public void initAndValidate() throws Exception {
        // what does this do and is it dangerous to call it or not to call it at the start or at the end??????
        super.initAndValidate();
        // do we need to initialize both species network and gene trees in this method???
    }

    @Override
    public void initStateNodes() throws Exception {
        // initialize population sizes to equal average branch length (equivalent to 2Ne = E[1/lambda])
        final NetworkNode speciesNetworkRoot = speciesNetworkInput.get().getRoot();

        double speciesNetworkLength = 0;
        for (final NetworkNode n : speciesNetworkRoot.getAllChildNodes()) {
            if(n.getLeftParent() != null)  speciesNetworkLength += n.getLeftLength();
            if(n.getRightParent() != null) speciesNetworkLength += n.getRightLength();
        }

        final int nSpeciesBranches = speciesNetworkInput.get().getBranchCount();
        final double averageBranchLength = speciesNetworkLength / (nSpeciesBranches - 1);

        final PopulationSizeModel populationModel = populationFunctionInput.get();
        populationModel.initPopSizes(averageBranchLength);

        // if( hasCalibrations ) { initWithCalibrations(); } else {
        final Method method = initMethod.get();
        switch( method ) {
            case POINT:
                fullInit();
                break;
            case ALL_RANDOM:
                randomInit();
                break;
        }
    }

    private double[] firstMeetings(final Tree gTree, final Map<String, Integer> tipName2Species, final int nSpecies) {
        final Node[] nodes = gTree.listNodesPostOrder(null, null);
        @SuppressWarnings("unchecked")
		final Set<Integer>[] tipsSpecies = new Set[nodes.length];
        for(int k = 0; k < tipsSpecies.length; ++k) {
            tipsSpecies[k] = new HashSet<>();
        }
        // d[i,j] = minimum height of node which has tips belonging to species i and j
        // d is is upper triangular
        final double[] dmin = new double[(nSpecies*(nSpecies-1))/2];
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
                final Set<Integer> u = new HashSet<>(sps[0]);
                u.retainAll(sps[1]);
                sps[0].removeAll(u);
                sps[1].removeAll(u);

                for (final Integer s1 : sps[0]) {
                    for (final Integer s2 : sps[1]) {
                        final int i = getDMindex(nSpecies, s1, s2);
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

    private double[] firstMeetings(final Network network, final Map<String, Integer> tipName2Species, final int nSpecies) {
        final NetworkNode[] nodes = network.getAllNodesAsArray();
        @SuppressWarnings("unchecked")
        final Set<Integer>[] tipsSpecies = new Set[nodes.length];
        for(int k = 0; k < tipsSpecies.length; ++k) {
            tipsSpecies[k] = new HashSet<>();
        }
        // d[i,j] = minimum height of node which has tips belonging to species i and j
        // d is is upper triangular
        final double[] dmin = new double[(nSpecies*(nSpecies-1))/2];
        Arrays.fill(dmin, Double.MAX_VALUE);

        for (final NetworkNode n : nodes) {
            if (n.isLeaf()) {
                tipsSpecies[n.getNr()].add(tipName2Species.get(n.getID()));
            } else {
                assert n.getChildCount() == 2;
                @SuppressWarnings("unchecked")
                final Set<Integer>[] sps = new Set[2];
                sps[0] = tipsSpecies[n.getLeftChild().getNr()];
                sps[1] = tipsSpecies[n.getRightChild().getNr()];
                final Set<Integer> u = new HashSet<>(sps[0]);
                u.retainAll(sps[1]);
                sps[0].removeAll(u);
                sps[1].removeAll(u);

                for (final Integer s1 : sps[0]) {
                    for (final Integer s2 : sps[1]) {
                        final int i = getDMindex(nSpecies, s1, s2);
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

    private int getDMindex(final int nSpecies, final int s1, final int s2) {
        final int mij = min(s1,s2);
        return (mij*(2*nSpecies-1 - mij))/2 + (abs(s1-s2)-1);
    }

    private void fullInit() throws Exception {
        // Build gene trees from  alignments

        final Function muInput = this.muInput.get();
        final double mu = (muInput != null) ? muInput.getArrayValue() : 1;

        final Network sNetwork = speciesNetworkInput.get();
        final TaxonSet species = sNetwork.taxonSetInput.get();
        final List<String> speciesNames = species.asStringList();
        final int nSpecies = speciesNames.size();

        final List<Tree> geneTrees = geneTreesInput.get();

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
        final Map<String, Integer> geneTips2Species = new HashMap<>();
        final List<Taxon> taxonSets = species.taxonsetInput.get();

        for(int k = 0; k < speciesNames.size(); ++k) {
            final Taxon nx = taxonSets.get(k);
            final List<Taxon> taxa = ((TaxonSet) nx).taxonsetInput.get();
            for( final Taxon n : taxa ) {
              geneTips2Species.put(n.getID(), k);
            }
        }
        final double[] dg = new double[(nSpecies*(nSpecies-1))/2];

        final double[][] genesDmins = new double[geneTrees.size()][];

        for(int ng = 0; ng < geneTrees.size(); ++ng) {
            final Tree gtree = geneTrees.get(ng);
            final double[] dmin = firstMeetings(gtree, geneTips2Species, nSpecies);
            genesDmins[ng] = dmin;

            for(int i = 0; i < dmin.length; ++i) {
                dg[i] += dmin[i];
                if (dmin[i] == Double.MAX_VALUE) {
                	// this happens when a gene tree has no taxa for some species-tree taxon.
                	// TODO: ensure that if this happens, there will always be an "infinite"
                	// distance between species-taxon 0 and the species-taxon with missing lineages,
                	// so i < nSpecies - 1.
                	// What if lineages for species-taxon 0 are missing? Then all entries will be 'infinite'.
                	String id = (i < nSpecies - 1? sNetwork.getLeafNodes().get(i+1).getID() : "unknown taxon");
                	if (i == 0) {
                		// test that all entries are 'infinite', which implies taxon 0 has lineages missing 
                		boolean b = true;
                		for (int k = 1; b && k < nSpecies - 1; k++) {
                			b = (dmin[k] == Double.MAX_VALUE);
                		}
                		if (b) {
                			// if all entries have 'infinite' distances, it is probably the first taxon that is at fault
                			id = sNetwork.getLeafNodes().get(0).getID();
                		}
                	}
                	throw new RuntimeException("Gene tree " + gtree.getID() + " has no lineages for species taxon " + id + " ");
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
                final int i = getDMindex(nSpecies, s1,s2);
                return dg[i];
            }
        };
        ctree.initByName("initial", sNetwork, "taxonset", species,"clusterType", "upgma", "distance", distance);

        final Map<String, Integer> sptips2SpeciesIndex = new HashMap<>();
        for(int i = 0; i < speciesNames.size(); ++i) {
            sptips2SpeciesIndex.put(speciesNames.get(i), i);
        }
        // initialize species network as a tree, so it should be okay to use the same function???
        final double[] spmin = firstMeetings(sNetwork, sptips2SpeciesIndex, nSpecies);

        for(int ng = 0; ng < geneTrees.size(); ++ng) {
            final double[] dmin = genesDmins[ng];
            boolean compatible = true;
            for(int i = 0; i < spmin.length; ++i) {
                if( dmin[i] <= spmin[i] ) {
                    compatible = false;
                    break;
                }
            }
            if(!compatible) {
                final Tree gtree = geneTrees.get(ng);
                final TaxonSet gtreeTaxa = gtree.m_taxonset.get();
                final Alignment alignment = gtreeTaxa.alignmentInput.get();
                final List<String> taxaNames = alignment.getTaxaNames();
                final int nTaxa =  taxaNames.size();
                // speedup
                final Map<Integer,Integer> g2s = new HashMap<>();
                for(int i = 0; i < nTaxa; ++i) {
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
                            final int i = getDMindex(nSpecies, s1,s2);
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

        final RealParameter lambda = birthRateInput.get();
        if(lambda != null) {
            final double rh = sNetwork.getRoot().getHeight();
            double l = 0;
            for(int i = 2; i < nSpecies+1; ++i) {
                l += 1./i;
            }
            lambda.setValue(l / rh);
        }

        final RealParameter nu = hybridRateInput.get();
        if(nu != null) {
            nu.setValue(0.001);
        }
    }

    private void randomInit() throws Exception {
        double lam = 1;
        final RealParameter lambda = birthRateInput.get();
        if(lambda != null) {
            lam = lambda.getArrayValue();
        }

        final Network sNetwork = speciesNetworkInput.get();
        final TaxonSet species = sNetwork.taxonSetInput.get();
        final int nSpecies = species.asStringList().size();
        double s = 0;
        for(int k = 2; k <= nSpecies; ++k) {
            s += 1.0/k;
        }
        final double rootHeight = (1/lam) * s;
        sNetwork.scale(rootHeight/sNetwork.getRoot().getHeight());

        // randomInitGeneTrees(rootHeight);
        final List<Tree> geneTrees = geneTreesInput.get();
        for (final Tree gtree : geneTrees) {
            gtree.makeCaterpillar(rootHeight, rootHeight/gtree.getInternalNodeCount(), true);

            /* figure out the mapping and write IntegerParameterList
            * -1 -> not passing the species network node
            * 0 -> left parent/branch
            * 1 -> right parent/branch
            */
            // IntegerParameterList mapping = gtree.getTreeMappingToNetwork();

        }
    }

    // private void initWithCalibrations() throws Exception {}

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        // if( hasCalibrations ) { stateNodes.add((Tree) calibratedYule.get().treeInput.get()); } else {
        stateNodes.add(speciesNetworkInput.get());

        for(final Tree gtree : geneTreesInput.get()) {
            stateNodes.add(gtree);
        }

        final RealParameter brate = birthRateInput.get();
        if(brate != null) {
            stateNodes.add(brate) ;
        }
    }
}