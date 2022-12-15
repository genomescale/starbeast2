package starbeast2.aimoperator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import starbeast2.NodeHeightComparator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@Description("Randomly selects true internal tree node (i.e. not the root) and move node height uniformly in interval "
		+ "restricted by the nodes parent and children.")
public class RankingAwareOperator extends TreeOperator {

	final public Input<RealParameter> migrationInput = new Input<>("m", "input of ne's of the species",
			Validate.REQUIRED);
	final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "input of ne's of the species",
			Validate.OPTIONAL);

	final static NodeHeightComparator nhc = new NodeHeightComparator();

	// empty constructor to facilitate construction by XML + initAndValidate
	public RankingAwareOperator() {
	}

	public RankingAwareOperator(Tree tree) {
		try {
			initByName(treeInput.getName(), tree);
		} catch (Exception e) {
			e.printStackTrace(); // To change body of catch statement use File | Settings | File Templates.
			throw new RuntimeException("Failed to construct Uniform Tree Operator.");
		}
	}

	@Override
	public void initAndValidate() {
	}

	public double treeProposal() {
		return 0.0;
	}

	/**
	 * change the parameter and return the hastings ratio.
	 *
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted *
	 */
	@Override
	public double proposal() {

		double HR = 0.0;
		final Tree tree = (Tree) InputUtil.get(treeInput, this);

		Node[] speciesNodesTmp = tree.getNodesAsArray();
		Node[] speciesNodes = new Node[speciesNodesTmp.length];
		// make a deep copy of every node
		for (int i = 0; i < speciesNodesTmp.length; i++)
			speciesNodes[i] = speciesNodesTmp[i].copy();

		// sort the species tree nodes by height
		Arrays.sort(speciesNodes, nhc);
		
		// make the move on the tree
		HR += treeProposal();

		Node[] newSpeciesNodesTmp = tree.getNodesAsArray();
		Node[] newSpeciesNodes = new Node[speciesNodesTmp.length];
		// make a deep copy of every node
		for (int i = 0; i < speciesNodesTmp.length; i++)
			newSpeciesNodes[i] = newSpeciesNodesTmp[i].copy();

		// sort the species tree nodes of the newly proposed by height
		Arrays.sort(newSpeciesNodes, nhc);

		// check if the ranked tree changed (both in topology and ranking)
		boolean same = true;
		for (int i = 0; i < newSpeciesNodes.length; i++)
			if (!nodeIsSame(newSpeciesNodes[i], speciesNodes[i]))
				same = false;
		
		// if the topology did not change, do nothing to the migration rates
		if (same)
			return HR;	

		// check which node numbers do not have an analouge in the new map
		boolean[] hasAnalogue = new boolean[newSpeciesNodes.length];
		for (int i = 0; i < newSpeciesNodes.length; i++)
			hasAnalogue[i] = nodeStillExists(newSpeciesNodes, speciesNodes[i]);

		
		// get the migration rate mapping before and after the move
		ArrayList<Integer[]> oldMigMap = getMigrationMap(speciesNodes);
		ArrayList<Integer[]> newMigMap = getMigrationMap(newSpeciesNodes);

		int[] mapping = new int[oldMigMap.size()];
		double[] newmigvals = new double[oldMigMap.size()];
		double[] newintvals = new double[oldMigMap.size()];
		

		for (int i = 0; i < oldMigMap.size(); i++)
			mapping[i] = getIndexOf(newMigMap, oldMigMap.get(i), hasAnalogue);
		

		// keeps track of the migration routes that have no one to one mapping
		List<Integer> removedVals = new ArrayList<>();
		// check which values do not exist
		for (int i = 0; i<mapping.length;i++)
			removedVals.add(i);
		
		for (int i = 0; i<mapping.length;i++)
			if (mapping[i]!=-1)
				removedVals.remove(removedVals.indexOf(mapping[i]));
				
		List<Integer> shuffeledVals = new ArrayList<>();
		while (removedVals.size() > 0) {
			int index = Randomizer.nextInt(removedVals.size());
			shuffeledVals.add(removedVals.get(index));
			removedVals.remove(index);
		}

		// shuffel the order of the parameter that appear and disappear
		int c = 0;
		for (int i = 0; i < oldMigMap.size(); i++) {
			int index = mapping[i];
			if (index == -1) {
				index = shuffeledVals.get(c);
				c++;
			}

			newmigvals[index] = migrationInput.get().getArrayValue(i);
			if (indicatorInput.get() != null)
				newintvals[index] = indicatorInput.get().getArrayValue(i);
		}

		for (int i = 0; i < newmigvals.length; i++) {
			migrationInput.get().setValue(i, newmigvals[i]);
			if (indicatorInput.get() != null)
				if (newintvals[i] > 0.5)
					indicatorInput.get().setValue(i, true);
				else
					indicatorInput.get().setValue(i, false);

		}

		

		return HR;
	}

	private boolean nodeIsSame(Node n1, Node n2) {
		List<Node> l1 = n1.getAllLeafNodes();
		List<Node> l2 = n2.getAllLeafNodes();

		if (l1.size() != l2.size())
			return false;

		List<Integer> nr1 = new ArrayList<>();
		List<Integer> nr2 = new ArrayList<>();

		for (int i = 0; i < l1.size(); i++) {
			nr1.add(l1.get(i).getNr());
			nr2.add(l2.get(i).getNr());
		}

		for (int i = 0; i < l1.size(); i++)
			if (nr1.indexOf(nr2.get(i)) == -1)
				return false;

		return true;
	}

	private boolean nodeStillExists(Node[] newNodes, Node n) {
		for (Node newN : newNodes) {
			if (nodeIsSame(newN, n))
				return true;
		}
		return false;
	}
	
	@SuppressWarnings({ "unchecked", "deprecation" })
	public ArrayList<Integer[]> getMigrationMap(Node[] speciesNodes) {

		double[] intervals = new double[speciesNodes.length];
		boolean[] isCoalescent = new boolean[speciesNodes.length];
		List<Node>[] lineagesAdded = new List[speciesNodes.length];
		List<Node>[] lineagesRemoved = new List[speciesNodes.length];

		ArrayList<Node> initAL = new ArrayList<>();
		intervals[0] = speciesNodes[0].getHeight();
		lineagesAdded[0] = new ArrayList<>();
		lineagesAdded[0].add(speciesNodes[0].copy());
		isCoalescent[0] = false;
		lineagesRemoved[0] = new ArrayList<>();
		lineagesRemoved[0].add(null);
		lineagesRemoved[0].add(null);
		for (int i = 1; i < speciesNodes.length; i++) {
			intervals[i] = speciesNodes[i].getHeight() - speciesNodes[i - 1].getHeight();
			lineagesAdded[i] = new ArrayList<>();
			lineagesAdded[i].add(speciesNodes[i].copy());
			if (!speciesNodes[i].isLeaf()) {
				isCoalescent[i] = true;
				lineagesRemoved[i] = new ArrayList<>();
				lineagesRemoved[i].add(speciesNodes[i].getLeft().copy());
				lineagesRemoved[i].add(speciesNodes[i].getRight().copy());
				lineagesRemoved[i].get(0).setParent(lineagesAdded[i].get(0));
				lineagesRemoved[i].get(1).setParent(lineagesAdded[i].get(0));
			} else {
				isCoalescent[i] = false;
				lineagesRemoved[i] = new ArrayList<>();
				lineagesRemoved[i].add(null);
				lineagesRemoved[i].add(null);
			}
		}

		ArrayList<Integer> activeStates = new ArrayList<>();

		ArrayList<Integer[]> migrationMap = new ArrayList<>();

		int a = 0;
		while (!isCoalescent[a]) {
			activeStates.add(lineagesAdded[a].get(0).getNr());
			a++;
		}
		// sorting ensures the correct order of migration rate elements
		Collections.sort(activeStates);

		for (int i = 0; i < activeStates.size(); i++) {
			for (int j = 0; j < activeStates.size(); j++) {
				if (i != j) {
					Integer[] migroute = new Integer[] { activeStates.get(i), activeStates.get(j) };
					migrationMap.add(migroute);

				}
			}
		}

		while (a < intervals.length) {
			int parent = lineagesAdded[a].get(0).getNr();
			activeStates.add(parent);
			int d1 = activeStates.indexOf(lineagesRemoved[a].get(0).getNr());
			activeStates.remove(d1);
			int d2 = activeStates.indexOf(lineagesRemoved[a].get(1).getNr());
			activeStates.remove(d2);
			for (int j = 0; j < activeStates.size() - 1; j++) {
				Integer[] migroute = new Integer[] { activeStates.get(j), parent };
				migrationMap.add(migroute);
			}
			for (int j = 0; j < activeStates.size() - 1; j++) {
				Integer[] migroute = new Integer[] { parent, activeStates.get(j) };
				migrationMap.add(migroute);
			}
			a++;
		}

		return migrationMap;

	}

	private int getIndexOf(ArrayList<Integer[]> maps, Integer[] route, boolean[] hasAnalogue) {
		if (!hasAnalogue[route[0]] || !hasAnalogue[route[1]])
			return -1;
		
		for (int i = 0; i < maps.size(); i++)
			if (maps.get(i)[0] == route[0] && maps.get(i)[1] == route[1])
				return i;

		return -1;
	}

}
