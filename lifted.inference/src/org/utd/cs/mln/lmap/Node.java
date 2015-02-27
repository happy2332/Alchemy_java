package org.utd.cs.mln.lmap;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import org.utd.cs.gm.utility.Pair;
import org.utd.cs.mln.alchemy.core.HyperCube;

public class Node {

	
	public ArrayList<ArrayList<HyperCube>> hyperCubesList = new ArrayList<ArrayList<HyperCube>>();
	public ArrayList<ArrayList<ArrayList<Integer>>> partialGroundings = new ArrayList<ArrayList<ArrayList<Integer>>>();
	//public ArrayList<Integer> groundedTermList = new ArrayList<Pair>();
	public ArrayList<Node> children = new ArrayList<Node>();
	public ArrayList<Integer> placeHolderList = new ArrayList<Integer>();
	public Node rootParent = null;
	//public Node nodeParent = new Node();
	public boolean isRoot = false;
	
	public Node(){}
	
	public Node(int maxPredId){
		for(int i = 0 ; i <= maxPredId ; i++){
			partialGroundings.add(new ArrayList<ArrayList<Integer>>());
		}
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
