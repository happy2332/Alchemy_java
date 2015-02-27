package org.utd.cs.mln.lmap;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import javax.print.attribute.Size2DSyntax;

import org.utd.cs.gm.utility.Pair;
import org.utd.cs.mln.alchemy.core.Atom;
import org.utd.cs.mln.alchemy.core.Evidence;
import org.utd.cs.mln.alchemy.core.HyperCube;
import org.utd.cs.mln.alchemy.core.MLN;
import org.utd.cs.mln.alchemy.core.PredicateNotFound;
import org.utd.cs.mln.alchemy.core.PredicateSymbol;
import org.utd.cs.mln.alchemy.core.Term;
import org.utd.cs.mln.alchemy.core.WClause;
import org.utd.cs.mln.alchemy.util.Parser;

public class LiftedPTP {

	public static int maxPredId = 0;
	public static double getPartition(MLN mln, boolean approximate){
		// base case : If every hyperCube of every clause
		if(isBaseCaseReached(mln) == true){
			return findBaseCaseZ(mln);
		}
		
		// Decomposition step
		ArrayList<Set<Pair>> predEquivalenceClasses = new ArrayList<Set<Pair>>();
		Decomposer d = new Decomposer();
		d.findEquivalenceClasses(mln, predEquivalenceClasses);
		ArrayList<MLN> decomposedMLNs = null;
		ArrayList<Integer> sizeOfSegments = new ArrayList<Integer>();
		for(int eqClassId = 0 ; eqClassId < predEquivalenceClasses.size() ; eqClassId++){
			if(d.isDecomposer(predEquivalenceClasses.get(eqClassId), mln)){
				System.out.println("Doing decomposition...");
				decomposedMLNs = d.initDecomposer(mln, predEquivalenceClasses.get(eqClassId), sizeOfSegments);
				if(decomposedMLNs != null){
					break;
				}
			}
		}
		 
		if(decomposedMLNs != null && decomposedMLNs.size() > 0){
			System.out.println("Decomposition Done...");
			double totalZ = 1;
			///System.out.println("Printing decomposed MLNs : ");
			for(int mlnId = 0 ; mlnId < decomposedMLNs.size() ; mlnId++){
				//System.out.println("MLN no. : " + mlnId+1);
				MLN decomposedMln = decomposedMLNs.get(mlnId);
				/*//
				for(WClause clause : decomposedMln.clauses){
					clause.print();
					System.out.println("HyperCubes : ");
					for(HyperCube hc : clause.hyperCubes){
						System.out.println(hc);
					}
				}*///
				totalZ *= Math.pow(getPartition(decomposedMln, approximate),sizeOfSegments.get(mlnId));
			}
			return totalZ;
		}
		///System.out.println("Doing Splitting...");
		return LiftedSplit_new.liftedSplit(mln, approximate);
	}
	
	// return (2^#groundings)*exp(sum(for each hypercube in each clause) #hyperCubegroundings*wt*num_copies
	private static double findBaseCaseZ(MLN mln) {
		double weightedSumClauses = 0.0;
		HashMap<Integer,Set<ArrayList<Integer>>> predsGroundings = new HashMap<Integer,Set<ArrayList<Integer>>>(); 
		for(WClause clause : mln.clauses){
			for(int childId = 0 ; childId < clause.root.children.size() ; childId++){
				int numSatisfiedHyperCubes = 0;
				for(Integer phId : clause.root.children.get(childId).placeHolderList){
					for(HyperCube hc : clause.root.hyperCubesList.get((int)phId)){
						if(hc.satisfied == false){
							continue;
						}
						ArrayList<ArrayList<Set<Integer>>> hyperCubeTuples = LiftedSplit.cartesianProd(hc.varConstants);
						for(Atom atom : clause.atoms){
							Set<ArrayList<Integer>> predGroundings = new HashSet<ArrayList<Integer>>();
							for(ArrayList<Set<Integer>> hyperCubeTuple : hyperCubeTuples){
								ArrayList<Integer> predGrounding = new ArrayList<Integer>();
								for(Term term : atom.terms){
									predGrounding.add(hyperCubeTuple.get(clause.terms.indexOf(term)).iterator().next());
								}
								predGroundings.add(predGrounding);
							}
							if(predsGroundings.containsKey(atom.symbol.id)){
								predsGroundings.get(atom.symbol.id).addAll(predGroundings);
							}
							else{
								predsGroundings.put(atom.symbol.id, predGroundings);
							}
						} // end of atom loop
						numSatisfiedHyperCubes += hyperCubeTuples.size()*hc.num_copies;
					}// end of hyperCube loop
				} // end of phloop
				weightedSumClauses += numSatisfiedHyperCubes * clause.weight.getValue();
			}
		} // end of clause loop
		int totalNumPredGroundings = 0;
		for(Integer predId : predsGroundings.keySet()){
			totalNumPredGroundings += predsGroundings.get(predId).size();
		}
		return Math.pow(2, totalNumPredGroundings)*Math.exp(weightedSumClauses);
	}
	
	// Base case is reached when every hyperCube of every clause is either satisfied or empty
	private static boolean isBaseCaseReached(MLN mln) {
		for(WClause clause : mln.clauses){
			for(int i = 0 ; i < clause.root.hyperCubesList.size() ; i++){
				for(HyperCube hc : clause.root.hyperCubesList.get(i)){
					if(hc.satisfied == false && !hc.isEmpty()){
						return false;
					}
				}
			}
		}
		return true;
	}
	/*
	public static ArrayList<MLN> baseSolver(MLN mln){
		
	}*/
	/**
	 * @param args
	 * @throws FileNotFoundException 
	 * @throws PredicateNotFound 
	 */
	public static void main(String[] args) throws FileNotFoundException, PredicateNotFound {
		long time = System.currentTimeMillis();
		MLN mln = new MLN();
		Parser parser = new Parser(mln);
		String filename = new String("smoke/smoke_mln.txt");
		//String filename = new String("PTP_data/random_pkb_mln2.txt");
		parser.parseInputMLNFile(filename);

		//ArrayList<Evidence> evidList = parser.parseInputEvidenceFile("PTP_data/random_evidence.txt");
		ArrayList<Evidence> evidList = parser.parseInputEvidenceFile("smoke/evidence.txt");
		MlnToHyperCube mlnToHyperCube = new MlnToHyperCube();
		HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap = mlnToHyperCube.createPredsHyperCube(evidList,mln);
	
		int origNumClauses = mln.clauses.size();
		boolean isNormal = false;
		for(int clauseId = 0 ; clauseId < origNumClauses ; clauseId++){
			mln.clauses.addAll(mlnToHyperCube.createClauseHyperCube(mln.clauses.get(clauseId), predsHyperCubeHashMap, isNormal));
		}
		for(int clauseId = origNumClauses-1 ; clauseId >= 0 ; clauseId--){
			mln.clauses.remove(clauseId);
		}
		LiftedPTP.maxPredId = mln.max_predicate_id;
		System.out.println("Time taken to create clauses in hyperCube form : " + (long)(System.currentTimeMillis() - time) + " ms");
		time = System.currentTimeMillis();
		/*
		mln.clauses.get(0).hyperCubes.get(0).varConstants.set(0,new HashSet<Integer>(Arrays.asList(0,1,2)));
		mln.clauses.get(0).hyperCubes.get(0).varConstants.set(1,new HashSet<Integer>(Arrays.asList(0)));*/
		//mln.clauses.get(0).hyperCubes.get(0).satisfied = true;
		//mln.clauses.get(1).hyperCubes.get(0).satisfied = true;
		//mln.clauses.get(1).hyperCubes.get(0).varConstants.set(0,new HashSet<Integer>(Arrays.asList(0)));*/
		/*
		HyperCube h = new HyperCube();
		h.varConstants.add(new HashSet<Integer>(Arrays.asList(1,2)));
		h.varConstants.add(new HashSet<Integer>(Arrays.asList(2)));
		mln.clauses.get(0).hyperCubes.add(new HyperCube(h));
		*
		*/
		setValidPredPos(mln);
		for(WClause clause : mln.clauses){
			createDummyLeaves(clause);
		}
		System.out.println("Printing initial set of hyperCubes for each clause");
		double Z = 0.0;
		boolean approximate = true;
		time = System.currentTimeMillis();
		if(approximate == true){
			int num_iter = 10;
			for(int iter = 0 ; iter < num_iter ; iter++){
				MLN tempMln = new MLN();
				for(WClause clause : mln.clauses){
					WClause newClause = MLN.create_new_clause(clause);
					tempMln.clauses.add(newClause);
				}
				double curZ = getPartition(tempMln, approximate);
				System.out.println("Iteration no. : " + (int)(iter+1));
				//System.out.println("curZ = " + curZ);
				Z += curZ;
			}
			Z = (Z/num_iter);
			System.out.println("Time taken for " + num_iter + " iterations : " + (double)(System.currentTimeMillis()-time)/1000.0 + " sec");
		}
		else{
			Z = getPartition(mln, approximate);
			System.out.println("Time taken for exact inference : " + (double)(System.currentTimeMillis()-time)/1000.0 + " sec");
		}
		
		System.out.println(Z);
	}

	private static void setValidPredPos(MLN mln) {
		for(int i  = 0  ; i <= mln.max_predicate_id ; i++){
			mln.validPredPos.add(new ArrayList<Boolean>());
		}
		for(WClause clause : mln.clauses){
			clause.print();
			for(Atom atom : clause.atoms){
				int predId = atom.symbol.id;
				if(mln.validPredPos.get(predId).size() == 0){
					for(int i = atom.terms.size() - 1 ; i >= 0 ; i--){
						mln.validPredPos.get(predId).add(true);
					}
				}
			}
			System.out.println("HyperCubes : ");
			clause.root.isRoot = true;
		}
		
	}

	private static void createDummyLeaves(WClause clause) {
		Node n = new Node(LiftedPTP.maxPredId);
		n.rootParent = clause.root;
		clause.root.children.add(n);
		for(int i = 0 ; i < clause.root.hyperCubesList.size() ; i++){
			n.placeHolderList.add(i);
		}
		HashMap<Integer,Integer> predOccurHashMap = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> predArityMap = new HashMap<Integer,Integer>();
		for(Atom atom : clause.atoms){
			if(!predArityMap.containsKey(atom.symbol.id)){
				predArityMap.put(atom.symbol.id, atom.terms.size());
			}
			if(predOccurHashMap.containsKey(atom.symbol.id)){
				predOccurHashMap.put(atom.symbol.id,predOccurHashMap.get(atom.symbol.id)+1);
			}
			else{
				predOccurHashMap.put(atom.symbol.id,0);
			}
		}
		for(Integer predId : predOccurHashMap.keySet()){
			for(int i = 0 ; i <= predOccurHashMap.get(predId) ; i++){
				ArrayList<Integer> entry = new ArrayList<Integer>();
				for(int j = 0 ; j < predArityMap.get(predId) ; j++){
					entry.add(-1);
				}
				n.partialGroundings.get(predId).add(entry);
			}
		}
		
	}

}
