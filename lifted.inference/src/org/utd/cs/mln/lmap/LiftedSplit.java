package org.utd.cs.mln.lmap;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.utd.cs.mln.alchemy.core.Atom;
import org.utd.cs.mln.alchemy.core.Evidence;
import org.utd.cs.mln.alchemy.core.HyperCube;
import org.utd.cs.mln.alchemy.core.MLN;
import org.utd.cs.mln.alchemy.core.PredicateNotFound;
import org.utd.cs.mln.alchemy.core.PredicateSymbol;
import org.utd.cs.mln.alchemy.core.Term;
import org.utd.cs.mln.alchemy.core.WClause;
import org.utd.cs.mln.alchemy.util.Parser;
import org.utd.cs.mln.alchemy.util.comb;

public class LiftedSplit {
	
	// Task : Given an MLN, it splits (in lifted manner) an atom and returns Z
	// Input : MLN mln
	// Output : A Double number indicating Z of input MLN
	
	public static double liftedSplit(MLN mln){
		
		// Find best pred to split on
		int predId = choosePredToSplitOn(mln);
		
		// See if predId is singleton or not, if it is singleton, then apply binomial, else split on all groundings
		boolean applyBinomial = false;
		if(isSingleton(predId)){
			applyBinomial = true;
		}
		
		// If applyBinomial is false, then basis segments of this predId are just all groundings, so just choose
		// any first grounding found in MLN, and split on that
		if(applyBinomial == false){
			// Ground the mln wrt predicate predId
			createGroundMlnPred(mln, predId);
			return nonBinomialLiftedSplitter(mln, predId);
		}
		// Else create Disjoint segments of unary predicate predId and apply binomial on them
		else{
			// We also need some helper data structures which will help in faster execution of binomial rule. These data structures are : 
			HashMap<HyperCube,Set<Integer>> hyperCubeClauseListHashMap = new HashMap<HyperCube,Set<Integer>>(); // For each segment we have of pred, we store set of clauseIds in which that segment appears
			ArrayList<HashMap<Integer,Integer>> clausesTermFreqCountList = new ArrayList<HashMap<Integer,Integer>>(); // For each clause, it stores frequency of each term in that clause. HashMap key : termId, value : frequency
			HashMap<Integer,ArrayList<Integer>> clausesBinomialTermIndices = new HashMap<Integer,ArrayList<Integer>>(); // For each clauseId, it stores list of terms which are arguments of predicate predId
			HashMap<Integer,ArrayList<Boolean>> clausesBinomialPredSigns = new HashMap<Integer,ArrayList<Boolean>>(); // For each clauseId, it stores list of signs for each occurrence of predId. Sign = True means predId came with negative sign, else positive.
			// Note that for each clause c, size(binomialTermIndices) = size(binomialPredSigns) because binomial is applied only on unary preds and we are not allowing S(x) | S(x) type, neither S(x) | !S(x) type
			
			// Now fill in above 4 data structures
			for(int clauseId = 0 ; clauseId < mln.clauses.size() ; clauseId++){
				WClause clause = mln.clauses.get(clauseId);
				System.out.println("Clause HyperCube : " + clause.hyperCubes);
				ArrayList<Integer> binomialClauseTermIndices = new ArrayList<Integer>();
				ArrayList<Boolean> binomialClausePredSigns = new ArrayList<Boolean>();
				HashMap<Integer,Integer> termFreqCount = new HashMap<Integer,Integer>();
				for(Atom atom : clause.atoms){
					for(Term term : atom.terms){
						if(termFreqCount.containsKey(clause.terms.indexOf(term))){
							termFreqCount.put(clause.terms.indexOf(term), termFreqCount.get(term)+1);
						}
						else{
							termFreqCount.put(clause.terms.indexOf(term), 1);
						}
					}
					if(predId == atom.symbol.id){
						binomialClauseTermIndices.add(clause.terms.indexOf(atom.terms.get(0))); // We need to get only first term because binomial atom is unary
						binomialClausePredSigns.add(clause.sign.get(clause.atoms.indexOf(atom)));
					}
				}
				clausesTermFreqCountList.add(termFreqCount);
				// If there is no binomial atom in this clause, then no need to go further
				if(binomialClauseTermIndices.size() == 0)
					continue;
				
				clausesBinomialTermIndices.put(clauseId, binomialClauseTermIndices);
				clausesBinomialPredSigns.put(clauseId, binomialClausePredSigns);
				System.out.println("termFreqCount : " + termFreqCount.toString());
				System.out.println("binomialClauseTermIndices : " + binomialClauseTermIndices.toString());
				System.out.println("binomialClausePredSigns : " + binomialClausePredSigns.toString());
				
				for(HyperCube hyperCube : clause.hyperCubes){
					for(Integer binomialTermIndex : binomialClauseTermIndices){
						HyperCube projectedHyperCube = new HyperCube(hyperCube.varConstants.get(binomialTermIndex));
						if(hyperCubeClauseListHashMap.containsKey(projectedHyperCube)){
							hyperCubeClauseListHashMap.get(projectedHyperCube).add(clauseId);
						}
						else{
							hyperCubeClauseListHashMap.put(projectedHyperCube, new HashSet<Integer>(Arrays.asList(clauseId)));
						}
					}
				}
			}
			
			/* Now we have to create disjoint segments of predId, which can be done using function createDisjointHyperCubes.
			 * That function takes 3 arguments : (1) list of hyperCubes, which we will store in data structure segments, 
			 * (2) segmentClauseList, which stores list of clauseIds in which each segment occurs,
			 * (3) hyperCubeClauseListHashMap, which will be updated according to disjoint segments
			 */
			ArrayList<HyperCube> segments = new ArrayList<HyperCube>();
			ArrayList<Set<Integer>> segmentClauseList = new ArrayList<Set<Integer>>();
			for(HyperCube hyperCube : hyperCubeClauseListHashMap.keySet()){
				segments.add(hyperCube);
				segmentClauseList.add(hyperCubeClauseListHashMap.get(hyperCube));
			}
			//System.out.println("HyperCubes list : " + segments);
			//System.out.println("ClauseList : " + segmentClauseList);
			//System.out.println("Final hyperCubes size : " + segments.size());
			Decomposer.createDisjointHyperCubes(segments, segmentClauseList, hyperCubeClauseListHashMap);
			
			// Finally, apply binomial rule on all these disjoint segments of predId. Also pass those 4 data structures we created earlier
			return binomialSplitter(mln, predId, hyperCubeClauseListHashMap, clausesTermFreqCountList, clausesBinomialTermIndices, clausesBinomialPredSigns);
		}
	}

	/* Task : Apply binomial rule on singleton predicate with id predId. From now on, this predicate will be referred as pred.
	 * Inputs : 
	 * (1) mln : MLN in which binomial rule is to be applied
	 * (2) predId : Id of singleton predicate on which binomial rule is to be applied
	 * (3) hyperCubeClauseListHashMap : For each segment of pred, stores set of clauseIds in which that segment appears
	 * (4) clausesTermFreqCountList : For each clause, stores frequency of each term appearing in that clause. HashMap key : termId, value : freq
	 * (5) clausesBinomialTermIndices : For a given clauseId, stores list of termIds which appear as arguments of pred in that clause  
	 * (6) clausesBinomialPredSigns : For a given clauseId, stores list of signs for each appearence of pred in that clause. Sign = true means pred came with negation in clause
	 * Output : Z of input mln. After applying binomial, this function calls liftedPartition, which then gives Z.
	 */
	private static double binomialSplitter(MLN mln, int predId,
			HashMap<HyperCube, Set<Integer>> hyperCubeClauseListHashMap,
			int segmentIndex,
			ArrayList<HashMap<Integer, Integer>> clausesTermFreqCountList,
			HashMap<Integer, ArrayList<Integer>> clausesBinomialTermIndices,
			HashMap<Integer, ArrayList<Boolean>> clausesBinomialPredSigns) {
		
		/* If hyperClauseListHashMap is empty, that means we have applied binomial on all segments, and thus we are done.
		 * Also, remove all occurrences of pred from all clauses now, and also remove all terms which were appearing only in pred.
		 * Also remove varConstants corresponding to those removed terms from each hyperCube 
		 */
		if(segmentIndex == hyperCubeClauseListHashMap.size()){
			for(WClause clause : mln.clauses){
				// First remove atoms with id = predId
				for(int atomId = clause.atoms.size()-1 ; atomId >= 0 ; atomId--){
					// If this atom is pred, then remove this atom
					if(clause.atoms.get(atomId).symbol.id == predId){
						clause.atoms.remove(atomId);
					}
				}
				
				// Now remove terms.
				// TODO : There must be a better method to know which terms to remove
				ArrayList<Integer> toRemoveTermIndices = new ArrayList<Integer>();
				if(clause.hyperCubes.size() > 0){
					HyperCube firstHyperCube = clause.hyperCubes.get(0);
					for(int termId =  firstHyperCube.varConstants.size() - 1; termId >= 0 ; termId--){
						if(firstHyperCube.varConstants.get(termId).isEmpty()){
							toRemoveTermIndices.add(termId);
						}
					}
					for(Integer termId : toRemoveTermIndices){
						clause.terms.remove(termId);
					}
					
					for(HyperCube hyperCube : clause.hyperCubes){
						for(Integer termId : toRemoveTermIndices){
							hyperCube.varConstants.remove(termId);
							
						}
					}
				}
			}
			// Finally call liftedPartition on this mln
			return LiftedPTP.getPartition(mln);
		} // end of base case i.e. when binomial has been applied on all segments
		
		/* Now we have to apply binomial on a segment. We can choose any segment, lets just choose first segment from hyperCubeClauseListHashMap.
		 * Also after choosing, delete that segment from this Map
		 * 
		 */
		HyperCube segment = null; // Segment on which we want to apply binomial
		Set<Integer> keyClauseList = null; // Set of clauseIds in which above segment appears
		int hyperCubeId = 0;
		for(HyperCube hyperCube : hyperCubeClauseListHashMap.keySet()){
			if(hyperCubeId > segmentIndex)
			{
				segment = new HyperCube(hyperCube);
				keyClauseList = new HashSet<Integer>(hyperCubeClauseListHashMap.get(hyperCube));
				segmentIndex = hyperCubeId;
				//hyperCubeClauseListHashMap.remove(hyperCube);
				break;
			}
			hyperCubeId++;
		}
		
		// Now when we apply binomial, we will get two things : binomialMLNs, which will be list of n+1 MLNs, where n is segment size, and binomialCoeff, which will be list of n+1 corresponding coefficients to be multiplied 
		ArrayList<MLN> binomialCasesMLNs = new ArrayList<MLN>();
		ArrayList<Double> binomialMLNCoeff = new ArrayList<Double>();
		binomialSplitterSegment(mln, predId, segment, keyClauseList, clausesTermFreqCountList, clausesBinomialTermIndices, clausesBinomialPredSigns, binomialCasesMLNs, binomialMLNCoeff);
		double Z = 0.0;
		for(int mlnId = 0 ; mlnId < binomialCasesMLNs.size() ; mlnId++){
			MLN binomialMln = binomialCasesMLNs.get(mlnId);
			Z += binomialMLNCoeff.get(mlnId)*binomialSplitter(binomialMln, predId, hyperCubeClauseListHashMap, segmentIndex, clausesTermFreqCountList, clausesBinomialTermIndices, clausesBinomialPredSigns);
		}
		return Z;
	}

	// TODO : Handle segments which could get lost due to tautology of hyperCube
	/**
	 * Task : Apply binomial rule on an MLN wrt to singleton predicate with id = predId, and a segment.
	 * @param parentMln
	 * @param predId
	 * @param segment
	 * @param keyClauseList
	 * @param clausesTermFreqCountList
	 * @param clausesBinomialTermIndices
	 * @param clausesBinomialPredSigns
	 * @param binomialCasesMLNs
	 * @param binomialMLNCoeff
	 */
	private static void binomialSplitterSegment(MLN parentMln, int predId,
			HyperCube segment, Set<Integer> keyClauseList,
			ArrayList<HashMap<Integer, Integer>> clausesTermFreqCountList,
			HashMap<Integer, ArrayList<Integer>> clausesBinomialTermIndices,
			HashMap<Integer, ArrayList<Boolean>> clausesBinomialPredSigns,
			ArrayList<MLN> binomialCasesMLNs, ArrayList<Double> binomialMLNCoeff) {

		// First, update the hyperCubes of each clause such that they become disjoint (or identical) wrt to param segment
		/* Extract first element of segment, we need only first element to check if param segment is same as other 
		 * segment of pred, because this param segment belongs to basis of pred's segments
		 */
		int segmentFirstElement = segment.varConstants.get(0).iterator().next();
		// Extract the set from hyperCube param segment. We know size of this hyperCube is 1 because it is unary pred
		Set<Integer> segmentSet = new HashSet<Integer>(segment.varConstants.get(0));
		
		// Create disjoint or identical segments wrt to param segment of all occurrences of a binomial atom in all clauses
		// For any two intersecting segments A and B, create B-A and A
		
		for(Integer clauseId : keyClauseList){
			System.out.println("In clauseId : " + clauseId);
			WClause clause = parentMln.clauses.get(clauseId);
			for(Integer termIndex : clausesBinomialTermIndices.get(clauseId)){
				int origNumHyperCubes = clause.hyperCubes.size();
				//System.out.println("In segment Index : " + termIndex);
				//System.out.println("Num of hyperCubes : " + origNumHyperCubes);
				for(int hcId = 0 ; hcId < origNumHyperCubes ; hcId++){
					if(clause.hyperCubes.get(hcId).varConstants.get(termIndex).contains(segmentFirstElement)){
						clause.hyperCubes.get(hcId).varConstants.get(termIndex).removeAll(segmentSet); // B-A
						if(!clause.hyperCubes.get(hcId).isEmpty())
							clause.hyperCubes.add(new HyperCube(clause.hyperCubes.get(hcId)));
						clause.hyperCubes.get(hcId).varConstants.set(termIndex, new HashSet<Integer>(segmentSet)); // A
					}
				}
				System.out.println("Now clause HyperCubes are : " + clause.hyperCubes);
			}
		}

		// Now loop over n+1 cases, where n is size of param segment
		int n = segmentSet.size();
		for(int k = 0 ; k <= n ; k++){
			MLN newMln = new MLN(); // create an empty mln
			double weightSatisfiedClauses = 0.0; // coeff to be multiplied
	
			// Go over all clauses. Clauses in which param segment doesn't appear are added as it is, and apply binomial on clauses in which param segment occurs
			for(int clauseId = 0 ; clauseId < parentMln.clauses.size() ; clauseId++){
				WClause newClause = MLN.create_new_clause(parentMln.clauses.get(clauseId));
				if(!keyClauseList.contains(clauseId)){
					newMln.clauses.add(newClause);
					continue;
				}
				// For current clause, extract freq count of each term
				HashMap<Integer,Integer> termFreqCount = clausesTermFreqCountList.get(clauseId);
				// For current clause, extract indices of terms which are arguments of pred
				ArrayList<Integer> binomialTermIndices = new ArrayList<Integer>(clausesBinomialTermIndices.get(clauseId));
				// For current clause, extract signs of all occurrences of pred in this clause. Note that size pf
				// binomialTermIndices and binomialPredSigns is same because for each occurrence of pred in a clause
				// different term appears i.e. no S(x) V S(x) type is allowed
				ArrayList<Boolean> binomialPredSigns = new ArrayList<Boolean>(clausesBinomialPredSigns.get(clauseId));
				
				
				// Go over each hyperCube of this clause and apply binomial. We are going from last hyperCube, because some hyperCube may get deleted
				for(int hcId = newClause.hyperCubes.size() - 1 ; hcId >= 0 ; hcId--){
					
					// For each occurrence of pred in clause, stores whether its term has same segment as param segment ?
					// So size(binomialTermIndices) = size(binomialPredSigns) = size(sameSegmentOccurs)					
					ArrayList<Boolean> sameSegmentOccurs = new ArrayList<Boolean>();
					// This stores the termIds which have same segments as param segment. This list will ease the calculation
					ArrayList<Integer> sameSegmentOccurTermList = new ArrayList<Integer>();
					// Now fill in sameSegmentOccurs
					for(Integer binomialTermIndex : binomialTermIndices){
						if(newUnSatHyperCube.varConstants.get((int)binomialTermIndex).contains(segmentFirstElement)){
							//numSameSegments++;
							sameSegmentOccurs.add(true);
							sameSegmentOccurTermList.add(binomialTermIndex);
						}
						else{
							sameSegmentOccurs.add(false);
						}
					}
					
					/*
					 *  Now we have to count number of true groundings of this hyperCube. It is easier to calculate number of false groundings of
					 *  a hyperCube.
					 *  NumFalseGroundings = {{for each term in sameSegment}prod(k if sign = true, n-k if sign = false}} * {for remaining terms}prod(segmentSize)
					 *  Then numTrueGroundings = prod(segmentSize of all terms) - numFalseGroundings
					 *  Also, we have to reduce domain of terms which occur in sameSegment but occur somewhere else also. 
					 */
					
					//int numFalseGroundings = 1;
					//int numTotalGroundings = 1;
					ArrayList<ArrayList<Set<Integer>>> sameSegmentsPartitions = new ArrayList<ArrayList<Set<Integer>>>();
					for(int i = 0 ; i < sameSegmentOccurs.size() ; i++){
						if(sameSegmentOccurs.get(i) == true){
							Set<Integer> falsePartiton = new HashSet<Integer>(segment.varConstants.get(0));
							Set<Integer> truePartiton = new HashSet<Integer>(segment.varConstants.get(0));
							if(binomialPredSigns.get(i) == false){
								for(int d = 0 ; d < k ; d++){
									falsePartiton.remove(falsePartiton.iterator().next());
								}
								truePartiton.removeAll(falsePartiton);
							}
							else{
								for(int d = 0 ; d < k ; d++){
									truePartiton.remove(truePartiton.iterator().next());
								}
								falsePartiton.removeAll(truePartiton);
							}
							ArrayList<Set<Integer>> partitions = new ArrayList<Set<Integer>>();
							partitions.add(falsePartiton);
							partitions.add(truePartiton);
							sameSegmentsPartitions.add(partitions);
						}
					}
					//HyperCube newUnSatHyperCube = newClause.hyperCubes.get(hcId); // Extract current hyperCube
					HyperCube curHyperCube = newClause.hyperCubes.get(hcId); // Extract current hyperCube
					ArrayList<HyperCube> newSatHyperCubes = new ArrayList<HyperCube>();
					int numHyperCubes = (int)Math.pow(2, sameSegmentOccurTermList.size());
					int numBits = sameSegmentOccurTermList.size();
					for(int c = 0 ; c < numHyperCubes ; c++){
						int temp = c;
						HyperCube h = new HyperCube(curHyperCube);
						for(int bitNum = 0 ; bitNum < numBits ; bitNum++){
							int termId = sameSegmentOccurTermList.get(bitNum);
							h.varConstants.set(termId, sameSegmentsPartitions.get(bitNum).get(temp%2));
							h.num_copies *= sameSegmentsPartitions.get(bitNum).get(temp%2).size();
							temp = temp/2;
						}
						if(c > 0){
							h.satisfied = true;
						}
						if(h.num_copies > 0){
							newSatHyperCubes.add(h);
						}
					} // end of c loop
					
					newClause.hyperCubes.remove(hcId);
//					for(int termId = 0 ; termId < newClause.terms.size() ; termId++){
//						//int segmentSize = newHyperCube.varConstants.get(termId).size();
//						if(sameSegmentOccurTermList.contains(termId)){
//							int binomialPredIndex = binomialTermIndices.indexOf(termId);
//							/*
//							if(binomialPredSigns.get(binomialPredIndex) == true){
//								numFalseGroundings *= k;
//							}
//							else{
//								numFalseGroundings *= (n-k);
//							}*/
//							// Now reduce the domain size appropriately
//							// It means if this term appeared here only, then empty this term
//							if(termFreqCount.get(termId) == 1){
//								newUnSatHyperCube.varConstants.get(termId).clear();
//								newSatHyperCube.varConstants.get(termId).clear();
//							}
//							else{
//								int domainReduceAmount = k;
//								if(binomialPredSigns.get(binomialPredIndex) == true){
//									domainReduceAmount = n-k;
//								}
//								newUnSatHyperCube.num_copies *= n-domainReduceAmount;
//								newSatHyperCube.num_copies *= domainReduceAmount;
//								for(int d = 0 ; d < domainReduceAmount ; d++){
//									Integer elemToRemove = newUnSatHyperCube.varConstants.get(termId).iterator().next();
//									newUnSatHyperCube.varConstants.get(termId).remove(elemToRemove);
//								}
//								for(int d = 0 ; d < n-domainReduceAmount ; d++){
//									Integer elemToRemove = newSatHyperCube.varConstants.get(termId).iterator().next();
//									newSatHyperCube.varConstants.get(termId).remove(elemToRemove);
//								}
//							}
//						}
//						/*
//						else{
//							// If this term was of pred and its segment size is 0, then this term has already been processed earlier
//							if(segmentSize == 0 && binomialTermIndices.contains(termId)){
//								segmentSize = 1;
//							}
//							numFalseGroundings *= segmentSize;
//						}
//						numTotalGroundings *= segmentSize;*/
//					} // termId loop ends here
					//int numTrueGroundings = numTotalGroundings - numFalseGroundings;
					//weightSatisfiedClauses += newClause.weight.getValue()*numTrueGroundings;
					/*
					if(numTrueGroundings == numTotalGroundings){
						newClause.hyperCubes.remove(hcId);
					}*/
					/*
					if(newUnSatHyperCube.num_copies == 0){
						newClause.hyperCubes.remove(hcId);
					}
					if(newSatHyperCube.num_copies > 0){
						newClause.hyperCubes.add(newSatHyperCube);
					}*/
					
				} // hyperCube loop ends here
				
			}// clauses loop ends here
			
			binomialCasesMLNs.add(newMln);
			binomialMLNCoeff.add(comb.findComb(n, k));
		}// end of k loop
	}

	private static double nonBinomialLiftedSplitter(MLN mln, int predId) {
		// TODO Auto-generated method stub
		return 0;
	}

	private static boolean isSingleton(int predId) {
		// TODO Auto-generated method stub
		return false;
	}

	private static int choosePredToSplitOn(MLN mln) {
		// TODO Auto-generated method stub
		return 0;
	}

	public static ArrayList<ArrayList<Set<Integer>>> cartesianProdVar(ArrayList<Set<Integer>> sets, Set<Integer> toGroundTermIdsSet){
		ArrayList<ArrayList<Set<Integer>>> result = new ArrayList<ArrayList<Set<Integer>>>();
		int numTuples = 1;
		for(Integer i : toGroundTermIdsSet){
			numTuples *= sets.get(i).size();
		}
		
		// add numTuples lists into result
		for(int i = 0 ; i < numTuples ; i++){
			result.add(new ArrayList<Set<Integer>>());
		}
		
		int numRepeats = numTuples;
		// Now fill each element of input sets
		for(int i = 0 ; i < sets.size() ; i++){
			// If this index's elements are to be remained as it is, then just copy them numTuples times 
			if(!toGroundTermIdsSet.contains(i)){
				for(int j = 0 ; j < numTuples ; j++){
					result.get(j).add(new HashSet<Integer>(sets.get(i)));
				}
			}
			// else create the set of single integers
			else{
				numRepeats = numRepeats/sets.get(i).size();
				int numFilledEntries = 0;
				while(numFilledEntries != numTuples){
					for(Integer elem : sets.get(i)){
						for(int j = 0 ; j < numRepeats ; j++){
							result.get(numFilledEntries).add(new HashSet<Integer>(Arrays.asList(elem)));
							numFilledEntries++;
						}
					}
				}
			}
		}
		return result;
	}
	
	// Create ground mln relative to predicate with id PredId
	public static void createGroundMlnPred(MLN mln, int predId){
		for(WClause clause : mln.clauses){
			// If predId appears in this clause, then only ground this clause, otherwise leave it as it is
			// Also, if it appears, then note the termIds which are arguments of this predId, so that we can
			// ground all hyperCubes on those termIds
			Set<Integer> toGroundTermIdsSet = new HashSet<Integer>();
			for(Atom atom : clause.atoms){
				if(atom.symbol.id == predId){
					for(Term term : atom.terms){
						toGroundTermIdsSet.add(clause.terms.indexOf(term));
					}
				}
			}
			// If there was predId in this clause i.e. if toGroundTermIdsSet is not empty, then ground each hyperCube
			// Append returned list of hyperCubes to hyperCubes list of this clause, also delete this hyperCube
			// So start from last hyperCube  so that indexing doesn't get damaged
			
			if(toGroundTermIdsSet.size() > 0){
				int origNumHyperCubes = clause.hyperCubes.size();
				for(int hcId = origNumHyperCubes - 1 ; hcId >= 0 ; hcId--){
					HyperCube hc = clause.hyperCubes.get(hcId);
					ArrayList<ArrayList<Set<Integer>>> groundedHyperCubes = cartesianProdVar(hc.varConstants, toGroundTermIdsSet);
					for(ArrayList<Set<Integer>> groundedVarConstants : groundedHyperCubes){
						HyperCube groundedHc = new HyperCube();
						groundedHc.varConstants = groundedVarConstants;
						clause.hyperCubes.add(groundedHc);
					}
					clause.hyperCubes.remove(hcId);
				}
			}
		}
	}
	
	/**
	 * @param args
	 * @throws FileNotFoundException 
	 * @throws PredicateNotFound 
	 */
	public static void main(String[] args) throws FileNotFoundException, PredicateNotFound {
		// TODO Auto-generated method stub
		MLN mln = new MLN();
		Parser parser = new Parser(mln);
		String filename = new String("smoke/smoke_mln.txt");
		parser.parseInputMLNFile(filename);

		ArrayList<Evidence> evidList = parser.parseInputEvidenceFile("smoke/evidence.txt");
		//ArrayList<Evidence> evidList = parser.parseInputEvidenceFile("entity_resolution/er-test-eclipse.db");
		MlnToHyperCube mlnToHyperCube = new MlnToHyperCube();
		HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap = mlnToHyperCube.createPredsHyperCube(evidList,mln);
	
		int origNumClauses = mln.clauses.size();
		for(int clauseId = 0 ; clauseId < origNumClauses ; clauseId++){
			mln.clauses.addAll(mlnToHyperCube.createClauseHyperCube(mln.clauses.get(clauseId), predsHyperCubeHashMap));
		}
		for(int clauseId = origNumClauses-1 ; clauseId >= 0 ; clauseId--){
			mln.clauses.remove(clauseId);
		}
		
		System.out.println("Printing initial set of hyperCubes for each clause");
		for(WClause clause : mln.clauses){
			clause.print();
			System.out.println("HyperCubes : ");
			for(HyperCube hc : clause.hyperCubes){
				System.out.println(hc);
			}
		}
		
		// Create ground mln wrt pred Smokes i.e. predId 0
		createGroundMlnPred(mln, 0);
		System.out.println("Printing final set of hyperCubes for each clause grounded wrt to S");
		for(WClause clause : mln.clauses){
			clause.print();
			System.out.println("HyperCubes : ");
			for(HyperCube hc : clause.hyperCubes){
				System.out.println(hc);
			}
		}
		
	}

}
