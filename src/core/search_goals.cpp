// ============================================================================
// ILOG Search Goals for Scheduling
// ============================================================================

#include <ilcp/cpext.h>
#include <cstring>

#include "search_goals.hpp"


using namespace std;

/**
 * Set times from permutation variables
 */
ILCGOAL2(SetTimesPermutation, IlcIntVarArray, permutation, IlcIntervalVarArray, acts) {
	for( int i = 0; i < permutation.getSize(); ++i ) {
		acts[permutation[i].getValue()].setStart(acts[permutation[i].getValue()].getStartMin());
	}
	return IlcGoalTrue(getCP());
}

/**
 * Wrapper for set times permutation
 **/
ILOCPGOALWRAPPER2(IloSetTimesPermutation, cp,
				  IloIntVarArray, _permutation,
				  IloIntervalVarArray, _tasks)
{
    IlcIntervalVarArray tasks(cp, _tasks.getSize());
    for( int i = 0; i < _tasks.getSize(); i++ )
        tasks[i] = cp.getInterval(_tasks[i]);
    return SetTimesPermutation(cp, cp.getIntVarArray(_permutation), tasks);
}


/**
 * Fix variable to its minimum value
 */
ILCGOAL1(FixIntVarMin, IlcIntVar, var) {
    var.setValue(var.getMin());
    return IlcGoalTrue(getCP());
}

/**
 * Wrapper to fix variable to its minimum value
 **/
ILOCPGOALWRAPPER1(IloFixIntVarMin, cp, IloIntVar, var)
{
    return FixIntVarMin(cp, cp.getIntVar(var));
}


/**
 * Set times of a sequence var
 */
ILCGOAL1(SetTimesSequence, IlcIntervalSequenceVar, seq) {
    // fix start times
    for( IlcIntervalVar act = seq.getEarliestInHead(); act.getImpl() != NULL; act = seq.getOneLaterInHead(act) ) {
    	act.setStart(act.getStartMin());
    }
    for( IlcIntervalVar act = seq.getEarliestInTail(); act.getImpl() != NULL; act = seq.getOneLaterInTail(act) ) {
    	act.setStart(act.getStartMin());
    }
	return IlcGoalTrue(getCP());
}

/**
 * Wrapper for set times sequence
 **/
ILOCPGOALWRAPPER1(IloSetTimesSequence, cp, IloIntervalSequenceVar, _seq)
{
    return SetTimesSequence(cp, cp.getIntervalSequence(_seq));
}


/**
 * Set precedences for SOP problems
 */
ILCGOAL3(SetPrecedences, IlcIntervalVarArray, tasks, IlcIntervalSequenceVar, seq, bool**, precedes) {
	for( int i = 0; i < tasks.getSize(); i++ ) {
		for( int j = 0; j < tasks.getSize(); j++ ) {
			if( precedes[i][j] ) {
				seq.setBefore(tasks[i], tasks[j]);
			}
		}
	}
    return IlcGoalTrue(getCP());
}

/**
 * Wrapper for SOP precedences
 */
ILOCPGOALWRAPPER3(IloSetPrecedences, cp, IloIntervalVarArray, _tasks, IloIntervalSequenceVar, _seq, bool**, precedes)
{
    IlcIntervalVarArray tasks(cp, _tasks.getSize());
    for( int i = 0; i < _tasks.getSize(); i++ )
        tasks[i] = cp.getInterval(_tasks[i]);
    return SetPrecedences(cp, tasks, cp.getIntervalSequence(_seq), precedes);
}


/**
 * Implementation goal for lexicographic search
 */
ILCGOAL2(LexSchedule, IlcIntervalVarArray, tasks, IlcIntervalSequenceVar, seq) {
    IloCP cp = getCP();
    for( int i = 0; i < tasks.getSize(); i++ ) {
        if( seq.isCandidateHead(tasks[i]) ) {
            return IlcAnd(seq.tryExtendHead(tasks[i]), LexSchedule(cp, tasks, seq));
        }
    }
    return IlcGoalTrue(cp);
}


/**
 * Set precedences for SOP problems
 */
ILCGOAL1(SetMinStartTime, IlcIntervalVarArray, tasks) {
	for( int i = 0; i < tasks.getSize(); i++ ) {
		tasks[i].setStartMax(tasks[i].getStartMin());
	}
	return IlcGoalTrue(getCP());
}

/**
 * Wrapper for set min start times
 */
ILOCPGOALWRAPPER1(IloSetMinStartTime, cp, IloIntervalVarArray, _tasks)
{
	IlcIntervalVarArray tasks(cp, _tasks.getSize());
	for( int i = 0; i < _tasks.getSize(); i++ )
		tasks[i] = cp.getInterval(_tasks[i]);
	return SetMinStartTime(cp, tasks);
}


/**
 * Wrapper for lexicographic search
 */
ILOCPGOALWRAPPER2(IloLexSchedule, cp, IloIntervalVarArray, _tasks, IloIntervalSequenceVar, _seq)
{
    IlcIntervalVarArray tasks(cp, _tasks.getSize());
    for( int i = 0; i < _tasks.getSize(); i++ )
        tasks[i] = cp.getInterval(_tasks[i]);
    return LexSchedule(cp, tasks, cp.getIntervalSequence(_seq));
}


/**
 * Goal for testing
 */
ILCGOAL1(IlcTest, char*, message) {
    cout << message << endl;
    return IlcGoalTrue(getCP());
}


/**
 * Sequencing based on activity belonging to shortest path in the MDD
 * TODO: little optimization: also save last fixed position
 */
ILCGOAL2(ShortestPathMDD, IlcIntVarArray, permutation, DisjunctivePropagator*, disj_propagator) {

  // search first unfixed activity in M
  int unfixed = -1;
  for( int i = 0; i < permutation.getSize() && unfixed == -1; ++i ) {
    if( !permutation[i].isFixed() ) {
      unfixed = i;
    }
  }
  if( unfixed == -1 ) {
    // activities are fixed: return
    return IlcGoalTrue(getCP());
  }

  // get MDD shortest path
  IloCP cp = getCP();
  int* optimum = disj_propagator->compute_shortest_path();
  return( IlcOr(IlcAnd(permutation[unfixed] == optimum[unfixed], ShortestPathMDD(cp, permutation, disj_propagator)),
                IlcAnd(permutation[unfixed] != optimum[unfixed], ShortestPathMDD(cp, permutation, disj_propagator))) );
}


/**
 * Wrapper: Sequencing based on activity belonging to shortest path in the MDD
 */
ILOCPGOALWRAPPER2(IloShortestPathMDD, cp, IloIntVarArray, permutation, DisjunctivePropagator*, disj_propagator) {
  return ShortestPathMDD(cp, cp.getIntVarArray(permutation), disj_propagator);
}





/**
 * Wrapper: Sequencing based on relaxed sum of setup times in MDD
 */
ILCGOAL2(SumSetupTimesMDD, IlcIntVarArray, permutation, DisjunctivePropagator*, disj_propagator) {

	IloCP cp = getCP();

	int pos, act;
	disj_propagator->get_activity_min_sum_setup_T(pos, act);
	if( pos == -1 ) {
		return IlcGoalTrue(cp);
	}

	//cout << "Trying " << pos << " --> " << act << endl;

	return( IlcOr(IlcAnd(permutation[pos] == act, SumSetupTimesMDD(cp, permutation, disj_propagator)),
			IlcAnd(permutation[pos] != act, SumSetupTimesMDD(cp, permutation, disj_propagator))) );
}



/**
 * Wrapper: Sequencing based on relaxed sum of setup times in MDD
 */
ILOCPGOALWRAPPER2(IloSumSetupTimesMDD, cp, IloIntVarArray, permutation, DisjunctivePropagator*, disj_propagator) {
  return SumSetupTimesMDD(cp, cp.getIntVarArray(permutation), disj_propagator);
}



/**
 * Fix node of the MDD, removing all other nodes of the corresponding layer
 */
ILCGOAL3(FixMDDNode, DisjunctivePropagator*, disj_propagator, int, layer, int, node) {
	IloCP cp = getCP();
	disj_propagator->fix_node(layer, node);
	return IlcGoalTrue(cp);
}

/**
 * Fix node of the MDD, removing all other nodes of the corresponding layer
 */
ILCGOAL3(RemoveMDDNode, DisjunctivePropagator*, disj_propagator, int, layer, int, node) {
	IloCP cp = getCP();
	disj_propagator->remove_node(layer, node);
	return IlcGoalTrue(cp);
}


/**
 * Wrapper: Sequencing based on relaxed sum of setup times in MDD
 */
ILCGOAL2(MDDNodeBranch, IlcIntVarArray, permutation, DisjunctivePropagator*, disj_propagator) {

	IloCP cp = getCP();

	int layer, node;
	disj_propagator->get_node_min_sum_setup_T(layer, node);
	if( layer == -1 ) {
		return IlcGoalTrue(cp);
	}

	//cout << "Trying " << pos << " --> " << act << endl;

	return( IlcOr(IlcAnd(FixMDDNode(cp, disj_propagator, layer, node), MDDNodeBranch(cp, permutation, disj_propagator)),
			IlcAnd(RemoveMDDNode(cp, disj_propagator, layer, node), MDDNodeBranch(cp, permutation, disj_propagator))) );
}



/**
 * Wrapper: Sequencing based on relaxed sum of setup times in MDD
 */
ILOCPGOALWRAPPER2(IloMDDNodeBranch, cp, IloIntVarArray, permutation, DisjunctivePropagator*, disj_propagator) {
  return MDDNodeBranch(cp, cp.getIntVarArray(permutation), disj_propagator);
}




/**
 * Goal for lexicographic search
 */
ILOCPGOALWRAPPER1(IloTest, cp, char*, message) {
	char* new_message = new(cp.getHeap()) char[256];
	strcpy(new_message, message);
    return IlcTest(cp, new_message);
}



