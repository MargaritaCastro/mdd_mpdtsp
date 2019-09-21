// ============================================================================
// ILOG Search Goals for Scheduling
// ============================================================================

#ifndef SEARCH_GOALS_HPP_
#define SEARCH_GOALS_HPP_

#include <ilcp/cpext.h>
#include "../constraints/disjunctive.hpp"

/**
 * Lexicographic sequencing of activities
 */
IloGoal IloLexSchedule(IloEnv env, IloIntervalVarArray activities, IloIntervalSequenceVar sequence);

/**
 * Set precedences for SOP problem
 */
IloGoal IloSetPrecedences(IloEnv env, IloIntervalVarArray activities, IloIntervalSequenceVar sequence, bool** precedes);


/**
 * Lexicographic sequencing of activities
 */
IloGoal IloLexSchedule(IloEnv env, IloIntervalVarArray activities, IloIntervalSequenceVar sequence);


/**
 * Set times of a disjunctive resource
 */
IloGoal IloSetTimesSequence(IloEnv env, IloIntervalSequenceVar seq);

/**
 * Set min start times for intervals
 */
IloGoal IloSetMinStartTime(IloEnv env, IloIntervalVarArray intervals);


/**
 * Set times of a disjunctive resource using permutation variables
 */
IloGoal IloSetTimesPermutation(IloEnv env, IloIntVarArray permutation, IloIntervalVarArray acts);


/**
 * Fix variable to its minimum value
 */
IloGoal IloFixIntVarMin(IloEnv env, IloIntVar var);

/**
 * Sequencing based on activity belonging to shortest path in the MDD
 */
IloGoal IloShortestPathMDD(IloEnv env, IloIntVarArray permutation, DisjunctivePropagator* disj_propagator);


/**
 * Sequencing based on sum of setup times
 */
IloGoal IloSumSetupTimesMDD(IloEnv env, IloIntVarArray permutation, DisjunctivePropagator* disj_propagator);


/**
 * Sequencing based on sum of setup times
 */
IloGoal IloMDDNodeBranch(IloEnv env, IloIntVarArray permutation, DisjunctivePropagator* disj_propagator);


#endif /* SEARCH_GOALS_HPP_ */
