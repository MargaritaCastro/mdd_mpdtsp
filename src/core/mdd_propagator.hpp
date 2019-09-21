//============================================================================
// MDD Constraint Definition for Propagators
//============================================================================


#ifndef MDD_CONSTRAINT_HPP_
#define MDD_CONSTRAINT_HPP_


#include <ilcp/cpext.h>
#include <vector>
#include "../mdd/mdd.hpp"

using namespace std;

/** Forward reference */
class IlcMDDPropagatorImp;

/**
 * Group for refinement
 */
typedef SimpleList<RevMDD::Arc*>* RefinementGroup;


/**
 * MDD Propagator - Definition Class
 */
class MDDPropagator {

  friend class IlcMDDPropagatorImp;
  friend class IlcMDDPropagatorSearch;

 public:

  virtual ~MDDPropagator() { }

  // Initialize constraint. Called only once during
  // initial propagation of handler.
  virtual void initialize() { }

  // Filter a layer during topdown
  virtual void filter_layer_topdown(int layer) = 0;

  // Filter a layer during bottomup
  virtual void filter_layer_bottomup(int layer) = 0;

  // Called when the incoming arcs of a node are modified
  virtual void setup_node_incoming(int layer, int node, RefinementGroup ref_group, bool new_node) = 0;

  // Called when the outgoing arcs of a node are modified
  virtual void setup_node_outgoing(int layer, int node, RefinementGroup ref_group, bool new_node) = 0;

  // Called at all MDD processing rounds before topdown pass starts
  virtual void initialize_topdown() { }

  // Update all lagragian multipliers that the MDD is using - Margaritas code
  virtual void set_all_lagr_multipliers() { }

  // Called at all MDD processing rounds before bottomup pass starts
  virtual void initialize_bottomup() { }

  // Split nodes to compose refinement groups during topdown
  virtual void partition_refinement_group(int layer, RefinementGroup* ref_groups, int &n_groups, int n_rounds) { }

  // Called when all variables of the MDD handler are fixed */
  virtual void all_vars_fixed() { }

  // Called when MDD propagation round ends. It cannot add or remove arcs */
  virtual void end_round() { }

  // Returns if propagator is in layer
  virtual bool in_layer(int layer) = 0;

  /** Print node states */
  virtual void print_states() { };


  // Extract model variables to solver variables. Called
  // once during CP Optimizer model extraction. You should use
  // current IloCP passed as parameter instead of propagator
  // IloCP, since handler is not initialized yet
  virtual void extract(IloCP cp, const IloCPConstraintI* ct) { }


  /*
   * --------------------------------------------------
   * MDD modification control
   * --------------------------------------------------
   */

//  /** Returns if node needs to be update according to its incoming states */
//  bool needs_update_incoming(int layer, int node);
//
//  /** Indicate that incoming arcs of node were modified */
//  void modified_in_node(int layer, int node);
//
//  /** Indicate that node was updated due to changes in its incoming arcs */
//  void reset_in_node(int layer, int node);
//
//  /** Update reversible modification incoming state */
//  void update_modified_in(int layer, int node);
//
//  /** Returns if node needs to be update according to its outgoing states */
//  bool needs_update_outgoing(int layer, int node);
//
//  /** Indicate that outgoing arcs of node were modified */
//  void modified_out_node(int layer, int node);
//
//  /** Indicate that node was updated due to changes in its outgoing arcs */
//  void reset_out_node(int layer, int node);
//
//  /** Update reversible modification outgoing state */
//  void update_modified_out(int layer, int node);
//
//  /** Update temporary state */
//  void set_modified_in(int layer, int node);
//
//  /** Update temporary state */
//  void set_modified_out(int layer, int node);


  /*
   * ------------------------------------------------
   * Utilities for constraints
   * ------------------------------------------------
   */

  RevMDD*  get_mdd();        /**< Get MDD reference */
  IloCP    get_CP();         /**< Get CP reference */
  void     fail();           /**< fail search node */

  void	 set_id(int _id);  /**< set propagator id */


  /**< Initialize internal parameters */
  void initialize_parameters(IloCP _cp, RevMDD* _mdd);


  /*
   * --------------------------------------------------
   * MDD modification control
   * --------------------------------------------------
   */

//  IlcRevBool***			incoming_modified_rev;			/**< rev: if node incoming arcs were modified (per layer/node) */
//  IlcRevBool***			outgoing_modified_rev;			/**< rev: if node outgoing arcs were modified (per layer/node) */
//
//  bool**					incoming_modified;				/**< if node incoming arcs were modified (per layer/node) */
//  bool**					outgoing_modified;				/**< if node outgoing arcs were modified (per layer/node) */



//  /** Get current timestamp */
//  int  get_timestamp();
//
//  /** Set current timestamp */
//  void set_timestamp(int _timestamp);

  /*
   * -------------------------------------
   * Statistics
   * -------------------------------------
   */

  Stats			stats;									/**< statistics */
  long int		num_changes_states;						/**< state updates that were effective */
  long int		num_state_updates;						/**< total number of state updates */


  int				prop_id;		  /**< propagator id */

private:

  IloCP           cp;               /**< CP base class */
  RevMDD*         mdd;              /**< MDD base class */

  int			  timestamp;		  /**< filter/refinement time process */

};



/**
 * Inline implementations
 */

/**
 * Fail search node
 */
inline void MDDPropagator::fail() {
  cp.fail();
}

/**
 * Get methods
 */

inline RevMDD* MDDPropagator::get_mdd() {
  return mdd;
}


inline IloCP MDDPropagator::get_CP() {
  return cp;
}

/**
 * Initialize internal parameters
 */
inline void MDDPropagator::initialize_parameters(IloCP _cp, RevMDD* _mdd) {
  cp = _cp;
  mdd = _mdd;

//  // initialize MDD modification control
//
//  incoming_modified = new( cp.getHeap() ) bool*[mdd->get_num_layers()];
//  incoming_modified_rev = new( cp.getHeap() ) IlcRevBool**[mdd->get_num_layers()];
//
//  outgoing_modified = new( cp.getHeap() ) bool*[mdd->get_num_layers()];
//  outgoing_modified_rev = new( cp.getHeap() ) IlcRevBool**[mdd->get_num_layers()];
//
//
//  for( int l = 0; l < mdd->get_num_layers(); l++ ) {
//    incoming_modified_rev[l] = new( cp.getHeap() ) IlcRevBool*[mdd->get_max_width()];
//    outgoing_modified_rev[l] = new( cp.getHeap() ) IlcRevBool*[mdd->get_max_width()];
//
//    incoming_modified[l] = new( cp.getHeap() ) bool[mdd->get_max_width()];
//    outgoing_modified[l] = new( cp.getHeap() ) bool[mdd->get_max_width()];
//
//    // initial state: all nodes need update
//    for( int w = 0; w < mdd->get_max_width(); w++ ) {
//      incoming_modified_rev[l][w] = new( cp.getHeap() ) IlcRevBool(cp, true);
//      outgoing_modified_rev[l][w] = new( cp.getHeap() ) IlcRevBool(cp, true);
//    }
//  }
//
//  // root only requires outgoing update
//  incoming_modified_rev[0][0]->setValue(cp, false);
//
//  // terminal only requires incoming update
//  outgoing_modified_rev[mdd->get_num_layers()-1][0]->setValue(cp, false);

}

/**
 * Set propagator id
 */
inline void MDDPropagator::set_id(int _prop_id) {
  prop_id = _prop_id;
}

///**
// * Time stamp control
// */
//inline int MDDPropagator::get_timestamp() {
//  return timestamp;
//}
//
//inline void MDDPropagator::set_timestamp(int _timestamp) {
//  timestamp = _timestamp;
//}


/**
 * MDD modification control
 */

///**
// * Returns if node needs to be update according to its incoming states
// **/
//inline bool MDDPropagator::needs_update_incoming(int layer, int node) {
//  return incoming_modified[layer][node];
//}
//
///**
// * Returns if node needs to be update according to its outgoing states
// **/
//inline bool MDDPropagator::needs_update_outgoing(int layer, int node) {
//  return outgoing_modified[layer][node];
//}
//
//
///**
// * Indicate that incoming arcs of node were modified
// **/
//inline void MDDPropagator::modified_in_node(int layer, int node) {
//  //cout << "\t\tNode " << layer << "," << node << " scheduled for incoming update - propagator " << prop_id << endl;
//  incoming_modified[layer][node] = true;
//}
//
///**
// * Indicate that node was updated due to changes in its incoming arcs
// **/
//inline void MDDPropagator::reset_in_node(int layer, int node) {
//  incoming_modified[layer][node] = false;
//}
//
///**
// * Indicate that outgoing arcs of node were modified
// **/
//inline void MDDPropagator::modified_out_node(int layer, int node) {
//  //cout << "\t\tNode " << layer << "," << node << " scheduled for outgoing update - propagator " << prop_id << endl;
//  outgoing_modified[layer][node] = true;
//}
//
///**
// * Indicate that node was updated due to changes in its outgoing arcs
// **/
//inline void MDDPropagator::reset_out_node(int layer, int node) {
//  outgoing_modified[layer][node] = false;
//}
//
///**
// * Update incoming state: rev
// **/
//inline void MDDPropagator::update_modified_out(int layer, int node) {
//  if( outgoing_modified[layer][node] != outgoing_modified_rev[layer][node]->getValue() ) {
//    outgoing_modified_rev[layer][node]->setValue(cp, outgoing_modified[layer][node]);
//  }
//}
//
///**
// * Update outgoing state: rev
// **/
//inline void MDDPropagator::update_modified_in(int layer, int node) {
//  if( incoming_modified[layer][node] != incoming_modified_rev[layer][node]->getValue() ) {
//    incoming_modified_rev[layer][node]->setValue(cp, incoming_modified[layer][node]);
//  }
//}
//
///**
// * Update incoming state
// **/
//inline void MDDPropagator::set_modified_in(int layer, int node) {
//  incoming_modified[layer][node] = incoming_modified_rev[layer][node]->getValue();
//}
//
///**
// * Update outgoing state
// **/
//inline void MDDPropagator::set_modified_out(int layer, int node) {
//  outgoing_modified[layer][node] = outgoing_modified_rev[layer][node]->getValue();
//}






#endif /* MDD_CONSTRAINT_HPP_ */
