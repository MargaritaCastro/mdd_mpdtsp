//=======================================================
// MDD Handler - Header File
//=======================================================


#ifndef MDD_HANDLER_HPP_
#define MDD_HANDLER_HPP_

#define MAX_TIME 	1000000				// max process timestamp
#define MAX_ROUNDS 2					// max. number of rounds per search node

#include <cstdlib>
#include <ilcp/cpext.h>
#include "../mdd/mdd.hpp"
#include "../util/stats.hpp"
#include "../util/util.hpp"
#include "mdd_propagator.hpp"


/**
 * MDD Handler Constraint
 * Process MDD and propagators
 */
IloConstraint IloMDDHandler(

    IloEnv                      env,              /**< ILOG env */
    IloIntVarArray              vars,             /**< variable array (in MDD order) */

    IloIntervalVarArray			interval_vars,	  /**< interval variables: this is required only for CP+MDD technique.
         	 	 	 	 	 	 	 	 	 	 	 	   If not empty, then the domains of 'vars' relate to
        												   the indices of the interval variable; e.g., if vars represents a sequence of
        												   intervals, then vars[i]=j indicates that the i-th position is occupied by interval
                                                           'interval_vars[j]' */

    IloIntervalSequenceVar		seq,			  /**< sequence var; this is required only for CP+MDD technique.
        												   If you want pure MDD propagation, you can pass an IloIntervalSequence
        												   object with empty handle. For CP+MDD, you have to pass a sequence variable
        												   containing all intervals to be fixed and they types must be consistent
        												   with the activity indices in the permutation 'vars' and interval_vars */

    IloArray<MDDPropagator*>    cons,             /**< MDD propagators */
    int                         max_width,        /**< max MDD width */
    //IloIntervalVarArray			commodity,
    const char*                name);            /**< propagator name */




/**
 * Typedef for propagator array
 */
typedef vector<MDDPropagator*> PropagatorArray;

/**
 * MDD Handler Implementation Class
 */
class IlcMDDHandlerImp : public IlcConstraintI {

public:

    /** Constructors */
    IlcMDDHandlerImp(IloCP cp,
    		IlcIntVarArray _vars, IlcIntervalVarArray _interval_vars, IlcIntervalSequenceVar _seq_var,
    		IloArray<MDDPropagator*> _cons, IlcIntExp* _next, int _max_width); //, IlcIntervalVarArray _commodity);

    /** Empty destructor */
    ~IlcMDDHandlerImp(){ };

    /** ILOG search control */
    virtual void propagate();                     /**< Propagation method */
    virtual void post();                          /**< Daemon setting */
    virtual IlcBool isViolated() const;           /**< Violation checking */


private:

    /** Update MDD arcs according to changes in variable domains */
    void synchronize_mdd_domains();

    /** Update domains according to changes in MDD arcs */
    void synchronize_domains_mdd();

    /** Process MDD in a top down fashion */
    void process_topdown();

    /** Process layer in a top down fashion */
    void filter_layer_topdown(int layer);

    /** Process MDD in a bottom-up fashion */
    void process_bottomup();

    /** Process layer in a top down fashion */
    void filter_layer_bottomup(int layer);

    /** Refine layer - top down */
    void refine_layer_topdown(int previous_layer);

    /** Refine layer - bottom up */
    void refine_layer_bottomup(int next_layer);

    /** Reserve memory */
    void reserve_memory();

    /** Remove nodes without incoming arcs */
    void clean_mdd_incoming();

    /** Remove nodes without outgoing arcs */
    void clean_mdd_outgoing();

    /** Remove dangling nodes  */
    void clean_mdd();

    /** Attributes */

    IlcIntVarArray              	   vars;           			/**< MDD variables */
    IlcIntervalVarArray			   	   interval_vars;			/**< interval variables */
    IlcIntervalSequenceVar			   seq_var;					/**< sequence variable */
    IloArray<MDDPropagator*>    	   propagators;    			/**< MDD propagators */
    vector< vector<MDDPropagator*> >   prop_layer;     			/**< MDD propagators for a particular layer */
    const int                   	   max_width;      			/**< maximum MDD width */

    RevMDD*                     	   mdd;            		    /**< MDD */
    IlcRevInt                   	   choice_points;  		    /**< last number of choice points considered */
    int                   	   		   rounds_node;  		    /**< number of rounds in a search node */
    IlcRevBool						   is_handler_active;		/**< if handler is active */

    int								   n_rounds;	   		    /**< number of rounds the process is called */


    /** MDD Processing */

    SimpleList<int>**           	   node_map;               /**< map from old to newly created nodes */

    RevMDD::Arc**               	   new_arcs_stack;         /**< new arcs to add in arc list */
    int                         	   size_arc_stack;         /**< size of arc stack */

    SimpleList<RevMDD::Arc*>*		   arcs_layer;         	   /**< arcs in a layer */
    int								   size_arcs_layer;		   /**< size of arcs in layer list */

    RefinementGroup*				   ref_groups;			   /**< refinement groups */
    int								   n_ref_groups;		   /**< current number of ref. groups */

    bool*							   node_taken;			   /**< auxiliary to indicate that a node was taken */
    SimpleList<int>*				   available_nodes;		   /**< nodes that are available */
    int*							   node_to_group;		   /**< map from nodes to groups */

    int**							   out_arcs_node;		   /**< number of out arcs per layer/node */
    int**							   in_arcs_node;		   /**< number of in arcs per layer/node */

    bool*							   node_removed;		   /**< if MDD node has been removed */
    int*							   arc_group;			   /**< arc group */

    bool*							   is_new_group;		   /**< if group of arcs will yield a new node */

    /**
     * Modification control
     **/

    bool*							   modified_in;			  /**< if incoming arcs of node were modified */
    bool*							   modified_out;	      /**< if outgoing arcs of node were modified */

    int								   timestamp;			  /**< modification date */

    int*							   arc_counter_A; 	 	  /**< how many incoming/outgoing arcs of each node */
    int*							   arc_counter_B;   	  /**< how many incoming/outgoing arcs of each node */

    IlcRevInt						   last_head_position;	  /**< last stored position in head-tail list */
    IlcRevInt						   first_tail_position;	  /**< first tail position in head-tail list */

    IlcIntExp*						   next_acts;			  /**< expression indicating the next activity in the sequence */
   // IlcIntervalVarArray			   	   commodity;			
};


/**
 * ----------------------------------------------------------------------
 * Inline implementations for IlcMDDHandlerImp
 * ----------------------------------------------------------------------
 */


/**
 * --------------------------------------------
 * Comparators
 * --------------------------------------------
 */
struct ArcGroupComparator {
	int* arc_group;
	ArcGroupComparator(int* _arc_group) : arc_group(_arc_group) { }
    inline bool operator()(RevMDD::Arc* arcA, RevMDD::Arc* arcB) {
        return( arc_group[arcA->get_layer_id()] < arc_group[arcB->get_layer_id()] );
    }
};


/**
 * ILC Constraint wrappers
 */

inline IlcConstraint IlcMDDHandler(
        IloCP cp,
        IlcIntVarArray _vars,
        IlcIntervalVarArray _interval_vars,
        IlcIntervalSequenceVar _seq_var,
        IloArray<MDDPropagator*> _cons,
        IlcIntExp* _next,
        int _max_width) //,IlcIntervalVarArray	_commodity)
{
    return new(cp.getHeap()) IlcMDDHandlerImp(cp, _vars, _interval_vars, _seq_var, _cons, _next, _max_width);//, _commodity);
}


/**
 * ILOG call for violation checking is disabled, since this will be done
 * during propagation
 */
inline IlcBool IlcMDDHandlerImp::isViolated() const {
    return IlcFalse;
}


/**
 * Process MDD in a top down fashion
 */
inline void IlcMDDHandlerImp::process_topdown() {

    // reset node mapping for first layer
    node_map[0]->clear();
    node_map[0]->add(0);

    // no nodes were removed at first
	memset( node_removed, false, sizeof(bool)*mdd->get_max_width() );

	// initialize top-down search for all propagators
	for( int p = 0; p < propagators.getSize(); p++ ) {
		propagators[p]->initialize_topdown();
    if( n_rounds <= 5) {
      //cout << "Updating lagrangian multipliers in round " <<  n_rounds<< endl;
      propagators[p]->set_all_lagr_multipliers();
    }
	}

	// compile each layer separately
	for( int l = 0; l < mdd->get_num_layers()-1; l++ ) {

		//cout << "\n****** [TopDown] Layer " << l << " ******" << endl;

//		// initialize process timestamp
//		timestamp = (timestamp+1) % MAX_TIME;
//		for( vector<MDDPropagator*>::iterator p = prop_layer[l].begin();
//				p != prop_layer[l].end();
//				p++ )
//		{
//			(*p)->set_timestamp(timestamp);
//		}

		// filter arcs top-down
		filter_layer_topdown(l);

		// refine next layer
		refine_layer_topdown(l);

		//cout << "\n****** [TopDown] End layer " << l << " ******" << endl;
	}
}


/**
 * Filter layer in a top down fashion
 */
inline void IlcMDDHandlerImp::filter_layer_topdown(int layer) {

	IloCP cp = getCP();

	// =========================================================================
	// The previous actions that affect this layer are only performed now.
	// Namely, we iterate through the MDD arcs and delete the arcs whose
	// source nodes were removed in previous filtering iterations.
	// In addition, we add the outgoing arcs of the nodes that were created
	// due to refinement in the previous layer
	// =========================================================================

	RevMDD::Arc* previous = mdd->arcs_layer[layer].get_header();

	RevMDD::Arc* new_arc = NULL;   // new arcs created due to refinement of previous layer
	size_arc_stack = 0;
	int target_node;

	for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
			arc != NULL;
			arc = mdd->arcs_layer[layer].get_next(arc) )
	{
		if( node_map[arc->source_node]->is_empty() ) {

			//cout << "Removing arc (" << layer << "," << arc->get_source() << "," << arc->get_val() << ")" << endl;

			// remove arc if its source was removed
			previous->next.setValue(cp, arc->next.getValue());

		} else {

			// create outgoing arcs of refined nodes
			for( int i = 0; i < node_map[arc->get_source()]->size; i++ ) {
				if( node_map[arc->get_source()]->elements[i] != arc->get_source() ) {

					// obtain new arc and set its target
					new_arc = mdd->get_arc(layer, node_map[arc->get_source()]->elements[i], arc->val);
					target_node = arc->get_target();
					if( new_arc->get_target_rev() != target_node ) {
						new_arc->set_target_rev(cp, target_node);
					}

					// update new arc target node
					new_arc->target_node = target_node;

					// add to new arc stack
					new_arcs_stack[size_arc_stack++] = new_arc;
				}
			}
			previous = arc;
		}
	}

	// add arcs of the stack to the overall list of arcs
	if( size_arc_stack > 0 ) {
		for( int i = 0; i < size_arc_stack; i++ ) {
			assert( new_arcs_stack[i] != NULL );
			previous = mdd->arcs_layer[layer].extend(previous, new_arcs_stack[i]);
		}
		mdd->arcs_layer[layer].close(previous);
	}

	// if layer is empty, problem is infeasible
	if( mdd->arcs_layer[layer].is_empty() ) {
		cp.fail();
	}

	// =============================================================================
	// Filtering is applied to ensure that the layer arcs are consistent w.r.t. to
	// all constraints before applying refinement
	// =============================================================================

	for( unsigned int p = 0; p < prop_layer[layer].size(); p++ ) {

		prop_layer[layer][p]->filter_layer_topdown(layer);

		// if layer is empty, problem is infeasible
		if( mdd->arcs_layer[layer].is_empty() ) {
			//cout << "layer is empty" << endl;
			cp.fail();
		}
	}
}


/**
 * Process MDD in a bottom-up fashion
 */
inline void IlcMDDHandlerImp::process_bottomup() {

    // reset out nodes for terminal node
	out_arcs_node[mdd->get_num_layers()-1][0] = 1;

	// initialize top-down search for all propagators
	for( int p = 0; p < propagators.getSize(); p++ ) {
		propagators[p]->initialize_bottomup();
	}

	// compile each layer separately
	for( int l = mdd->get_num_layers()-2; l >= 0; l-- ) {

		//cout << "\n****** [BottomUp] Layer " << l << " ******" << endl;

//		// initialize process timestamp
//		timestamp = (timestamp+1) % MAX_TIME;
//		for( vector<MDDPropagator*>::iterator p = prop_layer[l].begin();
//				p != prop_layer[l].end();
//				p++ )
//		{
//			(*p)->set_timestamp(timestamp);
//		}

		// filter arcs top-down
    //cout << "Filtering layer..." << endl;
		filter_layer_bottomup(l);

		// refine next layer
    //cout << "Saving information..." << endl;
		refine_layer_bottomup(l);

		//cout << "\n****** [BottomUp] End layer " << l << " ******" << endl;
	}

	// remove dangling nodes
	clean_mdd_incoming();
}

/**
 * Filter layer in a bottom up fashion.
 */
inline void IlcMDDHandlerImp::filter_layer_bottomup(int layer) {

	IloCP cp = getCP();

	// =========================================================================
	// The previous actions that affect this layer are only performed now.
	// Namely, we iterate through the MDD arcs and delete the arcs whose
	// source nodes were removed in previous filtering iterations.
	// In addition, we add the outgoing arcs of the nodes that were created
	// due to refinement in the previous layer
	// =========================================================================

	/**
	 * 1. Obtain list of arcs from MDD object. Already removes the values
	 *  that are inconsistent with variable domains and nodes in the layer
	 */
	RevMDD::Arc* previous = mdd->arcs_layer[layer].get_header();

	for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
			arc != NULL;
			arc = mdd->arcs_layer[layer].get_next(arc) )
	{
		if( out_arcs_node[layer+1][arc->get_target()] == 0 ) {

			// remove arc if its target node has no outgoing arcs
			previous->next.setValue(cp, arc->next.getValue());

			// mark its source node as modified
			mdd->modified_out_node(layer, arc->get_source());

		} else {
			previous = arc;
		}
	}

	// if layer is empty, problem is infeasible
	if( mdd->arcs_layer[layer].is_empty() ) {
		cp.fail();
	}

	/**
	 * 2. Run each propagator to eliminate infeasible arcs from layer
	 */

	for( unsigned int p = 0; p < prop_layer[layer].size(); p++ ) {

		prop_layer[layer][p]->filter_layer_bottomup(layer);

		// if layer is empty, problem is infeasible
		if( mdd->arcs_layer[layer].is_empty() ) {
			cp.fail();
		}
	}
}


/**
 * Remove nodes without incoming arcs.
 */
inline void IlcMDDHandlerImp::clean_mdd_incoming() {
  IloCP cp = getCP();

  // reset in_arcs of root node
  in_arcs_node[0][0] = 1;

  // we keep a pointer to the previous arc considered for removal purposes
  RevMDD::Arc* previous;

  for( int l = 0; l < mdd->get_num_layers()-1; l++ ) {
    memset(in_arcs_node[l+1], 0, sizeof(int)*mdd->get_max_width());

    previous = mdd->arcs_layer[l].get_header();
    for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
        arc != NULL;
        arc = mdd->arcs_layer[l].get_next(arc) )
    {
      if( in_arcs_node[l][arc->get_source()] == 0 ) {
        // if source does not have incoming nodes, we remove the arc from the MDD
        previous->next.setValue(cp, arc->next.getValue());
      } else {
        // count number of incoming arcs of the next layer
        in_arcs_node[l+1][arc->get_target()]++;
        previous = arc;
      }
      // fail if resulting layer is empty
      if( mdd->arcs_layer[l].is_empty() ) {
        cp.fail();
      }
    }
  }
  assert( in_arcs_node[mdd->get_num_layers()-1][0] > 0 );
}

/**
 * Remove nodes without outgoing arcs.
 */
inline void IlcMDDHandlerImp::clean_mdd_outgoing() {

  IloCP cp = getCP();

  // reset out_arcs of terminal node
  out_arcs_node[mdd->get_num_layers()-1][0] = 1;

  RevMDD::Arc* previous;
  for (int l = mdd->get_num_layers()-2; l >= 0; --l) {
    memset(out_arcs_node[l], 0, sizeof(int)*mdd->get_max_width());

    previous = mdd->arcs_layer[l].get_header();
    for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
        arc != NULL;
        arc = mdd->arcs_layer[l].get_next(arc) )
    {
      if( out_arcs_node[l+1][arc->get_target()] == 0 ) {
        // if source does not have incoming nodes, we remove the arc from the MDD
        previous->next.setValue(cp, arc->next.getValue());
        //cout << "effective!" << endl;
      } else {
        // count number of incoming arcs of the next layer
        ++out_arcs_node[l][arc->get_source()];
        previous = arc;
      }
    }
    // fail if resulting layer is empty
    if( mdd->arcs_layer[l].is_empty() ) {
      cp.fail();
    }
  }
  assert( out_arcs_node[0][0] > 0 );
}






#endif /* MDD_HANDLER_HPP_ */
