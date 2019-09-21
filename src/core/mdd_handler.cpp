//=======================================================
// MDD handler - Implementation File
//=======================================================


#include <cstdlib>
#include "mdd_handler.hpp"
#include "../core/intset.hpp"
#include <bitset>

/*
 * -------------------------------------------------------
 * ILOG demons and wrappers
 * -------------------------------------------------------
 */

/**
 * Wrapper for ILO object
 */
ILOCPCONSTRAINTWRAPPER5(IloMDDHandler,
		cp,
		IloIntVarArray, _vars,
		IloIntervalVarArray, _interval_vars,
		IloIntervalSequenceVar, _seq_var,
		IloArray<MDDPropagator*>, _cons,
		int, _max_width) //,IloIntervalVarArray,_commodity)
		{
		  use(cp, _vars);

	// 

	// compose sequence var
	IlcIntervalSequenceVar seq_var;
	if( _seq_var.getImpl() != NULL ) {
		use(cp, _seq_var);
		seq_var = cp.getIntervalSequence(_seq_var);
	}

	// compose interval var
	IlcIntervalVarArray interval_vars(cp, _interval_vars.getSize());
	for( int i = 0; i < interval_vars.getSize(); i++ ) {
		use(cp, _interval_vars[i]);
		interval_vars[i] = cp.getInterval(_interval_vars[i]);
	}
	
	
	// compose commodity
	//IlcIntervalVarArray commodity(cp, _commodity.getSize());
	//for( int i = 0; i < commodity.getSize(); i++ ) {
	//	use(cp, _commodity[i]);
	//	commodity[i] = cp.getInterval(_commodity[i]);
	//}

	// compose 'next' variables in a sequence
	IlcIntExp* next_acts = NULL;
	if( _seq_var.getImpl() != NULL ) {
	  // !!!
		// next_acts = new( cp.getHeap() ) IlcIntExp[interval_vars.getSize()];
		// IloIntExpr next_act;
		// for( int i = 0; i < interval_vars.getSize(); ++i ) {
		// 	next_act = IloTypeOfNext(_seq_var, _interval_vars[i], -1);
		// 	//use(cp, next_act);
		// 	cout << "nha!" << endl;
		// 	next_acts[i] = cp.getIntExp( next_act );
		// }
	}

	// extract specific variables from constraints
	for (int i = 0; i < _cons.getSize(); ++i) {
	  _cons[i]->extract(cp, this);
	}

	return IlcMDDHandler(cp, cp.getIntVarArray(_vars), interval_vars, seq_var, _cons, next_acts, _max_width); //, commodity);
}


/**
 * Monitor changes of each variable
 */
void IlcMDDHandlerImp::post() {
	// sequence variable changes
	if( seq_var.getImpl() != NULL ) {
		seq_var.whenNotSequenced(this);
	}
	// permutation variables changes
	for( IlcInt i = 0; i < vars.getSize(); i++ ) {
		vars[i].whenDomain(this);
	}
	//  for( int i = 0; i < interval_vars.getSize(); ++i ) {
	//    interval_vars[i].whenIntervalDomain(this);
	//  }
}


/**
 * MDD Handler Constructor.
 * Initialize MDD and reserve node/arc memory.
 */
IlcMDDHandlerImp::IlcMDDHandlerImp(IloCP cp,
		IlcIntVarArray _vars,
		IlcIntervalVarArray _interval_vars,
		IlcIntervalSequenceVar _seq_var,
		IloArray<MDDPropagator*> _cons,
		IlcIntExp* next,
		int _max_width) //,IlcIntervalVarArray _commodity)
: IlcConstraintI(cp), vars(_vars),
  interval_vars(_interval_vars), seq_var(_seq_var),
  propagators(_cons), max_width(_max_width),
  next_acts(next) //, commodity(_commodity)
{
	//cout << "\n[MDD-Handler] Initializing...\n";

	n_rounds = 0;
	timestamp = 0;
	is_handler_active.setValue(cp, true);

	// create reversible MDD
	mdd = new RevMDD(cp, vars, max_width);

	// initialize filter internal parameters
	prop_layer.resize(mdd->get_num_layers());
	for( int c = 0; c < propagators.getSize(); c++ ) {
		propagators[c]->initialize_parameters(cp, mdd);
		// register propagators in their respective layers
		for( int l = 0; l < mdd->get_num_layers(); l++ ) {
			if( propagators[c]->in_layer(l) ) {
				prop_layer[l].push_back(propagators[c]);
			}
		}
	}

	// print maximum number of constraints in a layer
	int max_num_ct = 0;
	for( int l = 0; l < mdd->get_num_layers(); l++ ) {
		max_num_ct = MAX(max_num_ct, (int)prop_layer[l].size());
	}
	cout << "\n[MDD-Handler] Maximum number of constraints in a layer: " << max_num_ct << endl << endl;

	// initialize filter propagation
	for( int c = 0; c < propagators.getSize(); c++ ) {
		propagators[c]->set_id(c);
		propagators[c]->initialize();
	}

	// reserve memory
	reserve_memory();

	//cout << "[MDD-Handler] done.\n" << endl;
	choice_points.setValue(cp, -1);
}

//
//// Goal testing
//ILCGOAL1(PrintX, IlcInt, x){
//   IloCP cp = getCP();
//   cout << "PrintX: a goal with one data member" << std::endl;
//   cout << x << std::endl;
//   return 0;
// }

/**
 * Propagate permutation
 */
void IlcMDDHandlerImp::propagate() {


	// -------------------------------------------------------------
	// Check if handler is active.
	// -------------------------------------------------------------

	if( !is_handler_active.getValue() ) {
		return;
	}

	// -------------------------------------------------------------
	// Initialize propagation round by synchronizing MDD arcs with
	// variable domains
	// -------------------------------------------------------------

	n_rounds++;
	//getCP().addReversibleAction( PrintX(getCP(), n_rounds) );

	//	  cout << "\n\n[MDD-Handler] Propagation round " << n_rounds << endl;
	 // cout << "Variables: " << endl;
	  //for( int v = 0; v < vars.getSize(); v++ ) {
	 //   cout << vars[v] << " ";
	//  }
	 // cout << endl;
	  //cout << "Activities: " << endl;
	 // for( int v = 0; v < vars.getSize(); v++ ) {
	 //   cout << "\t" << interval_vars[v] << endl;
	 // }
	  //cout << endl;
	 //cout << "num_fails: " << getCP().getInfo(IloCP::NumberOfFails) << endl;


	//   Uncomment this line to export MDD in GML format. You can visualize the
	//   graph using a software called yEd, which automatically creates a layer-by-layer format
	//   for you
	//char filename[256];
	//sprintf(filename, "graphs/before-round-%d.gml", n_rounds);
	//mdd->export_to_gml(filename);


	//    char filename[256];
	//    if( n_rounds == 125 ) {
	//      sprintf(filename, "graphs/before-round-%d.gml", n_rounds);
	//      mdd->export_to_gml(filename);
	//    }


	// -------------------------------------------------------------
	// Synchronize MDD layers with variable domains
	// -------------------------------------------------------------
	synchronize_mdd_domains();


	// -------------------------------------------------------------
	// Check number of rounds in current search node
	// -------------------------------------------------------------
	IloCP cp = getCP();
	if( cp.getInfo(IloCP::NumberOfChoicePoints) != choice_points.getValue() ) {
		choice_points.setValue(cp, cp.getInfo(IloCP::NumberOfChoicePoints));
		rounds_node = 0;
	} else {
		if( rounds_node > MAX_ROUNDS ) {

			// if all variables are fixed, setup sum of distances
			bool is_fixed = true;
			for( int i = 0; i < vars.getSize() && is_fixed; i++ ) {
				is_fixed = vars[i].isFixed();
			}

			if( is_fixed ) {
				// advise constraints that all variables are fixed
				for( int p = 0; p < propagators.getSize(); p++ ) {
					propagators[p]->all_vars_fixed();
				}
				// disable handler
				is_handler_active.setValue(getCP(), false);
				//cout << "constraint is not active" << endl;
			}
			return;
		}
		rounds_node++;
	}


	// -------------------------------------------------------------
	// Top-down filtering and refinement
	// -------------------------------------------------------------
	//char filename[256];
	//cout << endl << "top down..." << endl;
	process_topdown();
	//cout << "\tdone top down!" << endl;
	//  if( n_rounds == 125 ) {
	    //sprintf(filename, "graphs/after-topdown-%d.gml", n_rounds);
	    //mdd->export_to_gml(filename);
	//  }


	// -------------------------------------------------------------
	// Bottom-up filtering and refinement
	// -------------------------------------------------------------

	//cout << endl << "bottom-up..." << endl;
	process_bottomup();
	//cout << endl << "\tend bottom-up!" << endl;

	//  if( n_rounds == 99 ) {
	  // sprintf(filename, "graphs/after-bottomup-%d.gml", n_rounds);
	  // mdd->export_to_gml(filename);
	//  }

	// -------------------------------------------------------------
	// Synchronize variable domains with MDD arcs
	// -------------------------------------------------------------

	synchronize_domains_mdd();


	//char filename[256];
	//sprintf(filename, "graphs/after-round-%d.gml", n_rounds);
	//mdd->export_to_gml(filename);
	
	
	// finish round
	for( int p = 0; p < propagators.getSize(); p++ ) {
		propagators[p]->end_round();
	}

	//exit(1); //to stop the mdd at the root node

	/*
	cout << "\n\n[Disjunctive] Round " << n_rounds << " ended" << endl;
	
  	cout << "Activities: " << endl;
	for( int v = 0; v < vars.getSize(); v++ ) {
		cout << "\t" << interval_vars[v] << endl;
	}
	
	cout << "After variables: " << endl;
	for( int v = 0; v < vars.getSize(); v++ ) {
	 	cout << vars[v] << " ";
	}
	cout << endl;

	cout << "num_fails: " << getCP().getInfo(IloCP::NumberOfFails) << endl;
	*/
	

	// -------------------------------------------------------------
	// Disable handler if all MDD variables are fixed
	// -------------------------------------------------------------

	// if all variables are fixed, setup sum of distances
	bool is_fixed = true;
	for( int i = 0; i < vars.getSize() && is_fixed; i++ ) {
		is_fixed = vars[i].isFixed();
	}

	if( is_fixed ) {

	//	cout << "all vars fixed" << endl;

		// advise constraints that all variables are fixed
		for( int p = 0; p < propagators.getSize(); p++ ) {
			propagators[p]->all_vars_fixed();
		}

		// disable handler
		is_handler_active.setValue(getCP(), false);
		//cout << "constraint is not active" << endl;
	}

	//cout << "\n\n[Disjunctive] Round " << n_rounds << " ended" << endl;
	// exit(1);

	if( n_rounds == 30) {
		//exit(1);
	}
}



/**
 * Reserve memory
 */
void IlcMDDHandlerImp::reserve_memory() {

	IloCP cp = getCP();

	int max_var_domain = vars.getMaxMax();

	// initialize node map
	node_map = new( cp.getHeap() ) SimpleList<int>*[mdd->get_max_width()];
	for( int w = 0; w < mdd->get_max_width(); w++ ) {
		node_map[w] = new( cp.getHeap() ) SimpleList<int>(cp, mdd->get_max_width());
	}

	// list to stack newly created arcs
	new_arcs_stack = new( cp.getHeap() ) RevMDD::Arc*[(max_var_domain+1) * mdd->get_max_width()];

	// arcs in a layer
	arcs_layer = new( cp.getHeap() ) SimpleList<RevMDD::Arc*>(cp, mdd->get_max_arc_id()+1);

	// allocate refinement groups
	ref_groups = new( cp.getHeap() ) RefinementGroup[mdd->get_max_width()];
	for( int w = 0; w < mdd->get_max_width(); w++ ) {
		ref_groups[w] = new( cp.getHeap() ) SimpleList<RevMDD::Arc*>(cp, mdd->get_max_arc_id()+1);
	}

	// pre-allocate node availability/map
	node_taken = new( cp.getHeap() ) bool[mdd->get_max_width()];
	available_nodes = new( cp.getHeap() ) SimpleList<int>(cp, mdd->get_max_width());
	node_to_group = new( cp.getHeap() ) int[mdd->get_max_width()];

	// in/out arcs per node and layer
	out_arcs_node = new( cp.getHeap() ) int*[mdd->get_num_layers()];
	in_arcs_node = new( cp.getHeap() ) int*[mdd->get_num_layers()];
	for( int i = 0; i < mdd->get_num_layers(); i++ ) {
		out_arcs_node[i] = new( cp.getHeap() ) int[mdd->get_max_width()];
		in_arcs_node[i] = new( cp.getHeap() ) int[mdd->get_max_width()];
	}


	// node removal
	node_removed = new( cp.getHeap() ) bool[mdd->get_max_width()];

	// arc groups group
	arc_group = new( cp.getHeap() ) int[mdd->get_max_arc_id()+1];

	// arc counters
	arc_counter_A = new( cp.getHeap() ) int[mdd->get_max_width()];
	arc_counter_B = new( cp.getHeap() ) int[mdd->get_max_width()];

	is_new_group = new( cp.getHeap() ) bool[mdd->get_max_width()];

	// other modifiers
	last_head_position.setValue(cp, 0);
	first_tail_position.setValue(cp, vars.getSize()-1);
}





/**
 * Refine layer within a top-down perspective
 */
void IlcMDDHandlerImp::refine_layer_topdown(int previous_layer) {

	//cout << "\n[MDD-Handler] Refinement - previous layer: " << previous_layer << endl << endl;

	/*
	 * -----------------------------------------------------------------------------------
	 * 1. Organize refinement groups. Each group is initially composed
	 * by arcs that have the same target node
	 * -----------------------------------------------------------------------------------
	 */

	// update nodes that belong to the layer
	memset( node_removed, true, sizeof(bool)*mdd->get_max_width() );

	// reset node map and refinement groups
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		node_map[g]->clear();
		ref_groups[g]->clear();
	}

	int target_node;
	int n_groups = 0;

	for( RevMDD::Arc* arc = mdd->arcs_layer[previous_layer].get_first();
			arc != NULL;
			arc = mdd->arcs_layer[previous_layer].get_next(arc) )
	{
		target_node = arc->get_target();

		// update number of groups and node map
		if( node_removed[target_node] ) {
			node_removed[target_node] = false;
			node_map[target_node]->add(target_node);
			n_groups++;
		}

		// add arc to its corresponding group
		ref_groups[target_node]->add(arc);
	}

	/**
	 * -----------------------------------------------------------------------------------
	 * 2. Perform refinement if maximum width is not met
	 * -----------------------------------------------------------------------------------
	 */
	if( previous_layer+1 < mdd->get_num_layers()-1 ) {
		//cout << "\nRefining layer" << previous_layer+1 << endl;
		for( unsigned int c = 0; n_groups < mdd->get_max_width() && c < prop_layer[previous_layer].size(); c++ ) {
			prop_layer[previous_layer][c]->partition_refinement_group(previous_layer, ref_groups, n_groups, n_rounds);
			assert( n_groups <= mdd->get_max_width() );
		}
	}


	/**
	 * -----------------------------------------------------------------------------------
	 * 3. Create nodes from groups
	 * -----------------------------------------------------------------------------------
	 */

	// first define groups that will correspond to new nodes
	memset(is_new_group, false, sizeof(bool)*mdd->get_max_width());

	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		is_new_group[g] = (!ref_groups[g]->is_empty()) && (g != ref_groups[g]->get_first()->get_target());
	}

	// define states for new groups first, since they can use the state information from
	// the original nodes before refinement
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		if( is_new_group[g] ) {
			//cout << "New group !! ==>  " << g << endl;
			for( PropagatorArray::iterator prop = prop_layer[previous_layer].begin();
					prop != prop_layer[previous_layer].end();
					prop++ )
			{
				(*prop)->setup_node_incoming(previous_layer+1, g, ref_groups[g], true);
			}

			// add group to node map
			node_map[ref_groups[g]->get_first()->get_target()]->add(g);

			// reset arc targets
			for( int i = 0; i < ref_groups[g]->size; i++ ) {
				ref_groups[g]->elements[i]->set_target_rev(getCP(), g);
				ref_groups[g]->elements[i]->target_node = g;
			}
		}
	}

	// define states for existing nodes
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		// only consider non-empty groups
		if( !ref_groups[g]->is_empty() && !is_new_group[g] ) {
			for( PropagatorArray::iterator prop = prop_layer[previous_layer].begin();
					prop != prop_layer[previous_layer].end();
					prop++ )
			{
				(*prop)->setup_node_incoming(previous_layer+1, g, ref_groups[g], false);
			}
		}
	}
}



/**
 * Refine layer within a bottom-up perspective
 */
void IlcMDDHandlerImp::refine_layer_bottomup(int layer) {

	//cout << "\n[MDD-Handler] Refinement - previous layer: " << previous_layer << endl << endl;

	/*
	 * -----------------------------------------------------------------------------------
	 * 1. Organize refinement groups. Each group is initially composed
	 * by arcs that have the same target node
	 * -----------------------------------------------------------------------------------
	 */

	// reset mdd counters
	memset( out_arcs_node[layer], 0, sizeof(int)*mdd->get_max_width() );

	// clean current group list
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		ref_groups[g]->clear();
	}

	// organize arcs according to groups
	for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
			arc != NULL;
			arc = mdd->arcs_layer[layer].get_next(arc) )
	{
		// increment arc counters
		out_arcs_node[layer][arc->get_source()]++;

		// add to respective group
		ref_groups[arc->get_source()]->add(arc);
	}

	/**
	 * -----------------------------------------------------------------------------------
	 * 2. Setup nodes according to outgoing arcs
	 * -----------------------------------------------------------------------------------
	 */
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		if( !ref_groups[g]->is_empty() ) {

			// we check each propagator separately
			// todo: this can be done in an incremental way
			for( PropagatorArray::iterator prop = prop_layer[layer].begin(); prop != prop_layer[layer].end(); ++prop )	{
				(*prop)->setup_node_outgoing(layer, g, ref_groups[g], false);
			}

			// outgoing update was completed
			mdd->reset_out_node(layer, g);
		}
	}
}

/**
 * Synchronize MDD according to changes in variable domains.
 * It also prepares MDD for top-down and bottom-up refinement.
 *
 * TODO: changes in the MDD can be saved so that filtering/refinement
 * is only applied to nodes that really require it
 */
void IlcMDDHandlerImp::synchronize_mdd_domains() {
	IloCP cp = getCP();

	//cout << "\n\t[MDD-Handler] Synchronizing domains..." << endl;

	// If sequence variable is set, assume integer variables correspond to the permutation
	// of the sequence and propagate candidate/tail list accordingly.
	if( seq_var.getImpl() != NULL ) {

		// ----------------------------------------------------------
		// Propagate Head and Tail of the sequence
		// ----------------------------------------------------------

		// propagate head of the sequence
		int head = 0;
		for( IlcIntervalSequenceVar::Iterator task(seq_var, IlcIntervalSequenceVar::Head); task.ok(); ++task ) {
			if( !vars[head].isFixed() ) {
				vars[head].setValue(seq_var.getType(*task));
			}
			head++;
		}

		// setup pointer
		last_head_position.setValue(cp, head);

		// propagate tail of the sequence
		int tail = vars.getSize()-1;
		for( IlcIntervalSequenceVar::Iterator task(seq_var, IlcIntervalSequenceVar::Tail); task.ok(); ++task ) {
			if( !vars[tail].isFixed() ) {
				vars[tail].setValue(seq_var.getType(*task));
			}
			tail--;
		}

		// setup pointer
		first_tail_position.setValue(cp, tail);

		// propagate candidate lists
		// note: doing this seems to be stronger than iterating using IlcIntervalSequenceVar::CandidateHead
		if( head < tail ) {

			// propagate head candidate list
			for( IlcIntExpIterator next_task(vars[head]); next_task.ok(); ++next_task) {
				if( !seq_var.isCandidateHead(interval_vars[*next_task]) ) {
					vars[head].removeValue(*next_task);
				}
			}

			// propagate tail candidate list
			for( IlcIntExpIterator previous_task(vars[tail]); previous_task.ok(); ++previous_task) {
				if( !seq_var.isCandidateTail(interval_vars[*previous_task]) ) {
					vars[tail].removeValue(*previous_task);
				}
			}
		}
	}

	// initialize root arc counter
	in_arcs_node[0][0] = 1;

	// Check if any arc has been removed from the MDD
	for( int l = 0; l < vars.getSize(); ++l ) {

		memset(in_arcs_node[l+1], 0, sizeof(int)*mdd->get_max_width());

		// iterate through the arcs of this layer
		RevMDD::Arc* previous = mdd->arcs_layer[l].get_header();
		for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
				arc != NULL;
				arc = mdd->arcs_layer[l].get_next(arc) )
		{
			// update arc target
			arc->target_node = arc->get_target_rev();

			if( !vars[l].isInDomain(arc->val) || in_arcs_node[l][arc->get_source()] == 0 ) {
				// remove value
				previous->next.setValue(cp, arc->next.getValue());
			} else {
				previous = arc;
				++in_arcs_node[l+1][arc->get_target()];
			}
		}
		// fail if resulting layer is empty
		if( mdd->arcs_layer[l].is_empty() ) {
			cp.fail();
		}
	}

	// At the end of this procedure, you might still have some nodes
	// connected to R that do not have any outgoing arcs. We remove
	// these nodes to prevent refinement issues
	clean_mdd_outgoing();
}



/**
 * Synchronize variables domains according to changes in the MDD
 */
void IlcMDDHandlerImp::synchronize_domains_mdd() {

	//cout << "\n\t[MDD-Handler] Synchronizing domains..." << endl;

	// update variable domains
	mdd->update_variable_domains();

	/*
	 * If sequence variable is set, assume integer variables correspond to the permutation
	 * of the sequence and propagate candidate/tail list accordingly.
	 */

	// TODO: it is possible to implement this more efficiently by storing the last variable
	// that was modified from the last node, and doing all the process from there

	if( seq_var.getImpl() != NULL) {


		// ---------------------------------------------------
		// Head update
		// ---------------------------------------------------

		/*
		 * Set sequence according to permutation variables,
		 */

		// get first unfixed variable that was observed last time
		int head = last_head_position.getValue();

		// first activity in the sequence
		if( head == 0 ) {
			if( vars[head].isFixed() ) {
				if( seq_var.getEarliestInHead().getImpl() == NULL ) {
					seq_var.extendHead(interval_vars[vars[head].getValue()]);
				}
				head++;
			}
		}

		// activities 2..n
		for( ; head < vars.getSize(); ++head ) {
			if( vars[head].isFixed() ) {
				seq_var.setPrevious(interval_vars[vars[head-1].getValue()], interval_vars[vars[head].getValue()]);
			} else {
				break;
			}
		}

		// if all variables are fixed, nothing else need to be done: return!
		if( head == vars.getSize() ) {
			return;
		}

		// update pointer
		last_head_position.setValue(getCP(), head);


		/*
		 * The variable 'head' now contains the position of the last unfixed variable.
		 * We can thus update the head candidate list
		 */
		// H!!
		//if( head > 0 )
		//getCP().add( vars[head] == next_acts[vars[head-1].getValue()] );



		// ---------------------------------------------------
		// Tail update
		// ---------------------------------------------------

		/*
		 * Set sequence according to permutation variables,
		 */

		// get first unfixed variable that was observed last time
		int tail = first_tail_position.getValue();

		// last activity in the sequence
		if( tail == vars.getSize()-1 ) {
			if( vars[tail].isFixed() ) {
				if( seq_var.getLatestInTail().getImpl() == NULL ) {
					seq_var.extendTail(interval_vars[vars[tail].getValue()]);
				}
				tail--;
			}
		}

		// activities n-1..1
		for( ; tail >= 0; --tail ) {
			if( vars[tail].isFixed() ) {
				seq_var.setPrevious(interval_vars[vars[tail].getValue()], interval_vars[vars[tail+1].getValue()]);
			} else {
				break;
			}
		}

		// if all variables are fixed, nothing else need to be done: return!
		if( tail < 0 ) {
			return;
		}

		// update pointer
		first_tail_position.setValue(getCP(), tail);

		/*
		 * TODO: The variable 'tail' now contains the position of the last unfixed variable.
		 * We can thus update the head candidate list
		 */
		//getCP().add( vars[head] == next_acts[vars[head-1].getValue()] );

	}
}


