/*
 * ===============================================
 * Disjunctive Propagator Implementation -- last version
 * ===============================================
 */

#include <string>
#include <cmath>
#include "disjunctive.hpp"
#include "../core/held_karp_relax.hpp"
//#include "../codeFiles.h"

#define EPS_INF 1000
#define CAP_FILT 1 	// 1 means that I use the capacity filtering rule. (1)YES, (0)NO
#define PRIORITY 0 	// 1 means that I use my refinament rule.  (1)YES, (0)NO


/**
 * Reserve node state/pre-processed info memory
 **/
void DisjunctivePropagator::reserve_memory() {

	cout << "\n[Disjunctive] Reserving memory..." << endl;

	IloCP cp = get_CP();
	RevMDD* mdd = get_mdd();

	// initialize MDD node states
	// TODO: it is possible to compute more accurate bounds on the start / end time for each layer
	states = new( cp.getHeap() ) DisjunctiveState**[mdd->get_num_layers()];
	for( int l = 0; l < mdd->get_num_layers(); l++ ) {
		states[l] = new( cp.getHeap() ) DisjunctiveState*[mdd->get_max_width()];
		for( int w = 0; w < mdd->get_max_width(); w++ ) {
			if( obj_type == Makespan ) {
				states[l][w] = new( cp.getHeap() ) DisjunctiveState(cp, acts.size(), 0, obj_var.getMax());
			} else {
				// TODO: compute a relaxation for max end time here
				states[l][w] = new( cp.getHeap() ) DisjunctiveState(cp, acts.size(), 0, INF);
			}
		}
	}

	// initialize node sketches
	sketches = new( cp.getHeap() ) DisjunctiveSketch*[mdd->get_max_arc_id()+1];
	for( int w = 0; w < mdd->get_max_arc_id()+1; w++ ) {
		sketches[w] = new( cp.getHeap() ) DisjunctiveSketch;
	}

	// initialize auxiliaries
	num_node_checked = new( cp.getHeap() ) int[mdd->get_max_width()];
	//nodes_checked = new( cp.getHeap() ) SparseSet()

	min_est = new( cp.getHeap() ) int[mdd->get_max_width()];
	min_size_implied = new( cp.getHeap() ) int[mdd->get_max_width()];

	node_infeasible = new( cp.getHeap() ) bool[mdd->get_max_width()];
	act_checked = new( cp.getHeap() ) bool[acts.size()];

	// allocate and compute pre-processed data
	min_setup_from = new( cp.getHeap() ) int[acts.size()];
	min_setup_to = new( cp.getHeap() ) int[acts.size()];

    max_setup_from = new( cp.getHeap() ) int[acts.size()];
    max_setup_to = new( cp.getHeap() ) int[acts.size()];


	// todo: possibly disconsider depot / terminal
	for( unsigned int i = 0; i < acts.size(); i++ ) {
		min_setup_from[i] = INF;
		min_setup_to[i] = INF;
		max_setup_from[i] = -1;
		max_setup_to[i] = -1;
		for( unsigned int j = 0; j < acts.size(); j++ ) {
			if( i != j ) {
				min_setup_from[i] = MIN(min_setup_from[i], setup[i][j]);
				min_setup_to[i] = MIN(min_setup_to[i], setup[j][i]);

				if( setup[i][j] != INF ) {
				  max_setup_from[i] = MAX(max_setup_from[i], setup[i][j]);
				}

                if( setup[j][i] != INF ) {
                  max_setup_to[i] = MAX(max_setup_to[i], setup[j][i]);
                }
			}
		}
	}

	has_previous = new( cp.getHeap() ) bool[acts.size()];
	memset(has_previous, false, acts.size()*sizeof(bool));

	has_next = new( cp.getHeap() ) bool[acts.size()];
	memset(has_next, false, acts.size()*sizeof(bool));

	previous_act = new( cp.getHeap() ) int[acts.size()];
	next_act = new( cp.getHeap() ) int[acts.size()];

	precedence_status = new( cp.getHeap() ) int*[acts.size()];
	for( unsigned int i = 0; i < acts.size(); ++i ) {
		precedence_status[i] = new( cp.getHeap() ) int[acts.size()];
		memset( precedence_status[i], P_NOTCOMPUTED, sizeof(int)*acts.size() );
	}

	mem_min_penalty = new( cp.getHeap() ) int*[mdd->get_num_layers()];
	mem_min_tardy = new( cp.getHeap() ) int*[mdd->get_num_layers()];
	for( int l = 0; l < mdd->get_num_layers(); ++l ) {
		mem_min_penalty[l] = new( cp.getHeap() ) int[mdd->get_max_arc_id()+1];
		mem_min_tardy[l] = new( cp.getHeap() ) int[mdd->get_max_arc_id()+1];
	}

	q_min_topdown = new( cp.getHeap() ) int*[mdd->get_num_layers()];
	q_max_topdown = new( cp.getHeap() ) int*[mdd->get_num_layers()];
	lagrangian_bound_topdown = new( cp.getHeap() ) double*[mdd->get_num_layers()];
	lagrangian_bound_topdown = new( cp.getHeap() ) double*[mdd->get_num_layers()];
	
	for( int l = 0; l < mdd->get_num_layers(); ++l ) {
		q_min_topdown[l] = new( cp.getHeap() ) int[mdd->get_max_arc_id()+1];
		q_max_topdown[l] = new( cp.getHeap() ) int[mdd->get_max_arc_id()+1];
		lagrangian_bound_topdown[l] = new( cp.getHeap() ) double[mdd->get_max_arc_id()+1];
	}
	
	lagr_cost_task 	= new( cp.getHeap() ) IlcRevFloat[mdd->get_num_layers()];
	lagr_cost_layer = new( cp.getHeap() ) IlcRevFloat[mdd->get_num_layers()];
	lagr_cost_pre 	= new( cp.getHeap() ) IlcRevFloat[mdd->get_num_layers()];

	lagr_cost_task_static = new( cp.getHeap() ) double[mdd->get_num_layers()];
	lagr_cost_layer_static = new( cp.getHeap() ) double[mdd->get_num_layers()];
	lagr_cost_pre_static = new( cp.getHeap() ) double[mdd->get_num_layers()];

	sp_lagr = new( cp.getHeap() ) int[mdd->get_num_layers()];
	sp_normal = new( cp.getHeap() ) int[mdd->get_num_layers()];

	temp_set.resize(acts.size(), false);
	aux_v.reserve(acts.size());
	may_be_refined = new( cp.getHeap() ) bool[mdd->get_max_width()];
	node_checked = new( cp.getHeap() ) bool[mdd->get_max_width()];

	exact_acts_array = new( cp.getHeap() ) int[acts.size()];
	exact_acts.reserve(acts.size());

	act_est = new( cp.getHeap() ) SetMinValueList(cp, acts.size());
	act_lft = new( cp.getHeap() ) SetMaxValueList(cp, acts.size());

	is_act_exact_layer = new( cp.getHeap() ) IlcRevBool*[mdd->get_num_layers()];
	for( int l = 0; l < mdd->get_num_layers(); ++l ) {
	  is_act_exact_layer[l] = new( cp.getHeap() ) IlcRevBool[acts.size()];
	  for( int i = 0; i < (int)acts.size(); ++i ) {
	    is_act_exact_layer[l][i].setValue(cp, false);
	  }
	}

    is_act_exact = new( cp.getHeap() ) IlcRevBool[mdd->get_num_layers()];
    for( int i = 0; i < (int)acts.size(); ++i ) {
      is_act_exact[i].setValue(cp, false);
    }

    valid_acts.reserve( acts.size() );


    // Refinement
    requires_refinement = new( cp.getHeap() ) IlcRevBool[mdd->get_num_layers()];
    for( int i = 1; i < mdd->get_num_layers()-1; ++i ) {
    	requires_refinement[i].setValue(cp, true);
    }
    requires_refinement[0].setValue(cp, false);
    requires_refinement[mdd->get_num_layers()-1].setValue(cp, false);

    node_ref_state = new( cp.getHeap() ) RefineState[mdd->get_max_width()];

	available_nodes.reserve(mdd->get_max_width());
	taken_nodes.reserve(mdd->get_max_width());

	cout << "[Disjunctive] done.\n\n";
}



/**
 * Set initial activity states
 * TODO: do this at every search node?
 **/
void DisjunctivePropagator::set_initial_states() {

	cout << "[Disjunctive] Warning: no initial state set" << endl;
	cout << endl;
	return;
//
//	assert( (int)acts.size() + 1 == get_mdd()->get_num_layers() );

//	IloCP cp = get_CP();
//	RevMDD* mdd = get_mdd();

	cout << endl;
	cout << "[Disjunctive] Setting initial activity states..." << endl;

	int num_acts = acts.size();

	vector<int> est_r1(num_acts);
	vector<int> lft_r1(num_acts);
	vector<int> process_times(num_acts);
	int min_release = INF;
	int max_deadline = -1;
	vector<SpanningTree::Edge> travel_times;

	for( int i = 0; i < num_acts; ++i ) {
		est_r1[i] = acts[i]->release + acts[i]->processing;
		lft_r1[i] = acts[i]->deadline - acts[i]->processing;
		process_times[i] = acts[i]->processing;
		min_release = MIN(min_release, acts[i]->release);
		max_deadline = MAX(max_deadline, acts[i]->deadline);

		for( int j = i+1; j < num_acts; ++j ) {
			SpanningTree::Edge edge;
			edge.u = i;
			edge.v = j;
			edge.val = MIN( setup[i][j], setup[j][i] );
			travel_times.push_back(edge);
		}
	}

	// compute spanning tree relaxation
	SpanningTree spt(num_acts, travel_times);
	vector<int> tree = spt.compute();

	// compute relaxed states
	sort(est_r1.begin(), est_r1.end());
	sort(lft_r1.begin(), lft_r1.end());
	sort(process_times.begin(), process_times.end());

	DisjunctiveState* state;

	// EST relaxation
	int est_r2 = min_release;
	for( int k = 0; k < num_acts; ++k) {

		// state is always the single node of the (k+1)-th layer
		state = get_node_state(k+1, 0);

		// update relaxations
		est_r2 += process_times[k];
		if( k > 0 ) {
			est_r2 += travel_times[tree[k-1]].val;
		}

		// compute relaxation
		int est = MAX(est_r1[k], est_r2);
		//cout << "EST - layer " << k+1 << " --> " << est << endl;
		state->earliest_start.setValue(get_CP(), est);
	}

	// LFT relaxation
	if( max_deadline < INF ) {
		int lft_r2 = max_deadline;
		for( int k = num_acts-1; k >= 0; --k) {

			// state is always the single node of the k-th layer
			state = get_node_state(k, 0);

			// update relaxations
			lft_r2 -= process_times[num_acts-1-k];

			// compute relaxation
			int lft = MIN(lft_r1[k], lft_r2);
			//cout << "LFT: layer " << k << " --> " << lft << endl;
			state->latest_finish.setValue(get_CP(), lft);
		}
	}

	cout << "[Disjunctive] done.\n\n";
	return;

//	cout << "Spanning tree: " << endl;
//	for( int i = 0; i < tree.size(); ++i ) {
//		cout << "\t" << tree[i] << " --> ";
//		cout << travel_times[tree[i]].u << "," << travel_times[tree[i]].v;
//		cout << "," << travel_times[tree[i]].val << endl;
//	}
//
//	exit(1);

//	//cout << "\tWarning: this initial state assumes all activities are subject to an initial setup cost" << endl;
//
//	// collect adjusted start and finish
//	vector<int> late_finish( acts.size() );
//	vector<int> early_start( acts.size() );
//	int max_deadline = -1;
//	for( size_t i = 0; i < acts.size(); ++i ) {
//
//		// todo: is it possible to consider setup times for calculating the early start as well?
//		early_start[i] = acts[i]->release + acts[i]->processing;
//		late_finish[i] = acts[i]->deadline - acts[i]->processing - min_setup_to[i];
//
//		max_deadline = MAX(max_deadline, acts[i]->deadline);
//	}
//
//	// Sort late finish vector. The resulting array is such that 'late_finish[i]'
//	// correspond to an upper bound of the the earliest finish time of an activity
//	// if at least i+1 activities are performed last on the resource. An analogous
//	// reasoning can be considered for the early start
//
//	sort( late_finish.rbegin(), late_finish.rend() );
//	sort( early_start.begin(), early_start.end() );
//
//	//	for( int i = 0; i < (int)acts.size(); i++ ) {
//	//		cout << i << " - processing: " << early_finish[i] << endl;
//	//	}
//
//	// terminal node has the maximum deadline
//	get_node_state(mdd->get_num_layers()-1, 0)->latest_finish.setValue(cp, max_deadline);
//
//	// late finish (assumes only one node per layer in the initial MDD)
//	for( size_t i = 0; i < acts.size(); i++ ) {
//		get_node_state(mdd->get_num_layers()-i-2, 0)->latest_finish.setValue(cp, late_finish[i]);
//		//cout << "layer: " << mdd->get_num_layers()-i-2 << " :: i=" << i << ":  late finish: " << late_finish[i] << endl;
//	}
//
//	// early start (assumes only one node per layer in the initial MDD)
//	for( size_t i = 0; i < acts.size(); i++ ) {
//		get_node_state(i+1, 0)->earliest_start.setValue(cp, early_start[i]);
//		cout << "layer: " << i+1 << " :: i=" << i << ":  early start: " << early_start[i] << endl;
//	}
//
//	cout << "[Disjunctive] done." << endl << endl;
}



/**
 * Filter layer in a top-down perspective
 **/
void DisjunctivePropagator::filter_layer_topdown(int layer) {

	//	cout << endl << "[Disjunctive] Topdown filtering of layer " << layer << endl;
	//	cout << "- obj. function: " << makespan << endl << endl;

	RevMDD* mdd = get_mdd();
	IloCP cp = get_CP();

//	// initialize infeasible node marker
//	memset(node_infeasible, false, sizeof(bool)*mdd->get_max_width());
//	num_infeas_node = 0;

	// pointer to the previous arc considered, for removal purposes
	RevMDD::Arc* previous = mdd->arcs_layer[layer].get_header();

	// remove infeasible arcs from layer
	for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
			arc != NULL;
			arc = mdd->arcs_layer[layer].get_next(arc) )
	{
	  if( /*node_infeasible[arc->get_source()] ||*/ is_incoming_infeasible(layer, arc) ) {

			// arc is infeasible: remove from layer
			previous->next.setValue(cp, arc->next.getValue());

		} else {
			previous = arc;
		}
	}

//	// perform one additional pass to remove arcs from infeasible nodes
//	if( num_infeas_node > 0 ) {
//		previous = mdd->arcs_layer[layer].get_header();
//		for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
//				arc != NULL;
//				arc = mdd->arcs_layer[layer].get_next(arc) )
//		{
//			if( node_infeasible[arc->get_source()] ) {
//
//				// arc is infeasible: remove from layer
//				previous->next.setValue(cp, arc->next.getValue());
//
//			} else {
//				previous = arc;
//			}
//		}
//	}
}


/**
 * Filter layer in a bottom up perspective
 **/
void DisjunctivePropagator::filter_layer_bottomup(int layer) {

	//cout << endl << "[Disjunctive] Bottomup filtering of layer " << layer << endl;

	RevMDD* mdd = get_mdd();
	IloCP cp = get_CP();

//	// initialize infeasible node marker
//	memset(node_infeasible, false, sizeof(bool)*mdd->get_max_width());
//	num_infeas_node = 0;

	// pointer to the previous arc considered, for removal purposes
	RevMDD::Arc* previous = mdd->arcs_layer[layer].get_header();

	// remove infeasible arcs from layer
	for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
			arc != NULL;
			arc = mdd->arcs_layer[layer].get_next(arc) )
	{
	  if( /*node_infeasible[arc->get_target()] ||*/ is_outgoing_infeasible(layer, arc) ) {
			// arc is infeasible: remove from layer
			previous->next.setValue(cp, arc->next.getValue());
		} else {
			previous = arc;
		}
	}

//	// perform one additional pass to remove arcs from infeasible nodes
//	if( num_infeas_node > 0 ) {
//		previous = mdd->arcs_layer[layer].get_header();
//		for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
//				arc != NULL;
//				arc = mdd->arcs_layer[layer].get_next(arc) )
//		{
//			if( node_infeasible[arc->get_target()] ) {
//
//				// arc is infeasible: remove from layer
//				previous->next.setValue(cp, arc->next.getValue());
//
//			} else {
//				previous = arc;
//			}
//		}
//	}
//	//cout << endl << "[Disjunctive] Bottomup filter end" << endl;
}


/**
 * Partition refinement group (refinement). Max MDD width must be observed
 **/
void DisjunctivePropagator::partition_refinement_group(
		int prev_layer,
		RefinementGroup *ref_groups,
		int &n_groups,
		int n_rounds)
{
	int layer = prev_layer+1;
	if( !requires_refinement[layer].getValue() ) {
		//cout << "\tNo refinement required" << endl;
		return;
	}

	// full extension refinement
	//partition_full_extension(prev_layer, ref_groups, n_groups);
	//if( n_groups == get_mdd()->get_max_width() || !requires_refinement[layer].getValue() ) {
	//	return;
	//}

	// partition_exact_capacity

	if(cap_refinement && n_rounds <= 4) {
		partition_exact_capacity(prev_layer, ref_groups, n_groups);
		//cout << "Using Capacity refinement..." << endl;
	}
	
	// alldiff-based refinement
	block_est_refinement = false;
	partition_exact_acts(prev_layer, ref_groups, n_groups);
	if( n_groups == get_mdd()->get_max_width() || block_est_refinement ) {
		return;
	}

	// EST-based refinement
	//partition_est(prev_layer, ref_groups, n_groups);
}


/**
 * Partition refinement group based on full extension
 **/
void DisjunctivePropagator::partition_full_extension(
		int layer,
		RefinementGroup *ref_groups,
		int &n_groups)
{
	RevMDD* mdd = get_mdd();

	// compute how many extra nodes we require to perform full
	// extension of the layer

	int extra_nodes = 0;

	available_nodes.clear();
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		if( ref_groups[g]->size == 0 ) {
			available_nodes.push_back(g);
		} else {
			extra_nodes += ref_groups[g]->size - 1;
		}
	}

	if( extra_nodes <= (int)available_nodes.size() ) {
		// perform full extension; layer does not require any more refinement
		for( int g = 0; g < mdd->get_max_width(); g++ ) {
			while( ref_groups[g]->size > 1 ) {
				assert( available_nodes.size() > 0 );

				// transfer arc to new node
				ref_groups[available_nodes.back()]->add( ref_groups[g]->pop_up() );

				// update available nodes and number of groups
				available_nodes.pop_back();
				n_groups++;
			}
		}
		// previous layer must not require refinement as well
		if( !requires_refinement[layer].getValue() ) {
			requires_refinement[layer+1].setValue(get_CP(), false);
		}
	}
}

/**
 * Partition refinement group based on capacity
 **/
void DisjunctivePropagator::partition_exact_capacity(
		int layer,
		RefinementGroup *ref_groups,
		int &n_groups)
{
	RevMDD* mdd = get_mdd();
	
	available_nodes.clear();
	// create comparator for groups according to earliest start time (EST)
	Arc_Q_Comparator group_comp(sketches);

	//cout << "== Layer "<< layer <<" ==" << endl;
	
	//Store available nodes and sort arcs in each refinement group by Q_min
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		if( ref_groups[g]->is_empty() ) {
			available_nodes.push_back(g);
		}
		else{
		
			if( ref_groups[g]->size > 1 ) {
				// sort group
				sort(ref_groups[g]->elements, ref_groups[g]->elements+ref_groups[g]->size, group_comp);
				
				//Check if the sorting is ok
				for(int i = 0; i < ref_groups[g]->size; ++i){
					//cout << "arc " << i << "\t Q_min = " << sketches[ref_groups[g]->elements[i]->get_layer_id()]->min_Q; 
					//cout << "\t Q_max = " << sketches[ref_groups[g]->elements[i]->get_layer_id()]->max_Q << endl;
				}

			}
		}
	}
	
	//Allocate arcs in available nodes
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		if( ref_groups[g]->size > 1){

			while (true){
				if( sketches[ref_groups[g]->get_first()->get_layer_id()]->min_Q ==
					sketches[ref_groups[g]->get_last()->get_layer_id()]->min_Q  || available_nodes.empty() ){  // no more different arc weight --> stop!
						break;
				}
				
				else{
					int aux_weight = sketches[ref_groups[g]->get_last()->get_layer_id()]->min_Q;
					//cout << "New node has weight =  " << aux_weight << " with arcs: ";
					
					while(sketches[ref_groups[g]->get_last()->get_layer_id()]->min_Q == aux_weight){
						//cout << ref_groups[g] -> size -1 << " ";
						ref_groups[available_nodes.back()]->add( ref_groups[g]->pop_up() );
						
					}
					//cout << endl;
					
					available_nodes.pop_back();
					//cout << "Number of available nodes = " << available_nodes.size() << endl;
				}
			}
		}
	}
}

/**
 * Partition refinement group based on exact activities
 **/
void DisjunctivePropagator::partition_exact_acts(
		int layer,
		RefinementGroup *ref_groups,
		int &n_groups)
{

	//  static int iter=0;
	//  cout << "\n\n\nIteration " << iter++ << endl;
	//  if( iter > 100 )
	//    exit(1);

	RevMDD* mdd = get_mdd();
	assert( n_groups < mdd->get_max_width() );

	//  cout << "\n\n\n------------------------" << endl;
	//  cout << "Layer " << layer << endl;
	//  cout << endl;

	// ------------------------------------------------
	// Initialization: Create comparator according to
	// of EST (with size of all-paths as tie break).
	// Also obtain available nodes
	// ------------------------------------------------

	// reset node availability
	available_nodes.clear();
	taken_nodes.clear();
	//memset(may_be_refined, false, sizeof(bool)*acts.size());

	// create comparator for groups according to earliest start time (EST)
	Arc_Est_Comparator group_comp(sketches);

	// sort each group according to lower bound
	for( int g = 0; g < mdd->get_max_width(); g++ ) {
		if( ref_groups[g]->is_empty() ) {
			// if group is empty, add to available nodes
			available_nodes.push_back(g);

		} else {

			taken_nodes.push_back(g);
			if( ref_groups[g]->size > 1 ) {
				// compute min_est
				min_est[g] = INF;
				for( int i = 0; i < ref_groups[g]->size; ++i ) {
					min_est[g] = MIN( min_est[g],
							sketches[ref_groups[g]->elements[i]->get_layer_id()]->earliest_start );
				}
				//may_be_refined[g] = true;
				//cout << "\n\tgroup may be refined: ";
				//print_group(g, ref_groups[g]);
			}
		}
	}

//	  cout << "Available nodes: ";
//	  for( int i = 0; i < (int)available_nodes.size(); ++i ) {
//	    cout << available_nodes[i] << " ";
//	  }
//	  cout << endl;

	// -------------------------------------
	// Refinement
	// TODO: queue can be pre-allocated
	// -------------------------------------

	// create priority queue of nodes that must be refined. Target nodes are
	// ordered according to minimum EST.
	// These values are stored dynamically in the vectors 'min_est' and 'size_implied'
	GroupESTComp node_comp(min_est);
	priority_queue<int, vector<int>, GroupESTComp> refine_set(node_comp);

	// marks if layer is exact with respect to all diff
	//bool is_exact_alldiff = true;

	// Main loop: iterate on activity that must be exact in layer. Priority is given to
	// order defined in 'exact_acts'
	for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {

		int act = *it;
		if( is_act_exact[act].getValue() || is_act_exact_layer[layer][act].getValue() ) {
			continue;
		}

		// TODO: take list of nodes from previous set
		memset(node_ref_state, Unknown, sizeof(RefineState)*mdd->get_max_width());

		bool has_no_refinable = false;

		// create refine set according to target nodes
		for( int g = 0; g < mdd->get_max_width(); ++g ) {

//			cout << "\nChecking eligibility of group " << g << " -- act: " << act << endl;
//			print_group(g, ref_groups[g]);

			RefineState state = is_refinable(layer, act, ref_groups[g]);
			if( state == Refinable ) {
				//cout << "\tgroup " << g << " is refinable" << endl;
				refine_set.push(g);
			} else if( state == NotRefinable ) {
				has_no_refinable = true;
			}
		}

		if( has_no_refinable ) {
			block_est_refinement = true;
		} else {
			if( refine_set.empty() ) {
				is_act_exact_layer[layer][act].setValue(get_CP(), true);
			}
		}
		
		if(available_nodes.empty()){
			block_est_refinement = true;
			return;
		}

		// keep refining until possible
		while( !refine_set.empty() ) {

			// take group with smallest EST from refine list
			int to_refine = refine_set.top();
			refine_set.pop();

			//cout << "\t\tselected " << to_refine << " - min_est: " << min_est[to_refine] << endl;

			// take an available node
			int new_group = available_nodes.back();
			available_nodes.pop_back();
			taken_nodes.push_back(new_group);

			// transfer arcs that do not contain 'act' in their implied set to new_group
			RevMDD::Arc* arc;
			for( int i = 0; i < ref_groups[to_refine]->size; ) {
				arc = ref_groups[to_refine]->elements[i];
				if( arc->get_val() != act && !get_node_state(layer, arc->get_source())->all_paths[act] ) {
					ref_groups[new_group]->add(arc);
					ref_groups[to_refine]->elements[i] = ref_groups[to_refine]->pop_up();
				} else {
					++i;
				}
			}
			assert( ref_groups[to_refine]->size > 0 );
			assert( ref_groups[new_group]->size > 0 );

			//cout << "sizeA: " << ref_groups[to_refine]->size  << endl;
			//cout << "sizeB: " << ref_groups[new_group]->size  << endl;

			if( ref_groups[to_refine]->size == 0 || ref_groups[new_group]->size == 0 ) {
				cout << "ERROR: refinement group has no arcs" << endl;
				exit(1);
			}

			//print_group(to_refine, ref_groups[to_refine]);
			//print_group(new_group, ref_groups[new_group]);

			// update number of groups
			n_groups++;

			// if width is met, refinement is concluded
			if( available_nodes.size() == 0 ) {
				block_est_refinement = true;
				return;
			}

			// update min_est of modified nodes
			if( ref_groups[new_group]->size > 1 ) {
				min_est[new_group] = INF;
				for( int i = 0; i < ref_groups[new_group]->size; ++i ) {
					min_est[new_group] = MIN( min_est[new_group],
							sketches[ref_groups[new_group]->elements[i]->get_layer_id()]->earliest_start );
				}
			}

			if( ref_groups[to_refine]->size > 1 ) {
				min_est[to_refine] = INF;
				for( int i = 0; i < ref_groups[to_refine]->size; ++i ) {
					min_est[to_refine] = MIN( min_est[to_refine],
							sketches[ref_groups[to_refine]->elements[i]->get_layer_id()]->earliest_start );
				}
			}
		}
	}
}


/**
 * Partition refinement group based exclusively on EST
 **/
void DisjunctivePropagator::partition_est(int layer, RefinementGroup *ref_groups, int &n_groups) {


	RevMDD* mdd = get_mdd();
	assert( n_groups < mdd->get_max_width( ));

	/**
	 * Check if it is possible to create one node for each arc sketch
	 * in the refinement group.
	 */

	// Count total number of sketches and the number of available nodes
	available_nodes.clear();
	int total_size = 0;
	for( int g = 0; g < mdd->get_max_width(); g++ ) {

		if( ref_groups[g]->is_empty() ) {
			// if group is empty, add to available nodes
			available_nodes.push_back(g);

		} else {
			// if group is not empty, count how many sketches we have to move
			if( ref_groups[g]->size > 1 ) {
				total_size += ref_groups[g]->size - 1;
			}
		}
	}

//	// If there is available space, create one node for each arc sketch
//	if( total_size <= (int) available_nodes.size() ) {
//		//cout << "Also full extension as well!!!" << endl;
//		for( int g = 0; g < mdd->get_max_width(); g++ ) {
//			while( ref_groups[g]->size > 1 ) {
//				assert( available_nodes.size() > 0 );
//
//				// transfer arc to new node
//				ref_groups[available_nodes.back()]->add( ref_groups[g]->pop_up() );
//
//				// update available nodes and number of groups
//				available_nodes.pop_back();
//				n_groups++;
//			}
//		}
//		return;
//	}

	// TODO: test this!
	//return;

	/*
	 * ===============================================
	 *  Create priority queue of groups to refine
	 * ===============================================
	 */

	RefGroup_Est_Comparator refgroup_comparator(min_est, min_size_implied);
	priority_queue<int, vector<int>, RefGroup_Est_Comparator> group_queue(refgroup_comparator);

	// reset node availability
	available_nodes.clear();

	// create comparator for groups according to earliest start time (EST)
	Arc_Est_Comparator group_comp(sketches);

	// sort each group according to lower bound
	for( int g = 0; g < mdd->get_max_width(); g++ ) {

		if( ref_groups[g]->is_empty() ) {

			// if group is empty, add to available nodes
			available_nodes.push_back(g);

		} else {

			// only consider groups with more than one arc, since they are the ones that can be
			// refined

			if( ref_groups[g]->size > 1 ) {

				// sort group
				sort(ref_groups[g]->elements, ref_groups[g]->elements+ref_groups[g]->size, group_comp);
				//min_dist[g] = get_min_dist_R(ref_groups[g]);
				min_est[g] = get_min_est(ref_groups[g]);
				min_size_implied[g] = get_min_all_paths(ref_groups[g]);

				// print
				//print_group(g, ref_groups[g]);

				// add group to queue if min_size of size implied are distinct
				if( is_est_distinct(ref_groups[g]) || is_min_all_paths_distinct(ref_groups[g]) ) {
					group_queue.push(g);

					//cout << "\t\tadded group to queue" << endl;

				}

			}
		}
	}

	/*
	 * =====================================================
	 *  Refine until width is met or no more partitions
	 *  are possible
	 * =====================================================
	 */

	int to_refine, new_group;
	while( n_groups < mdd->get_max_width() && !group_queue.empty() ) {

		// select group to refine
		to_refine = group_queue.top();
		group_queue.pop();

		//cout << "\n\tGroup to refine: " << to_refine << endl;

		// take new node
		assert( !available_nodes.empty() );
		new_group = available_nodes.back();
		available_nodes.pop_back();

		// increment number of existing groups
		n_groups++;

		// auxiliary
		RevMDD::Arc* arc;

		if( is_est_distinct(ref_groups[to_refine]) ) {

			/*
			 * Refinement by min_dist_R
			 */

			// transfer arcs
			while( sketches[ref_groups[to_refine]->get_last()->get_layer_id()]->earliest_start == min_est[to_refine] ) {

				// take transfer arc
				arc = ref_groups[to_refine]->pop_up();

				// mark source node as modified
				mdd->modified_out_node(layer, arc->get_source());

				// add to new group
				ref_groups[new_group]->add( arc );
			}

			assert( ref_groups[to_refine]->size > 0 );
			assert( ref_groups[new_group]->size > 0 );

		} else {

			/*
			 * Refinement by implied up
			 */

			// transfer arcs
			while( sketches[ref_groups[to_refine]->get_last()->get_layer_id()]->size_implied == min_size_implied[to_refine] ) {

				// take transfer arc
				arc = ref_groups[to_refine]->pop_up();

				// mark source node as modified
				mdd->modified_out_node(layer, arc->get_source());

				// add to new group
				ref_groups[new_group]->add( arc );
			}

			assert( ref_groups[to_refine]->size > 0 );
			assert( ref_groups[new_group]->size > 0 );
		}

		// sort new group
		sort(ref_groups[new_group]->elements, ref_groups[new_group]->elements+ref_groups[new_group]->size, group_comp);

		// update groups sketch info

		min_est[to_refine] = get_min_est(ref_groups[to_refine]);
		min_size_implied[to_refine] = get_min_all_paths(ref_groups[to_refine]);

		min_est[new_group] = get_min_est(ref_groups[new_group]);
		min_size_implied[new_group] = get_min_all_paths(ref_groups[new_group]);


		//		// print
		//		print_group(to_refine, ref_groups[to_refine]);
		//		print_group(new_group, ref_groups[new_group]);

		// add groups to queue if they can still be refined
		if( is_est_distinct(ref_groups[new_group]) || is_min_all_paths_distinct(ref_groups[new_group]) ) {
			group_queue.push(new_group);

			//			cout << "\t\tadded group " << new_group << " to queue" << endl;

		}

		//cout << endl;

		if( is_est_distinct(ref_groups[to_refine]) || is_min_all_paths_distinct(ref_groups[to_refine]) ) {
			group_queue.push(to_refine);

			//cout << "\t\tadded group " << to_refine << " to queue" << endl;

		}

		//cout << endl;
	}

	//cout << "refinement completed" << endl << endl;
	//exit(1);
}



/**
 * Set the state of a node according to its incoming arcs.
 **/
void DisjunctivePropagator::setup_node_incoming(int layer, int node, RefinementGroup ref_group, bool new_node) {

	RevMDD* mdd = get_mdd();
	IloCP cp = get_CP();

	// initialize state of the new node
	DisjunctiveState* state = get_node_state(layer, node);

	state->some_paths.reset();
	state->all_paths.set();
	state->in_est->clear();
	state->in_Q_min->clear();
	state->in_Q_max->clear();
	state->min_lagrangian_bound_R->clear();

	switch( obj_type ) {
	case SumSetupTimes:
		state->in_setup->clear();
		break;
	case TotalTardiness:
		state->in_tard->clear();
		break;
	case Latency:
		state->in_latency->clear();
		break;
	default:
		break;
	}

	// if new node, we initialize the disjunctive states with the original node states
	if( new_node ) {
		DisjunctiveState* original_state = get_node_state(layer, ref_group->get_first()->get_target());
		state->latest_finish.setValue(cp, original_state->latest_finish.getValue());

		switch( obj_type ) {
		case SumSetupTimes:
			state->min_sumsetup_T.setValue(cp, original_state->min_sumsetup_T.getValue());
			break;
		case TotalTardiness:
			state->max_tard.setValue(cp, original_state->max_tard.getValue());
			break;
		case Latency:
			state->max_latency.setValue(cp, original_state->max_latency.getValue());
			break;
		default:
			break;
		}
	}

	// initialize auxiliaries
	memset( num_node_checked, 0, sizeof(int)*mdd->get_max_width() );

	int min_earliest_start = INF;
	int min_sum_setup_R = INF;
	int min_tard_R = INF;
	int min_lat = INF;
	double sum_lagrangian = INF;

	RevMDD::Arc* arc;
	for( int i = 0; i < ref_group->size; i++ ) {
		arc = ref_group->elements[i];

		// update start time info
		min_earliest_start = MIN( min_earliest_start, sketches[arc->get_layer_id()]->earliest_start );

		// update list of incoming arcs
		state->in_est->add(arc->get_val(), sketches[arc->get_layer_id()]->earliest_start);

		// update capacities
		state->in_Q_min->add(arc->get_val(), sketches[arc->get_layer_id()]->min_Q);
		state->in_Q_max->add(arc->get_val(), sketches[arc->get_layer_id()]->max_Q);
		
		//update lagrangian bounds at each state
		state->min_lagrangian_bound_R->add(arc->get_val(), sketches[arc->get_layer_id()]->min_lagrangian_bound_R);

		sum_lagrangian = MIN(sum_lagrangian, sketches[arc->get_layer_id()]->min_lagrangian_bound_R);

		switch( obj_type ) {
		case SumSetupTimes:
			min_sum_setup_R = MIN( min_sum_setup_R, sketches[arc->get_layer_id()]->min_sumsetup_R );
			state->in_setup->add(arc->get_val(), sketches[arc->get_layer_id()]->min_sumsetup_R);
			break;
		case TotalTardiness:
			min_tard_R = MIN( min_tard_R, sketches[arc->get_layer_id()]->min_tard );
			state->in_tard->add(arc->get_val(), sketches[arc->get_layer_id()]->min_tard);
			break;
		case Latency:
			min_lat = MIN( min_lat, sketches[arc->get_layer_id()]->min_latency );
			state->in_latency->add(arc->get_val(), sketches[arc->get_layer_id()]->min_latency);
			break;
		default:
			break;
		}

		// some_paths update: value will always belong to state
		state->some_paths.set(arc->get_val(), true);

		// check how many times node was verified
		switch( num_node_checked[arc->get_source()] ) {

		case 0:
			// update 'some_paths' state
			state->some_paths |= get_node_state(layer-1, arc->get_source())->some_paths;

			// update 'all_paths' state
			if( state->all_paths[arc->get_val()] ) {
				state->all_paths &= get_node_state(layer-1, arc->get_source())->all_paths;
				state->all_paths.set(arc->get_val(), true);

			} else {
				state->all_paths &= get_node_state(layer-1, arc->get_source())->all_paths;
			}

			num_node_checked[arc->get_source()]++;
			break;

		case 1:
			// just intersect 'all_paths'
			state->all_paths &= get_node_state(layer-1, arc->get_source())->all_paths;
			num_node_checked[arc->get_source()]++;

			break;

		default:
			break;

		}

	}

	// verify if some_paths defines a hall set
	state->hall_set = ((int)state->some_paths.count() == layer);

	// set earliest start state
	state->earliest_start.setValue(cp, min_earliest_start);

	// set additional states
	switch( obj_type ) {
	case Makespan:
		if( layer == mdd->get_num_layers()-1 ) {
			//cout << "\tRelaxation: " << state->earliest_start.getValue() << " - makespan: " << obj_var << endl;
			if( state->earliest_start.getValue() > obj_var.getMin() ) {
				//exit(1);
				cp.add( obj_var >= state->earliest_start.getValue() );
			}
		}
		break;

	case SumSetupTimes:
		state->min_sumsetup_R.setValue(cp, min_sum_setup_R);
		if( layer == mdd->get_num_layers()-1 ) {
			/*
			cout << "[Top Down]";
			cout <<	"\tCurrent bound = " << obj_var.getMin();
			cout <<	"\tnormal bound = " << state->min_sumsetup_R.getValue();
			cout <<	"\tlagrangian bound = " << (int)sum_lagrangian << endl;
			*/
			if( state->min_sumsetup_R.getValue() > obj_var.getMin() ) {
				cp.add( obj_var >= state->min_sumsetup_R.getValue() );
				//cout <<"* Improved! (normal) --> New bound = " << state->min_sumsetup_R.getValue() << "\n";
			}
			if( sum_lagrangian >0.0 ) {
				if( (int) sum_lagrangian > obj_var.getMin() ){
					cp.add( obj_var >= (int) sum_lagrangian );
					//cout <<"* Improved! (lagrangian)--> New bound = " << (int) sum_lagrangian << "\n";
				}
			}
		}
		break;

	case TotalTardiness:
		state->min_tard.setValue(cp, min_tard_R);
		if( layer == mdd->get_num_layers()-1 ) {
			if( state->min_tard.getValue() > obj_var.getMin() ) {
				//cout << "\tTardiness update: " << state->min_tard.getValue() << endl;
				cp.add( obj_var >= state->min_tard.getValue() );
			} else {
				//cout << "\tCurrent tardiness: " << obj_var << endl;
			}
		}
		break;

	case Latency:
		state->min_latency.setValue(cp, min_lat);
		if( layer == mdd->get_num_layers()-1 ) {
			if( state->min_latency.getValue() > obj_var.getMin() ) {
				//cout << "\tLatency update: " << state->min_lat.getValue() << endl;
				cp.add( obj_var >= state->min_latency.getValue() );
			} else {
				//cout << "\tCurrent latency: " << obj_var << endl;
			}
		}
		break;

	default:
		break;
	}

	//cout << "Incoming state: " << endl;
	//print_state(layer, node, state);
}




/**
 * Set the state of a node according to its outgoing arcs.
 **/
void DisjunctivePropagator::setup_node_outgoing(int layer, int node, RefinementGroup ref_group, bool new_node) {

	RevMDD* mdd = get_mdd();
	IloCP cp = get_CP();

	// initialize state of the new node
	DisjunctiveState* state = get_node_state(layer, node);

	state->some_paths.reset();
	state->all_paths.set();
	state->out_lft->clear();
	state->out_Q_min->clear();
	state->out_Q_max->clear();
	state->min_lagrangian_bound_T->clear();

	switch( obj_type ) {
	case SumSetupTimes:
		state->out_setup->clear();
		break;
	case TotalTardiness:
		state->out_tard->clear();
		break;
	case Latency:
		state->out_latency->clear();
		break;
	default:
		break;
	}

	int max_latest_finish = -1;
	int min_sum_setup_T = INF;
	int max_tard = -1;
	int max_lat = -1;
	int min_tard_T = INF;
	double sum_lagrangian= INF;

	// initialize auxiliaries
	memset( num_node_checked, 0, sizeof(int)*mdd->get_max_width() );

	RevMDD::Arc* arc;
	for( int i = 0; i < ref_group->size; i++ ) {
		arc = ref_group->elements[i];

		// update start time info
		max_latest_finish = MAX( max_latest_finish, sketches[arc->get_layer_id()]->latest_finish );

		// add val to outgoing set of node
		state->out_lft->add(arc->get_val(), sketches[arc->get_layer_id()]->latest_finish);

		// update capacities
		state->out_Q_min->add(arc->get_val(), sketches[arc->get_layer_id()]->min_Q);
		state->out_Q_max->add(arc->get_val(), sketches[arc->get_layer_id()]->max_Q);
		
		//update lagrangian bound at each state
		state->min_lagrangian_bound_T->add(arc->get_val(), sketches[arc->get_layer_id()]->min_lagrangian_bound_T, arc);

		sum_lagrangian= MIN(sum_lagrangian, sketches[arc->get_layer_id()]->min_lagrangian_bound_T);

		switch( obj_type ) {
		case SumSetupTimes:
			min_sum_setup_T = MIN( min_sum_setup_T, sketches[arc->get_layer_id()]->min_sumsetup_T );
			state->out_setup->add(arc->get_val(), sketches[arc->get_layer_id()]->min_sumsetup_T, arc);
			break; 
		case TotalTardiness:
			max_tard = MAX( max_tard, sketches[arc->get_layer_id()]->max_tard );
			min_tard_T = MIN( min_tard_T, sketches[arc->get_layer_id()]->min_tard_T );
			state->out_tard->add(arc->get_val(), sketches[arc->get_layer_id()]->max_tard);
			break;
		case Latency:
			max_lat = MAX( max_lat, sketches[arc->get_layer_id()]->max_latency );
			state->out_latency->add(arc->get_val(), sketches[arc->get_layer_id()]->max_latency);
			break;
		default:
			break;
		}

		// some_paths update: value will always belong to state
		state->some_paths.set(arc->get_val(), true);

		// check how many times node was verified
		switch( num_node_checked[arc->get_target()] ) {

		case 0:

			// update 'some_paths' state
			state->some_paths |= get_node_state(layer+1, arc->get_target())->some_paths;

			// update 'all_paths' state
			if( state->all_paths[arc->get_val()] ) {
				state->all_paths &= get_node_state(layer+1, arc->get_target())->all_paths;
				state->all_paths.set(arc->get_val(), true);

			} else {
				state->all_paths &= get_node_state(layer+1, arc->get_target())->all_paths;
			}

			num_node_checked[arc->get_target()]++;
			break;

		case 1:

			// just intersect 'all_paths'
			state->all_paths &= get_node_state(layer+1, arc->get_target())->all_paths;
			num_node_checked[arc->get_target()]++;

			break;

		default:
			break;
		}

	}

	// verify if some_paths defines a hall set
	state->hall_set = ((int)state->some_paths.count() == ( mdd->get_num_layers() - (layer+1) ) );

	// define states
	state->latest_finish.setValue(cp, max_latest_finish);

	// set additional states
	switch( obj_type ) {

	case SumSetupTimes:
		state->min_sumsetup_T.setValue(cp, min_sum_setup_T);
		if( layer == 0 ) {
			/*
			cout << "[Bottom UP]";
			cout <<	"\tCurrent bound = " << obj_var.getMin();
			cout <<	"\tnormal bound = " << state->min_sumsetup_T.getValue();
			cout <<	"\tlagrangian bound = " << (int)sum_lagrangian << endl;
			*/
			if( state->min_sumsetup_T.getValue() > obj_var.getMin() ) {
				cp.add( obj_var >= state->min_sumsetup_T.getValue() );
				//cout <<"* Improved! (normal) --> New bound = " << state->min_sumsetup_T.getValue() << "\n";
			}
			
			if(sum_lagrangian >0.0){
				if( (int)sum_lagrangian > obj_var.getMin() ){
					//cout <<"* Improved! (lagrangian)--> New bound = " << (int) sum_lagrangian << " Old bound = " << obj_var.getMin() <<" \n";
					cp.add( obj_var >= (int)sum_lagrangian);
				}
			}

			//cout << "relaxation provided by M: " << state->min_sumsetup_T.getValue() << endl;
			//cout << obj_var << endl;
//			exit(1);
		}
		break;

	case TotalTardiness:
		state->max_tard.setValue(cp, max_tard);
		state->min_tard_T.setValue(cp, min_tard_T);
		break;

	case Latency:
		state->max_latency.setValue(cp, max_lat);
		break;

	default:
		break;
	}

	if(layer == 0){
		save_shortest_path();
		save_shortest_path_lagr();
	}

//	cout << "Outgoing state: " << " - var: " << permutation[layer] << endl;
//	print_state(layer, node, state);
}



/**
 * Check if an in arc is not necessarily infeasible. It also
 * sets the corresponding arc sketch.
 *
 **/
bool DisjunctivePropagator::is_incoming_infeasible(int layer, RevMDD::Arc* arc) {
	if( is_incoming_infeasible_alldiff(layer, arc) ) {
		return true;
	}
	if( is_incoming_infeasible_disjunctive(layer, arc) ) {
		return true;
	}
	#if CAP_FILT > 0
	if( is_incoming_infeasible_capacity(layer, arc) ) {
		//cout << "remove arc due to capacity constraint- incoming" << endl;
		return true;
	}
	#endif
	
	is_incoming_infeasible_lagrangian_bound(layer, arc);
	
	return false;
	
}


/**
 * Check if an in arc is not necessarily infeasible according to alldiff. It also
 * sets the corresponding arc sketch.
 *
 **/
inline bool DisjunctivePropagator::is_incoming_infeasible_alldiff(int layer, RevMDD::Arc* arc) {

	DisjunctiveState* source_state = get_node_state(layer, arc->get_source());

	// --------------------------------------------------
	// All diff conditions
	// --------------------------------------------------

	// all_path condition
	if( source_state->all_paths[arc->get_val()] ) {
		return true;
	}

	// some paths condition
	if( source_state->hall_set && source_state->some_paths[arc->get_val()] )	{
		return true;
	}

	return false;
}

/**
 * Check if an in arc is not necessarily infeasible according to  the capacity. It also
 * sets the corresponding arc sketch.
 **/
inline bool DisjunctivePropagator::is_incoming_infeasible_capacity(int layer, RevMDD::Arc* arc) {

	DisjunctiveState* source_state = get_node_state(layer, arc->get_source());
	DisjunctiveSketch* arc_sketch = sketches[arc->get_layer_id()];

	if( layer == 0){
		arc_sketch->min_Q = 0;
		arc_sketch->max_Q = 0;
		q_min_topdown[layer][arc->get_layer_id()] = 0;
		q_max_topdown[layer][arc->get_layer_id()] = 0;
		return false;
	}

	int q_min = INF;
	int q_max = (-1)*capacity;

	for( int i = 0; i < source_state->in_Q_min->size; ++i ) {
	    int prev_city = source_state->in_Q_min->list[i];
	    //if( prev_city != arc->get_val()) {
	    	int q_min_src = source_state->in_Q_min->value(prev_city);
	    	int q_max_src = source_state->in_Q_max->value(prev_city);

	    	q_min = MIN(q_min, q_min_src + weights[arc->get_val()]);
	    	q_max = MAX(q_max, q_max_src + weights[arc->get_val()]);
	   //}
	}

	//cout << "Arc at (" << layer << "," << arc->get_source() << ") - value=" << arc->get_val();
	//cout << " :: [" << q_min << "," << q_max << "]" << endl;

	q_min = MAX(0, q_min);
	q_max = MIN(capacity, q_max);

	// check feasibility ----------------------------------------------------
	if (q_min > capacity){
		//cout << "Infeasible solution: Q_min greater than the capacity" << endl;
	 	return true;
	}
	 
	if (q_max < 0){
		//cout << "Infeasible solution: Q_max smaller than cero" << endl;
	 	return true;
	 }
	// ----------------------------------------------------------------------

	arc_sketch->min_Q = q_min;
	arc_sketch->max_Q = q_max;

	q_min_topdown[layer][arc->get_layer_id()] = arc_sketch->min_Q;
	q_max_topdown[layer][arc->get_layer_id()] = arc_sketch->max_Q;

	return false;
}

/**
 * Check if an in arc is not necessarily infeasible according to disjunctive. It also
 * sets the corresponding arc sketch.
 *
 **/
bool DisjunctivePropagator::is_incoming_infeasible_disjunctive(int layer, RevMDD::Arc* arc) {

	DisjunctiveState* source_state = get_node_state(layer, arc->get_source());


	// --------------------------------------------------
	// Disjunctive conditions
	// --------------------------------------------------

	/**
	 * 1. Compute arc 'est'
	 */

	//int min_setup_prev = INF;
	int est_arc = INF, est_arc_tmp;
	valid_acts.clear();

	if( layer > 0 ) {
	  int prev_value;
	  for( int i = 0; i < source_state->in_est->size; ++i ) {
	    prev_value = source_state->in_est->list[i];
	    if( prev_value != arc->get_val()
	        &&
	        !precedes_activity(arc->get_val(), prev_value) )
	    {
	      est_arc_tmp = source_state->in_est->value(prev_value)
	          + setup[prev_value][arc->get_val()];
	      est_arc = MIN(est_arc, est_arc_tmp);

	      valid_acts.push_back(prev_value);
	    }
	  }
	} else {
	  est_arc = acts[arc->get_val()]->start_time.getStartMin();
	}

	// if no activity can precede the current one, arc is infeasible
	if( est_arc == INF ) {
	  return true;
	}

	// take into account release and processing times
	est_arc = MAX(est_arc, acts[arc->get_val()]->start_time.getStartMin())
	    + acts[arc->get_val()]->processing;

	// compute sketch
	DisjunctiveSketch* sketch = sketches[arc->get_layer_id()];
    DisjunctiveState* target_state = get_node_state(layer+1, arc->get_target());
	sketch->earliest_start = MAX(est_arc, target_state->earliest_start.getValue());


	/**
	 * 2. Compute arc 'lft' and update lft
	 */

	// compute latest finish time
	sketch->latest_finish = target_state->latest_finish.getValue();
	sketch->latest_finish = MIN(sketch->latest_finish, acts[arc->get_val()]->start_time.getEndMax());
	if( sketch->earliest_start > sketch->latest_finish ) {
	  return true;
	}


	/**
	 * 3. Objective function propagator
	 */
	switch( obj_type ) {

	case Makespan:
		sketch->latest_finish = MIN(sketch->latest_finish, obj_var.getMax());
		if( sketch->earliest_start > sketch->latest_finish ) {
		  return true;
		}
		break;


	case SumSetupTimes:
		if( layer > 0 ) {
			sketch->min_sumsetup_R = INF;
		} else {
			sketch->min_sumsetup_R = 0;
		}
		for( vector<int>::iterator it = valid_acts.begin();
				it != valid_acts.end();
				++it )
		{
			sketch->min_sumsetup_R = MIN(sketch->min_sumsetup_R,
					source_state->in_setup->value(*it)
					+ setup[*it][arc->get_val()]);
		}
		sketch->min_sumsetup_R = MAX(sketch->min_sumsetup_R, target_state->min_sumsetup_R.getValue());
		if( sketch->min_sumsetup_R + target_state->min_sumsetup_T.getValue()
				> obj_var.getMax() )
		{
	      return true;
	    }
		break;


	case TotalTardiness:
	{
		if( layer > 0 ) {
			sketch->min_tard = INF;
		} else {
			sketch->min_tard = 0;
		}
		for( vector<int>::iterator it = valid_acts.begin();
				it != valid_acts.end();
				++it )
		{
			sketch->min_tard = MIN(source_state->in_tard->value(*it), sketch->min_tard);
		}
		int local_penalty = MAX(0, sketch->earliest_start - acts[arc->get_val()]->duedate)
				* acts[arc->get_val()]->weight;

		sketch->min_tard += local_penalty;

		if( sketch->min_tard > MIN(target_state->max_tard.getValue(), obj_var.getMax()) ) {
			return true;
		}
		mem_min_penalty[layer][arc->get_layer_id()] = local_penalty;
		mem_min_tardy[layer][arc->get_layer_id()] = sketch->min_tard;
		break;
	}

	case MaxTardiness:
	{
		if( layer > 0 ) {
			sketch->min_tard = INF;
		} else {
			sketch->min_tard = 0;
		}
		for( vector<int>::iterator it = valid_acts.begin();
				it != valid_acts.end();
				++it )
		{
			sketch->min_tard = MIN(source_state->in_tard->value(*it), sketch->min_tard);
		}
		int local_penalty = MAX(0, sketch->earliest_start - acts[arc->get_val()]->duedate)
				* acts[arc->get_val()]->weight;

		sketch->min_tard = MAX( sketch->min_tard, local_penalty );

		if( sketch->min_tard > MIN(target_state->max_tard.getValue(), obj_var.getMax()) ) {
			return true;
		}
		mem_min_penalty[layer][arc->get_layer_id()] = local_penalty;
		mem_min_tardy[layer][arc->get_layer_id()] = sketch->min_tard;
		break;
	}

	case Latency:
	{
		if( layer > 0 ) {
			sketch->min_latency = INF;
		} else {
			sketch->min_latency = 0;
		}
		for( vector<int>::iterator it = valid_acts.begin();
				it != valid_acts.end();
				++it )
		{
			sketch->min_latency = MIN(source_state->in_latency->value(*it), sketch->min_latency);
		}
		int local_penalty = sketch->earliest_start * acts[arc->get_val()]->weight;
		sketch->min_latency += local_penalty;

		if( sketch->min_latency > MIN(target_state->max_latency.getValue(), obj_var.getMax()) ) {
			return true;
		}
		mem_min_penalty[layer][arc->get_layer_id()] = local_penalty;
		mem_min_tardy[layer][arc->get_layer_id()] = sketch->min_latency;
		break;
	}

	default:
		break;
	}

	/**
	 * Precedence information
	 */
	if( has_input_prec ) {

		// --------------------------------
		// Definitely precedes filtering
		// --------------------------------
		if( source_state->all_paths.intersects(act_precedes[arc->get_val()]) ) {
			return true;
		}

		// --------------------------------
		// Definitely succeeds filtering
		// --------------------------------
		temp_set.set();
		temp_set -= source_state->some_paths;
		if( temp_set.intersects(act_succeeds[arc->get_val()]) ) {
			return true;
		}
		//exit(1);
	}

	// remaining sketch parameters
	act_est->add(arc->get_val(), est_arc - acts[arc->get_val()]->processing);

	// commented out due to partition_exact_acts
	//sketch->size_implied = source_state->all_paths.count()+1;
	sketch->size_implied = 0;

	return false;
}

/**
 * Check if an in arc is not necessarily infeasible according to cost of alldiff. It also
 * sets the corresponding arc sketch.
 		- Notice that this method only computes the cost and is never going to return that the arc is infeasible.
 **/
bool DisjunctivePropagator::is_incoming_infeasible_lagrangian_bound(int layer, RevMDD::Arc* arc){

	DisjunctiveState* source_state = get_node_state(layer, arc->get_source());
	DisjunctiveSketch* sketch = sketches[arc->get_layer_id()];
	
	// Compute the cost of the current arc
	if (layer > 0){
		sketch->min_lagrangian_bound_R = INF;
		int prev_value;
		
		for( int i = 0; i < source_state->in_est->size; ++i ) {
			prev_value = source_state->in_est->list[i];
		
			if( prev_value != arc->get_val() && !precedes_activity(arc->get_val(), prev_value) ){
                sketch->min_lagrangian_bound_R = MIN(sketch->min_lagrangian_bound_R,
                                                     source_state->min_lagrangian_bound_R->value(prev_value) + setup[prev_value][arc->get_val()]
                                                     + lagr_cost_task[arc->get_val()].getValue()  + lagr_cost_layer[layer].getValue()*weights[arc->get_val()]
                                                     + layer*lagr_cost_pre[arc->get_val()].getValue() );
			}
		}
	}
	else {
		sketch->min_lagrangian_bound_R = lagr_cost_task[arc->get_val()].getValue() + lagr_cost_layer[0].getValue() + lagr_cost_pre[0].getValue();
	}

	
	//cout << "Arc at (" << layer << "," << arc->get_source() << ") - with value=" << arc->get_val();
	//cout << " :: " << sketch->min_lagrangian_bound_R << endl;
	
	lagrangian_bound_topdown[layer][arc->get_layer_id()] = sketch->min_lagrangian_bound_R; // this is used in the bottom up filtering

	return false;
}


/**
 * Check if an in arc is not necessarily infeasible
 **/
bool DisjunctivePropagator::is_outgoing_infeasible(int layer, RevMDD::Arc* arc) {
	if( is_outgoing_infeasible_alldiff(layer, arc) ) {
		return true;
	}
	if( is_outgoing_infeasible_disjunctive(layer, arc) ) {
		return true;
	}
	#if CAP_FILT > 0
	if( is_outgoing_infeasible_capacity(layer, arc) ) {
		//cout << "remove arc due to capacity constraint -outgoing" << endl;
		return true;
	}
	#endif
	
	if( is_outgoing_infeasible_lagrangian_bound(layer, arc) ){
		return true;
	} 
	
	return false;
}

/**
 * Check if an in arc is not necessarily infeasible according to the capacity
 **/
inline bool DisjunctivePropagator::is_outgoing_infeasible_capacity(int layer, RevMDD::Arc* arc) {

	DisjunctiveState* target_state = get_node_state(layer+1, arc->get_target());
	DisjunctiveSketch* arc_sketch = sketches[arc->get_layer_id()];

	if (layer >= get_mdd()->get_num_layers() - 2) {
		arc_sketch->min_Q = 0;
		arc_sketch->max_Q = 0;
		return false;
	}
	
	//cout << "layer = " << layer << endl;

	int q_min = INF;
	int q_max = (-1)*INF;

	for( int i = 0; i < target_state->out_Q_min->size; ++i ) {
	    int prev_city = target_state->out_Q_min->list[i];
	    if( prev_city != arc->get_val()) {
	    	int q_min_trg = target_state->out_Q_min->value(prev_city);
	    	int q_max_trg = target_state->out_Q_max->value(prev_city);

	    	q_min = MIN(q_min, q_min_trg + weights[arc->get_val()]);
	    	q_max = MAX(q_max, q_max_trg + weights[arc->get_val()]);
	    	
	    }
	}

	//cout << "Arc at (" << layer << "," << arc->get_source() << ") - value=" << arc->get_val();
	//cout << " :: [" << q_min << "," << q_max << "]" << endl;
	
	q_min = MAX(-capacity, q_min);
	q_max = MIN(0, q_max);
	
	// check feasibility ----------------------------------------------------
	if (q_min > 0 )	{
		//cout << "Infeasible solution: Q_min greater than 0" << endl;
		return true;
	}
	if (q_max < -capacity){
		//cout << "Infeasible solution: Q_max smaller than  -capacity" << endl;
		return true;
	}
	
	
	
	if (q_min + q_min_topdown[layer][arc->get_layer_id()] - weights[arc->get_val()] > 0){
	 	//cout << "Infeasible solution: path Qmin > 0 " << endl;
	 	return true;
	}
	if (q_max + q_max_topdown[layer][arc->get_layer_id()] - weights[arc->get_val()] < 0){
	 	//cout << "Infeasible solution: path Qmax < 0 " << endl;
	 	return true;
	}
	// ----------------------------------------------------------------------

	arc_sketch->min_Q = q_min;
	arc_sketch->max_Q = q_max;

	return false;
}

/**
 * Check if an in arc is not necessarily infeasible according to alldiff
 **/
inline bool DisjunctivePropagator::is_outgoing_infeasible_alldiff(int layer, RevMDD::Arc* arc) {

	DisjunctiveState* target_state = get_node_state(layer+1, arc->get_target());

	// --------------------------------------------------
	// All diff conditions
	// --------------------------------------------------

	// all_path condition
	if( target_state->all_paths[arc->get_val()] ) {
		return true;
	}

	// some paths condition
	if( target_state->hall_set && target_state->some_paths[arc->get_val()] )	{
		return true;
	}

	return false;
}


/**
 * Check if an in arc is not necessarily infeasible accordding to disjunctive
 **/
bool DisjunctivePropagator::is_outgoing_infeasible_disjunctive(int layer, RevMDD::Arc* arc) {

	DisjunctiveState* target_state = get_node_state(layer+1, arc->get_target());

	// --------------------------------------------------
	// Disjunctive conditions
	// --------------------------------------------------


	/**
	 * 1. Compute arc ltf
	 */
    DisjunctiveState* source_state = get_node_state(layer, arc->get_source());

	// compute latest finish time of arc
	int lft_arc = -1, lft_arc_tmp;
	valid_acts.clear();

	if( layer < get_mdd()->get_num_layers()-2 ) {
	  int next_value;
	  for( int i = 0; i < target_state->out_lft->size; ++i ) {
	    next_value = target_state->out_lft->list[i];
	    if( next_value != arc->get_val() &&
	        !precedes_activity(next_value, arc->get_val()) )
	    {
	      lft_arc_tmp = target_state->out_lft->value(next_value)
                  - setup[arc->get_val()][next_value];
	      lft_arc = MAX(lft_arc, lft_arc_tmp);

	      valid_acts.push_back(next_value);
	    }
	  }
	} else {
	  lft_arc = acts[arc->get_val()]->start_time.getEndMax();
	}

	// if not activity can succeed the current one, arc is infeasible
	if( lft_arc == -1 ) {
	  return true;
	}

	// take into account deadline and processing times
	lft_arc = MIN(lft_arc, acts[arc->get_val()]->start_time.getEndMax())
    				- acts[arc->get_val()]->processing;

	lft_arc = MIN(lft_arc, source_state->latest_finish.getValue());


	DisjunctiveSketch* sketch = sketches[arc->get_layer_id()];
	sketch->latest_finish = lft_arc;


	/**
	 * 2. Check if arc has support
	 */

	// depending on the objective function, we can update the
	// latest finish time of the activity
	switch( obj_type ) {

	case Makespan:
		lft_arc = MIN(lft_arc, obj_var.getMax());
		break;

	case TotalTardiness:
		if( layer < get_mdd()->get_num_layers()-2 ) {
			sketch->max_tard = -1;
		} else {
			sketch->max_tard = obj_var.getMax();
		}

		// obtain maximum tardiness from outgoing nodes
		for( vector<int>::iterator it = valid_acts.begin();
				it != valid_acts.end();
				++it )
		{
			sketch->max_tard = MAX(sketch->max_tard, target_state->out_tard->value(*it));
		}
		if( mem_min_tardy[layer][arc->get_layer_id()] > sketch->max_tard ) {
			return true;
		}

		// search for minimum tardiness from incoming nodes
		if( sketch->max_tard < INF ) {

			// TODO: I think this was already computed before... you just need
			// to extract the minimum local penalty
			int min_tardy_R = INF;
			if( layer > 0 ) {
				int prev_value;
				for( int i = 0; i < source_state->in_tard->size; ++i ) {
					prev_value = source_state->in_tard->list[i];
					if( prev_value != arc->get_val() &&
							!precedes_activity(arc->get_val(), prev_value) )
					{
						min_tardy_R = MIN(min_tardy_R, source_state->in_tard->value(prev_value));
					}
				}
			} else {
				min_tardy_R = 0;
			}
			// compute maximum allowed lft
			if( acts[arc->get_val()]->weight > 0 ) {
				int max_lft = ((sketch->max_tard - min_tardy_R)/acts[arc->get_val()]->weight+1)
							+ acts[arc->get_val()]->duedate
							- acts[arc->get_val()]->processing;
				lft_arc = MIN(lft_arc, max_lft);
				if( lft_arc < 0 ) {
					cerr << "ERROR :: lft" << endl;
					exit(1);
				}
			}
			// adjust maximum tardiness to arc
			sketch->max_tard -= mem_min_penalty[layer][arc->get_layer_id()];
		}

		// compute sum of tardiness from the route
		sketch->min_tard_T = mem_min_penalty[layer][arc->get_layer_id()] + target_state->min_tard_T.getValue();
		break;


	case MaxTardiness:
		sketch->max_tard = obj_var.getMax();
		if( sketch->max_tard < INF - EPS_INF ) {
			if( acts[arc->get_val()]->weight > 0 ) {
				int max_lft = (sketch->max_tard/acts[arc->get_val()]->weight+1)	+ acts[arc->get_val()]->duedate;
				if( lft_arc < 0 ) {
					cerr << "ERROR :: lft: " << lft_arc << " - max_lft: " << max_lft;
					cerr << " - obj: " << sketch->max_tard << " - INF: " << INF << endl;
					exit(1);
				}
				lft_arc = MIN(lft_arc, max_lft);
			}
		}
		break;


	case Latency:
		if( layer < get_mdd()->get_num_layers()-2 ) {
			sketch->max_latency = -1;
		} else {
			sketch->max_latency = obj_var.getMax();
		}
		// obtain maximum tardiness from outgoing nodes
		for( vector<int>::iterator it = valid_acts.begin();
				it != valid_acts.end();
				++it )
		{
			sketch->max_latency = MAX(sketch->max_latency, target_state->out_latency->value(*it));
		}
		if( mem_min_tardy[layer][arc->get_layer_id()] > sketch->max_latency ) {
			return true;
		}

		// search for minimum tardiness from incoming nodes

		if( sketch->max_tard < INF ) {

			// TODO: I think this was already computed before... you just need
			// to extract the minimum local penalty
			int min_late_R = INF;
			if( layer > 0 ) {
				int prev_value;
				for( int i = 0; i < source_state->in_latency->size; ++i ) {
					prev_value = source_state->in_latency->list[i];
					if( prev_value != arc->get_val() &&
							!precedes_activity(arc->get_val(), prev_value) )
					{
						min_late_R = MIN(min_late_R, source_state->in_latency->value(prev_value));
					}
				}
			} else {
				min_late_R = 0;
			}
			// compute maximum allowed lft
			if( acts[arc->get_val()]->weight > 0 ) {
				int max_lft = ((sketch->max_latency - min_late_R)/acts[arc->get_val()]->weight+1)
						- acts[arc->get_val()]->processing;
				lft_arc = MIN(lft_arc, max_lft);
			}
			// adjust maximum tardiness to arc
			sketch->max_latency -= mem_min_penalty[layer][arc->get_layer_id()];
		}
		break;


	default:
		break;
	}

	// search for a support
	bool has_support = false;
	if( layer > 0 ) {
		int prev_value;
		for( int i = 0; i < source_state->in_est->size && !has_support; ++i ) {
			prev_value = source_state->in_est->list[i];
			if( prev_value != arc->get_val() &&
					!precedes_activity(arc->get_val(), prev_value) )
			{
				has_support =
						(source_state->in_est->value(prev_value)
								+ setup[prev_value][arc->get_val()]) <= lft_arc;
			}
		}
	} else {
		has_support = (lft_arc >= 0);
	}
	if( !has_support ) {
		return true;
	}


	/**
	 * 3. Objective function propagator
	 */
	RevMDD::Arc* next_arc = NULL;

	switch( obj_type ) {

	case SumSetupTimes:
		if( layer < get_mdd()->get_num_layers()-2 ) {
			sketch->min_sumsetup_T = INF;
		} else {
			sketch->min_sumsetup_T = 0;
		}

		
		for( vector<int>::iterator it = valid_acts.begin();
				it != valid_acts.end();
				++it )
		{
			int min_setup_temp = target_state->out_setup->value(*it) + setup[arc->get_val()][*it];

			if(sketch->min_sumsetup_T > min_setup_temp ){
				sketch->min_sumsetup_T = min_setup_temp;
				next_arc = (RevMDD::Arc*)target_state->out_setup->pointer(*it);
			}

			// 	sketch->min_sumsetup_T = MIN(sketch->min_sumsetup_T,
			// 			target_state->out_setup->value(*it)
			// 			+ setup[arc->get_val()][*it]);
		}
		arc->p_hook = (void*) next_arc;

		sketch->min_sumsetup_T = MAX(sketch->min_sumsetup_T, source_state->min_sumsetup_T.getValue());
		if( source_state->min_sumsetup_R.getValue() + sketch->min_sumsetup_T
				> obj_var.getMax() )
		{
	      return true;
	    }
		break;

	default:
		break;
	}

	/**
	 * Precedence information
	 */
	if( has_input_prec ) {

		// --------------------------------
		// Definitely succeeds filtering
		// --------------------------------
		if( target_state->all_paths.intersects(act_succeeds[arc->get_val()]) ) {
			return true;
		}

		// --------------------------------
		// Definitely preceeds filtering
		// --------------------------------
		temp_set.set();
		temp_set -= target_state->some_paths;
		if( temp_set.intersects(act_precedes[arc->get_val()]) ) {
			return true;
		}
	}



    act_lft->add(arc->get_val(), lft_arc+acts[arc->get_val()]->processing);
	return false;
}

bool DisjunctivePropagator::is_outgoing_infeasible_lagrangian_bound(int layer, RevMDD::Arc* arc){

	DisjunctiveState* target_state = get_node_state(layer+1, arc->get_target());
	DisjunctiveSketch* sketch = sketches[arc->get_layer_id()];
	RevMDD::Arc* next_arc = NULL;
	
	// Compute the cost of the current arc
	if (layer < get_mdd()->get_num_layers()-2) {
		sketch->min_lagrangian_bound_T = INF;
		int next_value;
		
		for( int i = 0; i < target_state->min_lagrangian_bound_T->size; ++i) {
			next_value = target_state->min_lagrangian_bound_T->list[i];
		
			if( next_value != arc->get_val() && !precedes_activity(next_value, arc->get_val()) ){
				double min_temp = target_state->min_lagrangian_bound_T->value(next_value) + setup[arc->get_val()][next_value];

				if( min_temp < sketch->min_lagrangian_bound_T  ){
					sketch->min_lagrangian_bound_T = min_temp;
					next_arc = (RevMDD::Arc*)target_state->min_lagrangian_bound_T->pointer(next_value);
				}

				//sketch->min_lagrangian_bound_T = MIN(sketch->min_lagrangian_bound_T, 
				//						 		 target_state->min_lagrangian_bound_T->value(next_value) + setup[arc->get_val()][next_value] );
			}
		}
		//for computing the shortest path
		arc->p_hook = (void*) next_arc;
	}
	else {
		sketch->min_lagrangian_bound_T = 0.0;
	}
	
	//Infeasibility rule
	if (sketch->min_lagrangian_bound_T > 0.0) {
		if ( (int) (sketch->min_lagrangian_bound_T + lagrangian_bound_topdown[layer][arc->get_layer_id()]) > obj_var.getMax() ){
		
		/*
		cout << "Arc was found infeasible due to lagrangian bound\n";
		cout << "Top-down = " << lagrangian_bound_topdown[layer][arc->get_layer_id()];
		cout << "\tBottom-Up = " << sketch->min_lagrangian_bound_T;
		cout << "\tCurrent Bound = " << obj_var.getMax() << endl;
		*/
		return true;
		}
	}
	
	if(layer > 0){
        sketch->min_lagrangian_bound_T += lagr_cost_task[arc->get_val()].getValue() + lagr_cost_layer[layer].getValue()*weights[arc->get_val()] + layer*lagr_cost_pre[arc->get_val()].getValue();  // add the lagrangian multiplier cost
	}else{
		sketch->min_lagrangian_bound_T += lagr_cost_task[arc->get_val()].getValue()+ lagr_cost_layer[layer].getValue() + lagr_cost_pre[arc->get_val()].getValue();  // add the lagrangian multiplier cost
	}
	
	// if (layer < get_mdd()->get_num_layers()-2) {
	// 	cout << "Arc at (" << layer << "," << arc->get_source() << ") - with value = " << arc->get_val();
	// 	cout << " ::  lagr bound = " << sketch->min_lagrangian_bound_T << endl;
	// 	cout << "\t:: hook to arc at (" << ((RevMDD::Arc*)arc->p_hook)->get_layer_id() << "," << ((RevMDD::Arc*)arc->p_hook)->get_source() << ") - with value = " << ((RevMDD::Arc*)arc->p_hook)->get_val() << endl;
	// 	cout << "\t:: hook to arc at (" <<  next_arc->get_layer_id() << "," << next_arc->get_source() << ") - with value = " << next_arc->get_val() << endl;
	// }

	return false;
}


// Compute makespan and update corresponding variable
// It is a simple breadth-first search that queues
// arcs instead of nodes

// TODO(priority!) If optimum solution is feasible, fix makespan
// and setup !!

// TODO !!!! should replace to recompute all objective
// functions, as there might be dangling nodes... !!!

// TODO do it more efficiently! perhaps it is
// really necessary to organize arcs by source node

// TODO after computing this, if we update the est
// of the nodes we do not need to do them again in
// the next top-down... perhaps it would be better
// to always start from bottom up
void DisjunctivePropagator::update_start_end_times() {

	// update activity start/end times
	for( int act = 0; act < (int)acts.size(); ++act ) {
		if( !act_est->contains(act) || !act_lft->contains(act) ) {
			get_mdd()->export_to_gml("graphs/error.gml");
			cout << "ERROR: MDD does not contain activity " << act << endl;
			exit(1);
		}
		if( act_est->value(act) > acts[act]->start_time.getStartMin() ) {
			//cout << "\tUpdate start times: " << acts[act]->start_time;
			acts[act]->start_time.setStartMin(act_est->value(act));
			//cout << " - to: " << acts[act]->start_time << endl;
		}
		if( act_lft->value(act) < acts[act]->start_time.getEndMax() ) {
			//cout << "\tUpdate end times: " << acts[act]->start_time;
			acts[act]->start_time.setEndMax(act_lft->value(act));
			//cout << " - to: " << acts[act]->start_time << endl;
		}
	}

    // MAYBE CAN LEAD TO IMPROVEMENTS!!!
	//for (int i = 0; i < (int)acts.size(); ++i) {
	  //if (permutation[i].isFixed()) {
	    //cout << "setting start time of " << acts[permutation[i].getValue()]->start_time << " to: " << acts[permutation[i].getValue()]->start_time.getStartMin() << endl;
	    //acts[permutation[i].getValue()]->start_time.setStart(
//								 acts[permutation[i].getValue()]->start_time.getStartMin());
	  //} else {
	   // break;
	  //}
//	}

  return;


  // verify shortest path
  RevMDD* mdd = get_mdd();

  
  // clean nodes
  for( int l = 0; l < mdd->get_num_layers()-1; ++l ) {
    for( int w = 0; w < mdd->get_max_width(); ++w ) {
      states[l+1][w]->in_est->clear();
      states[l+1][w]->est_tmp = INF;
    }
  }
  root_state()->est_tmp = 0;
  root_state()->in_est->clear();

  // compute shortest path
  DisjunctiveState *source_state, *target_state;
  RevMDD::Arc* prev_arc = NULL;
  int est_arc, est_arc_tmp, est_arc_relaxed, prev_value;
  for( int l = 0; l < mdd->get_num_layers()-1; ++l ) {
    for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
        arc != NULL;
        arc = mdd->arcs_layer[l].get_next(arc) )
    {
      source_state = states[l][arc->get_source()];
      target_state = states[l+1][arc->get_target()];

      est_arc = INF;
      if( l == 0 ) {
        est_arc = 0;
        prev_arc = NULL;
      } else {
        int prev_value;
        for( int i = 0; i < source_state->in_est->size; ++i ) {
          prev_value = source_state->in_est->list[i];
          if( prev_value != arc->get_val()
              &&
              !precedes_activity(arc->get_val(), prev_value) )
          {
            est_arc_tmp = source_state->in_est->value(prev_value)
                            + setup[prev_value][arc->get_val()];
            if( est_arc_tmp < est_arc ) {
              est_arc = est_arc_tmp;
              prev_arc = (RevMDD::Arc*)source_state->in_est->pointer(prev_value);
            }
          }
        }
        arc->p_hook = (void*)prev_arc;
      }
      if( est_arc != INF ) {
        target_state->in_est->add(arc->get_val(), est_arc, arc);
      }
    }
  }
  terminal_state()->est_tmp = INF;
  void* last = NULL;
  for( int i = 0; i < terminal_state()->in_est->size; ++i ) {
    if( terminal_state()->est_tmp > terminal_state()->in_est->value(terminal_state()->in_est->list[i]) ) {
      terminal_state()->est_tmp = terminal_state()->in_est->value(terminal_state()->in_est->list[i]);
      last = terminal_state()->in_est->pointer(terminal_state()->in_est->list[i]);
    }
  }
  // compose optimum solution
  int* optimum = terminal_state()->optimum_sol;
  int l = permutation.getSize()-1;
  while( last != NULL ) {
    terminal_state()->optimum_sol[l--] = ((RevMDD::Arc*)last)->get_val();
    last = ((RevMDD::Arc*)last)->p_hook;
  }



  // check if the solution is feasible
  memset(act_checked, false, sizeof(bool)*acts.size());

  int completion = acts[optimum[0]]->release + acts[optimum[0]]->processing; // assumed to be feasible
  act_checked[optimum[0]] = true;

  bool feasible = true;
  int sum_setups = 0;
  for( int i = 1; i < (int)acts.size(); ++i ) {
    if( act_checked[optimum[i]] ) {
      feasible = false;
      //cout << "\trepeated activity: " << optimum[i] << endl;
      break;
    }

    act_checked[optimum[i]] = true;
    completion = MAX(completion + setup[optimum[i-1]][optimum[i]], acts[optimum[i]]->release)
        + acts[optimum[i]]->processing;

    sum_setups += setup[optimum[i-1]][optimum[i]];

//    cout << "\t\t" << optimum[i-1] << " --> " << optimum[i] << ": " << setup[optimum[i-1]][optimum[i]] << endl;
//    cout << "\t\t\trelease: " << acts[optimum[i]]->release << endl;
//    cout << "\t\t\tcompletion: " << completion << endl;
//    cout << "\t\t\trelease: " << acts[optimum[i]]->release << endl;
//    cout << "\t\t\tprocessing: " << acts[optimum[i]]->release << endl;


    if( completion > acts[optimum[i]]->deadline ) {
      feasible = false;
      //cout << "\tcompletion time: " << optimum[i] << " - " << completion << " - " << acts[optimum[i]]->deadline << endl;
      break;
    }
  }

//  cout << "Permutation: " << endl;
//  for( int i = 0; i < permutation.getSize(); ++i ) {
//    cout << "\t" << permutation[i] << endl;
//  }
//  cout << "optimum sol: ";
//  for( int i = 0; i < permutation.getSize(); ++i ) {
//    cout << terminal_state()->optimum_sol[i] << " ";
//  }
//  cout << endl;
//  cout << "\tvalue: " << terminal_state()->est_tmp << " - sum setup: " << sum_setups << endl;

  if( !feasible ) {
    //cout << "\tinfeasible" << endl;
  } else {
    for( int i = 0; i < permutation.getSize(); ++i ) {
      permutation[i].setValue(optimum[i]);
    }
    //cout << "\t*feasible - sum setups: " << sum_setups << endl;
  }
  //cout << endl;

  //exit(1);

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
  return;

  //mdd->export_to_gml("graphs/test.gml");


  // reset arc states
  root_state()->in_est->clear();
  for( int l = 0; l < mdd->get_num_layers()-1; ++l ) {
    for( int w = 0; w < mdd->get_max_width(); ++w ) {
      states[l+1][w]->in_est->clear();
      states[l+1][w]->est_tmp = INF;
    }
  }
  act_est->clear();
  act_lft->clear();

  // arcs in layer 0: minimum start time correspond to their release dates
  for( RevMDD::Arc* arc = mdd->arcs_layer[0].get_first();
      arc != NULL;
      arc = mdd->arcs_layer[0].get_next(arc) )
  {
    states[1][arc->get_target()]->in_est->add(arc->get_val(),
                  acts[arc->get_val()]->start_time.getStartMin()+acts[arc->get_val()]->processing);
    act_est->add(arc->get_val(), acts[arc->get_val()]->start_time.getStartMin());
    //act_lft->add(arc->get_val(), acts[arc->get_val()]->start_time.getStartMin()+acts[arc->get_val()]->processing);
  }

  // remaining arcs: compute latest start time
  for( int l = 1; l < mdd->get_num_layers()-1; ++l ) {
    // organize arcs according to groups
    for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
        arc != NULL;
        arc = mdd->arcs_layer[l].get_next(arc) )
    {
      source_state = states[l][arc->get_source()];
      target_state = states[l+1][arc->get_target()];

      assert( source_state->in_est->size > 0 );
      est_arc = INF; est_arc_relaxed = INF;
      for( int i = 0; i < source_state->in_est->size; ++i ) {
        prev_value = source_state->in_est->list[i];

        est_arc_tmp = source_state->in_est->value(prev_value)
                + setup[prev_value][arc->get_val()];
        est_arc_relaxed = MIN(est_arc_relaxed, est_arc_tmp);

        if( prev_value != arc->get_val()
            &&
            !precedes_activity(arc->get_val(), prev_value) )
        {
          est_arc = MIN(est_arc, est_arc_tmp);
        }
      }

      // try to be as exact as possible... but we are only interested
      // in a relaxation here
      if( est_arc == INF ) {
        est_arc = est_arc_relaxed;
      }

      // take into account release and processing times
      est_arc = MAX(est_arc, acts[arc->get_val()]->start_time.getStartMin());
      act_est->add(arc->get_val(), est_arc);

      est_arc += acts[arc->get_val()]->processing;
      //act_lft->add(arc->get_val(), est_arc);

      target_state->in_est->add(arc->get_val(), est_arc);
      target_state->est_tmp = MIN(target_state->est_tmp, est_arc);
    }
  }

  //cout << "makespan: " << terminal_state()->est_tmp << endl;

//  // update makespan
//  if (makespan.getMin() < terminal_state()->est_tmp) {
//    get_CP().add( makespan >= terminal_state()->est_tmp );
//  }

//  // update activity start/end times
//  for( int act = 0; act < (int)acts.size(); ++act ) {
//    if( !act_est->contains(act) || !act_lft->contains(act) ) {
//      cout << "ERROR: MDD does not contain activity " << act << endl;
//      exit(1);
//    }
//    cout << "activity: " << act << endl;
//    if( act_est->value(act) > acts[act]->start_time.getStartMin() ) {
//      cout << "EST: before " << acts[act]->start_time.getStartMin();
//      cout << " - after: " << act_est->value(act) << endl;
//      acts[act]->start_time.setStartMin(act_est->value(act));
////      cout << "\tEST now:" << acts[act]->start_time.getStartMin() << endl;
////      cout << endl;
//    }
//    if( act_lft->value(act) < acts[act]->start_time.getEndMax() ) {
//      cout << "LFT: before " << acts[act]->start_time.getEndMax();
//      cout << " - after: " << act_lft->value(act) << endl;
//      acts[act]->start_time.setEndMax(act_lft->value(act));
//    }
//  }

  // makespan relaxation value is now in the terminal node
//  std::cout << "makespan: " << terminal_state()->est_tmp << std::endl;
//  cout << "optimum sol: ";
//  for( int i = 0; i < (int)acts.size(); ++i ) {
//    cout << terminal_state()->optimum_sol[i] << " ";
//  }
//  cout << endl;
//  exit(1);
}


/**
 * Compute shortest path in M, which corresponds to the best
 * solution of the MDD
 */
int* DisjunctivePropagator::compute_shortest_path(int& value) {

  // verify shortest path
  RevMDD* mdd = get_mdd();

  // clean nodes
  for( int l = 0; l < mdd->get_num_layers()-1; ++l ) {
    for( int w = 0; w < mdd->get_max_width(); ++w ) {
      states[l+1][w]->in_est->clear();
      states[l+1][w]->est_tmp = INF;
    }
  }
  root_state()->est_tmp = 0;
  root_state()->in_est->clear();


  // compute shortest path

  DisjunctiveState *source_state, *target_state;
  RevMDD::Arc* prev_arc = NULL;

  // first layer
  for( RevMDD::Arc* arc = mdd->arcs_layer[0].get_first();
      arc != NULL;
      arc = mdd->arcs_layer[0].get_next(arc) )
  {
    states[1][arc->get_target()]->in_est->add(arc->get_val(), 0, arc);
    arc->p_hook = NULL;
  }

  // remaining layers
  int est_arc, est_arc_tmp;
  for( int l = 1; l < mdd->get_num_layers()-1; ++l ) {
    for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
        arc != NULL;
        arc = mdd->arcs_layer[l].get_next(arc) )
    {
      source_state = states[l][arc->get_source()];
      target_state = states[l+1][arc->get_target()];

      est_arc = INF;
      int prev_value;
      for( int i = 0; i < source_state->in_est->size; ++i ) {
        prev_value = source_state->in_est->list[i];
        if( prev_value != arc->get_val()
            &&
            !precedes_activity(arc->get_val(), prev_value) )
        {
          est_arc_tmp = source_state->in_est->value(prev_value)
                                + setup[prev_value][arc->get_val()];
          if( est_arc_tmp < est_arc ) {
            est_arc = est_arc_tmp;
            prev_arc = (RevMDD::Arc*)source_state->in_est->pointer(prev_value);
          }
        }
      }
      arc->p_hook = (void*)prev_arc;

      if( est_arc != INF ) {
        target_state->in_est->add(arc->get_val(), est_arc, arc);
      }

      //cout << "Arc in layer " << l << " with value " << arc->get_val() << " has distance: " << est_arc << endl;
    }
  }

  // compute shortest path value and store it in terminal node
  terminal_state()->est_tmp = INF;
  void* last = NULL;
  for( int i = 0; i < terminal_state()->in_est->size; ++i ) {
  	//cout << "Arc going to the terminal node: " << terminal_state()->in_est->value(terminal_state()->in_est->list[i]) << endl;
    if( terminal_state()->est_tmp > terminal_state()->in_est->value(terminal_state()->in_est->list[i]) ) {
      terminal_state()->est_tmp = terminal_state()->in_est->value(terminal_state()->in_est->list[i]);
      last = terminal_state()->in_est->pointer(terminal_state()->in_est->list[i]);
    }
  }

  if( terminal_state()->est_tmp == INF ) {
    cout << "ERROR: MDD is not feasible" << endl;
    exit(1);
  }

  // compose optimum solution
  int* optimum = terminal_state()->optimum_sol;
  int l = permutation.getSize()-1;
  while( last != NULL ) {
    optimum[l--] = ((RevMDD::Arc*)last)->get_val();
    last = ((RevMDD::Arc*)last)->p_hook;
  }
  value = terminal_state()->est_tmp;
  //cout << "Optimum value: " << terminal_state()->est_tmp << "\n";

 // if(check_feasible_solution(optimum)){
 //	fix_solution(optimum);
 //}

  return optimum;
}


int* DisjunctivePropagator::compute_shortest_path(double& value, double* cost_task, double* cost_layer, double* cost_pre, int nMDD) {

  // verify shortest path
  RevMDD* mdd = get_mdd();


  // clean nodes
  for( int l = 0; l < mdd->get_num_layers()-1; ++l ) {
    for( int w = 0; w < mdd->get_max_width(); ++w ) {
      states[l+1][w]->in_est_double->clear();
      states[l+1][w]->est_tmp_double = INF;
    }
  }

  root_state()->est_tmp_double = 0.0;
  root_state()->in_est_double->clear();
  
  // compute shortest path

  DisjunctiveState *source_state, *target_state;
  RevMDD::Arc* prev_arc = NULL;

  // first layer
  for( RevMDD::Arc* arc = mdd->arcs_layer[0].get_first();
      arc != NULL;
      arc = mdd->arcs_layer[0].get_next(arc) )
  {
    states[1][arc->get_target()]->in_est_double->add(arc->get_val(), cost_task[arc->get_val()] + cost_layer[0] + cost_pre[0], arc);
    arc->p_hook = NULL;
  }

  // remaining layers
  double est_arc, est_arc_tmp;
  for( int l = 1; l < mdd->get_num_layers()-1; ++l ) {
    for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
        arc != NULL;
        arc = mdd->arcs_layer[l].get_next(arc) )
    {
      source_state = states[l][arc->get_source()];
      target_state = states[l+1][arc->get_target()];

      est_arc = INF;
      int prev_value;
      for( int i = 0; i < source_state->in_est_double->size; ++i ) {
        prev_value = source_state->in_est_double->list[i];
        if( prev_value != arc->get_val()
            &&
            !precedes_activity(arc->get_val(), prev_value) )
        {
          est_arc_tmp = source_state->in_est_double->value(prev_value)
                                + setup[prev_value][arc->get_val()];
          if( est_arc_tmp < est_arc ) {
            est_arc = est_arc_tmp;
            prev_arc = (RevMDD::Arc*)source_state->in_est_double->pointer(prev_value);
          }
        }
      }
      arc->p_hook = (void*)prev_arc;
      if(nMDD == 1){
      	est_arc += weights[arc->get_val()]*cost_layer[l] + cost_task[arc->get_val()] + l*cost_pre[arc->get_val()];
      }
      else{ // case where we are using two MDDs
      	est_arc += arc->get_val()*cost_layer[l] + cost_task[arc->get_val()];
      }

      if( est_arc != INF ) {
        target_state->in_est_double->add(arc->get_val(), est_arc, arc);
      }
      //cout << "Arc in layer " << l << " with value " << arc->get_val() << " has distance: " << est_arc << endl;
    }
  }

  // compute shortest path value and store it in terminal node
  terminal_state()->est_tmp_double = INF;
  void* last = NULL;
  for( int i = 0; i < terminal_state()->in_est_double->size; ++i ) {
  	//cout << "Arc going to the terminal node: " << terminal_state()->in_est->value(terminal_state()->in_est->list[i]) << endl;
    if( terminal_state()->est_tmp_double > terminal_state()->in_est_double->value(terminal_state()->in_est_double->list[i]) ) {
      terminal_state()->est_tmp_double = terminal_state()->in_est_double->value(terminal_state()->in_est_double->list[i]);
      last = terminal_state()->in_est_double->pointer(terminal_state()->in_est_double->list[i]);
    }
  }

  if( terminal_state()->est_tmp_double == INF ) {
    //cout << "ERROR: MDD is not feasible" << endl;
    //exit(1);
  }

  // compose optimum solution
  int* optimum = terminal_state()->optimum_sol;
  int l = permutation.getSize()-1;
  while( last != NULL ) {
    optimum[l--] = ((RevMDD::Arc*)last)->get_val();
    last = ((RevMDD::Arc*)last)->p_hook;
  }

  value += terminal_state()->est_tmp_double;
  //cout << "Optimum value: " << terminal_state()->est_tmp_double << "\n";

 // if(check_feasible_solution(optimum)){
 //	fix_solution(optimum);
 //}

  return optimum;
}

/* Compute shortest path inside the MDD */
void  DisjunctivePropagator::save_shortest_path(){
	
	// compute shortest path value and store it in terminal node
  	root_state()->est_tmp = INF;
  	void* first_arc = NULL;
  	for( int i = 0; i < root_state()->out_setup->size; ++i ) {
  	    	if( root_state()->est_tmp > root_state()->out_setup->value(root_state()->out_setup->list[i]) ) {
      		root_state()->est_tmp = root_state()->out_setup->value(root_state()->out_setup->list[i]);
      		first_arc = root_state()->out_setup->pointer(root_state()->out_setup->list[i]);
    	}
  	}

	val_sp_normal = root_state()->est_tmp;
	//cout << "Shortest Path extraction.  \t value = " << val_sp_normal << endl;

	// compose optimum shortest path
  	sp_normal = root_state()->optimum_sol;
	for (int i = 0; i < permutation.getSize()-1; i ++) {
    	sp_normal[i] = ((RevMDD::Arc*)first_arc)->get_val();
    	first_arc = ((RevMDD::Arc*)first_arc)->p_hook;

		if (first_arc == NULL) cout << "pointer is null" << endl;
  	}
	sp_normal[permutation.getSize()-1] = permutation.getSize()-1;
	  
	//Printing
	// cout << "Path = : ";
	// for (int i = 0; i < permutation.getSize(); i ++) {
    // 	cout << sp_normal[i] << " ";
  	// }
	// cout << endl;

}

void DisjunctivePropagator::save_shortest_path_lagr(){

	// compute shortest path value and store it in terminal node
  	root_state()->est_tmp_double = INF;
  	void* first_arc = NULL;
  	for( int i = 0; i < root_state()->min_lagrangian_bound_T->size; ++i ) {
  	//cout << "Arc going to the terminal node: " << terminal_state()->in_est->value(terminal_state()->in_est->list[i]) << endl;
    	if( root_state()->est_tmp_double > root_state()->min_lagrangian_bound_T->value(root_state()->min_lagrangian_bound_T->list[i]) ) {
      		root_state()->est_tmp_double = root_state()->min_lagrangian_bound_T->value(root_state()->min_lagrangian_bound_T->list[i]);
      		first_arc = root_state()->min_lagrangian_bound_T->pointer(root_state()->min_lagrangian_bound_T->list[i]);
    	}
  	}

	val_sp_lagr = root_state()->est_tmp_double;

	//cout << "Lagragian Shortest Path extraction.  \t value = " << val_sp_lagr << endl;

	// compose optimum shortest path
  	sp_lagr = root_state()->optimum_sol;
	for (int i = 0; i < permutation.getSize()-1; i ++) {
    	sp_lagr[i] = ((RevMDD::Arc*)first_arc)->get_val();
    	first_arc = ((RevMDD::Arc*)first_arc)->p_hook;

		if (first_arc == NULL) cout << "pointer is null" << endl;
  	}
	sp_lagr[permutation.getSize()-1] = permutation.getSize()-1;
	  
	//Printing
	// cout << "Path = : ";
	// for (int i = 0; i < permutation.getSize(); i ++) {
    // 	cout << sp_lagr[i] << " ";
  	// }
	// cout << endl;
}

int* DisjunctivePropagator::get_shortest_path(int & value){
	value = val_sp_normal;
	return sp_normal;
}

int* DisjunctivePropagator::get_shortest_path_lagr(double & value){
	value = val_sp_lagr;
	return sp_lagr;
}

int* DisjunctivePropagator::compute_shortest_path() {

  // verify shortest path
  RevMDD* mdd = get_mdd();

  // clean nodes
  for( int l = 0; l < mdd->get_num_layers()-1; ++l ) {
    for( int w = 0; w < mdd->get_max_width(); ++w ) {
      states[l+1][w]->in_est->clear();
      states[l+1][w]->est_tmp = INF;
    }
  }
  root_state()->est_tmp = 0;
  root_state()->in_est->clear();


  // compute shortest path

  DisjunctiveState *source_state, *target_state;
  RevMDD::Arc* prev_arc = NULL;

  // first layer
  for( RevMDD::Arc* arc = mdd->arcs_layer[0].get_first();
      arc != NULL;
      arc = mdd->arcs_layer[0].get_next(arc) )
  {
    states[1][arc->get_target()]->in_est->add(arc->get_val(), 0, arc);
    arc->p_hook = NULL;
  }

  // remaining layers
  int est_arc, est_arc_tmp;
  for( int l = 1; l < mdd->get_num_layers()-1; ++l ) {
    for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
        arc != NULL;
        arc = mdd->arcs_layer[l].get_next(arc) )
    {
      source_state = states[l][arc->get_source()];
      target_state = states[l+1][arc->get_target()];

      est_arc = INF;
      int prev_value;
      for( int i = 0; i < source_state->in_est->size; ++i ) {
        prev_value = source_state->in_est->list[i];
        if( prev_value != arc->get_val()
            &&
            !precedes_activity(arc->get_val(), prev_value) )
        {
          est_arc_tmp = source_state->in_est->value(prev_value)
                                + setup[prev_value][arc->get_val()];
          if( est_arc_tmp < est_arc ) {
            est_arc = est_arc_tmp;
            prev_arc = (RevMDD::Arc*)source_state->in_est->pointer(prev_value);
          }
        }
      }
      arc->p_hook = (void*)prev_arc;

      if( est_arc != INF ) {
        target_state->in_est->add(arc->get_val(), est_arc, arc);
      }

      //cout << "Arc in layer " << l << " with value " << arc->get_val() << " has distance: " << est_arc << endl;
    }
  }

  // compute shortest path value and store it in terminal node
  terminal_state()->est_tmp = INF;
  void* last = NULL;
  for( int i = 0; i < terminal_state()->in_est->size; ++i ) {
  	//cout << "Arc going to the terminal node: " << terminal_state()->in_est->value(terminal_state()->in_est->list[i]) << endl;
    if( terminal_state()->est_tmp > terminal_state()->in_est->value(terminal_state()->in_est->list[i]) ) {
      terminal_state()->est_tmp = terminal_state()->in_est->value(terminal_state()->in_est->list[i]);
      last = terminal_state()->in_est->pointer(terminal_state()->in_est->list[i]);
    }
  }

  if( terminal_state()->est_tmp == INF ) {
    cout << "ERROR: MDD is not feasible" << endl;
    exit(1);
  }

  // compose optimum solution
  int* optimum = terminal_state()->optimum_sol;
  int l = permutation.getSize()-1;
  while( last != NULL ) {
    optimum[l--] = ((RevMDD::Arc*)last)->get_val();
    last = ((RevMDD::Arc*)last)->p_hook;
  }
  //cout << "Optimum value: " << terminal_state()->est_tmp << "\n";

 // if(check_feasible_solution(optimum)){
 //	fix_solution(optimum);
 //}

  return optimum;
}


/**
 * Simple edge for maximum spanning tree computation
 */
struct EdgeSPT {
  int v;
  int u;
  int weight;

//  bool operator<(Edge &a) {
//    return a.weight > weight;
//  }
};

struct EdgeSPTComp {
  bool operator()(EdgeSPT a, EdgeSPT b) {
    return a.weight > b.weight;
  }
};


///**
// * Compute activities that have priority to be exact in any solution
// */
//void DisjunctivePropagator::compute_exact_acts() {
//
//  // collect activities that are not fixed
//  int num_acts = acts.size();
//
//  // create activity edges for maximum spanning tree
//  vector<EdgeSPT> edges;
//  for( int i = 0; i < num_acts-1; ++i ) {
//    for( int j = i+1; j < num_acts; ++j ) {
//      EdgeSPT edge;
//      edge.u = i;
//      edge.v = j;
//      edge.weight = MIN(setup[i][j]+acts[j]->processing, setup[j][i]+acts[i]->processing);
//      edges.push_back(edge);
//    }
//  }
//  sort(edges.begin(), edges.end(), EdgeSPTComp());
//
//  memset(act_checked, false, sizeof(bool)*num_acts);
//
//  for( int i = 0; i < (int)edges.size(); ++i ) {
//    if( !act_checked[edges[i].u] ) {
//      exact_acts.push_back(edges[i].u);
//      act_checked[edges[i].u] = true;
//    }
//    if( !act_checked[edges[i].v] ) {
//      exact_acts.push_back(edges[i].v);
//      act_checked[edges[i].v] = true;
//    }
//    if( (int)exact_acts.size() == num_acts ) {
//      break;
//    }
//  }
//
//  cout << "Exact activity priority: ";
//  for( int i = 0; i < (int)exact_acts.size(); ++i ) {
//    cout << exact_acts[i] << " ";
//  }
//  cout << endl;
//
//  //exit(1);
//}



/**
 * Compute activities that have priority to be exact in any solution
 */
void DisjunctivePropagator::compute_exact_acts() {

	if( input_refinement_file ) {

		cout << "\n[Disjunctive] Exact activities priority (read from file): ";
		for( int i = 0; i < (int)exact_acts.size(); ++i ) {
			cout << exact_acts[i] << " ";
		}
		cout << endl << endl;


		return;
	}

	int num_acts = acts.size();

	if( random_exact_refinement ) {

		exact_acts.clear();
		for( int i = 0; i < num_acts; ++i ) {
			exact_acts.push_back(-1);
		}

		boost::random::uniform_int_distribution<> dist(0, num_acts-1);
		for( int act = 0; act < num_acts; ++act ) {
			int pos = dist( random_generator );
			while( exact_acts[pos] != -1 ) {
				pos = dist( random_generator );
			}
			exact_acts[pos] = act;
		}

		ofstream order_file("ref-order.txt");
		for( int i = 0; i < (int)exact_acts.size(); ++i ) {
			order_file << exact_acts[i] << " ";
		}
		order_file << "\n";
		order_file.close();


		cout << "\n[Disjunctive] Exact activities priority (random): ";
		for( int i = 0; i < (int)exact_acts.size(); ++i ) {
			cout << exact_acts[i] << " ";
		}
		cout << endl << endl;

		return;
	}


	switch( obj_type ) {

	case Makespan:
	{
		// collect activities that are not fixed
		memset(act_checked, false, sizeof(bool)*num_acts);
		for( int i = 0; i < permutation.getSize(); ++i ) {
			if( permutation[i].isFixed() ) {
				act_checked[permutation[i].getValue()] = true;
			}
		}

		exact_acts.clear();
		for( int i = 0; i < num_acts; ++i ) {
			if( !act_checked[i] ) {
				exact_acts.push_back(i);
			}
		}
		if( exact_acts.size() <= 2 )
			return;

		// first activity has the maximum deadline, using the largest processing time
		// as tie-breaking
		int max_release = -1;
		int max_proc = -1;
		int first = INF;
		for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
			int i = *it;
			if( (max_release < acts[i]->start_time.getStartMin() + acts[i]->processing )
					|| (max_release == acts[i]->start_time.getStartMin() + acts[i]->processing && max_proc < acts[i]->processing) )
			{
				first = it - exact_acts.begin();
				max_release = acts[i]->start_time.getStartMin() + acts[i]->processing;
				max_proc = acts[i]->processing;
			}
		}

		// put activity in the first position
		int tmp = exact_acts[0];
		exact_acts[0] = exact_acts[first];
		exact_acts[first] = tmp;

		// sort remaining activities by largest (processing time+min_setup_from)
		aux_v.resize(num_acts);
		for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
			int i = *it;
			aux_v[i] = (acts[i]->start_time.getStartMin() + acts[i]->processing + max_setup_from[i]);
			//aux_v[i] = (acts[i]->start_time.getStartMin() + acts[i]->processing + min_setup_from[i]);
			//aux_v[i] = max_setup_from[i] + acts[i]->processing;
			if( acts[i]->weight > 0 ) {
				aux_v[i] *= acts[i]->weight;
			}
		}
		sort(exact_acts.begin()+1, exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
		//sort(exact_acts.begin(), exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
		//sort(exact_acts.rbegin(), exact_acts.rend());

		cout << "\n[Disjunctive] Exact activities priority: ";
		for( int i = 0; i < (int)exact_acts.size(); ++i ) {
			cout << exact_acts[i] << " ";
		}
		cout << endl << endl;
	}
	break;

	case SumSetupTimes:
	{
		exact_acts.clear();
		vector<bool> is_ordered( num_acts, false );

		// first activity has the maximum release date, using the largest processing time
		// as tie-breaking
		int max_release = -1;
		int max_proc = -1;
		int first = INF;
		for ( int i = 0; i < num_acts; ++i ) {
			if( (max_release < acts[i]->start_time.getStartMin() + acts[i]->processing )
					|| (max_release == acts[i]->start_time.getStartMin() + acts[i]->processing && max_proc < acts[i]->processing) )
			{
				first = i;
				max_release = acts[i]->start_time.getStartMin() + acts[i]->processing;
				max_proc = acts[i]->processing;
			}
		}

		// put activity in the first position
		exact_acts.push_back( first );
		is_ordered[ first ] = true;
			
		// add remaining activities according to sum of setup times
		while( (int)exact_acts.size() < num_acts ) {
			int max_sum = -1;
			int act = -1;
			for( int i = 0; i < num_acts; ++i ) {
				if( is_ordered[i] )
					continue;
				int sum = 0;
				for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
					if( input_prec[*it][i] ) {
						sum += setup[*it][i];
					} else {
						sum += setup[i][*it];
					}
				}
				if( sum > max_sum ) {
					max_sum = sum;
					act = i;
				}
			}
			cout << "act: " << act << " - max_sum: " << max_sum << endl;
			exact_acts.push_back( act );
			is_ordered[act] = true;
		}
		
		#if PRIORITY > 0
		
		vector<int> auxiliar(exact_acts);  //store the order given by the setup times
		exact_acts.clear();
		exact_acts.push_back(0);
		
		for(int i = 0; i < num_acts; ++i) {
			is_ordered[i] = false;
		}
		is_ordered[0] = true;
		
		while( (int)exact_acts.size() < num_acts ) {
			int max_cap = -100;
			int min_cap = 100;
			int index = 0;
			
			for( int i = 0; i < num_acts; ++i ) {
				if (is_ordered[i]) 
					continue;
				if (abs(weights[i]) < min_cap){
				//if (abs(weights[i]) > max_cap){
					max_cap = abs(weights[i]);
					min_cap = abs(weights[i]);
					index = i;
				//} else if (abs(weights[i]) == max_cap){ //break ties with the sum of setups
				} else if (abs(weights[i]) == min_cap){ //break ties with the sum of setups
					for(int j = 0; j < num_acts; j++){
						if (auxiliar[j] == i){  //the new one has more priority
							max_cap = abs(weights[i]);
							min_cap = abs(weights[i]);
							index = i;
							break;
						}
						else if (j == index) break;
					}
				}
			}
			exact_acts.push_back(index);
			is_ordered[index] = true;
		}
		
		cout << "\n[Disjunctive] Exact activities priority (old): ";
		for( int i = 0; i < (int)exact_acts.size(); ++i ) {
			cout << auxiliar[i] << " ";
		}
		#endif

		cout << "\n[Disjunctive] Exact activities priority: ";
		for( int i = 0; i < (int)exact_acts.size(); ++i ) {
			cout << exact_acts[i] << " ";
		}
		cout << endl << endl;
	}
		break;

	case TotalTardiness:
	case MaxTardiness:
	{
		// collect activities that are not fixed
		memset(act_checked, false, sizeof(bool)*num_acts);
		for( int i = 0; i < permutation.getSize(); ++i ) {
			if( permutation[i].isFixed() ) {
				act_checked[permutation[i].getValue()] = true;
			}
		}

		exact_acts.clear();
		for( int i = 0; i < num_acts; ++i ) {
			if( !act_checked[i] ) {
				exact_acts.push_back(i);
			}
		}
		if( exact_acts.size() <= 2 )
			return;

		// first activity has the maximum deadline, using the largest processing time
		// as tie-breaking
		int max_release = -1;
		int max_proc = -1;
		int first = INF;
		for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
			int i = *it;
			if( (max_release < acts[i]->start_time.getStartMin() + acts[i]->processing )
					|| (max_release == acts[i]->start_time.getStartMin() + acts[i]->processing && max_proc < acts[i]->processing) )
			{
				first = it - exact_acts.begin();
				max_release = acts[i]->start_time.getStartMin() + acts[i]->processing;
				max_proc = acts[i]->processing;
			}
		}

		// put activity in the first position
		int tmp = exact_acts[0];
		exact_acts[0] = exact_acts[first];
		exact_acts[first] = tmp;

		// sort remaining activities by largest (processing time+min_setup_from)
		aux_v.resize(num_acts);
		for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
			int i = *it;
			//aux_v[i] = (acts[i]->start_time.getStartMin() + acts[i]->processing + max_setup_from[i]);
			//aux_v[i] = (acts[i]->start_time.getStartMin() + acts[i]->processing + min_setup_from[i]);
			//aux_v[i] = (acts[i]->duedate - acts[i]->release - acts[i]->processing);
			aux_v[i] = -acts[i]->duedate;
//			if( acts[i]->weight > 0 ) {
//				aux_v[i] *= acts[i]->weight;
//			}
		}
		//sort(exact_acts.begin()+1, exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
		sort(exact_acts.begin(), exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
		//sort(exact_acts.rbegin(), exact_acts.rend());

		cout << "\n[Disjunctive] Exact activities priority: ";
		for( int i = 0; i < (int)exact_acts.size(); ++i ) {
			cout << exact_acts[i] << " ";
		}
		cout << endl << endl;

//		exact_acts.clear();
//			vector<bool> is_ordered( num_acts, false );
//
//			// first activity has the maximum release date, using the largest processing time
//			// as tie-breaking
//			int max_release = -1;
//			int max_proc = -1;
//			int first = INF;
//			for ( int i = 0; i < num_acts; ++i ) {
//				if( (max_release < acts[i]->start_time.getStartMin() + acts[i]->processing )
//						|| (max_release == acts[i]->start_time.getStartMin() + acts[i]->processing && max_proc < acts[i]->processing) )
//				{
//					first = i;
//					max_release = acts[i]->start_time.getStartMin() + acts[i]->processing;
//					max_proc = acts[i]->processing;
//				}
//			}
//
//			// put activity in the first position
//			exact_acts.push_back( first );
//			is_ordered[ first ] = true;
//
//			// add remaining activities according to sum of setup times
//			while( (int)exact_acts.size() < num_acts ) {
//				int max_sum = -1;
//				int act = -1;
//				for( int i = 0; i < num_acts; ++i ) {
//					if( is_ordered[i] )
//						continue;
//					int sum = 0;
//					for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
//						if( input_prec[*it][i] ) {
//							sum += setup[*it][i];
//						} else {
//							sum += setup[i][*it];
//						}
//					}
//					if( sum > max_sum ) {
//						max_sum = sum;
//						act = i;
//					}
//				}
//				cout << "act: " << act << " - max_sum: " << max_sum << endl;
//				exact_acts.push_back( act );
//				is_ordered[act] = true;
//			}
//
//			cout << "\n[Disjunctive] Exact activities priority: ";
//			for( int i = 0; i < (int)exact_acts.size(); ++i ) {
//				cout << exact_acts[i] << " ";
//			}
//			cout << endl << endl;
	}
		break;

	case Latency:
	{
		// collect activities that are not fixed
		memset(act_checked, false, sizeof(bool)*num_acts);
		for( int i = 0; i < permutation.getSize(); ++i ) {
			if( permutation[i].isFixed() ) {
				act_checked[permutation[i].getValue()] = true;
			}
		}

		exact_acts.clear();
		for( int i = 0; i < num_acts; ++i ) {
			if( !act_checked[i] ) {
				exact_acts.push_back(i);
			}
		}
		if( exact_acts.size() <= 2 )
			return;

		// first activity has the maximum deadline, using the largest processing time
		// as tie-breaking
		int max_release = -1;
		int max_proc = -1;
		int first = INF;
		for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
			int i = *it;
			if( (max_release < (acts[i]->start_time.getStartMin() + acts[i]->processing)*acts[i]->weight )
					|| (max_release == (acts[i]->start_time.getStartMin() + acts[i]->processing)*acts[i]->weight
							&& max_proc < acts[i]->processing) )
			{
				first = it - exact_acts.begin();
				max_release = acts[i]->weight*(acts[i]->start_time.getStartMin() + acts[i]->processing);
				max_proc = acts[i]->processing;
			}
		}

		// put activity in the first position
		int tmp = exact_acts[0];
		exact_acts[0] = exact_acts[first];
		exact_acts[first] = tmp;

		// sort remaining activities by largest (processing time+min_setup_from)
		aux_v.resize(num_acts);
		vector<int> v2(num_acts);
		for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
			int i = *it;
			aux_v[i] = (acts[i]->start_time.getStartMin() + acts[i]->processing + max_setup_from[i]);
			//aux_v[i] = (acts[i]->start_time.getStartMin() + acts[i]->processing + min_setup_from[i]);
			//aux_v[i] = max_setup_from[i] + acts[i]->processing;
//			aux_v[i] = -acts[i]->duedate;
			if( acts[i]->weight > 0 ) {
				aux_v[i] *= acts[i]->weight;
			}
			//aux_v[i] = -acts[i]->weight;
			//v2[i] = -(acts[i]->start_time.getStartMin() + acts[i]->processing + max_setup_from[i]);
		}
		sort(exact_acts.begin()+1, exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
		//sort(exact_acts.begin()+1, exact_acts.end(), IndexComparatorDescendingSTL2D(aux_v, v2));
		//sort(exact_acts.begin(), exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
		//sort(exact_acts.rbegin(), exact_acts.rend());

		cout << "\n[Disjunctive] Exact activities priority: ";
		for( int i = 0; i < (int)exact_acts.size(); ++i ) {
			cout << exact_acts[i] << " ";
		}
		cout << endl << endl;
	}
		break;

	default:
		cerr << "\nRefinement not specified!" << endl;
		exit(1);
		break;
	}



	//	exact_acts.clear();
//
//	int n_requests = (num_acts-2)/2;
//
//	for( int i = 0; i < n_requests; ++i ) {
//		exact_acts.push_back( 1 + i );
//	}
//
//	for( int i = 0; i < n_requests; ++i ) {
//		exact_acts.push_back( n_requests + 1 + i );
//	}
//
//
//	exact_acts.push_back( 0 );
//	exact_acts.push_back( num_acts - 1 );
//
//	cout << "\n[Disjunctive] Exact activities priority: ";
//	for( int i = 0; i < (int)exact_acts.size(); ++i ) {
//		cout << exact_acts[i] << " ";
//	}
//	cout << endl << endl;
//
//	return;
//
//
//
//	vector<bool> is_ordered(num_acts, false);
//	vector<int> act_order;
//
//	// select pair of activities with largest distance
//	int a1 = 0, a2 = 0;
//	for( int i = 0; i < num_acts; ++i ) {
//		for( int j = 0; j < num_acts-1; ++j ) {  // exclude last task
//			if( setup[i][j] > setup[a1][a2] ) {
//				a1 = i;
//				a2 = j;
//			}
//		}
//	}
//
//	cout << "largest pair: " << a1 << "," << a2 << " ==> " << setup[a1][a2] << endl;
//	is_ordered[a1] = true;
//	is_ordered[a2] = true;
//	act_order.push_back(a1);
//	act_order.push_back(a2);
//
////	while( (int)act_order.size() < num_acts ) {
////		int min_max = -1;
////		int act = -1;
////		for( int i = 0; i < num_acts; ++i ) {
////			if( is_ordered[i] )
////				continue;
////			int min = INF;
////			for( vector<int>::iterator it = act_order.begin(); it != act_order.end(); ++it ) {
////				if( input_prec[*it][i] ) {
////					min = MIN(min, setup[*it][i]);
////				} else if( input_prec[i][*it] ) {
////					min = MIN(min, setup[i][*it]);
////				} else {
////					min = MIN(min, setup[i][*it]);
////					min = MIN(min, setup[*it][i]);
////				}
////			}
////			if( min > min_max ) {
////				min_max = min;
////				act = i;
////			}
////		}
////		cout << "act: " << act << " - min_max: " << min_max << endl;
////		act_order.push_back(act);
////		is_ordered[act] = true;
////	}
//
//	while( (int)act_order.size() < num_acts ) {
//		int min_max = -1;
//		int act = -1;
//		for( int i = 0; i < num_acts; ++i ) {
//			if( is_ordered[i] )
//				continue;
//			int sum = 0;
//			for( vector<int>::iterator it = act_order.begin(); it != act_order.end(); ++it ) {
//				if( input_prec[*it][i] ) {
//					sum += setup[*it][i];
//				} else {
//					sum += setup[i][*it];
//				}
//			}
//			if( sum > min_max ) {
//				min_max = sum;
//				act = i;
//			}
//		}
//		cout << "act: " << act << " - min_max: " << min_max << endl;
//		act_order.push_back(act);
//		is_ordered[act] = true;
//	}
//
//
//	exact_acts = act_order;
//	cout << "\n[Disjunctive] Exact activities priority: ";
//	for( int i = 0; i < (int)exact_acts.size(); ++i ) {
//		cout << exact_acts[i] << " ";
//	}
//	cout << endl << endl;
//
//	return;
//
//////	exit(1);

	// ---------------------------------------------------------
//
//	// collect activities that are not fixed
//	memset(act_checked, false, sizeof(bool)*num_acts);
//	for( int i = 0; i < permutation.getSize(); ++i ) {
//		if( permutation[i].isFixed() ) {
//			act_checked[permutation[i].getValue()] = true;
//		}
//	}
//
//	exact_acts.clear();
//	for( int i = 0; i < num_acts; ++i ) {
//		if( !act_checked[i] ) {
//			exact_acts.push_back(i);
//		}
//	}
//	if( exact_acts.size() <= 2 )
//		return;
//
//	// first activity has the maximum deadline, using the largest processing time
//	// as tie-breaking
//	int max_release = -1;
//	int max_proc = -1;
//	int first = INF;
//	for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
//		int i = *it;
//		if( (max_release < acts[i]->start_time.getStartMin() + acts[i]->processing )
//				|| (max_release == acts[i]->start_time.getStartMin() + acts[i]->processing && max_proc < acts[i]->processing) )
//		{
//			first = it - exact_acts.begin();
//			max_release = acts[i]->start_time.getStartMin() + acts[i]->processing;
//			max_proc = acts[i]->processing;
//		}
//	}
//
//	// put activity in the first position
//	int tmp = exact_acts[0];
//	exact_acts[0] = exact_acts[first];
//	exact_acts[first] = tmp;
//
//	// sort remaining activities by largest (processing time+min_setup_from)
//	aux_v.resize(num_acts);
//	for( vector<int>::iterator it = exact_acts.begin(); it != exact_acts.end(); ++it ) {
//		int i = *it;
//		aux_v[i] = (acts[i]->start_time.getStartMin() + acts[i]->processing + max_setup_from[i]);
//		//aux_v[i] = (acts[i]->start_time.getStartMin() + acts[i]->processing + min_setup_from[i]);
//		//aux_v[i] = max_setup_from[i] + acts[i]->processing;
//		if( acts[i]->weight > 0 ) {
//			aux_v[i] *= acts[i]->weight;
//		}
//	}
//	sort(exact_acts.begin()+1, exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
//	//sort(exact_acts.begin(), exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
//	//sort(exact_acts.rbegin(), exact_acts.rend());
//
//	// H !!!!!!!!!!!!!!!!
//	//sort(exact_acts.begin(), exact_acts.end(), IndexComparatorDescendingSTL(aux_v));
//
//	cout << "\n[Disjunctive] Exact activities priority: ";
//	for( int i = 0; i < (int)exact_acts.size(); ++i ) {
//		cout << exact_acts[i] << " ";
//	}
//	cout << endl << endl;


	// ---------------------------------------------------------

	// for sum of setup times

}


/**
 * Get activity and position of the activity with best sum of setup times to T
 **/
void DisjunctivePropagator::get_activity_min_sum_setup_T(int& pos, int& act) {

	if( obj_type == TotalTardiness ) {
		RevMDD* mdd = get_mdd();
		//int prev_act = -1;
		int min_tard_T;
		DisjunctiveState *source_state, *target_state;

//		pos = 0;
//		while( mdd->arcs_layer[pos].is_unique() ) {
//			prev_act = mdd->arcs_layer[pos++].get_first()->get_val();
//		}

		//cout << "Layer " << pos << endl;
		min_tard_T = INF;

		if( pos > 0 ) {
			for( RevMDD::Arc* arc = mdd->arcs_layer[pos].get_first();
					arc != NULL;
					arc = mdd->arcs_layer[pos].get_next(arc) )
			{
				source_state = get_node_state(pos, arc->get_source());
				target_state = get_node_state(pos+1, arc->get_target());

				if( min_tard_T > source_state->min_tard.getValue() + target_state->min_tard_T.getValue() ) {
					act = arc->get_val();
					min_tard_T = source_state->min_tard.getValue() + target_state->min_tard_T.getValue();
					//cout << "\ttarget: " << arc->get_target() << " - val: " << arc->get_val() << " - min tard T: " << min_tard_T << endl;
				}
			}
		} else {
			for( RevMDD::Arc* arc = mdd->arcs_layer[pos].get_first();
					arc != NULL;
					arc = mdd->arcs_layer[pos].get_next(arc) )
			{
				source_state = get_node_state(pos, arc->get_source());
				target_state = get_node_state(pos+1, arc->get_target());
				if( source_state->min_tard.getValue() + target_state->min_tard_T.getValue() < min_tard_T ) {
					act = arc->get_val();
					min_tard_T = source_state->min_tard.getValue() + target_state->min_tard_T.getValue();
					//cout << "\ttarget2: " << arc->get_target() << " - min tard T: " << target_state->min_tard_T.getValue() << endl;
				}
			}
		}

		//exit(1);
		if( pos == mdd->get_num_layers()-1 ) {
			pos = -1;
		}
		return;
	}

	RevMDD* mdd = get_mdd();
	int prev_act = -1;
	int min_sum_setup_T;
	DisjunctiveState *source_state, *target_state;

	pos = 0;
	while( mdd->arcs_layer[pos].is_unique() ) {
		prev_act = mdd->arcs_layer[pos++].get_first()->get_val();
	}

	//cout << "Layer " << pos << endl;
	min_sum_setup_T = INF;

	if( pos > 0 ) {
		for( RevMDD::Arc* arc = mdd->arcs_layer[pos].get_first();
				arc != NULL;
				arc = mdd->arcs_layer[pos].get_next(arc) )
		{
			source_state = get_node_state(pos, arc->get_source());
			target_state = get_node_state(pos+1, arc->get_target());
			int relaxation_val = source_state->min_sumsetup_R.getValue()
					+ setup[prev_act][arc->get_val()]
					+ target_state->min_sumsetup_T.getValue();

			if( relaxation_val < min_sum_setup_T ) {
				act = arc->get_val();
				min_sum_setup_T = relaxation_val;
				//cout << "\ttarget: " << arc->get_target() << " - val: " << arc->get_val() << " - min sum setup: " << min_sum_setup_T << endl;
			}
		}
	} else {
		for( RevMDD::Arc* arc = mdd->arcs_layer[pos].get_first();
				arc != NULL;
				arc = mdd->arcs_layer[pos].get_next(arc) )
		{
			target_state = get_node_state(pos+1, arc->get_target());
			if( target_state->min_sumsetup_T.getValue() < min_sum_setup_T ) {
				act = arc->get_val();
				min_sum_setup_T = target_state->min_sumsetup_T.getValue();
				//cout << "\ttarget2: " << arc->get_target() << " - min sum setup: " << state->min_sumsetup_T.getValue() << endl;
			}
		}
	}

	//exit(1);
	if( pos == mdd->get_num_layers()-1 ) {
		pos = -1;
	}
	return;
}

/**
 * Get node with best sum of setup times to T
 **/
void DisjunctivePropagator::get_node_min_sum_setup_T(int& layer, int& node) {

	RevMDD* mdd = get_mdd();
	int prev_act = -1;
	int min_sum_setup_T;
	//int act;
	DisjunctiveState *source_state, *target_state;

	int pos = 0;
	while( mdd->arcs_layer[pos].is_unique() ) {
		prev_act = mdd->arcs_layer[pos++].get_first()->get_val();
	}

	//cout << "Layer " << pos << endl;
	min_sum_setup_T = INF;
	layer = -1;
	node = -1;

	if( layer > 0 ) {
		for( RevMDD::Arc* arc = mdd->arcs_layer[pos].get_first();
				arc != NULL;
				arc = mdd->arcs_layer[pos].get_next(arc) )
		{
			source_state = get_node_state(pos, arc->get_source());
			target_state = get_node_state(pos+1, arc->get_target());
			int relaxation_val = source_state->min_sumsetup_R.getValue()
					+ setup[prev_act][arc->get_val()]
					+ target_state->min_sumsetup_T.getValue();

			if( relaxation_val < min_sum_setup_T ) {
				//act = arc->get_val();
				min_sum_setup_T = relaxation_val;
				//cout << "\ttarget: " << arc->get_target() << " - val: " << arc->get_val() << " - min sum setup: " << min_sum_setup_T << endl;
				node = arc->get_target();
				layer = pos+1;
			}
		}
	} else {
		for( RevMDD::Arc* arc = mdd->arcs_layer[pos].get_first();
				arc != NULL;
				arc = mdd->arcs_layer[pos].get_next(arc) )
		{
			target_state = get_node_state(pos+1, arc->get_target());
			if( target_state->min_sumsetup_T.getValue() < min_sum_setup_T ) {
				//act = arc->get_val();
				min_sum_setup_T = target_state->min_sumsetup_T.getValue();
				node = arc->get_target();
				layer = pos+1;
				//cout << "\ttarget2: " << arc->get_target() << " - min sum setup: " << target_state->min_sumsetup_T.getValue() << endl;
			}
		}
	}

	//cout << "Selected node " << node << " in layer " << layer << endl;
	if( pos == mdd->get_num_layers()-1 ) {
		pos = -1;
//		layer = -1;
//		node = -1;
	}
	return;
}



/**
 * Returns if node is refinable with respect to a particular activity
 */
inline DisjunctivePropagator::RefineState DisjunctivePropagator::is_refinable(
		int layer,
		int act,
		RefinementGroup ref_groups)
{
	//cout << "\n[Disjunctive] is_refinable -- test " << endl;
//	if( ref_groups->size <= 1 )
//		return false;

	int node_with_act = 0;
	int node_no_act = 0;

	DisjunctiveState* source_state;
	for( int i = 0; i < ref_groups->size; ++i ) {
		RevMDD::Arc* arc = ref_groups->elements[i];

		if (arc->get_val() == act) {
			++node_with_act;

		} else {

			switch (node_ref_state[arc->get_source()]) {

			case Unknown:
				source_state = get_node_state(layer, arc->get_source());
				if( source_state->all_paths[act] ) {
					++node_with_act;
					node_ref_state[arc->get_source()] = Implied;
				} else if( !source_state->some_paths[act] ) {
					++node_no_act;
					node_ref_state[arc->get_source()] = NotImplied;
				} else {
					node_ref_state[arc->get_source()] = NotRefinable;
					return NotRefinable;
				}
				break;

			case Implied:
				++node_with_act;
				break;

			case NotImplied:
				++node_no_act;
				break;

			case NotRefinable:
				return NotRefinable;
				break;

			default:
				cout << "Error: node state is not known" << endl;
				exit(1);
			}
		}
	}

	if (node_with_act > 0 && node_no_act > 0) {
		return Refinable;
	} else if (node_no_act > 0) {
		return NotImplied;
	}
	return Implied;
}


///**
// * Compute activities that have priority to be exact in any solution
// * TODO: should be recomputed at each search node?
// */
//void DisjunctivePropagator::compute_exact_acts() {
//  int num_acts = acts.size();
//
//  // first activity has the maximum deadline, using the largest processing time
//  // as tie-breaking
//
//  int max_release = -1;
//  int max_proc = -1;
//  for( int i = 0; i < num_acts; ++i) {
//    if( (max_release < acts[i]->start_time.getStartMin())
//        || (max_release == acts[i]->start_time.getStartMin() && max_proc < acts[i]->processing) )
//    {
//      exact_acts_array[0] = i;
//      max_release = acts[i]->start_time.getStartMin();
//      max_proc = acts[i]->processing;
//    }
//  }
//
//  // compose remaining indices
//  for( int i = 1; i < num_acts; ++i ) {
//    exact_acts_array[i] = (i-1 < exact_acts_array[0] ? i-1 : i);
//  }
//
//  // sort remaining activities by largest (processing time+min_setup_from)
//  aux_v.resize(num_acts);
//  for( int i = 0; i < num_acts; ++i ) {
//    aux_v[i] = (acts[i]->release + acts[i]->processing + min_setup_from[i]);
//    //aux_v[i] = (acts[i]->release + acts[i]->processing + max_setup_from[i]);
//    //aux_v[i] = -acts[i]->processing - min_setup_from[i];
//    //aux_v[i] = (acts[i]->weight*acts[i]->processing);
//    //aux_v[i] = acts[i]->weight;
//    //aux_v[i] = -acts[i]->weight;
//  }
//  sort(exact_acts_array+1, exact_acts_array+num_acts, IndexComparatorDescendingSTL(aux_v));
//
//  cout << "\n[Disjunctive] Exact activities priority: ";
//  for( int i = 0; i < num_acts; ++i ) {
//    cout << exact_acts_array[i] << " ";
//  }
//  cout << endl << endl;
//}


/**
 * Apply Held-Karp relaxation
 */
void DisjunctivePropagator::apply_help_karp_relaxation() {

	//return;

	RevMDD* mdd = get_mdd();
	HeldKarpsRelaxation hk_relax( acts.size() );

	// initialize root state
	get_node_state(0, 0)->in_est->clear();

	// build graph
	for ( int l = 0; l < mdd->get_num_layers()-1; ++l ) {

		//cout << "Layer " << l << endl;

		// clean incoming list of the next layer
		for ( RevMDD::Arc* arc = mdd->arcs_layer[l+1].get_first();
				arc != NULL;
				arc = mdd->arcs_layer[l+1].get_next(arc) )
		{
			get_node_state( l+1, arc->get_source() )->in_list->clear();
		}

		// create arcs
		for( RevMDD::Arc* arc = mdd->arcs_layer[l].get_first();
				arc != NULL;
				arc = mdd->arcs_layer[l].get_next(arc) )
		{
			//cout << "\tArc: " << arc->get_val() << endl;

			// add arcs to graph
			DisjunctiveState* state = get_node_state( l, arc->get_source() );
			for ( int i = 0; i < state->in_list->size; ++i ) {
				int prev_value = state->in_list->elements[i];
				//cout << "\tChecking " << prev_value << " <--> " << arc->get_val() << endl;
				if ( prev_value != arc->get_val()
						&&
						!precedes_activity(arc->get_val(), prev_value) )
				{
					hk_relax.add_edge_min_weight(prev_value, arc->get_val(),
							setup[prev_value][arc->get_val()]);
				}
			}
			// add to target list
			get_node_state( l+1, arc->get_target() )->in_list->add( arc->get_val() );
		}
		//cout << endl;
	}

	// compute bound
	//hk_relax.print();
	int bound = hk_relax.compute( obj_var.getMax() );
//	cout << "Current MDD Lower bound: " << obj_var << endl;
//	cout << "Held-Karps bound: " << bound << endl;


	if( bound > obj_var.getMin() ) {
//		cout << "Improved!" << endl;
//		cout << "\tCurrent MDD Lower bound: " << obj_var << endl;
//		cout << "\tHeld-Karps bound: " << bound << endl;
		get_CP().add( obj_var >= bound );
	}
	//cout << endl;
	//exit(1);
}



/**
 * Called when all variables are fixed.
 * For this constraint, we just fix the costs
 **/
void DisjunctivePropagator::all_vars_fixed() {

	IloCP cp = get_CP();

	int start_time;
	int sum_setup_val = 0;
	int makespan_val = acts[permutation[0].getValue()]->release + acts[permutation[0].getValue()]->processing;
	int tardiness = 0;
	int max_tardiness = 0;
	int weighted_completion = acts[permutation[0].getValue()]->weight * makespan_val;

	// minimum time for first activity in the resource
	acts[permutation[0].getValue()]->start_time.setStartMin(
			MAX(acts[permutation[0].getValue()]->release, acts[permutation[0].getValue()]->start_time.getStartMin()));

	for( int i = 1; i < permutation.getSize(); i++ ) {

		//cout << "Activities: " << endl;
		//for (int j = 0; j < permutation.getSize(); ++j) {
		//    cout << "\tj=" << j << ": " << acts[permutation[j].getValue()]->start_time << endl;
		//}

		// minimum start time of activity
		start_time = MAX(acts[permutation[i].getValue()]->start_time.getStartMin(),
						 makespan_val + setup[acts[permutation[i-1].getValue()]->id][acts[permutation[i].getValue()]->id]);

	   // cout << "Setting the start time of " << permutation[i].getValue() << " to " << start_time << endl;
	    //cout << "\twhere " << acts[permutation[i].getValue()]->start_time << endl;

		acts[permutation[i].getValue()]->start_time.setStart(start_time);
		acts[permutation[i].getValue()]->start_time.setEnd(start_time+1);
		
		//cout << "Imposing precedence..." << endl;  

		// impose precedence relation with respect to previous activity
		IlcEndBeforeStart(acts[permutation[i-1].getValue()]->start_time,
							acts[permutation[i].getValue()]->start_time,
							setup[acts[permutation[i-1].getValue()]->id][acts[permutation[i].getValue()]->id]);

		// compute makespan / sum of setup times
		sum_setup_val += setup[acts[permutation[i-1].getValue()]->id][acts[permutation[i].getValue()]->id];
		makespan_val = start_time + acts[permutation[i].getValue()]->processing;


		// compute tardiness
		int tardiness_penalty = acts[permutation[i].getValue()]->weight *
			    MAX(0, makespan_val - acts[permutation[i].getValue()]->duedate);
		tardiness += tardiness_penalty;
		max_tardiness = MAX(max_tardiness, tardiness_penalty);

		weighted_completion += acts[permutation[i].getValue()]->weight * makespan_val;
	}

//	cout << "Permutation: " << endl;
//	for( int i = 0; i < permutation.getSize(); ++i ) {
//		cout << "\t" << permutation[i].getValue() << endl;
//	}
//	cout << "sum-setup: " << sum_setup_val << endl << endl;

	// update objective function
	switch( obj_type ) {

	case Makespan:
		cp.add( obj_var >= makespan_val );
		break;

	case SumSetupTimes:
		cp.add( obj_var == sum_setup_val );
		break;

	case TotalTardiness:
		cp.add( obj_var >= tardiness );
		break;

	case MaxTardiness:
		cp.add( obj_var >= max_tardiness );
		break;

	case Latency:
		cp.add( obj_var >= weighted_completion );
		break;

	default:
		break;
	}
}

// Checking if the shortest path is feasible and checking the solution

bool DisjunctivePropagator::check_feasible_solution(int* sequence) {

	// Print the solution
	//cout << "Optimum solution: ";
 	for( int i = 0; i < (int)acts.size(); ++i ) {
    	//cout << sequence[i] << " ";
  	}
  	//cout << endl;

  	// Check the allDifferent constraint
  	counter_task.resize((int)acts.size());

  	for(int i = 0; i < (int)acts.size(); i++){
  		counter_task[i] = 0;
  	}

  	for(int i = 0; i < (int)acts.size(); i++){
		counter_task[sequence[i]]++;

		if(counter_task[sequence[i]] > 1){
			//cout << "The solution violates the alldifferent constraint.\n";
			//cout << "The vehicle visits city " << sequence[i] << " twice \n"; 
			return false; 
		}
	}

	// Check capacity constraint
	int sum_cap(0);
	
	for (int i = 0; i < (int)acts.size(); i++){
		sum_cap += weights[sequence[i]];
		
		if (sum_cap > capacity || sum_cap < 0) {
			//cout << "The solution violates the capacity constraint.\n";
			//cout << "Vehicle's weight at city " << sequence[i] << " is: " << sum_cap << "\n";
		 	return false;
		}
	}

	//Check Precedence constraint
	for (int i = 0; i < (int)acts.size(); i++){
		for(int j = i+1; j < (int)acts.size(); j++){
			if ( precedes_activity(sequence[j], sequence[i]) ) {
				//cout <<"Location " << sequence[j] << " has to be preceded by location " << sequence[i] << endl;
				return false;
			}

		}
	}

	return true;
}


inline void DisjunctivePropagator::fix_solution(int* sequence){
	cout << "Fixing the solution!" << endl;
	for (int i = 0; i < (int)acts.size(); i++){
		permutation[i].setValue(sequence[i]);
	}


}

int DisjunctivePropagator::get_lagr_type(){
	return lagr_type;
}




