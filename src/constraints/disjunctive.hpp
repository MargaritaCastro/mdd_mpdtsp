//============================================================================
// Disjunctive MDD Propagator - Header File
//============================================================================


#ifndef DISJUNCTIVE_MDD_HPP_
#define DISJUNCTIVE_MDD_HPP_

#include <cassert>
#include <boost/dynamic_bitset.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <ilcp/cpext.h>

#include <queue>
#include <vector>
#include <list>

#include "../core/mdd_propagator.hpp"
#include "../util/util.hpp"
//#include "../codeFiles.h"


/**
 * Precedence Status Constants
 */
#define P_NOTCOMPUTED 		0
#define P_UNKNOWN          1
#define P_BEFORE 			2
#define P_AFTER			  	3

/**
 * Activity properties
 */
#define NO_DUEDATE         -1


using namespace std;


/**
 * Activity
 */
struct Activity {
	int					id;					/**< unique activity identification */
	int					release;			/**< original activity release time */
	int 				deadline;			/**< original activity deadline */
	int					processing;			/**< original activity processing time */
	int                 duedate;           /**< activity due-date */
	int                 weight;            /**< activity tardiness weight */

	IloIntervalVar		ilo_start_time;		/**< model interval variable for start time */
	IlcIntervalVar		start_time;			/**< interval variable for start time*/

	/** Setters */
	void set_processing(int p);
	void set_duedate(int d);
	void set_weight(int w);

	Activity(IloEnv env, int _id, int _release, int _deadline, int _processing);
	Activity(IloEnv env, int _id, int _release, int _deadline, int _processing, int _due_date, int _weight);
	Activity(IloEnv env, IloIntervalVar _ilo_act, int _id, int _release, int _deadline, int _processing);
};

/**
 * Array of activities
 */
typedef vector<Activity*> ActivityArray;


/**
 * Objective function
 */
enum DisjunctiveObjective {
	Makespan,
	SumSetupTimes,
	TotalTardiness,
	MaxTardiness,
	Latency,
	None
};


/**
 * Node state for disjunctive MDD --> here add the Q_min and Q_max
 */
struct DisjunctiveState {

	IlcRevInt					min_sumsetup_R;	  /**< min sum of setup times from root R up to this node */
	IlcRevInt					min_sumsetup_T;	  /**< min sum of setup times from terminal T up to this node */

    IlcRevInt                   min_tard;         /**< min sum of tardiness allowed at this node */
    IlcRevInt                   min_tard_T;       /**< min sum of tardiness from the terminal T */
    IlcRevInt                   max_tard;         /**< max sum of tardiness allowed at this node */

    IlcRevInt                   min_latency;      /**< min latency allowed at this node */
    IlcRevInt                   max_latency;      /**< max latency allowed at this node */

	IlcRevInt					earliest_start;	  /**< earliest start time of sequence starting from this node */
	IlcRevInt					latest_finish;	  /**< latest finish time of sequence starting from this node */

	boost::dynamic_bitset<>		all_paths;		  /**< values in all paths up to this node */
	boost::dynamic_bitset<>		some_paths;		  /**< values in some paths up to this node */

	RedundantList*				in_list;		  /**< list of incoming arcs */

	SetMinValueList*		    in_est;		      /**< incoming arc values w.r.t. to EST */
	SetMaxValueList*			out_lft;		  /**< outgoing arc values w.r.t. to LFT */

    SetMinValueList*            in_setup;         /**< incoming arc values w.r.t. to setup times */
    SetMinValueList*            out_setup;        /**< outgoing arc values w.r.t. to setup times */

	SetMinValueList*            in_tard;          /**< minimum tardiness of incoming arc values */
	SetMaxValueList*            out_tard;         /**< maximum tardiness of outgoing arc values */

	SetMinValueList*            in_latency;       /**< minimum latency of incoming arc values */
	SetMaxValueList*            out_latency;      /**< maximum latency of outgoing arc values */

	SetMinValueList*			in_Q_min;		  /**< minimum capcity of incoming arc values */
	SetMaxValueList*			in_Q_max;		  /**< maximum capcity of incoming arc values */

	SetMinValueList*			out_Q_min;		  /**< minimum capcity of outgoing arc values */
	SetMaxValueList*			out_Q_max;		  /**< maximum capcity of outgoing arc values */
	
	SetMinValueListDouble*		min_lagrangian_bound_R;	   /**< min cost using lagrangian multiplier for the all diff constraint from root R up to this node */
	SetMinValueListDouble*		min_lagrangian_bound_T;	   /**< min cost using lagrangian multiplier for the all diff constraint from terminal T up to this node */

	bool						hall_set;	      /**< if some_paths state defines a hall set */

	int                         est_tmp;          /**< temporary to compute earliest start time */
	int*                        optimum_sol;      /**< optimum partial solution from root node */

	SetMinValueListDouble*		in_est_double;		/**< incoming arc values w.r.t. to EST */
	double                      est_tmp_double;		/**< temporary to compute earliest start time */

	/** Constructor */
	DisjunctiveState(IloCP cp,  int num_acts, int _earliest_start, int _latest_finish);
};


/**
 * Relaxed stated information for a potential node to be created in a layer.
 */
struct DisjunctiveSketch {

	RevMDD::Arc*    arc;               /**< arc the sketch refers to */

	int				min_sumsetup_R;	   /**< min sum of setup times from root R up to this node */
	int				min_sumsetup_T;	   /**< min sum of setup times from terminal T up to this node */

    int             min_tard;          /**< min sum of tardiness from root R up to this node */
    int             min_tard_T;        /**< min sum of tardiness from terminal T up to this node */
    int             max_tard;          /**< max sum of tardiness allowed at this node */

    int             min_latency;       /**< min latency from root R up to this node */
    int             max_latency;       /**< max latency allowed at this node */

	int				earliest_start;	   /**< earliest start time of sequence starting from this node */
	int				latest_finish;	   /**< latest finish time of sequence starting from this node */

	int             size_implied;      /**< size of all_paths state */

	int 			min_Q;			   /**< minimum capcity */
	int 			max_Q;			   /**< maximum capcity */
	
	double			min_lagrangian_bound_R;	   /**< min cost using lagrangian multiplier for the all diff constraint from root R up to this node */
	double			min_lagrangian_bound_T;	   /**< min cost using lagrangian multiplier for the all diff constraint from terminal T up to this node */
};

/**
 * Disjunctive MDD Propagator
 */
class DisjunctivePropagator : public MDDPropagator {

public:

	/** Constructor: Propagator with objective function */
	DisjunctivePropagator(
			vector<Activity*>		&_acts,				/**< activities to be scheduled */
			IloIntervalSequenceVar	_ilo_sequence,		/**< sequence of activities variable */
			IloIntVarArray			_ilo_permutation,	/**< permutation variables */
			int**					_setup,				/**< setup time between activities */
			DisjunctiveObjective    _obj_type,		    /**< objective function type */
			IloIntVar				_obj_var,			/**< objective function variable */
			IloIntArray             _weights,
			IloInt                  _capacity
	);

	/** Constructor: Propagator with no objective function */
	DisjunctivePropagator(
			vector<Activity*>		&_acts,				/**< activities to be scheduled */
			IloIntervalSequenceVar	_ilo_sequence,		/**< sequence of activities variable */
			IloIntVarArray			_ilo_permutation,	/**< permutation variables */
			int**					_setup				/**< setup time between activities */
	);



	/** Destructor */
	virtual ~DisjunctivePropagator() { }


	/*
	 * -------------------------------------
	 * Input method
	 * -------------------------------------
	 */

	/** Set input precedence */
	void set_precedence(int act_indexA, int act_indexB);
	void set_precedences(bool** _input_precedences);
	
	/** Set lagrangian costs */
	void set_cost_task(double* cost, bool debug);
	void set_cost_layer(double* cost, bool debug);
	void set_cost_pre(double* cost, bool debug);
	void set_all_lagr_multipliers();
	void set_lagr_type(int type);
	int get_lagr_type();
	
	/** Set Capacity refinement */
	void set_cap_ref(bool type);


	/*
	 * -------------------------------------
	 * Propagator Base Methods
	 * -------------------------------------
	 */

	/** Initialize constraint */
	void initialize();

	/** Initialize top-down pass */
	void initialize_topdown();

	/** Initialize bottom-up pass */
	void initialize_bottomup();

	/** Filter layer in a top-down perspective */
	void filter_layer_topdown(int layer);

	/** Filter layer in a bottom up perspective */
	void filter_layer_bottomup(int layer);

	/** Partition refinement group (refinement). Max MDD width must be observed */
	void partition_refinement_group(int layer, RefinementGroup* ref_groups, int &n_groups, int n_rounds);

	/** Setup the state of a node according to its incoming arcs. */
	void setup_node_incoming (int layer, int node, RefinementGroup ref_group, bool new_node);

	/** Setup the state of a node according to its outgoing arcs. */
	void setup_node_outgoing(int layer, int node, RefinementGroup ref_group, bool new_node);

	/** Called when all variables are fixed */
	void all_vars_fixed();

	/** If constraint should be processed in a particular layer */
	bool in_layer(int layer);

	/** Print node states */
	void print_states();

	/** Called when MDD propagation round ended */
	void end_round();

	/** Get activity and position of the activity with best sum of setup times to T */
	void get_activity_min_sum_setup_T(int& pos, int& act);

	/** Get node with best sum of setup times to T */
	void get_node_min_sum_setup_T(int& layer, int& node);

    /** Compute shortest path corresponding to optimal solution in M */
    int* compute_shortest_path();
    int* compute_shortest_path(int& value);
    //int* compute_shortest_path(double& value, double* cost_task, double* cost_layer);
    int* compute_shortest_path(double& value, double* cost_task, double* cost_layer, double* cost_pre, int nMDD = 1);

	/** Shortest path computation inside the MDD */
	void save_shortest_path_lagr();
	void save_shortest_path();
	int* get_shortest_path_lagr(double & value);
	int* get_shortest_path(int & value);

    /** Remove node from the MDD */
    void remove_node(int layer, int node);

    /** Fix node in the MDD, removing all other nodes in the same layer */
    void fix_node(int layer, int node);

    /** Set activity selection for refinement as random */
    void set_random_exact_refinement();

    /** Refinement activity selection order */
    void set_exact_act_order(vector<int> &order);

  /** 
   * Extract model variables to solver variables. Called
   * once during CP Optimizer model extraction. You should use
   * current IloCP passed as parameter instead of propagator
   * IloCP, since handler is not initialized yet 
   **/
  void extract(IloCP cp, const IloCPConstraintI* ct);

private:

    /**< State of a node w.r.t. alldiff refinement */
	enum RefineState {
		Unknown = 0,
		Implied,
		NotImplied,
		NotRefinable,
		Refinable
	};


	/** Check if an in arc is not necessarily infeasible **/
	bool is_incoming_infeasible(int layer, RevMDD::Arc* arc);

	/** Check if an out arc is not necessarily infeasible **/
	bool is_outgoing_infeasible(int layer, RevMDD::Arc* arc);

	/** Check if an out arc is not necessarily infeasible according to all diff **/
	bool is_incoming_infeasible_alldiff(int layer, RevMDD::Arc* arc);

	/** Check if an out arc is not necessarily infeasible according to all diff **/
	bool is_outgoing_infeasible_alldiff(int layer, RevMDD::Arc* arc);

	/** Check if an out arc is not necessarily infeasible according to disjunctive **/
	bool is_incoming_infeasible_disjunctive(int layer, RevMDD::Arc* arc);

	/** Check if an out arc is not necessarily infeasible according to disjunctive **/
	bool is_outgoing_infeasible_disjunctive(int layer, RevMDD::Arc* arc);

	/** Check if an in arc is not necessarily infeasible according to the capacity **/
	bool is_incoming_infeasible_capacity(int layer, RevMDD::Arc* arc);

	/** Check if an in arc is not necessarily infeasible according to the capacity  **/
	bool is_outgoing_infeasible_capacity(int layer, RevMDD::Arc* arc);
	
	/** Check if an in arc is not necessarily infeasible according to the alldiff cost **/
	bool is_incoming_infeasible_lagrangian_bound(int layer, RevMDD::Arc* arc);
	
	/** Check if an in arc is not necessarily infeasible according to the alldiff cost  **/
	bool is_outgoing_infeasible_lagrangian_bound(int layer, RevMDD::Arc* arc);

	/** Obtain node state */
	DisjunctiveState* get_node_state(int layer, int index);

	/** Reserve node state/pre-processed memory */
	void reserve_memory();

	/** Partition refinement group based on capacity */
    void partition_exact_capacity(int layer, RefinementGroup* ref_groups, int &n_groups);

    /** Partition refinement group based on full extension */
    void partition_full_extension(int layer, RefinementGroup* ref_groups, int &n_groups);

    /** Partition refinement group based on exact activities */
    void partition_exact_acts(int layer, RefinementGroup* ref_groups, int &n_groups);

    /** Partition refinement group based exclusively on EST */
    void partition_est(int layer, RefinementGroup* ref_groups, int &n_groups);

	/** Get refinement group min setup */
	int	get_min_dist_R(RefinementGroup ref_group);

	/** Get refinement group earliest start time */
	int	get_min_est(RefinementGroup ref_group);

	/** Get refinement group min size implied up */
	int	get_min_all_paths(RefinementGroup ref_group);

	/** Return if min_dist_R are distinct in a group */
	bool is_min_dist_R_distinct(RefinementGroup ref_group);

	/** Return if earliest start times are distinct in a group */
	bool is_est_distinct(RefinementGroup ref_group);

	/** Return if min_size_implied_up are distinct in a group */
	bool is_min_all_paths_distinct(RefinementGroup ref_group);

	/** Print refinement group */
	void print_group(int id, RefinementGroup group);

	/** Get minimum setup time from an activity in the given list */
	int get_min_setup_from(int act, SetMinValueList* out_set);

	/** Get minimum setup time to an activity in the given list */
	int get_min_setup_to(int act, SetMinValueList* in_set);

	/** Print state */
	void print_state(int layer, int node, DisjunctiveState* state);

	/** Return if an activity precedes another one */
	bool precedes_activity(int actA, int actB);

	/** Set initial activity states */
	void set_initial_states();

	/** Return root state */
	DisjunctiveState* root_state();

    /** Return terminal state */
    DisjunctiveState* terminal_state();

    /** Update activity start/end times */
    void update_start_end_times();

    /** Compute activities that have priority to be exact in M */
    void compute_exact_acts();

    /** Returns if activity 'act' is relaxed in node group, i.e.
     * act \in avail_down_g and act \not \in all_down_g */
    bool is_activity_relaxed(int layer, int act, RefinementGroup ref_groups);

    /** Returns if a node is refinable with respect to activity 'act' */
    RefineState is_refinable(int layer, int act, RefinementGroup ref_groups);

    /** Held-Karp Relaxation */
    void apply_help_karp_relaxation();

    /** Check if a path of the MDD is Feasible **/
    bool check_feasible_solution(int* solution);

    /** Fix a solution if the shortest path of the MDD is Feasible **/
    void fix_solution(int* sequence);


	/*
	 * -------------------------------------
	 * Disjunctive attributes
	 * -------------------------------------
	 */

	vector<Activity*>		acts;				/**< activities to be scheduled */
	IloIntervalSequenceVar	ilo_sequence;		/**< model variable for sequence of activities */
	IloIntVarArray			ilo_permutation;	/**< permutation model variables */
	int**					setup;				/**< setup time between activities */

	IlcIntervalSequenceVar	sequence;			/**< variable for sequence of activities */
	IlcIntVarArray			permutation;		/**< permutation variables */

	DisjunctiveObjective	obj_type;			/**< type of objective function to consider */
    IloIntVar				ilo_obj_var;		/**< objective function variable: model object */
	IlcIntVar				obj_var;			/**< objective function variable */

	bool 								has_input_prec;		/**< if precedence info is given as input */
	bool**								input_prec;			/**< precedences input matrix */
	vector< boost::dynamic_bitset<> >	act_precedes;		/**< precedence set of each activity */
	vector< boost::dynamic_bitset<> >	act_succeeds;		/**< sucessor set of each activity */

	bool					random_exact_refinement;		/**< random exact refinement */
	bool					input_refinement_file;			/**< refinement order read from a file */
	bool					cap_refinement;					/**< indicates whether to use the capacity refinement or not */

	IloIntArray             weights;			/**< absolute weight at each location */
	IloInt                  capacity;			/**< maximum vehicle's capacity */
	
	int						lagr_type;			/**< type of lagragian multipliers that we are using:  0 == nothing; 1 == allDiff cost; 3 == capacity cost; 4 == cap+AllDiff; 5 == precedence; 6 == Prec+AllDiff */
	IlcRevFloat*			lagr_cost_task;		/**< lagrangian multipliers for each task */
	IlcRevFloat*			lagr_cost_layer;	/**< lagrangian multipliers for each layer of the MDD */
	IlcRevFloat*			lagr_cost_pre;		/**< lagrangian multipliers related to the precedence constraint */

	double*					lagr_cost_task_static;
	double*					lagr_cost_layer_static;
	double*					lagr_cost_pre_static;
	
	int*					sp_lagr;
	int*					sp_normal;
	double					val_sp_lagr;
	int						val_sp_normal;					



    //	IloIntVar				ilo_sum_setup;		/**< model variable for sum of setup times */
    //	IloIntVar				ilo_makespan;		/**< model variable for makespan */
    //  IloIntVar               ilo_total_tard;     /**< model variable for total tardiness */
    //	IlcIntVar				sum_setup;			/**< variable for sum of setup times */
    //	IlcIntVar				makespan;			/**< variable for makespan */
    //	IlcIntVar               total_tardiness;    /**< variable for total tardiness */


	/*
	 * -------------------------------------
	 * Propagation attributes
	 * -------------------------------------
	 */
	DisjunctiveState***		states;				/**< node states */
	DisjunctiveSketch**		sketches;			/**< node sketches */

	vector<int>				valid_acts;			/**< incoming/outgoing activities to be considered for filter */

	bool*					node_infeasible;    /**< mark if a node is infeasible */
	int						num_infeas_node;	/**< number of infeasible nodes in layer */
	int*					num_node_checked;	/**< number of times node was checked */

	bool					has_input_precedes;	/**< if precedences are given in input */
	bool**					precedes;			/**< if activity 'i' must precede activity 'j' (given by input) */

	int**					precedence_status;	/**< status of precedence relation between 'i' and 'j' */

	int**					mem_min_tardy;		/**< memoized minimum tardiness */
	int**					mem_min_penalty;	/**< memoized minimum delay penalty */
	
	int**					q_min_topdown;	    /**< q min from topdown code*/
	int**					q_max_topdown;	    /**< q max from topdown code*/
	
	double**				lagrangian_bound_topdown; /**< cost of the lagrangian relaxation topdown */

	SetMinValueList*        act_est;            /**< activity EST inferred from the MDD */
	SetMaxValueList*        act_lft;            /**< activity LFT inferred from the MDD */

	boost::dynamic_bitset<>	temp_set;			/**< temporary set for general operations */
    vector<int>             aux_v;              /**< temporary pre-allocated vector for general operations */


	/*
	 * -------------------------------------
	 * Refinement attributes
	 * -------------------------------------
	 */

	vector<int>				available_nodes;		/**< available nodes for refinement */
	vector<int>				taken_nodes;			/**< nodes that are already occupied */

    IlcRevBool*				requires_refinement;	/**< if layer requires refinement */
	int*					min_est;				/**< auxiliary: min est */
	int*					min_size_implied;		/**< auxiliary: min size of implied  */

	int*                    exact_acts_array;   	/**< activities that have priority to be exact in the MDD */
	vector<int>             exact_acts;         	/**< activities that have priority to be exact in the MDD */

	IlcRevBool**            is_act_exact_layer; 	/**< if activity 'i' is exact in layer 'l', indexed by (l,i) */
	IlcRevBool*             is_act_exact;       	/**< if activity 'i' is exact in all MDD */

	bool*                   may_be_refined;     	/**< groups that still may be refined */
	bool*                   node_checked;       	/**< if node was checked */
	bool*                   act_checked;        	/**< if activity was checked */

	RefineState*			node_ref_state;			/**< refinement state of a node */
	bool					block_est_refinement;	/**< if EST refinement should be blocked */

	/*
	 * -------------------------------------
	 * Preprocessed attributes
	 * -------------------------------------
	 */
	int*					min_setup_from;	    /**< min setup from an activity */
	int*					min_setup_to;		/**< min setup to an activity */

    int*                    max_setup_from;     /**< max setup from an activity */
    int*                    max_setup_to;       /**< max setup to an activity */

	bool*					has_previous;		/**< if activity 'i' has immediate predecessor */
	int*					previous_act;		/**< immediate sucessor of activity 'i' */

	bool*					has_next;			/**< if activity 'i' has immediate sucessor */
	int*					next_act;			/**< immediate predecessor of activity 'i' */
	vector<int>				counter_task;		/**< used for shortest path >**/


	/*
	 * -------------------------------------
	 * Control attributes
	 * -------------------------------------
	 */
	boost::random::mt19937 	random_generator;	/**< random generator */
};




/**
 * -----------------------------------------------------------
 * Inline implementations
 * -----------------------------------------------------------
 */


/**
 * Activity constructor
 * Initialize interval variable for activity
 */
inline Activity::Activity(IloEnv env, int _id, int _release, int _deadline, int _processing)
: id(_id), release(_release), deadline(_deadline), processing(_processing),
  duedate(NO_DUEDATE), weight(0)
{
	char name[256];
	sprintf(name, "act[%d]", id);
	ilo_start_time = IloIntervalVar(env, processing, name);

	ilo_start_time.setStartMin(release);
	ilo_start_time.setEndMax(deadline);
}

/**
 * Activity constructor
 * Initialize interval variable for activity
 */
inline Activity::Activity(IloEnv env, int _id, int _release, int _deadline, int _processing, int _duedate, int _weight)
: id(_id), release(_release), deadline(_deadline), processing(_processing),
  duedate(_duedate), weight(_weight)
{
    char name[256];
    sprintf(name, "act[%d]", id);
    ilo_start_time = IloIntervalVar(env, processing, name);

    ilo_start_time.setStartMin(release);
    ilo_start_time.setEndMax(deadline);
}

/**
 * Activity constructor
 * Initialize interval variable for activity
 */
inline Activity::Activity(IloEnv env, IloIntervalVar _ilo_act, int _id, int _release, int _deadline, int _processing)
: id(_id), release(_release), deadline(_deadline), processing(_processing),
  duedate(NO_DUEDATE), weight(0)
{
	char name[256];
	sprintf(name, "act[%d]", id);
	ilo_start_time = _ilo_act;
}


/** Activity setters */
inline void Activity::set_processing(int p) {
  processing = p;
  ilo_start_time.setSizeMin(p);
  ilo_start_time.setSizeMax(p);
}

inline void Activity::set_duedate(int d) {
  duedate = d;
}

inline void Activity::set_weight(int w) {
  weight = w;
}




/**
 * Disjunctive state constructor
 **/
inline DisjunctiveState::DisjunctiveState(IloCP cp,  int num_acts, int _earliest_start, int _latest_finish) {

	min_sumsetup_R.setValue(cp, 0);
	min_sumsetup_T.setValue(cp, 0);

    min_tard.setValue(cp, 0);
    min_tard_T.setValue(cp, 0);
    max_tard.setValue(cp, INF);

    min_latency.setValue(cp, 0);
    max_latency.setValue(cp, INF);

	earliest_start.setValue(cp, _earliest_start);
	latest_finish.setValue(cp, _latest_finish);

	all_paths.resize(num_acts, false);
	some_paths.resize(num_acts, false);

	
	in_est = new( cp.getHeap() ) SetMinValueList(cp, num_acts);
	out_lft = new( cp.getHeap() ) SetMaxValueList(cp, num_acts);

    in_setup = new( cp.getHeap() ) SetMinValueList(cp, num_acts);
    out_setup = new( cp.getHeap() ) SetMinValueList(cp, num_acts);

	in_tard = new( cp.getHeap() ) SetMinValueList(cp, num_acts);
	out_tard = new( cp.getHeap() ) SetMaxValueList(cp, num_acts);

	in_Q_min = new( cp.getHeap() ) SetMinValueList(cp, num_acts);
	in_Q_max = new( cp.getHeap() ) SetMaxValueList(cp, num_acts);
	out_Q_min = new( cp.getHeap() ) SetMinValueList(cp, num_acts);
	out_Q_max = new( cp.getHeap() ) SetMaxValueList(cp, num_acts);
	
	min_lagrangian_bound_R = new( cp.getHeap() ) SetMinValueListDouble(cp, num_acts);
	min_lagrangian_bound_T = new( cp.getHeap() ) SetMinValueListDouble(cp, num_acts);
	
	in_est_double = new( cp.getHeap() ) SetMinValueListDouble(cp, num_acts);

	in_latency = new( cp.getHeap() ) SetMinValueList(cp, num_acts);
	out_latency = new( cp.getHeap() ) SetMaxValueList(cp, num_acts);

	in_list = new( cp.getHeap() ) RedundantList(cp, num_acts);

	optimum_sol = new( cp.getHeap() ) int[num_acts];
	memset(optimum_sol, 0, sizeof(int)*num_acts);
}


/**
 * Disjunctive constructor: objective function +
 **/
inline DisjunctivePropagator::DisjunctivePropagator(
		ActivityArray			&_acts,				/**< activities to be scheduled */
		IloIntervalSequenceVar	_ilo_sequence,		/**< sequence of activities variable */
		IloIntVarArray			_ilo_permutation,	/**< permutation variables */
		int**					_setup,				/**< setup time between activities */
		DisjunctiveObjective    _obj_type,          /**< objective function type */
		IloIntVar				_obj_var,          	/**< objective variable */
		IloIntArray             _weights,			/**< absolute weight at each location */
		IloInt                  _capacity			/**< maximum vehicle's capacity */

) : acts(_acts), ilo_sequence(_ilo_sequence), ilo_permutation(_ilo_permutation),
	setup(_setup), obj_type(_obj_type), ilo_obj_var(_obj_var), weights(_weights),
	capacity(_capacity)
{
	// allocate activity precedence
	has_input_prec = false;
	int num_acts = acts.size();
	act_precedes.resize(num_acts);
	act_succeeds.resize(num_acts);
	for( int i = 0; i < num_acts; ++i ) {
		act_precedes[i].resize(num_acts, false);
		act_succeeds[i].resize(num_acts, false);
	}
	input_prec = new bool*[num_acts];
	for( int i = 0; i < num_acts; ++i ) {
		input_prec[i] = new bool[num_acts];
		memset(input_prec[i], false, sizeof(bool)*num_acts);
	}

	random_exact_refinement = false; // TODO: make it a parameter
	input_refinement_file = false;
    lagr_type = 0;
    cap_refinement= 0;
}



/**
 * Disjunctive constructor: no objective function
 **/
inline DisjunctivePropagator::DisjunctivePropagator(
		ActivityArray			&_acts,				/**< activities to be scheduled */
		IloIntervalSequenceVar	_ilo_sequence,		/**< sequence of activities variable */
		IloIntVarArray			_ilo_permutation,	/**< permutation variables */
		int**					_setup				/**< setup time between activities */

) : acts(_acts), ilo_sequence(_ilo_sequence), ilo_permutation(_ilo_permutation),
	setup(_setup), obj_type(None)
{
	has_input_prec = false;
	// allocate activity precedence
	int num_acts = acts.size();
	act_precedes.resize(num_acts);
	act_succeeds.resize(num_acts);
	for( int i = 0; i < num_acts; ++i ) {
		act_precedes[i].resize(num_acts, false);
		act_succeeds[i].resize(num_acts, false);
	}
	input_prec = new bool*[num_acts];
	for( int i = 0; i < num_acts; ++i ) {
		input_prec[i] = new bool[num_acts];
		memset(input_prec[i], false, sizeof(bool)*num_acts);
	}
	random_exact_refinement = false;
	input_refinement_file = false;
    lagr_type = 0;
    cap_refinement= 0;
}



/**
 * Initialization: extract variables
 */
inline void DisjunctivePropagator::initialize() {

	cout << "[Disjunctive] Initializing..." << endl;

	IloCP cp = get_CP();

	// extract activities
	cout << "\textracting activities..." << endl;
	for( ActivityArray::iterator act = acts.begin(); act != acts.end(); ++act ) {
		(*act)->start_time = cp.getInterval((*act)->ilo_start_time);
	}

	// extract sequence variables
	cout << "\textracting sequencing vars..." << endl;
	permutation = cp.getIntVarArray(ilo_permutation);
	if( ilo_sequence.getImpl() != NULL ) {
		sequence = cp.getIntervalSequence(ilo_sequence);
	}

	// extract objective function variable
	if( obj_type != None ) {
	  cout << "\textracting obj function..." << endl;
	  obj_var = cp.getIntVar(ilo_obj_var);
	}

	// extract precedeces
	has_input_precedes = false;

	// initialize random generator
	const unsigned int seed = std::time( 0 );
	random_generator.seed( seed );
	cout << "\tRandom seed: " << seed << endl;

	// reserve pre-processing/node state memory
	reserve_memory();

	// set initial node states
	set_initial_states();

	// compute initial exact set
	compute_exact_acts();
	
	// by default we dont use the capacity refinement
	//cap_refinement = false;

	cout << "[Disjunctive] done" << endl;
}



/**
 * Obtain node state
 */
inline DisjunctiveState* DisjunctivePropagator::get_node_state(int layer, int index) {
	assert( layer >= 0 && layer < get_mdd()->get_num_layers() );
	assert( index >= 0 && index < get_mdd()->get_max_width() );
	return states[layer][index];
}



/**
 * Initialize top-down pass
 **/
inline void DisjunctivePropagator::initialize_topdown() {

	//	cout << endl;
	//	cout << " ------------------------------------------------------------- ";
	//	cout << endl;
	//	for( ActivityArray::iterator act = acts.begin(); act != acts.end(); ++act ) {
	//		cout << (*act)->start_time << endl;
	//	}

    // used to compute activity time bounds inferred by the MDD
    act_est->clear();
    act_lft->clear();

	// initialize root node states related to alldiff constraint
	states[0][0]->all_paths.reset();
	states[0][0]->some_paths.reset();

	states[0][0]->hall_set = false;
	states[0][0]->in_est->clear();

	// update precedence relations
	if( sequence.getImpl() != NULL ) {
		for( unsigned int i = 0; i < acts.size(); ++i ) {
			memset( precedence_status[i], P_NOTCOMPUTED, sizeof(int)*acts.size() );
		}
//		for( unsigned int i = 0; i < acts.size(); ++i ) {
//			for( unsigned int j = 0; j < acts.size(); ++j ) {
//				if( i != j ) {
//					precedes[i][j] = sequence.isBefore(acts[i]->start_time, acts[j]->start_time);
//				}
//			}
//		}
	}

	// initialize activities that are exact
	IloCP cp = get_CP();
	for( int i = 0; i < (int)acts.size(); ++i ) {
	  if( permutation[i].isBound() && !is_act_exact[permutation[i].getValue()] ) {
	    int v = permutation[i].getValue();
	    is_act_exact[v].setValue(cp, true);
//	    for( int l = 0; l < get_mdd()->get_num_layers(); ++l ) {
//	      is_act_exact_layer[l][v].setValue(cp, true);
//	    }
	  }
	}

	//cout << "**total tardiness: " << total_tardiness << endl;

	// update activities that must be exact
	//compute_exact_acts();
}

/**
 * Initialize bottom-up pass
 **/
inline void DisjunctivePropagator::initialize_bottomup() {

	// initialize terminal node states related to alldiff constraint
	int last_layer = get_mdd()->get_num_layers()-1;

	states[last_layer][0]->all_paths.reset();
	states[last_layer][0]->some_paths.reset();

	states[last_layer][0]->hall_set = false;
	states[last_layer][0]->out_lft->clear();
}



/**
 * If a constraint should be processed in a particular layer
 */
inline bool DisjunctivePropagator::in_layer(int layer) {
	return true;
}

/**
 * Print node states
 **/
inline void DisjunctivePropagator::print_states() {

}

/**
 * Get refinement group min distance
 **/
inline int DisjunctivePropagator::get_min_dist_R(RefinementGroup ref_group) {
	return sketches[ref_group->get_last()->get_layer_id()]->min_sumsetup_R;
}

/**
 * Get refinement group min earliest start time
 **/
inline int DisjunctivePropagator::get_min_est(RefinementGroup ref_group) {
	return sketches[ref_group->get_last()->get_layer_id()]->earliest_start;
}


/**
 * Get refinement group min size implied up
 **/
inline int DisjunctivePropagator::get_min_all_paths(RefinementGroup ref_group) {
	return sketches[ref_group->get_last()->get_layer_id()]->size_implied;
}


/**
 * Return if min_dist_R are distinct in a group
 **/
inline bool DisjunctivePropagator::is_min_dist_R_distinct(RefinementGroup ref_group) {
	return sketches[ref_group->get_first()->get_layer_id()]->min_sumsetup_R
			!= sketches[ref_group->get_last()->get_layer_id()]->min_sumsetup_R;
}

/**
 * Return if earliest start times are distinct in a group
 **/
inline bool DisjunctivePropagator::is_est_distinct(RefinementGroup ref_group) {
	return sketches[ref_group->get_first()->get_layer_id()]->earliest_start
			!= sketches[ref_group->get_last()->get_layer_id()]->earliest_start;
}


/**
 * Return if min_size_implied_up are distinct in a group
 **/
inline bool DisjunctivePropagator::is_min_all_paths_distinct(RefinementGroup ref_group) {
	return sketches[ref_group->get_first()->get_layer_id()]->size_implied
			!= sketches[ref_group->get_last()->get_layer_id()]->size_implied;
}

/**
 * Set activity selection for refinement as random
 **/
inline void DisjunctivePropagator::set_random_exact_refinement() {
	random_exact_refinement = true;
}


/**
 * Print refinement group
 **/
inline void DisjunctivePropagator::print_group(int id, RefinementGroup group) {
	cout << "\tGroup " << id << endl;
	for( int i = 0; i < group->size; i++ ) {
		cout << "\t\tarc (";
		cout << group->elements[i]->get_source();
		cout << "," << group->elements[i]->get_val();
		cout << ") - min_dist_R = " << sketches[group->elements[i]->get_layer_id()]->min_sumsetup_R;
		cout << " - earliest start = " << sketches[group->elements[i]->get_layer_id()]->earliest_start;
		cout << " - size implied = " << sketches[group->elements[i]->get_layer_id()]->size_implied;
		cout << endl;
	}
	cout << "\tmin earliest start: " << min_est[id] << endl;
	cout << "\tmin size all: " << min_size_implied[id] << endl;
}




/**
 * Print state
 **/
inline void DisjunctivePropagator::print_state(int layer, int node, DisjunctiveState* state) {
	cout << "State of node " << layer << "," << node << endl;
	cout << "\tall_paths: { ";
	for( int val = state->all_paths.find_first(); val != (int)state->all_paths.npos; val = state->all_paths.find_next(val) ) {
		cout << val << " ";
	}
	cout << "}" << " - " << state->all_paths.count() << endl;

	cout << "\tsome_paths: { ";
	for( int val = state->some_paths.find_first(); val != (int)state->some_paths.npos; val = state->some_paths.find_next(val) ) {
		cout << val << " ";
	}
	cout << "}" << " - " << state->some_paths.count() << endl;

	cout << "\tin arcs: { ";
	for( int i = 0; i < state->in_est->size; i++ ) {
		cout << state->in_est->list[i] << " ";
	}
	cout << "}" << endl;

	cout << "\tout arcs: { ";
	for( int i = 0; i < state->out_lft->size; i++ ) {
		cout << state->out_lft->list[i] << " ";
	}
	cout << "}" << endl;


	cout << "\thall set: " << state->hall_set << endl;

	cout << "\tearliest start: " << state->earliest_start.getValue() << endl;
	cout << "\tlatest finish: " << state->latest_finish.getValue() << endl;

	cout << "\tmin sum R: " << state->min_sumsetup_R.getValue() << endl;
	cout << "\tmin sum T: " << state->min_sumsetup_T.getValue() << endl;

    cout << "\tmin tard: " << state->min_tard.getValue() << endl;
    cout << "\tmax tard: " << state->max_tard.getValue() << endl;
	cout << endl;
}


/**
 * Get minimum setup time from an activity in the given list. Notice it also takes into
 * account the alldiff constraint as well.
 **/
inline int DisjunctivePropagator::get_min_setup_from(int act, SetMinValueList* out_set) {

	if( out_set->size == 0 )
		return 0;

	int min_setup = INF;

	// notice that we are trying to find the minimum setup to activity 'act'
	// from the activities in the list; hence, we use the preprocessed data
	// min_setup_to

	// todo: remove input 'precedes' (include that in a different version)
	for( int i = 0; i < out_set->size && min_setup > min_setup_to[act]; i++ ) {
		if( out_set->list[i] != act && !precedes[act][out_set->list[i]] && (!precedes_activity(act, out_set->list[i])) ) {
			min_setup = MIN(min_setup, setup[out_set->list[i]][act] );
		}
	}

	//assert( min_setup != INF );
	return min_setup;
}

/**
 * Get minimum setup time from an activity in the given list. Notice it also takes into
 * account the alldiff constraint as well.
 **/
inline int DisjunctivePropagator::get_min_setup_to(int act, SetMinValueList* in_set) {

	if( in_set->size == 0 )
		return 0;

	int min_setup = INF;

	// notice that we are trying to find the minimum setup from activity 'act'
	// to the activities in the list; hence, we use the preprocessed data
	// min_setup_from

	// todo: remove input 'precedes' (include that in a different version)
	for( int i = 0; i < in_set->size && min_setup > min_setup_from[act]; i++ ) {
		if( in_set->list[i] != act && !precedes[in_set->list[i]][act] && !precedes_activity(in_set->list[i],act) ) {
			min_setup = MIN(min_setup, setup[act][in_set->list[i]] );
		}
	}

	return min_setup;
}

/**
 * Return if an activity precedes another one
 **/
inline bool DisjunctivePropagator::precedes_activity(int actA, int actB) {
//	return sequence.isBefore(acts[actA]->start_time, acts[actB]->start_time);
	if( precedence_status[actA][actB] == P_NOTCOMPUTED ) {
		if( sequence.getImpl() != NULL ) {
			if( input_prec[actA][actB] || sequence.isBefore(acts[actA]->start_time, acts[actB]->start_time) ) {
				precedence_status[actA][actB] = P_BEFORE;
				precedence_status[actB][actA] = P_AFTER;

			} else {
				precedence_status[actA][actB] = P_UNKNOWN;
			}
		} else {
			if( input_prec[actA][actB] || acts[actA]->start_time.getEndMax() < acts[actB]->start_time.getStartMin() ) {
				precedence_status[actA][actB] = P_BEFORE;
				precedence_status[actB][actA] = P_AFTER;

			} else {
				precedence_status[actA][actB] = P_UNKNOWN;
			}
		}
	}
	return( precedence_status[actA][actB] == P_BEFORE );
}


/** Return root state */
inline DisjunctiveState* DisjunctivePropagator::root_state() {
  return( states[0][0] );
}

/** Return terminal state */
inline DisjunctiveState* DisjunctivePropagator::terminal_state() {
  return( states[get_mdd()->get_num_layers()-1][0] );
}

/** Called when MDD propagation round ended. It cannot add or remove arcs */
inline void DisjunctivePropagator::end_round() {
  update_start_end_times();
  //apply_help_karp_relaxation();
  //compute_shortest_path();
}



/** 
 * Extract model variables to solver variables. Called
 * once during CP Optimizer model extraction. You should use
 * current IloCP passed as parameter instead of propagator
 * IloCP, since handler is not initialized yet 
 **/
inline void DisjunctivePropagator::extract(IloCP cp, const IloCPConstraintI* ct) {

  // // extract activities
  // for( ActivityArray::iterator act = acts.begin(); act != acts.end(); ++act ) {
  //   ct->use(cp, (*act)->ilo_start_time);
  //   (*act)->start_time = cp.getInterval((*act)->ilo_start_time);
  // }

  // // extract sequence variables
  // permutation = cp.getIntVarArray(ilo_permutation);
  // if( ilo_sequence.getImpl() != NULL ) {
  //   ct->use(cp, ilo_sequence);
  //   sequence = cp.getIntervalSequence(ilo_sequence);
  // }

  // extract objective function variable
  if( obj_type != None ) {
    ct->use(cp, ilo_obj_var);
    obj_var = cp.getIntVar(ilo_obj_var);
  }

  // // extract precedeces
  // has_input_precedes = false;
}





/**
 * Returns if activity 'act' is relaxed in node group, i.e.
 * act \in avail_down_g and act \not \in all_down_g. It computes
 * this indirectly by using the source states
 *
 * TODO: mark nodes that already contain it in all_paths to avoid
 * rechecking
 */
inline bool DisjunctivePropagator::is_activity_relaxed(int layer, int act, RefinementGroup ref_groups) {
  RevMDD::Arc* arc;
  DisjunctiveState* source_state;

  bool exists_in_all_paths = true;   // if activity exists in all all_paths
  bool exists_in_some_paths = false;  // if activity exists in some some_paths

  for( int i = 0; i < ref_groups->size; ++i ) {
    arc = ref_groups->elements[i];
    source_state = get_node_state(layer, arc->get_source());

    if( !exists_in_some_paths ) {
      exists_in_some_paths = (arc->get_val() == act
          || source_state->some_paths[act]);
    }

    if( exists_in_all_paths ) {
      exists_in_all_paths = (arc->get_val() == act
          || source_state->all_paths[act]);
    }

    if( exists_in_some_paths && !exists_in_all_paths ) {
      return true;
    }
  }
  return false;
}

/** Set input precedence */
inline void DisjunctivePropagator::set_precedence(int act_indexA, int act_indexB) {
	has_input_prec = true;
	act_precedes[act_indexA][act_indexB] = true;
	act_succeeds[act_indexB][act_indexA] = true;
	input_prec[act_indexA][act_indexB] = true;
}

/** Set input precedences */
inline void DisjunctivePropagator::set_precedences(bool** _input_precedences) {
	has_input_prec = true;
	for( int i = 0; i < (int)acts.size(); ++i ) {
		for( int j = 0; j < (int)acts.size(); ++j ) {
			input_prec[i][j] = _input_precedences;
			if( i != j && input_prec[i][j] ) {
				act_precedes[i][j] = true;
				act_succeeds[j][i] = true;
			}
		}
	}
}

/** Set cost from the lagrangian multipliers */
inline void DisjunctivePropagator::set_cost_task(double* cost, bool debug){
	lagr_cost_task_static = cost;
	IloCP cp = get_CP();
    
    if(debug){
        cout << "Using Lagrangian - tour cost:";
        for (int i = 0;  i < (int)acts.size(); ++i ){
            lagr_cost_task[i].setValue(cp, cost[i]) ;
            cout << lagr_cost_task[i].getValue() << " ";
        }
        cout << endl;
    }
}

/** Set cost from the lagrangian multipliers */
inline void DisjunctivePropagator::set_cost_layer(double* cost, bool debug){
	lagr_cost_layer_static = cost;
	IloCP cp = get_CP();
    
    if(debug){
        cout << "Using Lagrangian - capacity cost:";
        for (int i = 0;  i < (int)acts.size(); ++i ){
            lagr_cost_layer[i].setValue(cp, cost[i]) ;
            cout << lagr_cost_layer[i].getValue() << " ";
        }
        cout << endl;
    }
}

/** Set cost from the lagrangianmultipliers */
inline void DisjunctivePropagator::set_cost_pre(double* cost, bool debug){
	lagr_cost_pre_static = cost;
	IloCP cp = get_CP();
    
    if(debug){
        cout << "Using Lagrangian - precedence cost:";
        for (int i = 0;  i < (int)acts.size(); ++i ){
            lagr_cost_pre[i].setValue(cp, cost[i]) ;
            cout << lagr_cost_pre[i].getValue() << " ";
        }
        cout << endl;
    }
}

/** Set all the lagrangian multipliers */

inline void DisjunctivePropagator::set_all_lagr_multipliers(){ 
	IloCP cp = get_CP();
	for (int i = 0;  i < (int)acts.size(); ++i ){
		lagr_cost_task[i].setValue(cp, lagr_cost_task_static[i]) ;
		lagr_cost_layer[i].setValue(cp, lagr_cost_layer_static[i]) ;
		lagr_cost_pre[i].setValue(cp, lagr_cost_pre_static[i]) ;
	}

}

/** Set type of lagrangian multipliers */
inline void DisjunctivePropagator::set_lagr_type(int type){
	if(type > 6 || type <0) cout << "Error in the lagragian type" << endl;
	lagr_type = type;
}

/** Set capacity refinement */
inline void DisjunctivePropagator::set_cap_ref(bool type){
	cap_refinement = type;
}

/**
 * Refinement activity selection order
 **/
inline void DisjunctivePropagator::set_exact_act_order(vector<int> &order) {
	assert( order.size() == acts.size() );
	exact_acts = order;
	input_refinement_file = true;
}


/**
 * Remove node from the MDD
 **/
inline void DisjunctivePropagator::remove_node(int layer, int node) {
	//IloCP cp = get_CP();
	RevMDD* mdd = get_mdd();
	//cout << "Removing node " << layer << "," << node << endl;
	if( layer == 0 || layer == mdd->get_num_layers() - 1 ) {
		cout << "Error: cannot remove source or target node" << endl;
		exit(1);
	}
	temp_set.reset();

	// previous layer
	//RevMDD::Arc* previous = mdd->arcs_layer[layer-1].get_header();
	for( RevMDD::Arc* arc = mdd->arcs_layer[layer-1].get_first();
			arc != NULL;
			arc = mdd->arcs_layer[layer-1].get_next(arc) )
	{
		if( arc->get_target() == node ) {
			//previous->next.setValue(cp, arc->next.getValue());
		} else {
			//previous = arc;
			temp_set.set(arc->get_val(), true);
		}
	}
	bool some_removed = false;
	for( IlcIntExpIterator val(permutation[layer-1]); val.ok(); ++val ) {
		if( !temp_set[*val] ) {
			some_removed = true;
			permutation[layer-1].removeValue(*val);
		}
	}
	if( !some_removed ) {
		cout << "Error: permutation variable was not changed (node removal)" << endl;
		exit(1);
	}

//	// current layer
//	//previous = mdd->arcs_layer[layer].get_header();
//	for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
//			arc != NULL;
//			arc = mdd->arcs_layer[layer].get_next(arc) )
//	{
//		if( arc->get_source() == node ) {
//			//previous->next.setValue(cp, arc->next.getValue());
//		} else {
//			//previous = arc;
//			temp_set.set(arc->get_val(), true);
//		}
//	}
	// TODO: remove nodes without incoming arcs
	//cout << "\tdone." << endl;
}

/**
 * Fix node in the MDD, removing all other nodes in the same layer
 **/
inline void DisjunctivePropagator::fix_node(int layer, int node) {
	//IloCP cp = get_CP();
	RevMDD* mdd = get_mdd();

//	cout << "Fixing node " << layer << "," << node << " - var: " << permutation[layer-1] << endl;
//	mdd->export_to_gml("graphs/before-fixing.gml");

	if( layer == 0 || layer == mdd->get_num_layers() - 1 ) {
		cout << "Error: cannot remove source or target node" << endl;
		exit(1);
	}
	temp_set.reset();

	// previous layer
	//RevMDD::Arc* previous = mdd->arcs_layer[layer-1].get_header();
	for( RevMDD::Arc* arc = mdd->arcs_layer[layer-1].get_first();
			arc != NULL;
			arc = mdd->arcs_layer[layer-1].get_next(arc) )
	{
		if( arc->get_target() != node ) {
			//previous->next.setValue(cp, arc->next.getValue());
		} else {
			//previous = arc;
			temp_set.set(arc->get_val(), true);
		}
	}
	bool some_removed = false;
	for( IlcIntExpIterator val(permutation[layer-1]); val.ok(); ++val ) {
		if( !temp_set[*val] ) {
			some_removed = true;
			permutation[layer-1].removeValue(*val);
		}
	}

//	cout << "After fixing " << endl;
//	mdd->export_to_gml("graphs/after-fixing.gml");

	if( !some_removed ) {
		cout << "Error: permutation variable was not changed" << endl;
		exit(1);
	}

//	// current layer
//	//previous = mdd->arcs_layer[layer].get_header();
//	for( RevMDD::Arc* arc = mdd->arcs_layer[layer].get_first();
//			arc != NULL;
//			arc = mdd->arcs_layer[layer].get_next(arc) )
//	{
//		if( arc->get_source() != node ) {
//			//previous->next.setValue(cp, arc->next.getValue());
//		} else {
//			//previous = arc;
//			temp_set.set(arc->get_val(), true);
//		}
//	}
	//cout << "\tdone." << endl;
}





/**
 * --------------------------------------------
 * Comparators
 * --------------------------------------------
 */

/**
 * Comparator by ascending minimum distance to R. Tie breaking
 * is done via size of implied up
 */
struct RefGroup_MinDistR_Comparator {

	int* min_dist;		// group lower bounds
	int* size_implied;	// group upper bounds

	RefGroup_MinDistR_Comparator(int *_min_dist, int* _size_implied)
	: min_dist(_min_dist), size_implied(_size_implied) { }

	inline bool operator()(int gA, int gB) {
		if( min_dist[gA] == min_dist[gB] )
			return size_implied[gA] > size_implied[gB];
		return min_dist[gA] < min_dist[gB]; // PERHAPS THIS IS INVERTED ?!?!?!?
	}
};

/**
 * Comparator by ascending minimum earliest start time. Tie breaking
 * is done via size of implied up
 */
struct RefGroup_Est_Comparator {

	int* min_dist;		// group lower bounds
	int* size_implied;	// group upper bounds

	RefGroup_Est_Comparator(int *_min_dist, int* _size_implied)
	: min_dist(_min_dist), size_implied(_size_implied) { }

	inline bool operator()(int gA, int gB) {
		if( min_dist[gA] == min_dist[gB] )
			return size_implied[gA] > size_implied[gB];
		return min_dist[gA] > min_dist[gB];
	}
};

/**
 * Comparator by ascending minimum earliest start time.
 */
struct GroupESTComp {
    int* min_est;      // group lower bounds
    GroupESTComp(int *_min_est) : min_est(_min_est) { }
    inline bool operator()(int gA, int gB) {
      return( min_est[gA] > min_est[gB] );
    }
};



/**
 * Arc comparator by **descending** minimum distance to R. Tie breaking
 * is done via size of implied up
 */
struct Arc_MinDistR_Comparator {

	DisjunctiveSketch** 	sketches;

	Arc_MinDistR_Comparator(DisjunctiveSketch** _sketches) : sketches(_sketches) { }

	inline bool operator()(RevMDD::Arc* arcA, RevMDD::Arc* arcB) {
		if( sketches[arcA->get_layer_id()]->min_sumsetup_R == sketches[arcB->get_layer_id()]->min_sumsetup_R ) {
			return sketches[arcA->get_layer_id()]->size_implied < sketches[arcB->get_layer_id()]->size_implied;
		}
		return sketches[arcA->get_layer_id()]->min_sumsetup_R > sketches[arcB->get_layer_id()]->min_sumsetup_R;
	}
};


/**
 * Arc comparator by earliest start time (EST). Tie breaking
 * is done via size of implied up
 */
struct Arc_Est_Comparator {

	DisjunctiveSketch** 	sketches;

	Arc_Est_Comparator(DisjunctiveSketch** _sketches) : sketches(_sketches) { }

	inline bool operator()(RevMDD::Arc* arcA, RevMDD::Arc* arcB) {
		if( sketches[arcA->get_layer_id()]->earliest_start == sketches[arcB->get_layer_id()]->earliest_start ) {
			return sketches[arcA->get_layer_id()]->size_implied < sketches[arcB->get_layer_id()]->size_implied;
		}
		return sketches[arcA->get_layer_id()]->earliest_start > sketches[arcB->get_layer_id()]->earliest_start;
	}
};

/**
 * Arc comparator by minimum weight (Q_min). Tie breaking
 * is done via size of implied up
 */
struct Arc_Q_Comparator {

	DisjunctiveSketch** 	sketches;

	Arc_Q_Comparator(DisjunctiveSketch** _sketches) : sketches(_sketches) { }

	inline bool operator()(RevMDD::Arc* arcA, RevMDD::Arc* arcB) {
		if( sketches[arcA->get_layer_id()]->min_Q == sketches[arcB->get_layer_id()]->min_Q ) {
			return sketches[arcA->get_layer_id()]->size_implied < sketches[arcB->get_layer_id()]->size_implied;
		}
		return sketches[arcA->get_layer_id()]->min_Q > sketches[arcB->get_layer_id()]->min_Q;
	}
};



#endif // DISJUNCTIVE_MDD_HPP

