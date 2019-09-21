/*
 * ------------------------------------------------------------
 * Header - MDD Definitions
 * ------------------------------------------------------------
 */

#ifndef MDD_HPP_
#define MDD_HPP_

#include <cassert>
#include <ilcp/cp.h>
#include <ilcp/cpext.h>
#include <vector>

#include "../core/intset.hpp"
#include "../core/linkedlist.hpp"
#include "../util/util.hpp"

using namespace std;


/**
 * Reversible MDD data structure
 */
class RevMDD {

public:

	/** Constructor */
	RevMDD(IloCP _cp, IlcIntVarArray _vars, int _max_width, bool create_path = true);

	/**
	 * Arc data structure
	 */
	struct Arc : LinkedListElem {

		friend class RevMDD;
		friend class IlcMDDHandlerImp;

		/** Attributes */
		int             hook;               /**< general purpose hook */
		void*           p_hook;             /**< general purpose hook pointer */

		/** General */
		int     get_source()                      { return source_node; }
		int     get_target()                      { return target_node; }
		int     get_layer_id()                    { return layer_id; }
		int     get_val()                      	{ return val; }

	private:

		/** Private attributes */
		int				layer_id;			/**< layer unique integer identifier (sequential) */
		int             val;                /**< value this arc represents */
		int             source_node;        /**< source node */
		int				target_node;		/**< target node (constant) */
		IlcRevInt       target_node_id;     /**< target node id (reversible structure) */

		/** Methods */
		int     get_target_rev()                      { return target_node_id.getValue(); }
		void    set_target_rev(IloCP cp, int node_id) { target_node_id.setValue(cp, node_id); }

		/** Constructor */
		Arc(IloCP _cp, int _layer_id, int _val, int _source_node);
	};


	/*
	 * --------------------------------------------------
	 * Methods for MDD handling
	 * --------------------------------------------------
	 */

	/** Get arc of a particular type */
	Arc* get_arc(int layer, int source_node_id, int val);

	/** Get arc of a particular layer id */
	Arc* get_arc(int layer, int layer_id);


	/** Export to file  using a top-down pass */
	void export_to_gml(const char* filename);


	/*
	 * --------------------------------------------------
	 * Methods for propagation engine
	 * --------------------------------------------------
	 */

	/** Update variable domains according to changes in arcs */
	void update_variable_domains();


	/*
	 * --------------------------------------------------
	 * General
	 * --------------------------------------------------
	 */

	/** MDD info */
	const int       get_max_width();
	const int       get_width(int layer);
	const int       get_num_layers();
	const int		get_max_arc_id();


	/*
	 * --------------------------------------------------
	 * Attributes for MDD handling
	 * --------------------------------------------------
	 */

	// ** Arcs in a layer ***
	RevLinkedList<Arc>      *arcs_layer;       /**< active arcs in a layer */


	/*
	 * --------------------------------------------------
	 * MDD modification control
	 * --------------------------------------------------
	 */

	/** Returns if node needs to be update according to its incoming states */
	bool needs_update_incoming(int layer, int node);

	/** Indicate that incoming arcs of node were modified */
	void modified_in_node(int layer, int node);

	/** Indicate that node was updated due to changes in its incoming arcs */
	void reset_in_node(int layer, int node);

	/** Returns if node needs to be update according to its outgoing states */
	bool needs_update_outgoing(int layer, int node);

	/** Indicate that outgoing arcs of node were modified */
	void modified_out_node(int layer, int node);

	/** Indicate that node was updated due to changes in its outgoing arcs */
	void reset_out_node(int layer, int node);

	/** Updates reversible structure for incoming arcs */
	void update_modified_out(int layer, int node);

	/** Updates reversible structure for outgoing arcs */
	void update_modified_in(int layer, int node);

	/** Update temporary state */
	void set_modified_in(int layer, int node);

	/** Update temporary state */
	void set_modified_out(int layer, int node);



private:

	/*
	 * --------------------------------------------------
	 * Internal MDD control
	 * --------------------------------------------------
	 */

	IloCP                   cp;             /**< solver reference */
	IlcIntVarArray          vars;           /**< problem variables */
	const int               max_width;      /**< MDD max width */
	const int               num_layers;     /**< total layers */

	Arc****                 arc_buffer;     /**< arc buffer: layer/source-node/val */
	Arc***                  arc_layer;      /**< arc buffer: layer/layer_id */


	IntSet                  aux_set;        /**< aux set for domain computation */

	int                     min_dom_val;    /**< minimum domain value */
	int                     max_dom_val;    /**< maximum domain value */

	int*					max_arc_id;		/**< max arc id per layer (min is always 0) */

	bool*                   has_value;      /**< if layer has some arc with a given value */

	/*
	 * --------------------------------------------------
	 * MDD modification control
	 * --------------------------------------------------
	 */

	IlcRevBool***			incoming_modified_rev;		/**< rev: if node incoming arcs were modified (per layer/node) */
	IlcRevBool***			outgoing_modified_rev;		/**< rev: if node outgoing arcs were modified (per layer/node) */

	bool**					incoming_modified;			/**< if node incoming arcs were modified (per layer/node) */
	bool**					outgoing_modified;			/**< if node outgoing arcs were modified (per layer/node) */


	/** Allocate MDD memory */
	void reserve_memory();

	/** Construct initial path MDD */
	void construct_path();

};

/*
 * ------------------------------------------------
 * Inline implementations
 * ------------------------------------------------
 */

/**
 * Get arc of a particular type
 */
inline RevMDD::Arc* RevMDD::get_arc(int layer, int source_node_id, int val) {
	assert( layer >= 0 && layer < num_layers-1 );
	assert( source_node_id >= 0 && source_node_id < max_width );
	assert( val >= min_dom_val && val <= max_dom_val );
	return arc_buffer[layer][source_node_id][val];
}

/**
 * Get arc of a particular layer id
 */
inline RevMDD::Arc* RevMDD::get_arc(int layer, int layer_id) {
	assert( layer >= 0 && layer < num_layers-1 );
	assert( layer_id >= 0 && layer_id <= max_arc_id[layer] );
	return arc_layer[layer][layer_id];
}



/**
 * Arc Constructor
 */
inline RevMDD::Arc::Arc(IloCP _cp, int _layer_id, int _val, int _source_node)
: LinkedListElem(_cp), layer_id(_layer_id), val(_val), source_node(_source_node)
{
}


/**
 * General
 */
inline const int RevMDD::get_max_width() {
	return max_width;
}

inline const int RevMDD::get_num_layers() {
	return num_layers;
}

inline const int RevMDD::get_max_arc_id() {
	return max_arc_id[0];
}


/**
 * MDD modification control
 */

/**
 * Returns if node needs to be update according to its incoming states
 **/
inline bool RevMDD::needs_update_incoming(int layer, int node) {
	return incoming_modified[layer][node];
}

/**
 * Returns if node needs to be update according to its outgoing states
 **/
inline bool RevMDD::needs_update_outgoing(int layer, int node) {
	return outgoing_modified[layer][node];
}


/**
 * Indicate that incoming arcs of node were modified
 **/
inline void RevMDD::modified_in_node(int layer, int node) {
	//cout << "\t\tNode " << layer << "," << node << " scheduled for incoming update - global" << endl;
	incoming_modified[layer][node] = true;
}

/**
 * Indicate that node was updated due to changes in its incoming arcs
 **/
inline void RevMDD::reset_in_node(int layer, int node) {
	incoming_modified[layer][node] = false;
}

/**
 * Indicate that outgoing arcs of node were modified
 **/
inline void RevMDD::modified_out_node(int layer, int node) {
	//cout << "\t\tNode " << layer << "," << node << " scheduled for outgoing update - global" << endl;
	outgoing_modified[layer][node] = true;
}

/**
 * Indicate that node was updated due to changes in its outgoing arcs
 **/
inline void RevMDD::reset_out_node(int layer, int node) {
	outgoing_modified[layer][node] = false;
}

/**
 * Update incoming state
 **/
inline void RevMDD::update_modified_out(int layer, int node) {
	if( outgoing_modified[layer][node] != outgoing_modified_rev[layer][node]->getValue() ) {
		outgoing_modified_rev[layer][node]->setValue(cp, outgoing_modified[layer][node]);
	}
}

/**
 * Update outgoing state
 **/
inline void RevMDD::update_modified_in(int layer, int node) {
	if( incoming_modified[layer][node] != incoming_modified_rev[layer][node]->getValue() ) {
		incoming_modified_rev[layer][node]->setValue(cp, incoming_modified[layer][node]);
	}
}


/**
 * Update incoming state
 **/
inline void RevMDD::set_modified_in(int layer, int node) {
	incoming_modified[layer][node] = incoming_modified_rev[layer][node]->getValue();
}

/**
 * Update outgoing state
 **/
inline void RevMDD::set_modified_out(int layer, int node) {
	outgoing_modified[layer][node] = outgoing_modified_rev[layer][node]->getValue();
}


/**
 * Arc - output stream function
 */
inline std::ostream& operator<<(std::ostream &os, RevMDD::Arc &arc) {
	os << arc.get_source() << " --->(" << arc.get_val() << ") " << arc.get_target();
	return os;
}

/**
 * Arc comparator according to target node
 */
struct ArcTargetNodeComparator {
	bool operator()(RevMDD::Arc* a, RevMDD::Arc* b) {
		return a->get_target() < b->get_target();
	}
};





#endif /* MDD_HPP_ */

