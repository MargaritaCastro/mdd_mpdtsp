//============================================================================
// General definitions & utilities
//============================================================================


#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <algorithm>
#include <cassert>
#include <cstring>
#include <ilcp/cp.h>
#include <ilcp/cpext.h>
#include <vector>
#include "stats.hpp"
#include <limits>


using namespace std;


/**
 * -------------------------------------------------------------
 * Constants
 * -------------------------------------------------------------
 */
const int INF = numeric_limits<int>::max();        /**< infinite value for our purposes*/

/**
 * -------------------------------------------------------------
 * Macros
 * -------------------------------------------------------------
 */
#define MIN(_a_,_b_)              ((_a_ < _b_) ? _a_ : _b_)
#define MAX(_a_,_b_)              ((_a_ > _b_) ? _a_ : _b_)
#define FLOOR(_a_)                ( (int)_a_ )
#define CEIL(_a_)                 ( ((int)_a_ == _a_) ? _a_ : (((int)_a_)+1) )


/**
 * -------------------------------------------------------------
 * Global types
 * -------------------------------------------------------------
 */
typedef vector< vector< vector< vector<int> > > > v_int4;
typedef vector< vector< vector<int> > > v_int3;
typedef vector< vector<int> > v_int2;
typedef vector<int> v_int;



/**
 * -------------------------------------------
 * Axiliary Comparators
 * -------------------------------------------
 */
template <class T>
struct IndexComparatorDescending {
	int* v;
	IndexComparatorDescending(int* _v) : v(_v) { }
	bool operator()(const int a, const int b) {
		return( v[a] < v[b] );
	}
};

struct IndexComparatorDescendingSTL {
	vector<int> &v;
	IndexComparatorDescendingSTL(vector<int> &_v) : v(_v) { }
	bool operator()(const int a, const int b) {
		return( v[a] > v[b] );
	}
};

struct IndexComparatorDescendingSTL2D {
	vector<int> &v1;
	vector<int> &v2;
	IndexComparatorDescendingSTL2D(vector<int> &_v1, vector<int> &_v2)
	: v1(_v1), v2(_v2) { }
	bool operator()(const int a, const int b) {
		if( v1[a] == v1[b] )
			return( v2[a] > v2[b] );
		return v1[a] > v2[b];
	}
};



/**
 * -------------------------------------------------------------
 * Global utility functions
 * -------------------------------------------------------------
 */

/**
 * Perform inplace union of two vectors in array 'a'
 */
inline void inplace_union(std::vector<int>& a, const std::vector<int>& b)
{
	int mid = a.size(); //Store the end of first sorted range

	//First copy the second sorted range into the destination vector
	std::copy(b.begin(), b.end(), std::back_inserter(a));

	//Then perform the in place merge on the two sub-sorted ranges.
	std::inplace_merge(a.begin(), a.begin() + mid, a.end());

	//Remove duplicate elements from the sorted vector
	a.erase(std::unique(a.begin(), a.end()), a.end());
}

/**
 * -------------------------------------------------------------
 * Global objects
 * -------------------------------------------------------------
 */


/**
 * Set that maintains pairs (value_1,value_2), where no
 * two elements have the same element 'value_1'. For a
 * particular 'value_1', it keeps the smallest 'value_2' added
 */

struct SetMinValueListDouble {
	bool*       has_element;    /**< if list contains element */
	double*     list;           /**< list of elements in the set */
	double*     min_val;        /**< value of the i-th element */
	double*     aux_val;        /**< aux value of the i-th element */
	int         size;           /**< list current size */
	int         max_size;       /**< maximum list size */
	void**      pointer_list;   /**< auxiliary pointer list */

	/**
	 * Constructor
	 */
	SetMinValueListDouble(IloCP cp, int _max_size) : size(0), max_size(_max_size) {
		has_element = new( cp.getHeap() ) bool[max_size];

		list = new( cp.getHeap() ) double[max_size];
		min_val = new( cp.getHeap() ) double[max_size];
		aux_val = new( cp.getHeap() ) double[max_size];
		pointer_list = new( cp.getHeap() ) void*[max_size];
		memset(has_element, false, sizeof(bool)*max_size);
	}
	
	
	/** Add element to list */
	void add(int id, double val) {
		assert( id >= 0 && id < max_size );
		if( !has_element[id] ) {
			assert( size < max_size );
			has_element[id] = true;
			min_val[id] = val;
			list[size++] = id;
		} else {
			min_val[id] = MIN(min_val[id], val);
		}
	}

	/** Add element to list (with a pointer) */
	void add(int value_1, double value_2, void* pointer) {
		assert( value_1 >= 0 && value_1 < max_size );
		if( !has_element[value_1] ) {
			assert( size < max_size );
			has_element[value_1] = true;
			min_val[value_1] = value_2;
			pointer_list[value_1] = pointer;
			list[size++] = value_1;
		} else {
			if( min_val[value_1] > value_2 ) {
				min_val[value_1] = value_2;
				pointer_list[value_1] = pointer;
			}
		}
	}

	/** Check if element is in the list */
	bool contains(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		return( has_element[value_1] );
	}

	double value(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		assert( has_element[value_1] );
		return( min_val[value_1] );
	}

	/** Returns element pointer */
	void* pointer(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		assert( has_element[value_1] );
		return( pointer_list[value_1] );
	}

	/** Clear list */
	void clear() {
		size = 0;
		memset(has_element, false, sizeof(bool)*max_size);
	}

};

struct SetMinValueList {

	bool*       has_element;    /**< if list contains element */
	int*        list;           /**< list of elements in the set */
	int*        min_val;        /**< value of the i-th element */
	int*        aux_val;        /**< aux value of the i-th element */
	int         size;           /**< list current size */
	int         max_size;       /**< maximum list size */
	void**      pointer_list;   /**< auxiliary pointer list */

	/**
	 * Constructor
	 */
	SetMinValueList(IloCP cp, int _max_size) : size(0), max_size(_max_size) {
		has_element = new( cp.getHeap() ) bool[max_size];

		list = new( cp.getHeap() ) int[max_size];
		min_val = new( cp.getHeap() ) int[max_size];
		aux_val = new( cp.getHeap() ) int[max_size];
		pointer_list = new( cp.getHeap() ) void*[max_size];
		memset(has_element, false, sizeof(bool)*max_size);
	}

	/** Add element to list */
	void add(int id, int val) {
		assert( id >= 0 && id < max_size );
		if( !has_element[id] ) {
			assert( size < max_size );
			has_element[id] = true;
			min_val[id] = val;
			list[size++] = id;
		} else {
			min_val[id] = MIN(min_val[id], val);
		}
	}


	/** Add element to list (with auxiliary) */
	void add(int id, int val, int aux) {
		assert( id >= 0 && id < max_size );
		if( !has_element[id] ) {
			assert( size < max_size );
			has_element[id] = true;
			min_val[id] = val;
			aux_val[id] = aux;
			list[size++] = id;
		} else {
			if( min_val[id] > val ) {
				min_val[id] = val;
				aux_val[id] = aux;
			} else if( min_val[id] == val ) {
				aux_val[id] = MIN(aux, aux_val[id]);
			}
		}
	}


	/** Add element to list (with a pointer) */
	void add(int value_1, int value_2, void* pointer) {
		assert( value_1 >= 0 && value_1 < max_size );
		if( !has_element[value_1] ) {
			assert( size < max_size );
			has_element[value_1] = true;
			min_val[value_1] = value_2;
			pointer_list[value_1] = pointer;
			list[size++] = value_1;
		} else {
			if( min_val[value_1] > value_2 ) {
				min_val[value_1] = value_2;
				pointer_list[value_1] = pointer;
			}
		}
	}

	/** Check if element is in the list */
	bool contains(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		return( has_element[value_1] );
	}

	/** Returns corresponding pair of the element */
	int value(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		assert( has_element[value_1] );
		return( min_val[value_1] );
	}

	/** Returns corresponding auxiliary value of the element */
	int aux(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		assert( has_element[value_1] );
		return( aux_val[value_1] );
	}


	/** Returns element pointer */
	void* pointer(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		assert( has_element[value_1] );
		return( pointer_list[value_1] );
	}


	/** Clear list */
	void clear() {
		size = 0;
		memset(has_element, false, sizeof(bool)*max_size);
	}
};

/**
 * Set that maintains pairs (value_1,value_2), where no
 * two elements have the same element 'value_1'. For a
 * particular 'value_1', it keeps the smallest 'value_2' added
 */
struct SetMaxValueList {

	bool*       has_element;    /**< if list contains element */
	int*        list;           /**< list of elements in the set */
	int*        max_val;         /**< value of the i-th element */
	int         size;           /**< list current size */
	int         max_size;       /**< maximum list size */

	/**
	 * Constructor
	 */
	SetMaxValueList(IloCP cp, int _max_size) : size(0), max_size(_max_size) {
		has_element = new( cp.getHeap() ) bool[max_size];

		list = new( cp.getHeap() ) int[max_size];
		max_val = new( cp.getHeap() ) int[max_size];
		memset(has_element, false, sizeof(bool)*max_size);
	}

	/** Add element to list */
	void add(int value_1, int value_2) {
		assert( value_1 >= 0 && value_1 < max_size );
		if( !has_element[value_1] ) {
			assert( size < max_size );
			has_element[value_1] = true;
			max_val[value_1] = value_2;
			list[size++] = value_1;
		} else {
			max_val[value_1] = MAX(max_val[value_1], value_2);
		}
	}

	/** Check if element is in the list */
	bool contains(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		return( has_element[value_1] );
	}

	/** Returns corresponding pair of the element */
	int value(int value_1) {
		assert( value_1 >= 0 && value_1 < max_size );
		assert( has_element[value_1] );
		return( max_val[value_1] );
	}

	/** Clear list */
	void clear() {
		size = 0;
		memset(has_element, false, sizeof(bool)*max_size);
	}
};


/**
 * List that keeps a duplicate bitset and a list for easy iterate
 */
struct RedundantList {

	bool*       has_element;    /**< if list contain element */
	int*        elements;       /**< list elements */
	int         size;           /**< list size */
	int         max_size;       /**< maximum list size */

	/**
	 * Constructor
	 */
	RedundantList(IloCP cp, int _max_size) : size(0), max_size(_max_size) {
		has_element = new( cp.getHeap() ) bool[max_size];
		elements = new( cp.getHeap() ) int[max_size];
		memset(has_element, false, sizeof(bool)*max_size);
	}

	/**
	 * Add element to list
	 */
	void add(int elem) {
		assert( elem >= 0 && elem < max_size );
		if( !has_element[elem] ) {
			assert( size < max_size );
			has_element[elem] = true;
			elements[size++] = elem;
		}
	}

	/**
	 * Check if element is in the list
	 */
	bool contains(int elem) {
		assert( size < max_size );
		assert( elem >= 0 && elem < max_size );
		return( has_element[elem] );
	}

	/**
	 * Clear list
	 */
	void clear() {
		size = 0;
		memset(has_element, false, sizeof(bool)*max_size);
	}

};


/**
 * Sparse set data structure
 *
 * Efficient operations for add(), remove(), contains(), clear(),
 * and forall()
 */
struct SparseSet {

	int*   sparse;    /**< sparse elements */
	int*   dense;     /**< dense elements */
	int    members;   /**< number of elements */

	int    max_size;  /**< maximum list size */

	/**
	 * Constructor
	 */
	SparseSet(IloCP cp, int _max_size) : members(0), max_size(_max_size) {
		dense = new( cp.getHeap() ) int[max_size];
		sparse = new( cp.getHeap() ) int[max_size];
		memset( dense, 0, sizeof(int)*max_size );
		memset( sparse, 0, sizeof(int)*max_size );
	}

	/**
	 * Add element to set (check if it is redundant)
	 */
	inline void add(int k) {
		assert ( k >= 0 && k < max_size );
	    if ( contains(k) )
	    	return;
	    sparse[k] = members;
	    dense[members] = k;
	    members++;
	}

	/**
	 * Add element to set (do not check if it is redundant)
	 */
	inline void add_not_redundant(int k) {
		assert ( k >= 0 && k < max_size );
	    sparse[k] = members;
	    dense[members] = k;
	    members++;
	}


	/**
	 * Check if element is in the list
	 */
	inline bool contains(int k) {
		assert ( k >= 0 && k < max_size);
		return ( sparse[k] < members && dense[sparse[k]] == k );
	}

	/**
	 * Clear list
	 */
	inline void clear() {
		members = 0;
	}

};




/**
 * Stream output function
 */
inline std::ostream& operator<<(std::ostream &os, RedundantList &list) {
	os << "[ ";
	for( int j = 0; j < list.size; j++ ) {
		cout << list.elements[j] << " ";
		assert( list.has_element[list.elements[j]] );
	}
	os << "]";
	return os;
}


/**
 * Simple list of elements
 */
template <class T>
struct SimpleList {

	T*        	elements;       /**< list elements */
	int         size;           /**< list size */
	int         max_size;       /**< maximum list size */

	/**
	 * Constructor
	 */
	SimpleList(IloCP cp, int _max_size) : size(0), max_size(_max_size) {
		elements = new( cp.getHeap() ) T[max_size];
		//elements = new int[max_size];
	}

	/**
	 * Add element to list
	 */
	inline void add(T elem) {
		assert( size < max_size );
		elements[size++] = elem;
	}

	/**
	 * Clear list
	 */
	inline void clear() {
		size = 0;
	}

	/**
	 * Return first element
	 */
	inline T get_first() {
		assert( size > 0 );
		return elements[0];
	}


	/**
	 * Return last element
	 */
	inline T get_last() {
		return elements[size-1];
	}

	/**
	 * Pop-up last element
	 */
	inline T pop_up() {
		return elements[--size];
	}

	/**
	 * Return if the list is empty
	 */
	inline bool is_empty() {
		return( size == 0 );
	}
};


/**
 * Simple representation of an integer set by booleans
 */
struct BoolSet {

	int         max_size;        /**< maximum list size */
	bool*       is_in;           /**< if element is in the set */


	/**
	 * Empty constructor
	 */
	BoolSet() : max_size(0) {

	}

	/*
	 * Initialization
	 */
	void initialize(IloCP cp, int _max_size) {
		max_size = _max_size;
		is_in = new( cp.getHeap() ) bool[max_size];
		memset(is_in, false, sizeof(bool)*max_size);
	}

	/**
	 * Add element to list
	 */
	void add(int elem) {
		assert( elem >= 0 && elem < max_size );
		is_in[elem] = true;
	}

	/**
	 * Check if an element is in the list
	 */
	bool contains(int elem) {
		assert( elem >= 0 && elem < max_size );
		return( is_in[elem] );
	}

	/**
	 * Intersects with another bool set
	 */
	void intersects_with(BoolSet &set) {
		assert( max_size == set.max_size );
		for( int i = 0; i < max_size; i++ )
			is_in[i] &= set.is_in[i];
	}


	/**
	 * Union with another bool set
	 */
	void union_with(BoolSet &set) {
		assert( max_size == set.max_size );
		for( int i = 0; i < max_size; i++ )
			is_in[i] |= set.is_in[i];
	}

	/**
	 * Get size of the list (DO NOT USE THIS OFTEN!!!)
	 */
	int get_size() {
		int size = 0;
		for( int i = 0; i < max_size; i++ )
			size += is_in[i] ? 1 : 0;
		return size;
	}


	/**
	 * Clear list
	 */
	void clear() {
		memset(is_in, false, sizeof(bool)*max_size);
	}
};



/**
 * Class for spanning tree computation on an
 * undirected graph
 */
class SpanningTree {
public:

	/* Graph edge */
	struct Edge {
		int u;
		int v;
		int val;
	};

	/** Constructor */
	SpanningTree(int _num_vertices, vector<Edge> &_edges)
	: num_v(_num_vertices), edges(_edges) { }

	/** Compute the spanning tree */
	vector<int> compute();

	/** Edge comparator according to indices */
	bool operator()(int e1, int e2) {
		return edges[e1].val < edges[e2].val;
	}

private:
	/** Union / find procedures */
	int find(int index, vector<int> &group) {
		return group[index] == index ? index : find(group[index], group);
	}
	void unionf(int indexA, int indexB, vector<int> &group) {
		group[find(indexA,group)] = indexB;
	}

	int  		  num_v;			// number of vertices
	vector<Edge>  &edges;         	// graph edges
	vector<int>  spanning_tree;  	// indices of the spanning tree edges
	vector<int>  edges_id;  	 	// edge indices
};




#endif /* UTIL_HPP_ */
