/*
 * Held-Karp relaxation
 *
 */

#ifndef HELD_KARP_RELAX_HPP_
#define HELD_KARP_RELAX_HPP_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace boost;


// ------------------------------------------------------
// Types
// ------------------------------------------------------

// Node potentials and edge weights
struct potential_t {
	typedef boost::vertex_property_tag kind;
};
typedef property< potential_t, int >	NodePotential;
typedef property< edge_weight_t, double >	EdgeWeight;


// Graph is an adjacent list where both vertex sets and
// edges sets are represented using std::vectors.
typedef adjacency_list < vecS, vecS, undirectedS,
		NodePotential, EdgeWeight > Graph;


// Edges and vertices of G
typedef graph_traits < Graph >::edge_descriptor 	Edge;
typedef graph_traits < Graph >::vertex_descriptor 	Vertex;
typedef std::pair<int, int> E;

// ----------------------------------------

class HeldKarpsRelaxation {

public:

	// Constructor
	HeldKarpsRelaxation( int num_vertices ) {
		g = Graph( num_vertices );
		weightmap = get(edge_weight, g);
	}

	// Add an edge with particular weight
	void add_edge_min_weight(int u, int v, int w) {
		Edge e;
		bool exists;
		tie( e, exists ) = edge( u, v, g );
		if ( !exists ) {
			e = add_edge(u, v, g).first;
			weightmap[e] = w;
			//cout << "\tadded edge " << u << " <--> " << v << endl;
		} else {
			int curr_w = weightmap[e];
			if ( curr_w > w ) {
				weightmap[e] = w;
			}
		}
	}

	// Compute Held-Karp relaxation
	int compute( int U );

	// Print file in TSP LIB format
	void print();

private:
	Graph g;
	property_map< Graph, edge_weight_t >::type weightmap;
};





#endif /* HELD_KARP_RELAX_HPP_ */
