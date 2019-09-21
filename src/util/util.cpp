/*
 * Utilities - implementation function
 */

#include "util.hpp"
#include <algorithm>

using namespace std;

/**
 * Compute spanning tree of a graph
 */
vector<int> SpanningTree::compute() {

	// Sort edges
	int max_vertex = -1;
	edges_id.resize(edges.size());
	for( int i = 0; i < (int)edges.size(); ++i ) {
		edges_id[i] = i;
		max_vertex = MAX(edges[i].u, edges[i].v);
	}
	sort(edges_id.begin(), edges_id.end(), (*this));

//	cout << "Sort order: ";
//	for( int i = 0; i < (int)edges.size(); ++i ) {
//		cout << edges_id[i] << " ";
//	}
//	cout << endl;

	// Create vertex group
	vector<int> group(max_vertex+1);
	for( int i = 0; i < (int)group.size(); ++i ) {
		group[i] = i;
	}

	// Create spanning tree
	spanning_tree.clear();
	int c_edge = 0; // current edge
	while( (int)spanning_tree.size() < num_v-1 ) {
		if( c_edge >= (int)edges.size() ) {
			break;
		}
		Edge& edge = edges[edges_id[c_edge]];
		if( find(edge.u, group) != find(edge.v, group) ) {
			spanning_tree.push_back(edges_id[c_edge]);
			unionf(edge.u, edge.v, group);
		}
		c_edge++;
	}

	return spanning_tree;
}
