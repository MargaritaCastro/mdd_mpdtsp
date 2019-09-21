/*
 * held_karp_relax.cpp
 *
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include "held_karp_relax.hpp"
#include "../util/util.hpp"

using namespace std;

inline int myceil( double v ) {
	return( (v - (double)((int)v)) >= 0.2 ? (int)(v+1.0) : (int)v );
}


/**
 * Compute HeldKarps bound
 */
int HeldKarpsRelaxation::compute( int U ) {

//	cout << "HeldKarps Relaxation" << endl;
//	cout << "\ttarget value: " << U << endl;

	int num_v = num_vertices(g);
	double* potentials = new double[ num_v ];
	int* degree = new int[ num_v ];

	// create original edge weight map
	map< Edge, double > original_weights;
	graph_traits<Graph>::edge_iterator ei, eiend;
	for ( tie(ei, eiend) = edges(g); ei != eiend; ++ei ) {
		original_weights[*ei] = weightmap[*ei];
	}

	//int opt[] = {0, 16, 9, 19, 17, 18, 10, 5, 15, 1, 11, 12, 6, 13, 7, 2, 4, 8, 20, 3, 14, 21};
	//int opt[] = {11, 10, 9, 7, 5, 6, 8, 2, 0, 4, 1, 3};



//	int total_w2 = 0;
//	for( int i = 0; i < 21; ++i ) {
//		Edge e; bool exists;
//		tie(e,exists) = edge(v[i], v[i+1], g);
//		cout << exists;
//		total_w2 += original_weights[e];
//	}
//	cout << endl << "total_w2: " << total_w2 << endl;

	// -----------------------------------
	// Initialization
	// -----------------------------------

	memset( potentials, 0.0, sizeof(double)*num_v );

	int num_iterations = 20;
	//cout << "\tnum iterations: " << num_iterations << endl;

	double h_star = -1;
	double h;
	double t_small = 0.001;
	double alpha = 2;
	double beta = 0.5;
	double potential_adj = 0;

	// Algorithm
	for( int c = 0; c < 20; ++c ) {
		//cout << "alpha: " << alpha << endl;
		for ( int k = 0; k < num_iterations; ++k ) {

			// update edge costs
			graph_traits<Graph>::edge_iterator ei, eiend;
			for (tie(ei, eiend) = edges(g); ei != eiend; ++ei) {
				weightmap[*ei] = original_weights[*ei]
				                   - potentials[ source(*ei, g) ]
				                   - potentials[ target(*ei, g) ];
//				cout << "\t(" << source(*ei,g) << "," << target(*ei,g) << ") " << original_weights[*ei] << " - "
//						<< potentials[ source(*ei, g) ] << " - "
//						<< potentials[ target(*ei, g) ] << " --> " << weightmap[*ei] << endl;
			}
//			cout << endl;

			// compute potentials
			//cout << "potentials: " << endl;
			potential_adj = 0;
			for( int v2 = 0; v2 < num_v; ++v2 ) {
				potential_adj += potentials[v2];
				//cout << "\t(" << v2 << "): " << potentials[v2] << endl;
			}
			potential_adj *= 2.0;

//			cout << endl;
//			cout << "total potential: " << potential_adj << endl;

			// get the two largest potentials, corresponding to node v_1
			double max_1 = -1, max_2 = -1;
			int v_max_1 = -1, v_max_2 = -1;
			for ( int v = 0; v < num_v; ++v ) {

				if ( max_1 == -1 ) {
					v_max_1 = v;
					max_1 = potentials[v];
				} else if ( max_2 == -1 ) {
					v_max_2 = v;
					max_2 = potentials[v];
				} else {
					if ( potentials[v] > max_1 ) {
						max_2 = max_1;
						v_max_2 = v_max_1;
						max_1 = potentials[v];
						v_max_1 = v;
					} else if( potentials[v] > max_2 ) {
						max_2 = potentials[v];
						v_max_2 = v;
					}
				}
			}

			//cout << "\t(" << v_max_1 << "), max_1: " << max_1 << " (" << v_max_2 << "), max_2: " << max_2 << endl;
			assert( v_max_1 != -1 );
			assert( v_max_2 != -1 );

//			int v_max_1 = 0, v_max_2 = num_v-1;
//			double max_1 = potentials[v_max_1], max_2 = potentials[v_max_2];
//			from_edges -= max_1;
//			from_edges -= max_2;

//			cout << "Optimum: " << endl;
//			double total_w2 = 0;
//			for( int i = 0; i < 11; ++i ) {
//				Edge e; bool exists;
//				tie(e,exists) = edge(opt[i], opt[i+1], g);
//
//				if ( !exists ) {
//					cout << "Error: edge " << opt[i] << "," << opt[i+1] << " does not exist." << endl;
//					exit(1);
//				}
//
//				cout << exists;
//				total_w2 += weightmap[e];
//
//				cout << "\t(" << source(e,g) << "," << target(e,g) << ") " << original_weights[e] << " - "
//						<< potentials[ source(e, g) ] << " - "
//						<< potentials[ target(e, g) ] << " --> " << weightmap[e] << endl;
//			}
//
//			cout << "\n\tsum edges: " << total_w2 << endl;
//			cout << "\tpotential_adj: " << potential_adj << endl;
//			cout << "\ttotal: " << (total_w2 + potential_adj - max_1 - max_2) << endl;
//			cout << "\tliq: " << 378.0 - (total_w2 + potential_adj - max_1 - max_2) << endl;
//			cout << endl;

			// compute optimum spanning tree
			vector < Edge > spanning_tree;
			kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

			// compute weight of 1-tree and degree
			double total_w = 0;
			memset( degree, 0, sizeof(int)*num_v );
			for (std::vector < Edge >::iterator ei = spanning_tree.begin();
					ei != spanning_tree.end();
					++ei)
			{
				++degree[source(*ei, g)];
				++degree[target(*ei, g)];
				total_w += weightmap[*ei];
//				   std::cout << source(*ei, g) << " <--> " << target(*ei, g)
//				      << " with weight of " << weightmap[*ei]
//				      << std::endl;
			}
			// update degrees for vertices connected to v1
			++degree[ v_max_1 ];
			++degree[ v_max_2 ];


//			cout << "\ttotal_w: " << total_w << " - potential: " << potential_adj << endl;

			h = (double)total_w + potential_adj - max_1 - max_2;
			if ( h > h_star ) {
				h_star = h;
			}
			//cout << "Iteration " << k << ": H=" << h << " - BestH=" << h_star << endl;

			// check if tour was found
			bool all_degree_2 = true;
			for ( int v2 = 0; v2 < num_v && all_degree_2; ++v2 ) {
				all_degree_2 = ( degree[v2] == 2 );
			}
			if ( all_degree_2 ) {
				//cout << "\tTour was found!" << endl;
				return myceil(h_star);
			}

			int adj = 0;
			for ( int v = 0; v < num_v; ++v ) {
				adj += ( 2 - degree[v] ) * ( 2 - degree[v] );
			}
//			cout << "\tdenominator: " << adj << endl;

			double t = alpha * ((double)U - h) / (double)adj;
//			cout << "\tt=" << t << endl;
			if ( t <= t_small ) {
//				//exit(1);
				return myceil(h_star);
			}

			// replace potentials
			for ( int v = 0; v < num_v; ++v ) {
				//cout << "\t" << v << " - potential: " << potentials[v] << " --> ";
				potentials[v] += t * (double)( 2 - degree[v] );
				//cout << potentials[v] << " - degree: " << degree[v] << endl;
			}
		}
		alpha = beta * alpha;
	}

	//cout << endl << endl;


	delete[] potentials;
	delete[] degree;

	//exit(1);

	return myceil(h_star);
}

// Print file in TSP LIB format
void HeldKarpsRelaxation::print() {

	// print graph
	static int fileindex = 0;

	char filename[256];
	sprintf(filename, "tsp/tsp_%d.tsp", fileindex);

	cout << "Index: " << fileindex << endl;
	++fileindex;

	int num_v = num_vertices(g);

	ofstream graph(filename);

	graph << "NAME: " << filename << endl;
	graph << "TYPE: TSP" << endl;
	graph << "COMMENT: automatically generated" << endl;
	graph << "DIMENSION: " << num_v+1 << endl;
	graph << "EDGE_WEIGHT_TYPE: EXPLICIT" << endl;
	graph << "EDGE_WEIGHT_FORMAT: FULL_MATRIX" << endl;
	graph << "EDGE_WEIGHT_SECTION" << endl;

	// root node
	graph << "1000000" << endl;
	for ( int i = 0; i < num_v; ++i ) {
		graph << "0 ";
	}
	graph << endl;

	for ( int i = 0; i < num_v; ++i ) {
		graph << "0 ";
		for ( int j = 0; j < num_v; ++j ) {
			if ( i == j ) {
				graph << "1000000 ";
				continue;
			}
			Edge e; bool exists;
			tie(e, exists) = edge(i,j,g);
			if( exists ) {
				graph << weightmap[e] << " ";
				//cout << "(" << source(e,g) << "," << target(e,g) << ") = " << weightmap[e] << endl;
			} else {
				graph << "1000000" << " ";
			}
		}
		graph << endl;
	}

	graph.close();
}




