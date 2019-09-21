/*
 * ------------------------------------------------------------
 * MDD Data Structure Implementations
 * ------------------------------------------------------------
 */

#include "mdd.hpp"


/**
 * Constructor for Reversible MDD
 *
 * Variable ordering is imposed by the order presented in 'vars' argument
 */
RevMDD::RevMDD(IloCP _cp, IlcIntVarArray _vars, int _max_width, bool create_path)
    : cp(_cp), vars(_vars), max_width(_max_width), num_layers(_vars.getSize()+1)
{
    reserve_memory();
    if( create_path )
        construct_path();
}


/**
 * Allocate MDD memory
 */
void RevMDD::reserve_memory() {

    cout << "\t[MDD] Reserving memory..." << endl;

    cout << "\t\tArc buffers...\n";

    min_dom_val = vars.getMinMin();
    max_dom_val = vars.getMaxMax();

    arc_buffer = new( cp.getHeap() ) Arc***[num_layers];
    arc_layer  = new( cp.getHeap() ) Arc**[num_layers];
    max_arc_id = new( cp.getHeap() ) int[num_layers];

    for( int l = 0; l < num_layers; l++ ) {
    	int layer_id = 0;
        arc_buffer[l] = new( cp.getHeap() ) Arc**[max_width];
        for( int w = 0; w < max_width; w++ ) {
            arc_buffer[l][w] = new( cp.getHeap() ) Arc*[max_dom_val+1];
            for( int v = min_dom_val; v <= max_dom_val; v++ ) {
                arc_buffer[l][w][v] = new Arc(cp, layer_id++, v, w);
            }
        }
        max_arc_id[l] = layer_id;
    }


    for( int l = 0; l < num_layers; l++ ) {
    	int layer_id = 0;
    	arc_layer[l] = new( cp.getHeap() ) Arc*[max_arc_id[l]];
        for( int w = 0; w < max_width; w++ ) {
            for( int v = min_dom_val; v <= max_dom_val; v++ ) {
                arc_layer[l][layer_id] = arc_buffer[l][w][v];
                layer_id++;
            }
        }
    }


    cout << "\t\tArc layer list... " << endl;
    arcs_layer = new( cp.getHeap() ) RevLinkedList<Arc>[num_layers];
    for( int l = 0; l < num_layers; l++ ) {
        Arc* header = new Arc(cp, -1, -1, -1);
        arcs_layer[l].initialize(cp, header, num_layers * max_width * max_width);
    }


    // initialize MDD modification control
    incoming_modified = new( cp.getHeap() ) bool*[num_layers];
    incoming_modified_rev = new( cp.getHeap() ) IlcRevBool**[num_layers];

    outgoing_modified = new( cp.getHeap() ) bool*[num_layers];
    outgoing_modified_rev = new( cp.getHeap() ) IlcRevBool**[num_layers];


    for( int l = 0; l < num_layers; l++ ) {
    	incoming_modified_rev[l] = new( cp.getHeap() ) IlcRevBool*[max_width];
    	outgoing_modified_rev[l] = new( cp.getHeap() ) IlcRevBool*[max_width];

    	incoming_modified[l] = new( cp.getHeap() ) bool[max_width];
    	outgoing_modified[l] = new( cp.getHeap() ) bool[max_width];

    	// initial state: all nodes need update
    	for( int w = 0; w < max_width; w++ ) {
    		incoming_modified_rev[l][w] = new( cp.getHeap() ) IlcRevBool(cp, true);
    		outgoing_modified_rev[l][w] = new( cp.getHeap() ) IlcRevBool(cp, true);
    	}
    }

    // root only requires outgoing update
    incoming_modified_rev[0][0]->setValue(cp, false);

    // terminal only requires incoming update
    outgoing_modified_rev[num_layers-1][0]->setValue(cp, false);

    // other initializations
    aux_set.resize(0, max_dom_val, false);
    has_value = new( cp.getHeap() ) bool[max_dom_val+1];

//    cout << "\t[MDD] done." << endl;
}


/**
 * Construct initial path MDD
 */
void RevMDD::construct_path() {

    cout << "\t[MDD] Constructing path MDD...\n";

    // just add arcs corresponding to variables domains in active arc list
    Arc* previous;

    for( int l = 0; l < num_layers-1; l++ ) {

        previous = arcs_layer[l].get_header();

        // add each var domain value
        for( IlcIntExpIterator iter(vars[l]); iter.ok(); ++iter ) {
            previous = arcs_layer[l].extend(previous, get_arc(l, 0, *iter));
            previous->set_target_rev(cp, 0);
        }

        arcs_layer[l].close(previous);
    }

    cout << "\t[MDD] Exporting to graphs/initial_mdd.gml...\n";
    export_to_gml("graphs/initial_mdd.gml");

//    cout << "\t[MDD] done.\n";
}


/**
 * Update variable domains according to changes in arcs
 * TODO: only update layers that had some change
 */
void RevMDD::update_variable_domains() {

    for( int l = 0; l < num_layers-1; l++ ) {

        // get the union of arc values
        memset(has_value, false, sizeof(bool)*(max_dom_val+1));
        for( Arc* arc = arcs_layer[l].get_first(); arc != NULL; arc = arcs_layer[l].get_next(arc) ) {
            has_value[arc->val] = true;
        }

        // remove values that do not belong to arc domain
        for( IlcIntExpIterator val(vars[l]); val.ok(); ++val ) {
            if( !has_value[*val] ) {
            	//cout << "\t\t[MDD] Removing val " << *val << " from variable " << l << endl;
            	vars[l].removeValue(*val);
            }
        }
    }
}


/**
 * Export to GML file
 */
void RevMDD::export_to_gml(const char* filename) {

    ofstream mdd_file(filename);
    mdd_file << "graph [\n";
    //mdd_file << "comment \"automatically generated MDD graph \"\n";
    mdd_file << "\t directed 1\n";
    mdd_file << "\t hierarchic 1\n";


    // --- first create nodes for each layer ---

    // create remaining nodes
    bool** node_created = new bool*[num_layers];
    for( int i = 0; i < num_layers; i++ ) {
        node_created[i] = new bool[max_width];
        memset(node_created[i], false, sizeof(bool)*max_width);
    }


    for( int l = 0; l < num_layers-1; l++ ) {
        for( Arc* arc = arcs_layer[l].get_first(); arc != NULL; arc = arcs_layer[l].get_next(arc) ) {

        	// fix arc target
        	arc->target_node = arc->target_node_id.getValue();

            if( !node_created[l][arc->get_source()] ) {
                node_created[l][arc->get_source()] = true;
                mdd_file << "node [ \n";
                mdd_file << "\t id " << l*max_width + arc->get_source() << "\n";
                mdd_file << "\t label \"" << l << "," << arc->get_source() << "\"\n";

                mdd_file << "\t graphics [ \n";
                mdd_file << "\t\t type \"ellipse\"\n";
                mdd_file << "\t\t hasFill 0 \n";
                mdd_file << "\t\t ] \n";
                mdd_file << "\t ]\n";
            }

            if( !node_created[l+1][arc->get_target()] ) {
                node_created[l+1][arc->get_target()] = true;
                mdd_file << "node [ \n";
                mdd_file << "\t id " << (l+1)*max_width + arc->get_target() << "\n";
                mdd_file << "\t label \"" << (l+1) << "," << arc->get_target() << "\"\n";

                mdd_file << "\t graphics [ \n";
                mdd_file << "\t\t type \"ellipse\"\n";
                mdd_file << "\t\t hasFill 0 \n";
                mdd_file << "\t\t ] \n";
                mdd_file << "\t ]\n";
            }
        }
    }
    for( int i = 0; i < num_layers; i++ ) {
        delete[] node_created[i];
    }
    delete[] node_created;


    // --- now create arcs ---
    for( int l = 0; l < num_layers; l++ ) {
        for( Arc* arc = arcs_layer[l].get_first(); arc != NULL; arc = arcs_layer[l].get_next(arc) ) {
            mdd_file << "edge [ \n";
            mdd_file << "\t source " << l*max_width + arc->get_source() << "\n";
            mdd_file << "\t target " << (l+1)*max_width + arc->get_target() << "\n";
            mdd_file << "\t label \"" << arc->val << "\"\n";
            mdd_file << "\t ]\n";
        }
    }

//    // create arcs
//    for( IlcInt l = 0; l < num_layers; l++ ) {
//        //cout << "layer " << l << " - active_nodes: " << active_nodes[l] << endl;
//        Node* node = nodes_layer[l].get_first();
//        while( node != NULL ) {
//
//            ArcLL* arc = node->out_arcs.get_first();
//            while( arc != NULL ) {
//
//                mdd_file << "edge [ \n";
//                mdd_file << "\t source " << l*max_width + node->id << "\n";
//                mdd_file << "\t target " << (l+1)*max_width + arc->node_id << "\n";
//                mdd_file << "\t label \"" << *(arc->get_domain()) << "\"\n";
//                mdd_file << "\t ]\n";
//
//                arc = node->out_arcs.get_next(arc);
//            }
//
//            node = nodes_layer[l].get_next(node);
//        }
//    }
//
    mdd_file << "]"; // graph [
    mdd_file.close();
}
