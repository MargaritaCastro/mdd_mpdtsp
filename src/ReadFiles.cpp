#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>

using namespace std;


// --------------
//  ReadFile Class 1
// --------------

void ReadFileClass1(string fileName, int& NbCities, int& Capacity, int& NbCommodities, vector<int>& dd, vector<int>& commodity) {

	ifstream myfile;
	string words = (fileName);
	myfile.open(words.c_str());

	if (myfile.is_open()) { 
		// ignore words
		myfile >> words;
		while(words.compare("DIMENSION:") != 0 ){
      			myfile >> words;
    		}

		// read the number of cities
		myfile >> NbCities; // this include the 2 dummy cities
        
        // ignore words
		while(words.compare("EDGE_WEIGHT_SECTION:") != 0 ){
      			myfile >> words;
    	}        
		
		// read the matrix of distances
		dd.resize(NbCities*NbCities);
		
		for (int i = 0; i < NbCities; i ++) {
			for (int j = 0; j < NbCities; j ++) {
				myfile >> dd[i*NbCities + j];
			}
		}
		
		// change -1 in the distance matrix
		//for (int i = 0; i < NbCities; i ++) {
			//for (int j = 0; j < NbCities; j ++) {
				//if(dd[i*NbCities + j] < 0){
					//dd[i*NbCities + j] = dd[NbCities-1]*9;
				//}
			//}
		//}
		dd[NbCities-1] = -1;
		
		// read capacity and number of comodities
		myfile >> words >> Capacity;
		myfile >> words >> NbCommodities;
		myfile >> words;
		
		// read the commodity matrix
		commodity.resize(NbCities*NbCommodities);

		for (int i = 0; i < NbCities; i ++) {
			for (int j = 0; j < NbCommodities; j ++) {
				if (j == 0){
					myfile >> words;
				}
				myfile >> commodity[i*NbCommodities + j];
			}
		}
		myfile.close();
	}
	else
        cout << "Error: Couldn't open the file  of class 1 :" << fileName << endl;
}

// ----------------------
//  ReadFile Class 2 & 3
// ----------------------

void ReadFileClass23(string fileName, int& NbCities, int& Capacity, int& NbCommodities, vector<int>& dd, vector<int>& commodity) {

	ifstream myfile;
	string words = (fileName);
	myfile.open(words.c_str());
	vector<double> coord; //coordinates for each cities

	if (myfile.is_open()) { 
		// ignore words
		myfile >> words;
		while(words.compare("DIMENSION:") != 0 ){
      			myfile >> words;
    		}
		
		// read the number of cities, capacity and comodities
		myfile >> NbCities; // this include the 2 dummy cities
		myfile >> words >> Capacity;
		myfile >> words >> words >>  words >> NbCommodities;
		
		dd.resize(NbCities*NbCities);
		commodity.resize(NbCities*NbCommodities);

		double aux(0.0);
        
	    // read the matrix of coordinates
      	myfile >> words; 

		for (int i = 0; i < NbCities; i ++) {
			for (int j = 0; j < 2; j ++) {
				if (j == 0){
					myfile >> words;
				}
				myfile >> aux;
				coord.push_back(aux); 
			}
		}
		
		//generate distance matrix
		for (int i = 0; i < NbCities; i ++) {
			for (int j = 0; j < NbCities; j ++) {
				if (i == j){
					dd[i*NbCities + j] = 0;
				} else{
					aux = (coord[i*2 + 0] - coord[j*2 + 0])*(coord[i*2 + 0] - coord[j*2 + 0])  
							+ (coord[i*2 + 1] - coord[j*2 + 1])*(coord[i*2 + 1] - coord[j*2 + 1]) ;
					aux = sqrt(aux);
					dd[i*NbCities + j] = (int)round(aux);
				}
			}
		}
		
		// read the commodity matrix
		myfile >> words;
		for (int i = 0; i < NbCities; i ++) {
			for (int j = 0; j < NbCommodities; j ++) {
				if (j == 0){
					myfile >> words;
				}
				myfile >> commodity[i*NbCommodities + j];
			}
		}
		myfile.close();
	}
	else {
        cout << "Error: Couldn't open the file: " << fileName << endl;
	}
}

void ReadFileTSPPD(string fileName, int& NbCities, int& Capacity, int& NbCommodities, vector<int>& dd, vector<int>& commodity) {


	ifstream input;
	string path = (fileName);
	input.open(path.c_str());

	if( !input.is_open() ) {
		cerr << "Error: could not open file " << fileName << endl;
		return;
	}

	// num of activities
	input >> NbCities; 	// includes the depot only once
	NbCommodities =trunc((NbCities-1)/2);	//Number of commodities is always half of the number of cities
	NbCities++;
	Capacity = NbCommodities;
	
	dd.resize(NbCities*NbCities, 0);
	commodity.resize(NbCities*NbCommodities, 0);

	// read coordinates and precedences
	int x,y;
	int index, kind, rel;
	vector< pair<int,int> > coords;

	// depot coordinates
	input >> index;
	input >> x;
	input >> y;
	coords.push_back( pair<int,int>(x,y) );

	// requests
	int aux_item=0;
	for( int i = 0; i < NbCities -2; ++i ) {
		input >> index;
		input >> x; input >> y;
		input >> kind; input >> rel;
		if( kind == 0 ) {
			// delivery
			index = index - 1;  //delivery city
			rel = rel - 1;		//pickup city
			commodity[index*NbCommodities + aux_item] = -1;
			commodity[rel*NbCommodities + aux_item] = 1;
			aux_item++;
		}
		coords.push_back( pair<int,int>(x,y) );
	}

	// terminal depot
	coords.push_back( coords[0] );
	
	//generate distance matrix
	double aux(0.0);
	for (int i = 0; i < NbCities; i ++) {
		for (int j = 0; j < NbCities; j ++) {
			if (i == j){
				dd[i*NbCities + j] = 0;
			} else{
				aux = (coords[i].first-coords[j].first)*(coords[i].first-coords[j].first)
				+ (coords[i].second-coords[j].second)*(coords[i].second-coords[j].second);
				aux = sqrt(aux);
				dd[i*NbCities + j] = (int)round(aux);
			}
		}
	}

	input.close();
}


