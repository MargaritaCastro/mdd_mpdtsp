// -------------------------------------------------------------- -*- C++ -*-
// File: Check.cpp
//-------------
// Code that check if the solution computed by the solver is actualy a feasible solution
//-------------

#include <iostream>
#include <vector>

using namespace std;
typedef vector<vector<int> > vector2;

bool CheckObjective(vector<int> sequence, vector<int> dist, int value ){

	int obj = 0;
	int n = sequence.size();
	
	for (size_t i = 0; i < sequence.size()-1; i++){
		obj += dist[sequence[i]*n + sequence[i+1]]; 
	}


	return true;
}

bool CheckAlldiff(vector<int> sequence){

	for(size_t i = 0; i < (sequence.size()- 1); i++){
		for(size_t j = i+1; j < sequence.size(); j++){
			if(sequence[i] == sequence[j]){
				cout << "The solution violates the alldifferent constraint.\n";
				cout << "The vehicle visits city " << sequence[i] << " twice \n"; 
				return false; 
			}
		}
	}

	return true;
}

bool CheckPrecedence(vector<int> sequence, vector2 pairs){
	
	vector<int> inverse(sequence.size());
	
	for(size_t i = 0; i < sequence.size(); i++){
		inverse[sequence[i]] = i;
	}
	
	for(size_t k = 0; k < pairs.size(); k++){
		if( inverse[pairs[k][0]] >  inverse[pairs[k][1]] ){
			cout << "The solution violates a precedence constraint.\n";
			cout << "The vehicle visits city " << pairs[k][1] << " before visiting city " << pairs[k][0] << "\n"; 
			return false; 
		} 
	}


	return true;
}

bool CheckCapacity(vector<int> sequence, vector<int> weight, int C){

	int sum_cap(0);
	
	for (size_t i = 0; i < sequence.size(); i++){
		sum_cap += weight[sequence[i]];
		
		if (sum_cap > C || sum_cap < 0) {
			cout << "The solution violates the capacity constraint.\n";
			cout << "Vehicle's weight at city " << sequence[i] << " is: " << sum_cap << "\n";
		 	return false;
		}
	}


	return true;
}
