// -------------------------------------------------------------- -*- C++ -*-
// File: lagrangian.cpp
//-------------
// Code for computting the larangian multipliyers
//-------------

#include <iostream>
#include <string>
#include <vector>
#include <algorithm> 						// for random_shuffle, min, max and sort
#include <ctime>							// measure time
#include <ilcplex/ilocplex.h>				// solve cplex
#include "constraints/disjunctive.hpp"		// For the MDD disjunctive propagator

using namespace std;
typedef vector<vector<double> > vector2_double;
typedef vector<vector<int> > vector2;

//Functions
void subGradientAllDiff( vector<double>& grad, int* sol);
void subGradientCapacity( vector<double>& grad, int* sol, int* weight, int C);
void subGradientPrecedence( vector<double>& grad, int* sol, vector2& precedence, int n);

void CostAllDiff(vector<double>& u_new, double* cost);
void CostCapacity(vector<double>& u_new, double* cost, int C);
void CostCapacity2(vector<double>& u_new, double* cost, int C);
void CostPrecedence (vector<double>& u_new, double* cost, vector2& precedence, int n);

double BundleMethod(vector<double>& theta, vector2_double& g, vector2_double& u, vector<double>& u_new, int UB, double t, int type, int nPre);
double CuttingMethod(vector<double>& theta, vector2_double& g, vector2_double& u, vector<double>& u_new, int UB, int type) ;


// Lagrangian Structure
struct Lag_sol {
	double LB_original;
	double LB_lagr;
	double iterations;
	double time;
	double solver_time;
	double sp_time;
};

struct Parameters {
	int* weight;				// city weight
	int C;						// Capacity
	vector2 precedence;			// precedence matrix
};

bool debug_l = false;

//===========================================================================================
// Main
//===========================================================================================

void computeLagrangian(DisjunctivePropagator* disj_prop, int n, int UB, double* cost_task, double* cost_layer, double* cost_pre, Lag_sol& lagr, int type, Parameters& param, bool debug) {

	
    debug_l = debug;
    
	// Parameters
	vector<double> theta(0);
	vector2_double g(0);
	int k = 0; 				// number of commodities
	int nMulti = n; 		// number of multipliers
	
	// See the type of the lagrangean:  0 = all diff; 1 = capacity; 2 = capacity from both sides; 3 = alldiff + capacity; 4 = precedence; 5 = alldiff + precedence
	if(type == 2) 		nMulti = 2*n;		// number of multipliers
	else if(type == 3) 	nMulti = 3*n;		// number of multipliers
	else if(type == 4){
		k = param.precedence.size();
		nMulti = k;
	}
	else if(type == 5){
		k = param.precedence.size();
		nMulti = n + k;
	}

	//First iteration
	vector<double> u_new(nMulti, 0.0);
	vector2_double multi(0);
	multi.push_back(u_new);
	int i = 0;
	double value = 0.0;
	double dual = 0.0;
	double diff = 1.0;
	double t = 1.0;
    
    //Auxilary vectors for computing gradiants
    vector<double> grad(nMulti);
    vector<double> g_aux(n);
    vector<double> g_aux2(k);

	while(true){
		//-- Compute shortest path
		double value_old = value;
		int* sequence = NULL;
		value = 0.0;

		clock_t start_shortest_path = clock();
		sequence = disj_prop->compute_shortest_path(value, cost_task, cost_layer, cost_pre);
		clock_t end_shortest_path = clock();
		
		lagr.sp_time += double(end_shortest_path - start_shortest_path) / CLOCKS_PER_SEC;
		theta.push_back(value);
		

		if(debug_l){
			cout << "\n #== Iteration " << i << "\n \n";
			cout << "Shortest path value: " << value << endl;
			cout << "Shortest path: ";
			for (int j = 0; j < n; j++){
				cout << sequence[j] <<" ";
			}
			cout << endl;
		}

		//-- Compute subgradient
        grad.assign(nMulti, 0);
		
		if (type == 0) {
			subGradientAllDiff(grad, sequence);
		}
		else if(type == 1){
			subGradientCapacity(grad, sequence, param.weight, param.C);
		}
		else if(type == 2){
			g_aux.assign(n, 0);
			subGradientCapacity(g_aux, sequence, param.weight, param.C); // for the <= C constraint
			for(int i = 0; i < n; i++) grad[i] = g_aux[i];
			
			subGradientCapacity(g_aux, sequence, param.weight, 0);	// for the >= 0 constraint
			for(int i = n; i < nMulti; i++) grad[i] = g_aux[i-n];
			
		}
		else if(type == 3){	// AllDiff + capacity constraint
			g_aux.assign(n, 0);
            
			subGradientAllDiff(g_aux, sequence);
			for(int j = 0; j < n; j++) grad[j] = g_aux[j];	//for the allDiff constraint
			
			subGradientCapacity(g_aux, sequence, param.weight, param.C); // for the <= C constraint
			for(int j = n; j < 2*n; j++) grad[j] = g_aux[j-n];
			
			subGradientCapacity(g_aux, sequence, param.weight, 0);	// for the >= 0 constraint
			for(int j = 2*n; j < nMulti; j++) grad[j] = g_aux[j-2*n];
		}
		else if (type == 4){	// Precedence Constraint 
			subGradientPrecedence(grad, sequence, param.precedence, n);
		}
		else if (type == 5){	// AllDiff + Precedence
            g_aux.assign(n, 0);
            g_aux2.assign(n, 0);
			
			subGradientAllDiff(g_aux, sequence);
			for(int j = 0; j < n; j++) grad[j] = g_aux[j];
			
			subGradientPrecedence(g_aux2, sequence, param.precedence, n);
			for(int j = 0; j < k; j++){
				grad[n + j] = g_aux2[j];
			}

		}
		g.push_back(grad);
		
		if(debug_l) {
			cout << "Final Subgradient: ";
			for(int j = 0; j < (int) grad.size(); j++){
				cout << grad[j] << " ";
			}
			cout << endl;
		}

		//-- Solve the cutting plane problem
		double dual_old = dual;
		
		clock_t start_cplex = clock();
		//double dual = CuttingMethod(theta, g, multi, u_new, UB, type);
		dual = BundleMethod(theta, g, multi, u_new, UB, t,  type, k);
		clock_t end_cplex = clock();
		
		lagr.solver_time += double(end_cplex - start_cplex) / CLOCKS_PER_SEC;
		
		multi.push_back(u_new);

		if(debug_l){
			cout << "Iteration = " << i;
			cout << "\tprimal = " << value;
			cout << "\tdual=" << dual << '\n';
            cout << "Lambda = ";
            for(int i=0; i < u_new.size(); i++){
                cout << u_new[i] << " ";
            }
            cout << '\n';
		}
		
		
		//-- Compute cost of each task
		if (type == 0) { 
			CostAllDiff(u_new, cost_task);
		} 
		else if (type == 1) {
			CostCapacity(u_new,  cost_layer, param.C); // one-side capacity
		}
		else if (type == 2){
			CostCapacity2(u_new,  cost_layer, param.C);	//two side capacity
		}
		else if (type == 3){
			vector<double> u_task(u_new.begin(), u_new.begin() + n);
			vector<double> u_layer(u_new.begin() + n, u_new.end());
			
			CostAllDiff(u_task, cost_task);
			CostCapacity2(u_layer,  cost_layer, param.C);
		}
		else if (type == 4){
			CostPrecedence(u_new, cost_pre, param.precedence, n);
		
		}
		else if (type == 5){
			vector<double> u_task(u_new.begin(), u_new.begin() + n);
			vector<double> u_pre(u_new.begin() + n, u_new.end());
		
			CostAllDiff(u_task, cost_task);
			CostPrecedence(u_pre, cost_pre, param.precedence, n);
		}

		//-- Stopping criteria
		if(i > 0){
			//if( abs(value - theta[i-1]) < 0.1){	// This condition is not working!. Is not giving me the optimal multipliers
			if(abs(value - dual) < 0.5 ){ 					//Stop if primal is equal to dual
				if(abs(value - value_old) < 0.01 || abs(dual - dual_old) < 0.01 ){			//Also compare the objective value to make sure it converge
					lagr.LB_lagr = value;	
					lagr.iterations = i;
					//cout<< "End of lagrangian. Primal= " << value << " and Dual = " << dual << "\n"; 
					break;
				}
			}
			//if(abs(value - dual) < 0.1 ) break;   
			if (i >= 400){
				lagr.LB_lagr = value;	
				lagr.iterations = i;
				cout << "Possible Problem with lagrangian = " << type + 1 << "\n";
				cout << "Iteration " << i <<": Primal= " << value << " and Dual = " << dual << "\n"; 
				break;
			}
			//-- Update step for the bundle method  --> it makes it much more faster
			double t_aux = 2*t*(1- (value - value_old)/diff);
			t = min(max(max(t_aux, t/10), 0.01),1.0);
			//cout << "Aux step = " << t_aux << endl;
			//cout << "New step = " << t << endl;
			
		}
		else{
			lagr.LB_original = value;
		}
		
		diff = dual - value;
		i++;
	}
}


//===========================================================================================
// STOPPING CRITERIA
//===========================================================================================

bool stopping_criteria( double primal, double dual, double gap, double dual_old){
	
	if( abs(primal - dual)/abs(primal) < gap || abs(dual- dual_old) < 0.1){ 
		return true;
	}
	return false;
}

//===========================================================================================
// SUBGRADIENT
//===========================================================================================

void subGradientAllDiff( vector<double>& grad, int* sol){ // each cell of the subgradiente represent the number of times a solution is choosen minus 1  

	for(int i = 0; i < (int) grad.size(); i++){
		grad[i] = -1;
	}

	for(int i = 0; i < (int) grad.size(); i++){
		grad[sol[i]] ++;
	}

	if(debug_l){
		cout << "Subgradient AllDiff: ";
		for(int i = 0; i < (int) grad.size(); i++){
			cout << grad[i] << " ";
		}
		cout << endl;
	}
}

void subGradientCapacity( vector<double>& grad, int* sol, int* weight, int C){


	//Adds the weight of the solution
	grad[0] = weight[sol[0]];
	for(int i = 1; i < (int) grad.size(); i++){
		grad[i] = grad[i-1] + weight[sol[i]];
	}
	
	// Add the - capacity term to everything
	for(int i = 0; i < (int) grad.size(); i++){
		grad[i] -= C;
	}
	
	if(debug_l){
		cout << "Subgradient (C = " << C <<") : ";
		for(int i = 0; i < (int) grad.size(); i++){
			cout << grad[i] << " ";
		}
		cout << endl;
	}

}

void subGradientPrecedence( vector<double>& grad, int* sol, vector2& precedence, int n) {
	
	// Compute gradiente by conunting the distances
	for(int k = 0; k < (int)grad.size(); k++){
		for(int i = 0; i < n; i++){
			if (sol[i] == precedence[k][0]) grad[k] += i;
			else if (sol[i] == precedence[k][1]) grad[k] -= i;
		}
	}
	
	//Add one to all the gradientes
	for(int i = 0; i < (int) grad.size(); i++){
		grad[i] += 1;
	}
	
	if(debug_l){
		cout << "Subgradient Precedence : ";
		for(int k = 0; k < (int) grad.size(); k++){
			cout << grad[k] << " ";
		}
		cout << endl;
	}
	
}


//===========================================================================================
// COMPUTE COST
//===========================================================================================

void CostAllDiff(vector<double>& u_new, double* cost){

	//Cost for each possible value
	for(int i = 0; i < (int) u_new.size(); i++) cost[i] = u_new[i];
	
	//Add constant at the root node
	for(int i = 0; i < (int) u_new.size(); i++) cost[0] -= u_new[i];

	//Print cost
    if(debug_l){
        cout << "Cost tour :";
        for(int i = 0; i < (int) u_new.size(); i++) cout << cost[i] <<" ";
	
        cout << endl;
    }
}

void CostCapacity(vector<double>& u_new, double* cost, int C){

	int n = (int) u_new.size();
	cost[n-1] = u_new[n-1];
	cost[0] = 0.0;
	
	//Additive cost
	for(int i = (n-2); i > 0; i--) cost[i] = cost[i+1] + u_new[i];
	
	//Cost at the root node
	for(int i = 0; i < n; i++) cost[0] -= u_new[i]*C;
	
	//Print cost
    if(debug_l){
        cout << "New cap(>=C) :";
        for(int i = 0; i < (int) u_new.size(); i++) cout << cost[i] <<" ";

        cout << endl;
    }
}

void CostCapacity2(vector<double>& u_new, double* cost, int C){

	int n = (int) u_new.size()/2;
	cost[n-1] = u_new[n-1] + u_new[(int) u_new.size() - 1];
	cost[0] = 0.0;
	
	//Additive cost
	for(int i = (n-2); i > 0; i--) cost[i] = cost[i+1] + u_new[i] + u_new[n+i] ;
	
	//Cost at the root node
	for(int i = 0; i < n; i++) cost[0] -= u_new[i]*C;

	//Print cost
    if(debug_l) {
        cout << "Cost cap :";
        for(int i = 0; i < n; i++) cout << cost[i] <<" ";
	
        cout << endl;
    }
}

void CostPrecedence (vector<double>& u_new, double* cost, vector2& precedence, int n ){

	int m = u_new.size(); // number of commodities
	
	for(int i = 0; i < n; i++)cost[i] = 0;

	for (int k = 0; k < m; k++){
		cost[precedence[k][0]] += u_new[k]; 	//pickup cost
		cost[precedence[k][1]] -= u_new[k]; 	//delivery cost
		cost[0] += u_new[k];
	}
	
	//Print cost
    if(debug_l){
        cout << "Cost pre :";
        for(int i = 0; i < n; i++) cout << cost[i] <<" ";
	
        cout << endl;
    }

}

//===========================================================================================
// SOLVING METHODS
//===========================================================================================

double BundleMethod(vector<double>& theta, vector2_double& g, vector2_double& u, vector<double>& u_new, int UB, double t, int type = 0, int nPre = 0) {

	IloEnv env;
	double r(0.0);
	
	try{
	
		//Initial parameters
		IloModel model(env);
		
		int K = g.size();
		int m = u[0].size();  //number of variables
		char varName[100];  
	
		//Define variable		:	z \forall i = 1, ..., m
		IloNumVarArray z(env, m, -IloInfinity, IloInfinity);
		sprintf(varName, "w"); 
		IloNumVar w(env, -IloInfinity, IloInfinity, varName);  // for the objective value
		
		for (int i = 0; i < m; i++){
			sprintf(varName, "z.(%d)", i); 
         	z[i].setName(varName);
         	if (type == 1) z[i].setLB(0.0);
         	else if (type == 2){
         		if(i < m/2) z[i].setLB(0.0);
         		else		z[i].setUB(0.0);
         	}
         	else if (type == 3){
         		if(i >= m/3 && i < 2*m/3)	z[i].setLB(0.0);
         		else if (i >= 2*m/3 )		z[i].setUB(0.0);
         	}
         	else if (type == 4) z[i].setLB(0.0);
         	else if (type ==5){
         		if(i >= m - nPre)	z[i].setLB(0.0);
         	}
		}
		
		//Create constraints
		for	(int k = 0; k < K; k++){
			IloExpr constr(env);
			double b(0.0);
			
			constr += w;
			b += theta[k];
			
			for(int i = 0; i < m; i++){
				constr -= g[k][i]*z[i];
				b -= g[k][i]*u[k][i];
			}
			
			model.add( constr <= b);
		}
		
		//Cuatradic term from the bundle method
		sprintf(varName, "obj_var"); 
		IloNumVar obj_var(env, -IloInfinity, UB, varName); 
		
		IloExpr obj(env);
		for(int i = 0; i < m; i++){
			obj -= t*(z[i] - u[K-1][i]) * (z[i]- u[K-1][i])/2;
		}
		model.add(obj_var <= w + obj);
		
		//Objective Function
		model.add( IloMaximize(env, obj_var) );
		
		//Solve IP problem
		IloCplex cplex(model); 
		cplex.setOut(env.getNullStream());
		cplex.setParam(IloCplex::Threads, 1);
		//cplex.exportModel("cutting.lp");
		
		if (cplex.solve()) {
			r = cplex.getObjValue();
			//cout << "Optimal r = " << r << endl;
			//cout << "u_new values: ";
			for(int i = 0; i < m; i++){
				if (cplex.isExtracted(z[i]) ) {
					u_new[i] = cplex.getValue(z[i]);
				}
				else u_new[i] = 0.0;

				//cout << u_new[i] << " ";
			}
			//cout << endl;
		} else{
			r = INF;
			cout << "No solution found" << endl;
		}
	
	}catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	
	env.end();
	
	return r;
}

double CuttingMethod(vector<double>& theta, vector2_double& g, vector2_double& u, vector<double>& u_new, int UB, int type = 0) {

	IloEnv env;
	double r(0.0);
	
	try{
	
		//Initial parameters
		IloModel model(env);
		
		int K = g.size();
		int m = u[0].size();  //number of variables
		char varName[100];  
	
		//Define variable		:	z \forall i = 1, ..., m
		IloNumVarArray z(env, m, -IloInfinity, IloInfinity);
		sprintf(varName, "w"); 
		IloNumVar w(env, -IloInfinity, UB, varName);  // for the objective value
		
		for (int i = 0; i < m; i++){
			sprintf(varName, "z.(%d)", i); 
         	z[i].setName(varName);
         	if (type == 1) z[i].setLB(0.0);
         	else if (type == 2){
         		if(i < m/2) z[i].setLB(0.0);
         		else		z[i].setUB(0.0);
         	}
         	else if (type == 3){
         		if(i >= m/3 && i < 2*m/3)	z[i].setLB(0.0);
         		else if (i >= 2*m/3 )		z[i].setUB(0.0);
         	}
		}
		
		//Create constraints
		for	(int k = 0; k < K; k++){
			IloExpr constr(env);
			double b(0.0);
			
			constr += w;
			b += theta[k];
			
			for(int i = 0; i < m; i++){
				constr -= g[k][i]*z[i];
				b -= g[k][i]*u[k][i];
			}
			
			model.add( constr <= b);
		}
		
		//Objective Function
		model.add( IloMaximize(env, w) );
		
		//Solve IP problem
		IloCplex cplex(model); 
		cplex.setOut(env.getNullStream());
		cplex.setParam(IloCplex::Threads, 1);
		//cplex.exportModel("cutting.lp");
		
		if (cplex.solve()) {
			r = cplex.getObjValue();
			//cout << "Optimal r = " << r << endl;
			cout << "u_new values: ";
			for(int i = 0; i < m; i++){
				if (cplex.isExtracted(z[i]) ) {
					u_new[i] = cplex.getValue(z[i]);
				}
				else u_new[i] = 0.0;

				cout << u_new[i] << " ";
			}
			cout << endl;
		} else{
			cout << "No solution found" << endl;
		}
	
	}catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	
	env.end();
	
	return r;
}

