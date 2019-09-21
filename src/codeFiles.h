#ifndef codeFiles_H
#define codeFiles_H

#include <string>
#include <vector>
#include "constraints/disjunctive.hpp"

using std::vector;
using std::string;
typedef vector<vector<int> > vector2;
typedef vector<vector<double> > vector2_double;

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
	int* weight;			//city weight
	int C;					// Capacity
	vector2 precedence;		// precedence matrix
};

//Lagrangian Multipliyers
void subGradientAllDiff( vector<double>& grad, int* sol);
void subGradientCapacity( vector<double>& grad, int* sol, int* weight, int C);
void subGradientPrecedence( vector<double>& grad, int* sol, vector2& precedence, int n);
void CostPrecedence (vector<double>& u_new, double* cost, vector2& precedence, int n );
void CostAllDiff(vector<double>& u_new, double* cost);
void CostCapacity(vector<double>& u_new, double* cost, int C);
void CostCapacity2(vector<double>& u_new, double* cost, int C);
double BundleMethod(vector<double>& theta, vector2_double& g, vector2_double& u, vector<double>& u_new, int UB, double t, int type, int nPre);
double CuttingMethod(vector<double>& theta, vector2_double& g, vector2_double& u, vector<double>& u_new, int UB, int type);
void computeLagrangian(DisjunctivePropagator* disj_prop, int n, int UB, double* cost_task, double* cost_layer, double* cost_pre, Lag_sol& lagr, int type, Parameters& param, bool debug);

//extra Lagragian multipliers functions
bool stopping_criteria( double primal, double dual, double gap, double dual_old);

//Read Files
void ReadFileClass1(string fileName, int& NbCities, int& Capacity, int& NbCommodities, vector<int>& dd, vector<int>& commodity);
void ReadFileClass23(string fileName, int& NbCities, int& Capacity, int& NbCommodities, vector<int>& dd, vector<int>& commodity);
void ReadFileTSPPD(string fileName, int& NbCities, int& Capacity, int& NbCommodities, vector<int>& dd, vector<int>& commodity);

// Checking feasibility
bool CheckObjective(vector<int> sequence, vector<int> dist, int value);
bool CheckAlldiff(vector<int> sequence);
bool CheckPrecedence(vector<int> sequence, vector2 pairs);
bool CheckCapacity(vector<int> sequence, vector<int> weight, int C);


#endif
