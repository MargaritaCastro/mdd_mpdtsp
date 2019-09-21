// -------------------------------------------------------------- -*- C++ -*-
// File: mpDTSP_MDD.cpp
// Input parameters: 
//		(1) Class Number
//		(2) Max width
//		(3) File name
//		(4) DFS where  0 == false; 1 == true
//		(5) Lagrangian 0 == nothing; 1 == allDiff cost; 3 == capacity cost; 4 == cap+AllDiff; 5 == precedence; 6 == Prec+AllDiff 
//		(6) Cap_refinement 0 == false; 1 == true
//		(7) second_MDD 0 == false; 1 == true
//		(8) Name of output file
//-------------
// CP model with MDDs for the one-to-one multicommodity pickup and delivery TSP
//-------------

#include <ilcp/cp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> 		//for random_shuffle, min, max and sort
#include <ctime>			//measure time

// MDD includes
#include "core/mdd_handler.hpp"
#include "core/search_goals.hpp"
#include "constraints/disjunctive.hpp"
#include "util/stats.hpp"
#include "util/util.hpp"

// My includes
#include "codeFiles.h"

/**
// -----------------------------------------------
// Comment out this block to use Frangioni's Lagrangian solver
//Lagragian include - Frangioni's code
#include "lagrangian/mddFiOrcl.hpp"
#include "lagrangian/SubGrad.h"
#include "lagrangian/PrimalDual.h"

//Bundle method includes - Frangioni's code
#include "bunddle/Bundle.h"
#include "bunddle/QPPnltMP.h"
// -----------------------------------------------
**/

using namespace std;

// macro parameters
#define TIME_LIMIT 7200 //time limit. 2 hours = 7200 sec. 30 min = 1800 sec.
#define WRITE 0 		// to decide if I write individual solutions for the problems or not

//TypeDef
typedef vector<vector<int> > vector2;
typedef vector<vector<double> > vector2_double;

// Solution structure
struct Solution {
	int 	value;			// value obtain by the solver
	double 	time;			// time of computation
	int 	branches;		// number of branches
	int 	fails;			// number of fails
	string	status;			// status of the solution ( eg. Feasible, Optimal, no Solution, Infeasible)
};

// For random 
int myrandom2 (int i) { return std::rand()%i;}
int rnd_seed = 0;

//Parameters of the model
bool 	DFS = false;
int 	MDDwidth = 0; 		// MDD's width. 0 = no MDD use
int 	classNumber;		// Set of problem instances. Can be 1, 2 or 3
string 	fileName;			// name of the file
int 	lagrangian = 0;
bool	cap_refinement = false;

Lag_sol lagr;

//Lagragian cost functions
double* cost_task;
double* cost_layer;
double* cost_pre;

bool f_buddle = false;  //Frangioni buddle method implementation
bool debug_lagrangian = false;
bool debug = false;


// Functions
int model (Solution& solution);
void PrintParameters( IloTransitionDistance distance, int n, IloIntArray2 pairs);
void PrintSolution(Solution& solution, int NbCities, int NbCommodities, int capacity, vector<int> stime = vector<int>(), vector<int> sequence = vector<int>());
void CheckSolution(vector<int> sequence, IloIntArray AbsWeight, IloIntArray2 pairs, int C, vector<int> dd, int obj);


//
// Shortest Path Search Strategy
//
ILCGOAL3(ShortestPathDD, int, pos, IlcIntVarArray, permutation, DisjunctivePropagator*, prop) {
    
	while (pos < permutation.getSize() && permutation[pos].isFixed()) {
        pos++;
    }

    // if last activity, fix remaining variables
    if (pos == permutation.getSize()) {
        //cout << "\tdone!" << endl;
        //return IlcGoalTrue(getCP());
        return IlcSimpleCompletionGoal(getCP());
    }

    int value = -1;
    double value_lagr = -1;
    int city = -1;
    int* shortest_path = prop->get_shortest_path(value);
    
    city = shortest_path[pos];

    if( prop->get_lagr_type() != 0){  // try the lagrangian shortest path
        int* shortest_path_lagr = prop->get_shortest_path_lagr(value_lagr);
        
        //cout << "original = " <<  value << "\t lagr = " << value_lagr << endl;
        if ( (int)value_lagr > value) city = shortest_path_lagr[pos];

    }
    
    return  IlcOr( IlcAnd(permutation[pos] == city,
                     ShortestPathDD(getCP(), pos+1, permutation, prop)),
              IlcAnd(permutation[pos] != city,
                     ShortestPathDD(getCP(), pos, permutation, prop)) );
}


//
// Wrapper: Sequencing based on activity belonging to shortest path in the MDD
//
ILOCPGOALWRAPPER2(IloShortestPathDD, cp, IloIntVarArray, permutation, 
		  DisjunctivePropagator*, params) 
{
  return ShortestPathDD(cp, 0, cp.getIntVarArray(permutation), params);
}


//-------------
// Main : read all input files 
//-------------
int main( int argc, char* argv[] ){ // parameters: classNumb, MDDWidth, fileName, DFS
	
	//Read input parameters
    classNumber = 1;
	fileName = argv[1];		//
	string outputFile = "";
	
	bool lagr_cap = false;
	bool lagr_pre = false;
	bool lagr_tour = false;

	//-----------------------------------------------------------------------------------
	//Read input paramertes
	int argcount = 1; 
	
	while(argcount < argc) {
		if(argv[argcount][0] == '-') {
			switch(argv[argcount][1]) {
                case 'c':
                {
                    classNumber = atoi(&(argv[argcount][2])); // 0 = TSPPD
                    
                    if (classNumber < 0 || classNumber >= 4){
                        cout << "Error - The code doesn't support class number " << classNumber << endl;
                        exit(1);
                    }
                    
                    break;
                }
				case 'd': // DFS
				{
					DFS = true;
					break;	
				}
                case 'f':
                {
                    f_buddle = true;        // Frangianni's buddle method implementation
                    break;
                }
				case 'l':  //lagrangian type
				{
					if(argv[argcount][2] == 'c'){
						lagr_cap = true;
						cout << "Using Lagrangian relaxation with capacity constraint" << endl;
					}
					else if(argv[argcount][2] == 't'){
						lagr_tour = true;
						cout << "Using Lagrangian relaxation with tour constraint" << endl;
					}
					else if(argv[argcount][2] == 'p'){
						lagr_pre = true;
						cout << "Using Lagrangian relaxation with precedence constraint" << endl;
					}
					break;
				} 
				case 'o':  //output file
				{
					outputFile = &(argv[argcount][2]);
					break;
				}
				case 'r': // refinement strategy
				{
					if(argv[argcount][2] == 'c'){
						cap_refinement = true;  //use the capacity refienement first
					}
					else {
						cap_refinement = false; // use only the the tour refinement
					}
					break;
				}
                case 'v': //flag for debugging
                {
                    debug = true;
                    
                    if(argv[argcount][2] == 'l'){
                        debug_lagrangian = true; // debug lagrangian code
                    }
                    
                    break;
                }
				case 'w':	// MDD width ... w0 --> CP model
				{
					if(argv[argcount][2] > 1)	MDDwidth = atoi(&(argv[argcount][2]));
               	 	else 						MDDwidth = 0;

					break;	
				}
				default:
				{ 
					cout << "No new parameters" << endl;
					break;
				}
			};
		
		}
		
		++argcount;
	}
	
	//Setting lagrangian parameters
	if 		( lagr_tour && !lagr_cap && !lagr_pre) lagrangian = 1;
	else if (!lagr_tour &&  lagr_cap && !lagr_pre) lagrangian = 3;
	else if ( lagr_tour &&  lagr_cap && !lagr_pre) lagrangian = 4;
	else if (!lagr_tour && !lagr_cap &&  lagr_pre) lagrangian = 5;
	else if ( lagr_tour && !lagr_cap &&  lagr_pre) lagrangian = 6;
	else if (!lagr_tour && !lagr_cap && !lagr_pre) lagrangian = 0;	// no lagrangian relaxation
	else{
		cout << "The lagragian combination is not currently supported" << endl;
		exit(1);
	}
	
	//Write running parameters
	cout << "\nRunning file: " << fileName  << endl;
	cout << "class number: "<< classNumber << endl;
	cout << "\n[Parameters]\n"; 
	if(MDDwidth <= 0){
		cout << "\tCP model\n";
	}else{
		cout << "\tMDD model\n";
		cout << "\twidth: " << MDDwidth << endl;
		cout << "\tLagrangian: " << lagrangian << endl;
		cout << "\tCapcity Refinement: " << cap_refinement << endl;
	}
	cout << "\tDFS: " << DFS << "\n\n";
	
	
	//-----------------------------------------------------------------------------------
	// Run model
	
	//Initialize lagrangian solution
	lagr.LB_original = 0.0;
	lagr.LB_lagr = 0.0;
	lagr.iterations = 0.0;
	lagr.time = 0.0;
	lagr.solver_time = 0.0;
	lagr.sp_time = 0.0;
	
	Solution sol;
    model(sol);
    
    //-----------------------------------------------------------------------------------
    //Write Solution summary of solution
    
    if( outputFile.size() >0 ){
    	fstream outputMaster;
    	
    	outputMaster.open(outputFile.c_str(), fstream::in | fstream::out | fstream::app);	
    	//outputMaster << "classNumb	fileName	time	value	fails	branches	status	MDDWidth" << endl;
    	outputMaster << classNumber			<< "	" 
    			 << fileName 			<< "	" 
    			 << sol.time 			<< "	" 
    			 << sol.value 			<< "	"
    			 << sol.fails 			<< "	"
    			 << sol.branches		<< "	" 
    			 << sol.status 			<< " 	"
    			 << MDDwidth 			<< " 	" 
    			 << lagrangian 			<< " 	"
    			 << lagr.iterations		<< "	" 
    			 << lagr.time			<< "	"
    			 << lagr.sp_time		<< "	"
    			 << lagr.solver_time	<< "	"
    			 << lagr.LB_original	<< " 	" 
    			 << lagr.LB_lagr		<< " 	" 
    			 << cap_refinement		<< endl;

		outputMaster.close();
    	
    }

	return 0;
}

//-------------
// CP model - with MDDs
//-------------

int model (Solution& solution) {
	// declare environment
	IloEnv env;
	try {
		// declare model	
		IloModel model(env);
		
		// Declare parameters
		int NbCities(0); 			//two dummy cities (depot)
		int NbCommodities(0); 		//number of comodities
		int capacity(0); 			// max capacity of the vehicle
		vector<int> dd(1); 			// distance matrix
		vector<int> comMatrix(1); 	// commodity matrix
		char varName[100];
		
		// Read File
		if (classNumber == 1){ //read files of class 1
			ReadFileClass1(fileName, NbCities, capacity, NbCommodities, dd, comMatrix);
		}
		else if (classNumber == 2 || classNumber == 3) { //read files of class 2 and 3
			ReadFileClass23(fileName, NbCities, capacity, NbCommodities, dd, comMatrix);
		}else{
			ReadFileTSPPD(fileName, NbCities, capacity, NbCommodities, dd, comMatrix);
		}

    	//Declare auxiliar parameters
		IloIntArray2 pairs(env, NbCommodities);	// pairs of source and sink for each commodity
		IloIntArray type(env, NbCities); 		// city type (auxiliar)
		IloIntArray AbsWeight(env, NbCities); 	// absolute weight in each city
		vector<int> weight(NbCommodities); 		// weight of each commodity
		vector<int> pickup(NbCities);			// number of pickups in each city
		vector<int> delivery(NbCities);			// number of deliveries in each city
		
		//Determinate weight in each node/city
		for (IloInt i = 0; i < NbCities; i++) {
			AbsWeight[i] = 0;
			for (IloInt j = 0; j < NbCommodities; j++) {
				AbsWeight[i] += comMatrix[i*NbCommodities + j];
				
				if (comMatrix[i*NbCommodities + j] > 0 ){
					 pickup[i] ++;
					 weight[j] = comMatrix[i*NbCommodities + j];
				}
				if (comMatrix[i*NbCommodities + j] < 0 ) delivery[i] += 1;
			}	
		}
		
		//Create matrix of pairs source-sink
		for (IloInt i = 0; i < NbCommodities; i++){
			pairs[i] = IloIntArray(env,2);
			for (IloInt j = 0; j < NbCities; j++) {
				if (comMatrix[i + j*NbCommodities] > 0){ 		//pickup node
					pairs[i][0] = j; 
				}
				else if (comMatrix[i + j*NbCommodities] < 0){	//delivery node
					pairs[i][1] = j;
				}
			}
		}
		
		//Create the distance matrix
		int dist_inf = 9999999;
		IloTransitionDistance distance(env, NbCities);
		for (IloInt i = 0; i < NbCities; ++i) {
			for (IloInt j = 0; j < NbCities; j++) {
				if (dd[i*NbCities + j] < 0){ //this only apply to class 1
					distance.setValue(i , j, dist_inf);
				}
				else {
					distance.setValue(i , j, dd[i*NbCities + j]);
				}
			}
		}
		
		//PrintParameters(distance, NbCities, pairs);

		//Declare variables
		IloIntervalVarArray city(env, NbCities);  				// interval variable representing each city
		
		//Create uperbounds
    	long int ub_obj = 0;
    	long int aux_max = 0;
    	for (int i = 0; i < NbCities; ++i) {
    		aux_max = 0;
      		for (int j = 0; j < NbCities; ++j) {
				if(aux_max < distance.getValue(i,j) && distance.getValue(i,j) != dist_inf) aux_max = distance.getValue(i,j);
      		}
      		ub_obj += aux_max;
   		}
   		
   		//Create objective variable
    	IloIntVar obj_var(env, 0, ub_obj, "ObjectiveVar");

    	//Create the city variable
		for (IloInt i = 0; i < NbCities; i++) {
			city[i] = IloIntervalVar(env, 1); 	//one time unit in each city
			type[i] = i; 						//all cities with a different type
			sprintf(varName, "city.(%d)", (int) i); 
         	city[i].setName(varName);
         	city[i].setEndMax(ub_obj + NbCities);	
		}
		
		//Create sequence variable
		IloIntervalSequenceVar tour(env, city, type); //sequence of cities
		
		//Create cumulative function
		IloCumulFunctionExpr currentCapacity (env);  			// current capacity of the vehicle
		
		// ---------------------------------------------------------------
		// CP constraints
		if(MDDwidth == 0){
			// Constraint: pickup before delivery
			for (IloInt j = 0; j < NbCommodities; j++){
 		  		model.add(IloBefore(env, tour, city[pairs[j][0]],   city[pairs[j][1]]));	
			}

			//Constraint: precedence for depot
			model.add(IloFirst(env, tour, city[0]));
			model.add(IloStartOf(city[0]) == 0);
			model.add(IloLast(env, tour, city[NbCities-1]));	
		
    		//Constraint: max capcity in each city
			for (IloInt i = 0; i < NbCities; i++) {
				if (AbsWeight[i] < 0) //to deal with negative weights
					currentCapacity -= IloStepAtStart(city[i], IloAbs(AbsWeight[i]));
				else if (AbsWeight[i] > 0)
					currentCapacity += IloStepAtStart(city[i], AbsWeight[i]);
			}	
			model.add(currentCapacity <= capacity);

    		//Constraint: NoOverlap with setup times
			model.add(IloNoOverlap(env, tour, distance, IloTrue));
        }
		
    	//Constraint; Relate objective values with the variable
		model.add(obj_var == (IloEndOf(city[NbCities -1]) - NbCities));

    	// ---------------------------------------------------------------
    	// MDD IMPLEMENTATION
	
		DisjunctivePropagator* disj_prop = NULL;
		IloIntVarArray permutation;
		
		if (MDDwidth > 0) {  //width = 0 --> run only the CPO
	
    		//create MDD activities
    		ActivityArray acts(NbCities);
    		for (int i = 0; i < NbCities; ++i) {
     			acts[i] = new(env) Activity(env,                // ILOG environment
					city[i],       		// IloIntervalVar associated with activity
					i,                  // index of activity in IloIntervalVarArray
					0,                  // release date
					ub_obj + NbCities,  // deadline (just an upper bound here)
					1                   // processing time (I assume "1" here)  
			    	);
    		}
    		
    		// create permutation variables
    		permutation = IloIntVarArray(env, NbCities);

    		permutation[0] = IloIntVar(env, 0, 0, "firstPerm");  // firt city
    		permutation[NbCities-1] = IloIntVar(env, NbCities-1, NbCities-1, "lastPerm"); //last city
    	
    		// remaining cities
    		for (int i = 1; i < NbCities-1; ++i) {
      			sprintf(varName, "permutation[%d]", i);
      			permutation[i] = IloIntVar(env, 1, NbCities-2, varName);
    		}

			// cities are pairwise distinct
    		//model.add( IloAllDiff(env, permutation) );
			
    		// create setup time matrix
    		int** setup = new(env) int*[NbCities];
    		for (int i = 0; i < NbCities; ++i) {
      			setup[i] = new(env) int[NbCities];
      			for (int j = 0; j < NbCities; ++j) {
					setup[i][j] = distance.getValue(i,j);
      			}
    		}
    
    		// create MDD disjunctive propagator
			disj_prop = new DisjunctivePropagator( acts,				// set of activities
				 tour,				// sequence var for disjunctive
				 permutation,		// permutation variables
				 setup,				// setup matrix
				 SumSetupTimes,		// type of objective function
				 obj_var,			// objective function variable
				 AbsWeight,			// total weight in each city
				 capacity			// maximum capacity of the vehicle
				 ); 
				 
			// array with all MDD propagators (in this case we only have one, the disjunctive)
			IloArray<MDDPropagator*> propagators(env);
    		propagators.add( disj_prop );

  			//Add commodity precedence constraints to the permutation variable
			for (IloInt i = 0; i < NbCommodities; i++){
      			disj_prop->set_precedence(pairs[i][0], pairs[i][1]);
			}
			
			//Set capacity refinement
			disj_prop->set_cap_ref(cap_refinement);
			
    		model.add( IloMDDHandler(env,          
			     permutation,		// permutation variables
			     city,				// set of interval vars
			     tour,				// IloIntervalSequence var
			     propagators,		// all propagators
			     MDDwidth,			// MDD width
			     "mdd-handler") );
		}
		// ---------------------------------------------------------------
		
		
    	//Add the objective function
		model.add(IloMinimize(env, obj_var));

    	//Create an instance of IloCP
		IloCP cp(model);
    	
    	if(MDDwidth > 0) cp.setParameter(IloCP::DefaultInferenceLevel, IloCP::Extended); //MDD parameter --> increace inference in each node
		
		/**
		//---------------------------------------------------------
        // Comment out this block to use Frangioni's Lagrangian solver
		// Lagrangian approach - using Frangioni's implementation bundle method implementation
		if (MDDwidth > 0 && lagrangian > 0 && f_buddle) {
            cp.propagate();
            disj_prop->set_lagr_type(lagrangian);
			
			cout << " == Lagrangian implementation == " << endl;
		
            //Create and initialize Fi object
            mddFiOrcl* Fi = new mddFiOrcl( disj_prop, NbCities, NbCommodities, lagrangian -1, debug_lagrangian);
            Fi->set_capacity(capacity);
            Fi->set_precedences(pairs);
            Fi->set_weight(AbsWeight);
            
            //Set upper bound of the function
            double LB = -ub_obj;
            Fi->SetLowerBound(LB);
			
            //--------------------------------------------------------
            //Bunddle method code
            const char *const ParFileBunddle = "ParValue.qp";
            
            ifstream ParFileB( ParFileBunddle );
            if( ! ParFileB.is_open() )cerr << "Warning: cannot open parameters file """ << ParFileBunddle << """" << endl;
            
            Bundle *s = new Bundle( &ParFileB );
            QPPenaltyMP *MP = new QPPenaltyMP( &ParFileB );
            
            Fi->SetFiTime();
            s->SetNDOTime();
            
            s->SetMPSolver( MP );
            s->SetFiOracle( Fi );
            Fi->SetNDOSolver( s );
                
            int lvl;
            DfltdSfInpt( &ParFileB , lvl , int( 0 ) );
            if( lvl ) s->SetNDOLog( &clog , char( lvl ) );
            
            //Solve the lagrangian dual problem
            Bundle::NDOStatus lagr_status;
            lagr_status = s->Solve();

            //Collect information about  the solution
            lagr.LB_lagr        = -s->ReadBestFiVal();
            lagr.LB_original    = Fi->GetOriginalBound();
            lagr.iterations     = s->NrIter();
                
            lagr.sp_time    = Fi->FiTime();
            lagr.time       = s->NDOTime();
            lagr.solver_time = lagr.time - lagr.sp_time;
                
            Index D;
            cIndex_Set I;
            cLMRow Lambda =  s->ReadBestSol(I, D);  //get set of optimal lambdas
            Fi->SetMDDMultipliers(Lambda);
                
            if (debug_lagrangian) {
                if( s->IsOptimal() ) cout << "\n(optimal value)" << endl;
                else cout << "\n(not provably optimal)" << endl;
                    
                cout << "Lambda* =" << endl;
                if( I ){
                    for( Index i = 0 ; ( i = *(I++) ) < Inf<Index>() ; )
                        cout << "[" << i << "]\t" << *(Lambda++) << endl;
                }
                else {
                    for( Index i = 0 ; i < s->GetNumVar() ; i++ )
                        cout << "[" << i << "]\t" << *(Lambda++) << endl;
                }
            }
            
            delete s;
            delete MP;
		}
        //---------------------------------------------------------
        **/
        
        // Lagrangian approach - our buddle method implementation
		if (MDDwidth > 0 && !f_buddle) {
			cost_task 	= new double[NbCities];
			cost_layer 	= new double[NbCities];
			cost_pre 	= new double[NbCities];
			double* cost_layer2 = new double[NbCities];
            
			for( int i = 0; i < NbCities; i++){
				cost_task[i]	= 0.0;
				cost_layer[i]	= 0.0;
				cost_pre[i] 	= 0.0;
				cost_layer2[i] 	= 0.0;
			}

			if (lagrangian > 0) {
				cp.propagate();
				
				disj_prop->set_cost_task(cost_task, false);
				disj_prop->set_cost_layer(cost_layer, false);
				disj_prop->set_cost_pre(cost_pre, false);
				disj_prop->set_lagr_type(lagrangian);
			
				// Parameters
				Parameters param;
				param.weight = new int[NbCities];
				param.C = capacity;

				// --- Inlcudes Capacity
				if( lagrangian >= 2){
					for(int i = 0; i < NbCities; i++){
						param.weight[i] = AbsWeight[i];
					}
				}
				// --- Include Precedence Constraints
				if (lagrangian >= 5){
					vector<int> aux_pre(2);
					
					for(int k = 0; k < NbCommodities; k++){
						aux_pre[0] = pairs[k][0];
						aux_pre[1] = pairs[k][1];
						
						param.precedence.push_back(aux_pre);
					}
				}
			
				// Compute optimal lagrangian multipliers at the root node
				clock_t start_lagr = clock();
                
                computeLagrangian(disj_prop, NbCities, ub_obj, cost_task, cost_layer, cost_pre, lagr, lagrangian - 1, param, debug_lagrangian);
				
                clock_t end_lagr = clock();
			
				lagr.time = double(end_lagr - start_lagr) / CLOCKS_PER_SEC;
				
                disj_prop->set_cost_task(cost_task, debug_lagrangian);
                disj_prop->set_cost_layer(cost_layer, debug_lagrangian);
                disj_prop->set_cost_pre(cost_pre, debug_lagrangian);
                

			}
		}
        
        if (lagrangian > 0) {
            cout <<"\n-- Summary Lagrangian --" <<endl;
            cout << " - Original Bound    = " << lagr.LB_original<< endl;
            cout << " - Lagrangian Bound  = " << lagr.LB_lagr << endl;
            cout << " - Number iterations = " << lagr.iterations << endl;
            cout << " - Total Lagrangian time = " << lagr.time << "\n\n";
        }
        
		// ---------------------------------------------------------

		bool is_solved = false;
        
        //Solver parameters
        double tol = 0.0001;
        cp.setParameter(IloCP::Workers, 1);
        cp.setParameter(IloCP::TimeLimit, TIME_LIMIT);
        cp.setParameter(IloCP::RelativeOptimalityTolerance, tol);
        
        //DFS parameters:
        if (DFS) {
            cp.setParameter(IloCP::SearchType, IloCP::DepthFirst);
            cp.setParameter(IloCP::FailureDirectedSearch, IloCP::Off);
            cp.setParameter(IloCP::Presolve, IloCP::Off);
        }
        //cp.setParameter( IloCP::ConflictRefinerOnVariables, IloCP::On);
        //cp.setParameter(IloCP::FailLimit, 30000);

        if (MDDwidth == 0) 	{
            //Running CP model - setting up search phase
            IloIntervalSequenceVarArray seqs(env);
            seqs.add(tour);
            IloSearchPhase search_seq(env, seqs);
            IloSearchPhaseArray search(env);
            search.add(search_seq);
            cp.setSearchPhases(search);
            
            is_solved = cp.solve();
        }
        else{
            //Running MDD approach
            is_solved = cp.solve(IloShortestPathDD(env, permutation, disj_prop));
        }

        
		if (is_solved) {
	 		solution.time =  cp.getTime();
	 		solution.fails = cp.getInfo(IloCP::NumberOfFails);
	 		solution.branches = cp.getInfo(IloCP::NumberOfBranches);
	 		solution.value = cp.getObjValue();
	 		
	 		if(cp.getStatus() == IloAlgorithm::Optimal)
	 			solution.status = "Optimal";
 			else if (cp.getStatus() == IloAlgorithm::Feasible)
 				solution.status = "Feasible";
	 		else 
 				solution.status = "Other";
	 		
            if(debug){
                vector<int> startTime(NbCities);
                vector<int> sequence(NbCities);
	 				
                // get startTime
                for (IloInt i = 0; i < NbCities; ++i) startTime[i] = cp.getStart(city[i]);
      		
                //get first in the sequence
                string name = "";
                IloIntervalVar aux(env,1);
                aux = cp.getFirst(tour);
                name = aux.getName();
      		
                name = name.substr(6,2);
                if (name[1]==')') name.erase(1);

                istringstream convert(name);
                convert >> sequence[0];
      		
                //get others in the sequence
                for (int i = 1; i < NbCities; ++i){
                    aux = cp.getNext(tour, aux);
                    name = aux.getName();
                    name = name.substr(6,2);
      			
                    if (name[1]==')') name.erase(1);

                    istringstream convert(name);
                    convert >> sequence[i];
                }
      		
                PrintSolution(solution, NbCities, NbCommodities, capacity, startTime, sequence);
                CheckSolution(sequence, AbsWeight, pairs, capacity, dd, solution.value);
            }
            
    	} else {
      		
      		solution.status = "No Solution";
      		solution.time = cp.getTime();
      		solution.fails = cp.getInfo(IloCP::NumberOfFails);
	 		solution.branches = cp.getInfo(IloCP::NumberOfBranches);
      		solution.value = 0;
      		
      		#if WRITE > 0	
      		PrintSolution(solution, NbCities, NbCommodities, capacity);
      		#endif	
    	}

  	} catch (IloException& ex) {
    	env.out() << "Error: " << ex << endl;
  	}
  	env.end();
  	return 0;
}


//===========================================================================================
//  Print and check solutions
//===========================================================================================

void PrintParameters( IloTransitionDistance distance, int n, IloIntArray2 pairs){

	int m = pairs.getSize();

	cout << "Distance Matrix: " << endl;
	for(IloInt i = 0; i < n; i++){
		for( IloInt j =0; j < n; j ++){
			cout << distance.getValue(i,j) << " ";
		}

		cout<< endl;
	}
	cout << "Commodities = " <<  m << endl; 
	cout << "Precedences:"  << endl;
	for(int j = 0; j < m; j++){
		cout << pairs[j][0] << " --> " << pairs[j][1] << endl;
	} 

}

void PrintSolution(Solution& solution, int NbCities, int NbCommodities, int capacity, vector<int> stime, vector<int> sequence) {

	ofstream outputFile;
	string output = "checking_solution.txt";
	outputFile.open (output.c_str());
	
	// Print initial parameters
	outputFile << "Class: " <<  classNumber << " Name: " << fileName << "\n";
	outputFile << "Time(sec): " <<  solution.time  << endl;
	outputFile << "Number of Fails: " <<  solution.fails  << endl;
	outputFile << "Number of Branches: " <<  solution.branches  << endl;
	outputFile << "Objective: " << solution.value << " Status: " <<  solution.status <<  endl;
	outputFile << "n(cities):" << NbCities << " m(commodity):" << NbCommodities << " Q(capacity):" << capacity << endl;
	
	//Print StartTimes of sequence
	if(solution.status != "No Solution"){
		outputFile << "Solution start times: \n"; 
		for(int i = 0; i < NbCities; i++){
			outputFile << "City " << sequence[i]  << " : [" << stime[sequence[i]]  << " , " << stime[sequence[i]] + 1 << "]\n";
		}
	}
	outputFile.close();
}

void CheckSolution(vector<int> sequence, IloIntArray AbsWeight, IloIntArray2 pairs, int C, vector<int> dd, int obj) {

	int n = sequence.size();
	int m = pairs.getSize();

	vector<int> weight(n);
	vector2 pd(m);
	
	for (int j = 0; j < m; j++){
		pd[j] = vector<int>(2);
		pd[j][0] = pairs[j][0]; 
		pd[j][1] = pairs[j][1]; 
	}
	
	for (int i = 0; i < n; i++){
		weight[i] = AbsWeight[i];
	}
	
	//Check all conditions
	if( !CheckAlldiff(sequence) ) return;
	if( !CheckPrecedence(sequence, pd) ) return;
	if( !CheckCapacity(sequence, weight, C) ) return;
	if( !CheckObjective(sequence, dd, obj) ) return;
	
	cout << "Solution is feasible! \n";
}
