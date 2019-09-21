

/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __mddFiOrcl
#define __mddFiOrcl /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"
#include "../constraints/disjunctive.hpp"		// For the MDD disjunctive propagator
#include <ilcp/cp.h> 							// only for the parameters used by cp optimizer


using namespace std;
typedef vector<vector<double> > vector2_double;
typedef vector<vector<int> > vector2;


class mddFiOrcl : public FiOracle {

	public:
		
		//Constructor
    mddFiOrcl( DisjunctivePropagator* disj_prop, int n, int k, int lagr_type, bool debug1=false, istream *iStrm = 0);
		
		//Set problem (mPDTSP) parameters
		void set_capacity(int C);
		void set_weight(IloIntArray& AbsWeight);
		void set_precedences(IloIntArray2& pairs);
		
		void set_debug(){
			debug = true;
		}
		
		//Solving the Lagrangian subproblem
		HpNum Fi(cIndex wFi);
		bool NewGi( cIndex wFi = Inf<Index>() );
		Index GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name = Inf<Index>() , cIndex strt = 0 , Index stp = Inf<Index>() );
		HpNum GetVal( cIndex Name = Inf<Index>() );
		
		//get and set  functions
		bool GetUC(cIndex i); 		// sign of the variable; true = unconstraint; false = nonegative
        LMNum GetUB(cIndex i);
        LMNum GetBndEps();
    
        HpNum GetLowerBound( cIndex wFi = Inf<Index>() );
        void SetLowerBound( HpNum lb);
    
        HpNum GetGlobalLipschitz(cIndex wFi = Inf<Index>());
    
        double GetOriginalBound(){
            if(!first_call) return first_bound;
            else            return -1;
        }
    
        //Set parameter to mdd
        void SetMDDMultipliers( cLMRow Lmbd );
		
		//Destructor
		virtual ~mddFiOrcl();
		
	
	protected:
	
	private:
	
		DisjunctivePropagator* mdd;
		bool debug;
        bool first_call;        //first call to Fi function
    
        //Lower bounds
        HpNum lower_bound;
        double sp_mdd;          //current value computed by mdd
        double first_bound;        //first bound computed by the approach
		
		//Parameter for different lagrangian relaxations
		const int type; 						// type of lagrangian relaxation
		bool* var_type;					// type of variable ( unbounded = true, bounded( >=0 ) = false) 
		
		//Parameters of the problem
		int nCities;				// number of locations
		int nCom;					// number of commodities
		int cap;					// maximum capacity of the problem
		int* weight;				// set of weights for each city 
		vector2 precedence;			// pair of precedences
		
		// Cost variables for the mdd
		double* cost_task;
		double* cost_cap;				// called cost_layer in the mdd
		double* cost_pre;
		
		//Optimal shortest path
		int* sequence;
    
        // Lipschitz constant
        bool lip_computed;
        HpNum Lip_const;
        HpNum Lip_tour;
        HpNum Lip_cap;
        HpNum Lip_pre;
		
		//--------- FUNCTIONS ---------//
		
		//Update variables of the problem depending on the lagrangian type
		void update_variables();
		
		// Compute costs to give to mdd
		void update_cost();
		void cost_tour(vector<double>& u_new);
		void cost_precedence(vector<double>& u_new);
		void cost_capacity(vector<double>& u_new);
		
		// Compute gradientes
		void gradient_tour(double* grad);
		void gradient_cap(double* grad);
		void gradient_pre(double* grad);
    
        //Check feasibility of solution
        bool check_solution();
    
        // Lipschitz constant computation
        void compute_lipschitz();
        void compute_lipschitz_tour();
        void compute_lipschitz_cap();
        void compute_lipschitz_pre();

};

#endif  /* mddFiOrcl.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File KnpsFiOrcl.h --------------------------*/
/*--------------------------------------------------------------------------*/




