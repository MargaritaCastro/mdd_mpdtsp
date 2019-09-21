

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "mddFiOrcl.hpp"
#include "OPTvect.h"
#include <stdlib.h>


// -------------------------------------------------------------------------------
//Constructor
mddFiOrcl::mddFiOrcl( DisjunctivePropagator* disj_prop, int n, int k, int lagr_type, bool debug1, istream *iStrm) : type(lagr_type) {

	//Initialize variables
	mdd = disj_prop;
	nCities = n;
    nCom = k;
	debug = debug1;
	
	//Set cost that we will give to mdd
	cost_task 	= new double[nCities];
	cost_cap 	= new double[nCities];
	cost_pre 	= new double[nCities];
	
	for( int i = 0; i < nCities; i++){
		cost_task[i]= 0.0;
		cost_cap[i]	= 0.0;
		cost_pre[i] = 0.0;
	}
        
    lip_computed = false;
	
	//Initialize Fi variables
	update_variables();
    first_call = true;
    
    //Initialize Lipschitz constants
    Lip_tour = 0.0;
    Lip_cap = 0.0;
    Lip_pre = 0.0;
    Lip_const = 0.0;
        
    if(debug) {
        cout << "\nStarting Fi Oracle " <<endl;
        cout << "\tLagrangian type = " << type << endl;
        cout << "\tNumber of variables = " << NumVar << endl;
    }
	
}


// -------------------------------------------------------------------------------
// Set parameters of the mPDTSp problem
void mddFiOrcl::set_capacity(int C) {
	cap = C;
}

void mddFiOrcl::set_weight(IloIntArray& AbsWeight) {

	weight = new int[nCities];

	for(int i = 0; i < nCities; i++){
		weight[i] = AbsWeight[i];
	}

}

void mddFiOrcl::set_precedences(IloIntArray2& pairs) {

	vector<int> aux_pre(2);

	for(int k = 0; k < nCom; k++){
		aux_pre[0] = pairs[k][0];
		aux_pre[1] = pairs[k][1];
						
		precedence.push_back(aux_pre);
	}

}

// -------------------------------------------------------------------------------
// Methods to solve the lagrangian subproblem

void mddFiOrcl::update_variables() {
    
    if(debug) cout << "Update variable types" << endl;
	
    if (type == 0)  NumVar = nCities;               // Case 0 : Tour constraint relaxation
    else if( type == 2) NumVar = nCities*2;         // Case 2 : Capacity constraint relaxation
    else if( type == 3) NumVar = nCities*3;         // Case 3 : Tour + Capacity constraint relaxation
    else if( type == 4) NumVar = nCom;              // Case 4 : Precedence constraint relaxation
    else if( type == 5) NumVar = nCities + nCom;    // Case 5 : Tour + precedence constraint relaxation
    else{
        cout << "We don't support this type for lagrangian relaxation";
        exit(1);
    }
    
    var_type = new bool[NumVar];
    
    if (type == 0 || type == 3 || type == 5 ){
        for (unsigned int i = 0; i < NumVar; i++){
          if (i < nCities) var_type[i] = true;      // equality constriant for tour
          else             var_type[i] = false;     // greater or equal constraint for the other constraints
        }
    }else {
        for (unsigned int i = 0; i < NumVar; i++) var_type[i] = false;
        
        if(debug) cout << "\t All variables are >= 0" << endl;
    }

}

HpNum mddFiOrcl::Fi( cIndex wFi ) {
	
    if( Fit ) Fit->Start();
    
	HpNum out = 0.0 ; //double

	//Option 1 : wFi == 0 --> constant value of the function
	//Option 2 : wFi == 1 --> shortest path computation (without the constant value)
	//Option 3 : wFi == Inf<Index>() --> both values sum together
	
    if(debug) cout << "\n==== NEW ITERATION ====" << endl;

	//Step 2: Constant part
    HpNum constant = 0;
	
	if (wFi == 0) {
		out += constant;
	}
	//Step 3: Shortest path computation
	else {
        if(debug) cout << " == Computing MDD shortest path ==" << endl;
        
        //Step 1: Compute the cost array for the mdd
        update_cost();
        
		sp_mdd = 0.0;
		
		sequence = mdd->compute_shortest_path(sp_mdd, cost_task, cost_cap, cost_pre);
        
        if(debug){
            cout << "\tSequence : ";
            for(int i=0; i < nCities; i++) cout << sequence[i] << " ";
            cout << "\n\tCost : " << sp_mdd << " constant = " << constant << endl;
        }
		
		if (wFi != Inf<Index>())	out = sp_mdd;
		else 						out = sp_mdd - constant; 		// Do not consider constant values;
        
        //out = sp_mdd;
        
        if (first_call) {
            first_bound = sp_mdd;
            first_call = false;
        }
	}
	
    if(Fit) Fit->Stop();
    
	return (-out);
}

// -------------------------------------------------------------------------------
// Subgradient functions

bool mddFiOrcl::NewGi( cIndex wFi){
    //TODO: change this to terminate
	return LHasChgd;
}

Index mddFiOrcl::GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name , cIndex strt , Index stp){

    //Timing the function
    if(Fit) Fit->Start();
    
    if(debug) cout << " == Computing Subgradients ==" << endl;

	//Case 1: constant part of the subgradient
	if( Name == MaxName ) {
        for (int i = 0; i < NumVar; i++) SubG[i] = 0;
	
	//Case 2: the full subgradient
	}else {
        
		if (type == 0 || type == 3 || type == 5) {
			double* grad_aux = new double[nCities];
			gradient_tour(grad_aux);
			for(int i = 0; i < nCities; i++) SubG[i] = -grad_aux[i];
            
            if(type == 3) {
                double* grad_aux2 = new double[nCities*2];
                gradient_cap(grad_aux2);
                for(int i = nCities; i < NumVar; i++) SubG[i] = -grad_aux2[i - nCities];
                
                delete[] grad_aux2;
            }
            else if(type == 5) {
                double* grad_aux2 = new double[nCom];
                gradient_pre(grad_aux2);
                for(int i = nCities; i < NumVar; i++) SubG[i] = -grad_aux2[i - nCities];
                
                delete[] grad_aux2;
            }
            delete[] grad_aux;
		}
        else{
            double* grad_aux = new double[NumVar];
            
            if (type == 2)  gradient_cap(grad_aux);
            else            gradient_pre(grad_aux);
            
            for(int i = 0; i < NumVar; i++) SubG[i] = -grad_aux[i];
            
            delete[] grad_aux;
        }
        
        //Check if solution is optimal
        bool optimal = false;
        
        if(optimal){
            lower_bound = -sp_mdd; //TODO
        }
	
	} 

   	SGBse = 0;
	LHasChgd = false; 
	
    if( Fit ) Fit->Stop();
	
	return NumVar; // number of variables that are non-zero
}

HpNum mddFiOrcl::GetVal( cIndex Name){
	if( Name < MaxName ) throw( NDOException( "KnpsFiOrcl::GetVal: Not implemented yet" ) );

 	return( 0 );
}

// -------------------------------------------------------------------------------
// Get functions... for the solver

bool mddFiOrcl::GetUC(cIndex i) {
	return var_type[i];
}

LMNum mddFiOrcl::GetUB(cIndex i) {
    //return (-lower_bound);
    return( Inf<LMNum>());
}

LMNum mddFiOrcl::GetBndEps(){
    return 0.5;
}

HpNum mddFiOrcl::GetLowerBound( cIndex wFi){
    return lower_bound;
}
void mddFiOrcl::SetLowerBound( HpNum lb){
    lower_bound = lb;
}

HpNum mddFiOrcl::GetGlobalLipschitz(cIndex wFi) {
    if (!lip_computed) compute_lipschitz();
    return Lip_const;
}

void mddFiOrcl::SetMDDMultipliers( cLMRow Lmbd ){
    
    Lambda = Lmbd;
    
    update_cost();
    
    //Set parameters inside mdd
    mdd->set_cost_task(cost_task, true);
    mdd->set_cost_layer(cost_cap, true);
    mdd->set_cost_pre(cost_pre, true);
    
}

// -------------------------------------------------------------------------------
// Compute cost to give to mdd

void  mddFiOrcl::update_cost() {
    
    if(debug) cout << " == Updating cost == " << endl;

	vector<double> u_task(nCities, 0.0);
	vector<double> u_cap(nCities*2, 0.0);
	vector<double> u_pre(nCom, 0.0);

	// Case 0 : Tour constraint relaxation 
	if(type == 0 || type == 3 || type == 5) {
		for (int i = 0; i < nCities; i++) u_task[i] = Lambda[i];
		
        if (type == 3) {
            for (int i = nCities; i < NumVar; i++) u_cap[i - nCities] = Lambda[i];

        }
        else if(type == 5 ){
            for (int i = nCities; i < NumVar; i++) u_pre[i - nCities] = Lambda[i];
        }
	
	} else {
        if (type == 2) {
            for (int i = 0; i < NumVar; i++) u_cap[i] = Lambda[i];
        }else{
            for (int i = 0; i < NumVar; i++) u_pre[i] = Lambda[i];
        }
    }
    
    //Print Lambda
    if(debug){
        cout << "\tLambda = ";
        for (int i=0; i < NumVar; i++){
            cout << Lambda[i] << " ";
        }
        cout << endl;
    }
    
    cost_tour(u_task);
    cost_capacity(u_cap);
    cost_precedence(u_pre);
}

void mddFiOrcl::cost_tour(vector<double>& u_new) {
	
	//Cost for each possible value
	for(int i = 0; i < nCities; i++)    cost_task[i] = u_new[i];

	//Add constant at the root node
	for(int i = 0; i < nCities; i++)    cost_task[0] -= u_new[i];
	
	//Print cost
	if(debug) {
		cout << "\tCost tour :";
		for(int i = 0; i < nCities; i++) cout << cost_task[i] <<" ";
		
		cout << endl;
	}

}


void mddFiOrcl::cost_precedence(vector<double>& u_new) {
	
	
	for(int i = 0; i < nCities; i++) cost_pre[i] = 0.0;
	
	for (int k = 0; k < nCom; k++){
		cost_pre[ precedence[k][0] ] += u_new[k]; 	// pickup cost
		cost_pre[ precedence[k][1] ] -= u_new[k]; 	// delivery cost
		cost_pre[0] += u_new[k];
	}
	
	//Print cost
	if(debug) {
		cout << "\tCost precedence :";
		for(int i = 0; i < nCities; i++) cout << cost_pre[i] << " ";

		cout << endl;
	}

}


void mddFiOrcl::cost_capacity(vector<double>& u_new) {
    
	int n = nCities;
	cost_cap[n-1] = u_new[n-1] - u_new[2*n - 1];
	cost_cap[0] = 0.0;
	
	//Additive cost
	for(int i = (n-2); i > 0; i--)  cost_cap[i] = cost_cap[i + 1] + u_new[i] - u_new[n + i];
	
	//Cost at the root node
	for(int i = 0; i < n; i++)  cost_cap[0] -= u_new[i]*cap;
	
	if(debug) {
		cout << "\tCost_cap :";
		for(int i = 0; i < n; i++) cout << cost_cap[i] <<" ";
		
		cout << endl;
	}

}

// -------------------------------------------------------------------------------
// Sub gradient computation according to last solution found

void mddFiOrcl::gradient_tour(double* grad) {

	for(int i = 0; i < nCities; i++) grad[i] = -1;

	for(int i = 0; i < nCities; i++) grad[sequence[i]]++;

	if(debug) {
		cout << "\tSubgradient Tour: ";
		for(int i = 0; i < nCities; i++) cout << grad[i] << " ";
		
		cout << endl;
	}
}

void mddFiOrcl::gradient_pre(double* grad) {

	//Add one to all the gradientes
	for(int k = 0; k < nCom; k++) grad[k] = 1;
    //for(int k = 0; k < nCom; k++) grad[k] = 0;

	// Compute gradiente by counting the distances
	for(int k = 0; k < nCom; k++){
		for(int i = 0; i < nCities; i++){
			if (sequence[i] == precedence[k][0]) 		grad[k] += i;       //pickup location
			else if (sequence[i] == precedence[k][1]) 	grad[k] -= i;       //delivery location
		}
	}
	
	if(debug) {
		cout << "\tSubgradient Precedence : ";
		for(int k = 0; k < nCom; k++){
			cout << grad[k] << " ";
		}
		cout << endl;
	}

}

void mddFiOrcl::gradient_cap(double* grad) {
	
	//Adds the weight of the solution
	grad[0] = weight[sequence[0]];          //constraint <= C
	grad[nCities] = -weight[sequence[0]];   //constraint >= 0
	for(int i = 1; i < nCities; i++){
		grad[i] = grad[i - 1] + weight[sequence[i]];						    // constraint <= C
		grad[i + nCities] = grad[i + nCities - 1] - weight[sequence[i]]; 	// constraint >= 0 (signs are inverse)
	}
	
	// Add the - capacity term to every component of the <= C constraint
	for(int i = 0; i < nCities; i++){
		grad[i] -= cap;					// only for constraint <= C
	}
	
	if(debug){
		cout << "\tSubgradient cap : ";
		int n = nCities*2;
		for(int i = 0; i < n; i++){
			cout << grad[i] << " ";
		}
		cout << endl;
	}
}

bool mddFiOrcl::check_solution(){
    //TODO
    
    bool optimal = true;
    
    //Check tour feasibility
    
    
    //Check precedence feasibility
    
    
    //Check capcity feaisbility
    
    
    return optimal;
}

// ----------- Lipschitz constant
void mddFiOrcl::compute_lipschitz(){
    
    // Compute the individual constants
    if (type == 0 || type == 3 || type == 5)    compute_lipschitz_tour();
    if (type == 2 || type == 3 )                compute_lipschitz_cap();
    if (type == 4 || type == 5 )                compute_lipschitz_pre();
    
    // Compute the general constant
    Lip_const = Lip_tour + Lip_cap + Lip_pre ;
    
    lip_computed = true;
}

void mddFiOrcl::compute_lipschitz_tour() {
    Lip_tour = (NumVar-3)*(NumVar-3) +(NumVar -3);
    
    if(debug) cout << "Lip_tour = " << Lip_tour << endl;
}
void mddFiOrcl::compute_lipschitz_cap() {
    
    //Step 1: compute maximum Delta weight
    int max_weight = 0;
    
    for(int i=0; i < nCities; i++){
        if( weight[i] > max_weight) max_weight = weight[i];
    }
    
    //Step 2: compute >=C part;
    Lip_cap = 0.0;
    
    for (int i= 0; i < nCities; i++) {
        if (i == 0){
            Lip_cap += cap*cap;
        }
        else if (i == nCities - 1) {
            Lip_cap += ( (i)*max_weight - cap)*( (i)*max_weight - cap);
        }
        else{
            Lip_cap += ( (i+1)*max_weight - cap)*( (i+1)*max_weight - cap);
        }
    }
    
    //Step 3: compute <=0 part;
    for (int i= 1; i < nCities; i++) {
        if (i == nCities - 1) {
            Lip_cap += ( i*max_weight)*( i*max_weight );
        }else{
            Lip_cap += ( (i+1)*max_weight)*( (i+1)*max_weight );
        }
    }
    
    
    if(debug) cout << "Lip_cap = " << Lip_cap << endl;
}
void mddFiOrcl::compute_lipschitz_pre() {
    // (n-2)*(n-1)/2 +1 per component
    double component = (NumVar-2)*(NumVar-1)/2 +1;
    Lip_pre = (NumVar-1)*(NumVar-1)*nCom;
    
    if(debug) cout << "Lip_pre = " << Lip_pre << endl;
}

// -------------------------------------------------------------------------------

mddFiOrcl::~mddFiOrcl() {
    //Delete variables
    delete [] var_type;
    delete [] weight;
    delete [] sequence;
    
    //Delete MDD cost
    if(false){
        delete [] cost_cap;
        delete [] cost_task;
        delete [] cost_pre;
    }
    
    
}

