#ifndef PROBLEM_H
#define PROBLEM_H

#include "common_definitions.h"
#include "cell.h"
#include "output_manager.h"
#include "tools.h"

class ProblemPoisson
{
	public:
    /* Constructor. */
    ProblemPoisson(); 
    /* Run this problem. */
    void run(); 
    	
	private:
         
    VectorXd x, b;
    SparseMatrix<double> A;
    vector<T> triplet_list;

    void pre_processing();
    void assemble_system();
    void solve();
    void output_results();	
    void post_processing();	


};

#endif