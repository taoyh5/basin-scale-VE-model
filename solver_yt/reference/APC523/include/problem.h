#ifndef PROBLEM_H
#define PROBLEM_H

#include "common_definitions.h"

#include "element.h"
#include "fe_value.h"
#include "cg_solver.h"
#include "operations.h"
#include "dof.h"
#include "print.h"
#include "bc_conditions.h"

class Problem
{
	public:
        /* Constructor. */
    	Problem(); 
        /* Run this problem. */
    	void run(); 
    	
	private:
         
        /* Using finite element method, we finally get a linear system Ax=b
           {global_matrix} is A
           {global_rhs} is b
           {solution} is x */
		matrix global_matrix; 
		vec global_rhs;       
		vec solution;

        /* {setup_system} initializes {global_matrix}, {global_rhs} and {solution}*/
    	void setup_system();

    	/* {assemble_system} is a very important function in our code!
            In this function, we loop over each finite element cell, then loop
            over each quadrature point, and then loop over each basis function to
            collection contributions. After each cell, we would distribute local
            contributions to {global_matrix} and {global_rhs}. This is how we
            get our linear system in a big picture.
            As mentioned earlier, this code is written in a generic way. If anyone wants
            to use this code to write DG-FEM method to solve another equation, he/she would 
            only need to change this {assemble_system} function! (Of course he/she would also
            need to rewrite the cg_solver if the system he/she gets is non-symmetric.) */
        void assemble_system();

        /* Solve the linear system using some linear solver. In this case, the congugate
           gradient method is used. */
    	void solve();

        /* Output {solution} to both terminal and a seperate .vtk file. */
    	void output_results();	

        /* Do L2 norm error estimation */
        void estimate_error();	


};

#endif