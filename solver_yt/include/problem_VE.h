#ifndef PROBLEM_VE_H
#define PROBLEM_VEH


#include "common_definitions.h"
#include "cell.h"
#include "output_manager.h"
#include "tools.h"

class ProblemVE
{
	public:
    /* Constructor. */
    ProblemVE(); 
    /* Run this problem. */
    void run(); 
    	
	private:
         
    VectorXd new_Pc_solution, new_Sc_solution, old_Sc_solution, old_Sc_solution, new_Pb_solution, new_Sb_solution, old_Sb_solution, old_Sb_solution, system_rhs, flag_vector;
    SparseMatrix<double> system_matrix;
    vector<T> triplet_list;

    double time;

    void run_one_step();

    void assemble_system_Pb();
    void solve_Pb();
    void compute_Pc();
    void compute_Sc();
    void compute_Sb();
    void output_results(int cycle);	
    void post_processing();


};



const double mu_c = 4.25e-5;

const double mu_b = 3e-4;

const double grav = 9.8;

const double rho_c = 710; 

const double rho_b = 1000;

const double rg_c = rho_c * grav;

const double rg_b = rho_b * grav;

const double sb_res = 0;

const double krc_star = 1;

const double formation_thickness = 50; 

const double porosity = 0.03; 

const double Pcap_entry = 0; 

#endif
