#ifndef DOF_H
#define DOF_H

#include "common_definitions.h"

/* Given which element we are sitting on and which dof number we have locally, return global dof number*/
int number_local_to_global(int ele_id, int local_dof_number);

/* Assign the global dof numbers to {dof_indices} for some element with id being {ele_id} */
void get_dof_indicies(int ele_id, int_vec& dof_indices);

/* Distribute contributions from {local_matrix} to {global_matrix}
   {dof_indices1} corresponds to the local index {i} while 
   {dof_indices2} corresponds to the local {j} */
void distribute_local_to_global(matrix &global_matrix, matrix &local_matrix,
		 int_vec &dof_indices1, int_vec& dof_indices2);

/* Distribute contributions from {local_rhs} to {global_rhs} */
void distribute_local_to_global(vec &global_rhs, vec &local_rhs,
		 int_vec &dof_indices);

#endif