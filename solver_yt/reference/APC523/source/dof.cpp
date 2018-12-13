#include "dof.h"

// All functions in this file serves as the dof handler. 
// They support the transition from local contributions. to global contributions.

int number_local_to_global(int ele_id, int local_dof_number)
{
	int global_dof_number = ele_id*DOF_PER_ELE + local_dof_number;
	return global_dof_number;
}


void get_dof_indicies(int ele_id, int_vec& dof_indices)
{
	for (int i = 0; i < dof_indices.size(); ++i)
	{
		dof_indices[i] = number_local_to_global(ele_id, i);
	}
}

void distribute_local_to_global(matrix &global_matrix, matrix &local_matrix,
		 int_vec &dof_indices1, int_vec& dof_indices2)
{
	for (int i = 0; i < DOF_PER_ELE; ++i)
	{
		for (int j = 0; j < DOF_PER_ELE; ++j)
		{
			global_matrix[dof_indices1[i]][dof_indices2[j]] += local_matrix[i][j];
		}
	}
}

void distribute_local_to_global(vec &global_rhs, vec &local_rhs,
		 int_vec &dof_indices)
{
	for (int i = 0; i < DOF_PER_ELE; ++i)
	{
		global_rhs[dof_indices[i]] += local_rhs[i];
	}
}