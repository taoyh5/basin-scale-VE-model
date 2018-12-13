#include "problem.h"


Problem::Problem()
{
	// Nothing happens here in constructor
}

void Problem::run()
{
    setup_system();
    assemble_system();
    solve();
    output_results();
    estimate_error();
}

void Problem::setup_system()
{
	for (int i = 0; i < TOTAL_DOF; ++i)
	{
		global_rhs.push_back(0);
		solution.push_back(0);
		vec vec_temp(TOTAL_DOF, 0);
		global_matrix.push_back(vec_temp);
	}
}

void Problem::assemble_system()
{
	// Loop over each finite elemnt cell.
	// {id} refers to the id of each finite element cell.
	for (int id = 0; id < ELE_NUMBER; ++id)
	{
		Element element(id);
		FeValue fe_value(&element);

		// {local_matrix} and {local_rhs} are containers for local contributions.
		// We need to distribute local contributions to {global_matrix} and {global_rhs}
		// later on.
		matrix local_matrix(DOF_PER_ELE, vec(DOF_PER_ELE, 0));
		vec local_rhs(DOF_PER_ELE, 0);

		// Obtain a map between local indices and global indices
		int_vec dof_indices(DOF_PER_ELE);
		get_dof_indicies(element.id, dof_indices);

		cout << "*************************************************" << endl;
		cout << "element id is " << element.id << "  " << endl;

		// Loop over each quadrature point
		for (int q = 0; q < QUAD_NUM*QUAD_NUM; ++q)
		{
			// Let {fe_value} know which quadrature point it is focusing on
			fe_value.get_quadrature(q);
			// Prepare for the value of f
			vec x(DIM, 0);
			x = fe_value.get_physical_crds(fe_value.xi, element.x_list);
			// Loop over each shape function
			for (int i = 0; i < DOF_PER_ELE; ++i)
			{
				// This contribution comes from the term "f" in the poisson equation
          		local_rhs[i] += fe_value.shape_value(i)*f_value(x)*fe_value.JxW();
          		// Loop over each shape function
				for (int j = 0; j < DOF_PER_ELE; ++j)
				{
					// This contribution comes from the term "grad(u) \cdot grad(v)"
					vec vec_temp1(DIM); vec_temp1 = fe_value.shape_gradient(i);
					vec vec_temp2(DIM); vec_temp2 = fe_value.shape_gradient(j);
					local_matrix[i][j] += vec_mul_vec(vec_temp1, vec_temp2)*fe_value.JxW();
				}
			}
		}

		// Loop over each face of each element. In our case, an element(square) has
		// four faces(sides). This is necessary because we are using DG method. 
		// Standard finite element method doesn't need to loop over interior faces.
		for (int face_id = 0; face_id < FACE_PER_ELE; ++face_id)
		{	
			ElementFace element_face(id, face_id);
			FeFaceValue fe_face_value(&element, &element_face);

			// It is important to distinguish what is a interior face from what is a boundary face!
			if(element_face.is_boundary)
			{
				// Now we are sitting on a boundary face

				vec normal_vector(DIM); 
				normal_vector = element_face.normal_vector;

				// Loop over each quadrature point
				// Notice that the quadrature points of a face are different from that of a cell
				for (int q = 0; q < QUAD_NUM; ++q)
				{
					fe_face_value.get_quadrature(q);

					// Prepare for the boundary conditions
					vec x(DIM, 0);
					x = fe_face_value.get_physical_crds(fe_face_value.xi, element_face.x_list);

					for (int i = 0; i < DOF_PER_ELE; ++i)
					{
						vec vec_temp1(DIM); vec_temp1 = fe_face_value.shape_gradient(i);
					
						// The penalty terms
    					local_rhs[i] += GAMMA/H*fe_face_value.shape_value(i)*
                   						Dirichlet_bc_value(x)*fe_face_value.JxW();
    					local_rhs[i] += -vec_mul_vec(vec_temp1, normal_vector)*
    									Dirichlet_bc_value(x)*fe_face_value.JxW();

						for (int j = 0; j < DOF_PER_ELE; ++j)
						{
							vec vec_temp2(DIM); vec_temp2 = fe_face_value.shape_gradient(j);
							
							// Contributions. For details, pls see the report of how DG method is formulated
							local_matrix[i][j] += GAMMA/H*fe_face_value.shape_value(i)*
                           						  fe_face_value.shape_value(j)*fe_face_value.JxW();
        					local_matrix[i][j] += -vec_mul_vec(vec_temp1, normal_vector)*      
                  								   fe_face_value.shape_value(j)*fe_face_value.JxW();
        					local_matrix[i][j] += -vec_mul_vec(vec_temp2, normal_vector)*      
                  								   fe_face_value.shape_value(i)*fe_face_value.JxW();
						}
					}
				}

			}
			else
			{
				// Now we are sitting on an interior face

				// Because each interior face will be looped twice because two elements may share a face,
				// we have to make a choice. We only want it once rather than twice.
				if(element_face.id > element_face.neigh_id)
				{

					// Only if the face belong to an element whose id is relatively larger, shall we 
					// count the contributions.

					// Initialize {element_neigh} as a neighbour element
					// Initialize the corresponding {fe_face_value_neigh}
					Element element_neigh(element_face.neigh_id);
					ElementFace element_face_neigh(element_face.neigh_id, element_face.neigh_face_id);
					FeFaceValue fe_face_value_neigh(&element_neigh, &element_face_neigh);

					// Create local containers
					matrix vi_ui_matrix(DOF_PER_ELE, vec(DOF_PER_ELE, 0));
					matrix vi_ue_matrix(DOF_PER_ELE, vec(DOF_PER_ELE, 0));
					matrix ve_ui_matrix(DOF_PER_ELE, vec(DOF_PER_ELE, 0));
					matrix ve_ue_matrix(DOF_PER_ELE, vec(DOF_PER_ELE, 0));

					int_vec this_dof_indices(DOF_PER_ELE);
					int_vec neigh_dof_indices(DOF_PER_ELE);

					get_dof_indicies(element_face.id, this_dof_indices);
					get_dof_indicies(element_face_neigh.id, neigh_dof_indices);

					vec normal_vector(DIM); 
					normal_vector = element_face.normal_vector;
					vec normal_vector_neigh(DIM);
					normal_vector_neigh = element_face_neigh.normal_vector;

					// The following contributions are made over the interior faces. Details pls
					// refer to the report.
 					for (int q = 0; q < QUAD_NUM; ++q)
					{
						fe_face_value.get_quadrature(q);
						fe_face_value_neigh.get_quadrature(q);
    
						for (int i = 0; i < DOF_PER_ELE; ++i)
						{
							vec vec_temp1(DIM); vec_temp1 = fe_face_value.shape_gradient(i);

      						for(int j = 0; j < DOF_PER_ELE; ++j)
      						{	
      							vec vec_temp2(DIM); vec_temp2 = fe_face_value.shape_gradient(j);
							     						 	
     						 	vi_ui_matrix[i][j] += (GAMMA/H*fe_face_value.shape_value(i)*
      						 						  fe_face_value.shape_value(j)*
      						 						  vec_mul_vec(normal_vector, normal_vector)-
													  0.5*vec_mul_vec(vec_temp1, normal_vector)*
      						 						  fe_face_value.shape_value(j)-
					 						 		  0.5*vec_mul_vec(vec_temp2, normal_vector)*
      						 						  fe_face_value.shape_value(i)
      						 						  )*fe_face_value.JxW();

      						}

      						for(int j = 0; j < DOF_PER_ELE; ++j)
      						{	
      							vec vec_temp2(DIM); vec_temp2 = fe_face_value_neigh.shape_gradient(j);

      						 	vi_ue_matrix[i][j] += (GAMMA/H*fe_face_value.shape_value(i)*
      						 						  fe_face_value_neigh.shape_value(j)*
      						 						  vec_mul_vec(normal_vector_neigh, normal_vector)-
													  0.5*vec_mul_vec(vec_temp1, normal_vector_neigh)*
      						 						  fe_face_value_neigh.shape_value(j)-
					 						 		  0.5*vec_mul_vec(vec_temp2, normal_vector)*
      						 						  fe_face_value.shape_value(i)
      						 						  )*fe_face_value.JxW();
      						}

      					}

						for (int i = 0; i < DOF_PER_ELE; ++i)
						{
							vec vec_temp1(DIM); vec_temp1 = fe_face_value_neigh.shape_gradient(i);

      						for(int j = 0; j < DOF_PER_ELE; ++j)
      						{	
      							vec vec_temp2(DIM); vec_temp2 = fe_face_value.shape_gradient(j);

      						 	ve_ui_matrix[i][j] += (GAMMA/H*fe_face_value_neigh.shape_value(i)*
      						 						  fe_face_value.shape_value(j)*
      						 						  vec_mul_vec(normal_vector, normal_vector_neigh)-
													  0.5*vec_mul_vec(vec_temp1, normal_vector)*
      						 						  fe_face_value.shape_value(j)-
					 						 		  0.5*vec_mul_vec(vec_temp2, normal_vector_neigh)*
      						 						  fe_face_value_neigh.shape_value(i)
      						 						  )*fe_face_value.JxW();
      						}

      						for(int j = 0; j < DOF_PER_ELE; ++j)
      						{	
      							vec vec_temp2(DIM); vec_temp2 = fe_face_value_neigh.shape_gradient(j);

      						 	ve_ue_matrix[i][j] += (GAMMA/H*fe_face_value_neigh.shape_value(i)*
      						 						  fe_face_value_neigh.shape_value(j)*
      						 						  vec_mul_vec(normal_vector_neigh, normal_vector_neigh)-
													  0.5*vec_mul_vec(vec_temp1, normal_vector_neigh)*
      						 						  fe_face_value_neigh.shape_value(j)-
					 						 		  0.5*vec_mul_vec(vec_temp2, normal_vector_neigh)*
      						 						  fe_face_value_neigh.shape_value(i)
      						 						  )*fe_face_value.JxW();
      						}

      					}

  					}
  					// Distribute our local contributions to {global_matrix}
  					distribute_local_to_global(global_matrix, vi_ui_matrix, this_dof_indices, this_dof_indices);
  					distribute_local_to_global(global_matrix, vi_ue_matrix, this_dof_indices, neigh_dof_indices);
  					distribute_local_to_global(global_matrix, ve_ui_matrix, neigh_dof_indices, this_dof_indices);
  					distribute_local_to_global(global_matrix, ve_ue_matrix, neigh_dof_indices, neigh_dof_indices);
				}
			}

		}
		cout << "*************************************************" << endl << endl;

		// Distribute our local contributions to {global_matrix} and {global_rhs}
		distribute_local_to_global(global_matrix, local_matrix, dof_indices, dof_indices);
		distribute_local_to_global(global_rhs, local_rhs, dof_indices);
	}


}

void Problem::solve()
{

	// Call {cg_solver} to solve the system. 
	solution = cg_solver(global_matrix, global_rhs);
}

void Problem::output_results()
{
	// Create a vtk file for output
	ofstream vtkstream("output/u.vtk");

	// Follow the format of vtk files and write our mesh and solution into this file.
	vtkstream << "# vtk DataFile Version 3.0" << endl;
	vtkstream << "#This file was generated by the IP_DG_poisson solver (author: T. Xue)" << endl;
	vtkstream << "ASCII" << endl;
	vtkstream << "DATASET UNSTRUCTURED_GRID" << endl << endl;
	vtkstream << "POINTS " << VERTICES_PER_ELE*ELE_NUMBER << " double" << endl;
	for (int id = 0; id < ELE_NUMBER; ++id)
	{
		Element element(id);
		for (int v = 0; v < VERTICES_PER_ELE; ++v)
		{
			vtkstream << element.vertices[v][0] << " "
			<< element.vertices[v][1] << " " << "0" << endl;
		}
	}
	vtkstream << endl << "CELLS " << ELE_NUMBER << " " << 5*ELE_NUMBER << endl;
	for (int id = 0; id < ELE_NUMBER; ++id)
	{
		vtkstream << 4 << "\t" << 0 + 4*id << "\t"
		<< 1 + 4*id << "\t" << 3 + 4*id << "\t"
		<< 2 + 4*id << endl;
	}
	vtkstream << endl << "CELL_TYPES " << ELE_NUMBER << endl;
	for (int id = 0; id < ELE_NUMBER; ++id)
	{
		vtkstream << " " << CELL_TYPE;
	}
	vtkstream << endl << "POINT_DATA" << " " << VERTICES_PER_ELE*ELE_NUMBER << endl;
	vtkstream << "SCALARS u double 1" << endl; 
	vtkstream << "LOOKUP_TABLE default" << endl;
	for (int id = 0; id < ELE_NUMBER; ++id)
	{
		for (int i = 0; i < VERTICES_PER_LINE; ++i)
    	{	
        	for (int j = 0; j < VERTICES_PER_LINE; ++j)
        	{
        		vtkstream << solution[id*DOF_PER_ELE + i*(DOF_PER_ELE - ORDER - 1) + j*ORDER]
        		<< " ";
        	}
        }
    }
    cout << endl;


    // Also, print the result to terminal to have a general look
	cout << "Solves Ax = y" << endl;
    cout << "Total DOF is " << TOTAL_DOF << endl;
    print( "\nx:", solution );
    // print( "\nCheck AX:", matrixTimesVector( global_matrix, solution ) );
}


void Problem::estimate_error()
{

	double u_h, u_e, diff;
	double error = 0;
	double area = pow(H/ORDER, 2);
	double error_inf = 0;

	double error1;
	double error2;

	double test_area = 0;

	for (int id = 0; id < ELE_NUMBER; ++id)
	{
		Element element(id);

		cout << "element " << id << endl;

		error1 = error;

		for (int i = 0; i < DOF_PER_ELE; ++i)
		{
			vec x(DIM, 0);
			x = element.x_list[i];
			u_e = u_exact(x);
			u_h = solution[id*DOF_PER_ELE + i];
			diff = fabs(u_e - u_h);

			error_inf = diff > error_inf ? diff : error_inf;


			if (i < ORDER + 1 || i > (ORDER + 1)*ORDER - 1)
			{
				if (i%(ORDER + 1) == 0 || (i + 1)%(ORDER + 1) == 0)
				{
					error += pow(diff, 2)*area/4.;
					test_area += area/4.;
				}
				else
				{
					error += pow(diff, 2)*area/2.;	
					test_area += area/2.;		
				}
			}
			else
			{
				if (i%(ORDER + 1) == 0 || (i + 1)%(ORDER + 1) == 0)
				{
					error += pow(diff, 2)*area/2.;
					test_area += area/2.;				
				}
				else
				{
					error += pow(diff, 2)*area/1.;	
					test_area += area/1.;		
				}				
			}
		}

		error2 = error;
		cout << " error on this element is " << error2 - error1 << endl << endl;
	}
	error = sqrt(error);

	cout << "Estimated L2 norm error is " << error << endl;
	cout << "Estimated max norm error is " << error_inf << endl;
	cout << "Test area is " << test_area << endl;

}




