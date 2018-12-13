#include "problem_poisson.h"


ProblemPoisson::ProblemPoisson()
:
x(CELL_NUMBER),
b(CELL_NUMBER),
A(CELL_NUMBER, CELL_NUMBER)
{}

void ProblemPoisson::run()
{
  pre_processing();
  assemble_system();
  solve();
  output_results();
  post_processing();
}

void ProblemPoisson::pre_processing()
{
  // Do we have anything to set up?
}

void ProblemPoisson::assemble_system()
{

  int counter = 0;
  // Iterate over cells
  for (int id = 0; id < CELL_NUMBER; ++id)
  {
    Cell cell(id);
    // Check out the syntax of C++ STL map if you are not familiar
    // It's like the dictionary in Python
    map<DIRECTION, int>::iterator it;
    // Iterate over faces
    for (it = cell.neighbour_ids_.begin(); it != cell.neighbour_ids_.end(); it++)
    {
      // Face contribution
      if (it->second == OUTSIDE_CELL_ID)
      {
        triplet_list.push_back(T(cell.id_, cell.id_, -2));
        double x_boundary_value = get_boundary_value_x(cell, it->first);
        b(cell.id_) += -2*x_boundary_value;
      }
      else
      {
        triplet_list.push_back(T(cell.id_, cell.id_, -1));
        triplet_list.push_back(T(cell.id_, it->second, 1));
      }
    }
    for (it = cell.neighbour_ids_.begin(); it != cell.neighbour_ids_.end(); it++)
    {
      // Purely for test purpose
      if (it->second == OUTSIDE_CELL_ID)
        counter++;   
    }
    // Volume contribution
    double rhs_f = rhs_function(cell.cell_center_);
    b(cell.id_) += -rhs_f*Area;
  }

  // To do: Jiarong -> You may want to change the following "assert" into an execption "try-catch" in the test code
  assert(("Bug! Wrong calculation of total number of boundary faces.", counter == 4*LINE_CELL_NUMBER));

  // This is a very efficient way of assembling the matrix A. The function is provided by the external library Eigen.
  A.setFromTriplets(triplet_list.begin(), triplet_list.end());

}


void ProblemPoisson::solve()
{

  // SimplicialCholesky<SparseMatrix<double>> chol(A);  // performs a Cholesky factorization of A
  // VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side

  SparseLU<SparseMatrix<double>>   solver;
  solver.analyzePattern(A); 
  solver.factorize(A); 
  x = solver.solve(b); 

  // cout << MatrixXd(A) << endl;
  // cout << x << endl;
  // cout << b << endl;
  // cout << x.maxCoeff() << endl;
}



void ProblemPoisson::output_results()
{
  OutputManager output_manager(0);
  output_manager.scalar_output(x, "u");
}



void ProblemPoisson::post_processing()
{
  double error = 0;
  for (int i = 0; i < CELL_NUMBER; ++i)
  {
    Cell cell(i);
    error += pow(boundary_function_x(cell.cell_center_) - x(i), 2)*Area; 
  }
  error = sqrt(error);

  // L2 norm should convege in second order w.r.t. mesh size
  cout << "L2 norm error is " << error << endl;
}





