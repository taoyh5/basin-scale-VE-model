#include "problem_flow.h"


ProblemFlow::ProblemFlow()
:
new_p_solution(CELL_NUMBER), 
new_S_solution(CELL_NUMBER), 
old_p_solution(CELL_NUMBER), 
old_S_solution(CELL_NUMBER), 
system_rhs(CELL_NUMBER),
flag_vector(CELL_NUMBER),
system_matrix(CELL_NUMBER, CELL_NUMBER),
time(0)
{
  // Do we need to initialize those vectors and matrix with all zero? It seems to be automatically done.
}


void ProblemFlow::run()
{

  // interpolate_S(old_S_solution);
  // interpolate_p(old_p_solution);

  interpolate_S(old_S_solution, time);
  interpolate_p(old_p_solution, time);

  cout << "Total number of dofs is " << CELL_NUMBER << endl;

  for (int i = 0; i < 100; ++i)
  {

    time += DELTA_T;

    triplet_list.clear();
    system_rhs.setZero();
    system_matrix.setZero();

    cout << "This is time step " << i << endl;
    run_one_step();
    if (i%10 == 0)
    {
      output_results(i/10);
    }

    // assert(("Numerical issue! Saturation is out of bound!",  new_p_solution.maxCoeff() < 2 &&  new_p_solution.maxCoeff() >  -1));
  }

  post_processing();
}


void ProblemFlow::run_one_step()
{
  assemble_system_p();
  solve_p();
  compute_S();
}


void ProblemFlow::assemble_system_p()
{

  int counter = 0;
  for (int id = 0; id < CELL_NUMBER; ++id)
  {
    Cell cell(id);

    map<DIRECTION, int>::iterator it;
    for (it = cell.neighbour_ids_.begin(); it != cell.neighbour_ids_.end(); it++)
    {
      // double K_component = compute_K(cell, it->first);
      // double saturation = get_saturation(cell, it->first, old_p_solution, old_S_solution);
      // double mobility = get_mobility(saturation);

      double saturation = get_saturation(cell, it->first, old_p_solution, old_S_solution, time);
      double kappa = get_kappa(saturation);

      // if (saturation > 0.9)
      // {
      //   // We can have a test here. When will saturation be greater than 0.9? Must be left boundary.
      //   cout << "(" << cell.face_centers_[it->first][0] << " , " << cell.face_centers_[it->first][1] << ")" << endl << endl; 
      // }

      // // Face contribution
      // if (it->second == OUTSIDE_CELL_ID)
      // {
      //   double p_boundary_value = get_boundary_value_p(cell, it->first);
      //   system_rhs(cell.id_) += 2*K_component*mobility*p_boundary_value;
      //   triplet_list.push_back(T(cell.id_, cell.id_, 2*K_component*mobility));
      // }
      // else
      // {
      //   triplet_list.push_back(T(cell.id_, cell.id_, 1*K_component*mobility));
      //   triplet_list.push_back(T(cell.id_, it->second, -1*K_component*mobility));
      // }



      // Face contribution
      if (it->second == OUTSIDE_CELL_ID)
      {
        double p_boundary_value = get_boundary_value_p(cell, it->first, time);
        system_rhs(cell.id_) += 2*kappa*p_boundary_value;
        triplet_list.push_back(T(cell.id_, cell.id_, 2*kappa));
      }
      else
      {
        triplet_list.push_back(T(cell.id_, cell.id_, 1*kappa));
        triplet_list.push_back(T(cell.id_, it->second, -1*kappa));
      }




    }

    // Volume contribution
    system_rhs(cell.id_) += 0*Area;
  }

  // This is a very efficient way of assembling the matrix A. The function is provided by the external library Eigen.
  system_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());

}


void ProblemFlow::solve_p()
{

  // SimplicialCholesky<SparseMatrix<double>> chol(A);  // performs a Cholesky factorization of A
  // VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side

  SparseLU<SparseMatrix<double>>   solver;
  solver.analyzePattern(system_matrix); 
  solver.factorize(system_matrix); 
  new_p_solution = solver.solve(system_rhs); 

  // cout << MatrixXd(system_matrix) << endl;
  // cout << new_p_solution << endl;
  // cout << system_rhs << endl;
  // cout << new_p_solution.maxCoeff() << endl;

  old_p_solution = new_p_solution;

}


void ProblemFlow::compute_S()
{

  new_S_solution.setZero();

  for (int id = 0; id < CELL_NUMBER; ++id)
  {
    Cell cell(id);
    map<DIRECTION, int>::iterator it;
    for (it = cell.neighbour_ids_.begin(); it != cell.neighbour_ids_.end(); it++)
    {

      // double K_component = compute_K(cell, it->first);
      // double saturation = get_saturation(cell, it->first, old_p_solution, old_S_solution);
      // double mobility = get_mobility(saturation);
      // double fraction = get_fraction(saturation);

      double wind;

      double saturation = get_saturation(cell, it->first, old_p_solution, old_S_solution, time);
      double kappa = get_kappa(saturation);


      bool is_upwind = determine_upwind(cell, it->first, old_p_solution, time);

      // Face contribution
      if (it->second == OUTSIDE_CELL_ID)
        wind = (get_boundary_value_p(cell, it->first, time) - old_p_solution(cell.id_))/(H/2.);
      else
        wind = (old_p_solution(it->second) - old_p_solution(cell.id_))/H;

      // new_S_solution(cell.id_) += fraction*mobility*K_component*wind*H*DELTA_T/Area; 

      new_S_solution(cell.id_) += saturation*kappa*wind*H*DELTA_T/Area; 

    }

    // Volume contribution
    new_S_solution(cell.id_) += old_S_solution(cell.id_);


  }
  old_S_solution = new_S_solution;

}


void ProblemFlow::output_results(int cycle)
{


  VectorXd exact_p_solution(CELL_NUMBER);
  interpolate_p(exact_p_solution, time);

  string p_name = "p";
  string S_name = "S";
  OutputManager output_manager(cycle);

  output_manager.scalar_output(old_p_solution, p_name);
  output_manager.scalar_output(new_S_solution, S_name);  
  output_manager.scalar_output(exact_p_solution, "p_exact");
 
}



void ProblemFlow::post_processing()
{
  double error = 0;
  for (int i = 0; i < CELL_NUMBER; ++i)
  {
    Cell cell(i);
    error += pow(boundary_function_S(cell.cell_center_, time) - new_S_solution(i), 2)*Area; 
  }
  error = sqrt(error);

  // L2 norm should convege in first order w.r.t. mesh size
  cout << "L2 norm error is " << error << endl;
}

