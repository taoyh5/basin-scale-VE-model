#include "tools.h"

// Rubbish code, change later


void interpolate(VectorXd &solution /* To do: There should be a function pointer as another argument here */)
{

  for (int id = 0; id < solution.size(); ++id)
  {
    Cell cell(id);
    // To fill up
  }

}


void interpolate_Sc(VectorXd &solution, double time /* This function is a rubbish version. Got to be deleted and replaced. */)
{
  for (int i = 0; i < solution.size(); ++i)
  {
    Cell cell(i);
    solution(i) = boundary_function_Sc(cell.cell_center_, time);
  } 
}



void interpolate_Pc(VectorXd &solution, double time/* This function is a rubbish version. Got to be deleted and replaced. */)
{
  for (int i = 0; i < solution.size(); ++i)
  {
    Cell cell(i);
    solution(i) = boundary_function_Pc(cell.cell_center_, time);
  } 
}



double compute_K(Cell &cell, DIRECTION direction)
{
  // 2/(1/K + 1/K)
  return 1e-12*formation_thickness;  

  
  /*  if (direction == WEST || direction == EAST)
    {
      return 1;   
    }
  else
    {
      return 1;
      }*/

  // if (direction == WEST || direction == EAST)
  // {
  //   srand(cell.id_);
  //   return (rand()%RONDOM_N + 1)/(double)RONDOM_N;   
  // }
  // else
  // {
  //   srand(cell.id_ + CELL_NUMBER);
  //   return (rand()%RONDOM_N + 1)/(double)RONDOM_N;
  // }
}




// Garbage code
bool determine_upwind(Cell &cell, DIRECTION direction, VectorXd &Pc_solution, double time)
{
  double Pc_self = Pc_solution(cell.id_);

  double Pc_other;
  if (cell.neighbour_ids_[direction] == OUTSIDE_CELL_ID)
    Pc_other = get_boundary_value_Pc(cell, direction, time);
  else
    Pc_other = Pc_solution(cell.neighbour_ids_[direction]);

  return Pc_self - Pc_other >= 0 ? false : true; 
}



double get_Sc(Cell &cell, DIRECTION direction, VectorXd &Pc_solution, VectorXd &Sc_solution, double time)
{
  bool is_upwind = determine_upwind(cell, direction, Pc_solution, time);

  if (is_upwind)
  {
    return Sc_solution(cell.id_);
  }
  else
  {
    if (cell.neighbour_ids_[direction] == OUTSIDE_CELL_ID)
      return get_boundary_value_Sc(cell, direction, time);
    else  
      return Sc_solution(cell.neighbour_ids_[direction]);
  }
}



double get_boundary_value_Pb(Cell &cell, DIRECTION direction, double time)
{
  assert(("Bug! This face should be a boundary face but it is not!", cell.neighbour_ids_[direction] == OUTSIDE_CELL_ID));
  return boundary_function_Pb(cell.face_centers_[direction], time);
}


double get_boundary_value_Sc(Cell &cell, DIRECTION direction, double time)
{
  assert(("Bug! This face should be a boundary face but it is not!", cell.neighbour_ids_[direction] == OUTSIDE_CELL_ID));
  return boundary_function_Sc(cell.face_centers_[direction], time);
}



/*double get_boundary_value_x(Cell &cell, DIRECTION direction)
{
  assert(("Bug! This face should be a boundary face but it is not!", cell.neighbour_ids_[direction] == OUTSIDE_CELL_ID));
  return boundary_function_x(cell.face_centers_[direction]);
  }*/



double boundary_function_Pb(vec &point, double time)
{
  return 0.2/M_PI*cos(M_PI*(point[0] + point[1] - 2*time)) + 0.5*(point[0] + point[1]);
}



double boundary_function_Sc(vec &point, double time)
{
  return sin(M_PI*(point[0] + point[1] - 2*time));
}



/*double boundary_function_x(vec &point)
{
  return sin(2*M_PI*(point[0] + point[1]));
  }*/

double rhs_function(vec &point)
{
  // Mannualy set the rhs function so that the home made analytical solution is satisfied in the poisson equation
  return 8*pow(M_PI, 2)*sin(2*M_PI*(point[0] + point[1]));
}

