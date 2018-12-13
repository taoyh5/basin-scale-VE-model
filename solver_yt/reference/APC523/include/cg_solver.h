#ifndef CG_SOLVER_H
#define CG_SOLVER_H

#include "common_definitions.h"
#include "operations.h"

/* The solver uses the conjugate gradient method to solve Ax=b and returns the value of x.
   It is the user's responsibility to check whether this matrix is symmetric. */
vec cg_solver(matrix &A,  vec &b);

#endif