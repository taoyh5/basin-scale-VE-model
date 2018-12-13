#ifndef BC_CONDITIONS_H
#define BC_CONDITIONS_H

#include "common_definitions.h"

/* Returns essential boundary conditions given a point {x} on boundary */
double Dirichlet_bc_value(vec &x);

/* Returns natural boundary conditions given a point {x} on boundary */
double Neumann_bc_value(vec &x);

/* We also put right hand side function f here */
double f_value(vec &x);

/* We also put exact solution here */
double u_exact(vec &x);

#endif