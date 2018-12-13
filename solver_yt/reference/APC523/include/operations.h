#ifndef OPERATIONS_H
#define OPERATIONS_H

#include "common_definitions.h"

/* Compute Ax */
vec matrix_mul_vector(matrix &A, vec &x);

/* Compute au + bv */
vec vector_combination( double a, vec &u, double b, vec &v );

/* Compute <x,y> */
double vec_mul_vec(vec &x, vec &y);

/* Compute ||x|| */
double norm(vec &x);

/* Compute det(A) */
double determinant(matrix &A);

/* Compute A^{-1} */
matrix inverse(matrix &B);

/* Compute A^{T} */
matrix transpose(matrix &B);

#endif