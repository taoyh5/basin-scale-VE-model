#ifndef PRINT_H
#define PRINT_H

#include "common_definitions.h"

/* Print to see if anything wrong. Only used for debugging */
void print_vector(vec &x);

/* Print to see if anything wrong. Only used for debugging */
void print_matrix(matrix &X);

/* Format printing for a vector */
void print(string title, vec &v);

/* Format printing for a matrix */
void print(string title, matrix &A);

#endif