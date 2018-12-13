#ifndef COMMON_DEFINITIONS_H
#define COMMON_DEFINITIONS_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

/* Some common constant variables are defined in this file. */

/* The order of the base functions */
const int ORDER = 1;

/* How many times we want to refine the mesh.
   For each time, one square element is divided into four. */
const int DIVISION = 1;

/* We want to play safe to choose how many quadrature we want to use when integrating*/
const int QUAD_NUM = ORDER + 1;

/* The dimention of the problem 
   The code is not written in a dimensional independent way so keep this unchanged.
   A dimensional independent code can be finished by the help of C++ function template */
const int DIM = 2;

/* Total number of elements*/
const int ELE_NUMBER = pow(4, DIVISION);

/* How many elements along on line */
const int LINE_DIVISION = pow(2, DIVISION);

/* Degree of freedom for each element */
const int DOF_PER_ELE = pow(ORDER + 1, 2);

/* Total number of dof */
const int TOTAL_DOF = DOF_PER_ELE*ELE_NUMBER;

/* A square element has 4 faces */
const int FACE_PER_ELE = 4;

/* A square element has 4 vertices */
const int VERTICES_PER_ELE = 4;

/* Along one line, there are 2 vertices for a square element */
const int VERTICES_PER_LINE = 2;

/* lower and upper limites of the lagrangian polynomial domain (aka, the parametric domain) */ 
const double XI_1 = -1;
const double XI_2 = 1;

/* {X_1} is the lower left point while {X_2} is the upper right point of the physical domain  */
const double X_1[DIM] = {-1, -1};
const double X_2[DIM] = {1, 1};

/* Mesh size */
const double H = (X_2[0] - X_1[0])/LINE_DIVISION;

/* Penalty factor */
const double GAMMA = 100;

/* Maximum order that one can use.
   To change this number, you need to add more Gauss quadrature points in the table.
   Warning: When the order is too high, certain issues may occur. 
   A uniformly spaced lagragian polynomial may not be suitable for constructing the base functions. */
const int MAX_ORDER = 5;

/* Used in the linear solver */
const double NEARZERO = 1.0e-10;
const double TOLERANCE = 1.0e-10;

/* Used in outputing VTK files */
const int CELL_TYPE = 9;

/* The quadrature table */
const double QUADRATURE_TABLE[MAX_ORDER + 1][MAX_ORDER + 1] = 
{
	{0, 0, 0, 0, 0 ,0},
	{-0.5773502691896257, 0.5773502691896257, 0, 0, 0, 0},
	{0.0000000000000000, -0.7745966692414834, 0.7745966692414834},
	{-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526, 0, 0},
	{0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640, 0},
	{0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521}
};

/* The weight table */
const double WEIGHT_TABLE[MAX_ORDER + 1][MAX_ORDER + 1] =
{
	{0, 0, 0, 0, 0, 0},
	{1.0000000000000000, 1.0000000000000000, 0, 0, 0, 0},
	{0.8888888888888888, 0.5555555555555556, 0.5555555555555556, 0, 0, 0},
	{0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538, 0, 0},
	{0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891,0.2369268850561891, 0},
	{0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704}
};

/* Pre-define some data type. Just some aliases. */
typedef vector<double> vec;
typedef vector<vec> matrix;
typedef vector<int> int_vec;

#endif