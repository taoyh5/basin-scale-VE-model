//
//
//                       _oo0oo_
//                      o8888888o
//                      88" . "88      |-------------------------------------|
//                      (| -_- |)  --> | You shall have no bug in this code. |
//                      0\  =  /0      |-------------------------------------|
//                    ___/`---'\___
//                  .' \\|     |// '.
//                 / \\|||  :  |||// \
//                / _||||| -:- |||||- \
//               |   | \\\  -  /// |   |
//               | \_|  ''\---/''  |_/ |
//               \  .-\__  '-'  ___/-. /
//             ___'. .'  /--.--\  `. .'___
//          ."" '<  `.___\_<|>_/___.' >' "".
//         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
//         \  \ `_.   \_ __\ /__ _/   .-` /  /
//     =====`-.____`.___ \_____/___.-`___.-'=====
//                       `=---='
//
//
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                Praying for no bug
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*----------------------------------------------------------------------------------/
| @Author: Tianju Xue														        |
| @Email: txue@princeton.edu										    	  	    |
| @Last modified: 04/07/2018														|
| @Brief: APC523 Final Project(A 2D DG-FEM solver for Poisson equation) 	  	    |
| @Description: This is a discontinuous Galerkin finite element solver written in   |
|				a generic way to solve the Poisson equation. The solver utilizes    |
|				the classic symmetric interior penalty Galerkin (SIPG) method.      |
|				The mesh can be refined uniformly and the order can be as high if   |
|				the quadrature points are given. The solver uses conjugate gradient |
|				method to solve the linear system obtained since our system matrix  |
|				A is symmetric. The output solution is in VTK format.               |													      |
/----------------------------------------------------------------------------------*/

/*

A few remarks about the comments for this code.

1. I have commeted both .cpp files and .h files. In .cpp files, I made in-text
comments while in .h files, I had a general comment for each individual function.

2. I use {} to denote an actual variable, function name, etc. 
eg. Suppose I there is "double a = 1;" somewhere. I would use {a} to denote this
variable "a" in the actual code.

*/

#include "common_definitions.h"
#include "problem.h"


#include <Eigen/Dense>

using namespace std;

int main()
{
	Problem poisson;
	poisson.run();
	return 0;
}
