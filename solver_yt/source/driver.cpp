//
//
//                       _oo0oo_
//                      o8888888o
//                      88" . "88      |-------------------------------------|
//                      (| -_- |)  --> | You shall have no bug in this code! |
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
//              Let's pray for no bug...
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*----------------------------------------------------------------------------------/
| @Authors: Tianju Xue(*), Yiheng Tao and Jiarong Wu                                |
| @Email: txue@princeton.edu(*)                                                     |
| @Last modified: 12/03/2018                                                        |
| @Brief: APC524 Final Project(A multi-phase flow solver)                           |
| @Description: To do                                                               |                            
/----------------------------------------------------------------------------------*/

/*

A few remarks about the comments for this code.

1. We have commeted both .cpp files and .h files. In .cpp files, we made in-text
comments while in .h files, we had a general comment for each individual function.

2. We use {} to denote an actual variable, function name, etc. 
eg. Suppose there is "double a = 1;" somewhere. We would use {a} to denote this
variable "a" in the actual code.

*/


#include "problem_poisson.h"
#include "problem_flow.h"
#include "cell.h"

int main()
{

  // ProblemPoisson poisson;
  // poisson.run();
  
  ProblemFlow two_phase_flow;
  two_phase_flow.run();
  
  return 0;

}