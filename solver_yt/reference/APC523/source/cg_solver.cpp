#include "cg_solver.h"

vec cg_solver(matrix &A, vec &b)
{
   int n = A.size();
   vec x(n, 0);

   vec r = b;
   vec p = r;
   int k = 0;

   cout << "Start solving the system with size to be " << n << endl;

   while (k < n)
   {
      // Store previous residual
      vec r_old = r;                                         
      vec Ap = matrix_mul_vector(A, p);

      double alpha = vec_mul_vec(r, r)/max(vec_mul_vec(p, Ap), NEARZERO);
      // Next estimate of solution
      x = vector_combination(1.0, x, alpha, p);
      // Residual            
      r = vector_combination(1.0, r, -alpha, Ap);           
      // Convergence test
      if (norm(r) < TOLERANCE) break;             

      double beta = vec_mul_vec(r, r)/max(vec_mul_vec(r_old, r_old), NEARZERO);
      // Next gradient
      p = vector_combination(1.0, r, beta, p);             
      k++;
   }

   return x;
}
