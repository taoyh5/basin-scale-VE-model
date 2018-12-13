#include "operations.h"

vec matrix_mul_vector(matrix &A, vec &x)
{
	vec y(x.size(), 0);
	for (int i = 0; i < x.size(); ++i)
	{
		y[i] = vec_mul_vec(A[i], x);
	}
	return y;
}

vec vector_combination( double a, vec &u, double b, vec &v )        
{
   int n = u.size();
   vec w( n );
   for ( int j = 0; j < n; j++ ) w[j] = a * u[j] + b * v[j];
   return w;
}


double vec_mul_vec(vec &x, vec &y)
{
	return inner_product( x.begin(), x.end(), y.begin(), 0.0 );
}

double norm(vec &x)
{
	return sqrt(vec_mul_vec(x, x));
}

double determinant(matrix &A)
{
	return (A[0][0]*A[1][1]-A[0][1]*A[1][0]);
}

matrix inverse(matrix &B)
{

	matrix A(DIM, vec(DIM, 0));
	
	double j = determinant(B);

	if (j == 0)
	{
		cout << "Warning: determinant is zero!" << endl;
	}

	A[0][0] = B[1][1]/j;
	A[0][1] = -B[0][1]/j;
	A[1][0] = -B[1][0]/j;
	A[1][1] = B[0][0]/j;

	return A;
}

matrix transpose(matrix &B)
{
	matrix A(DIM, vec(DIM, 0));

	A[0][0] = B[0][0];
	A[0][1] = B[1][0];
	A[1][0] = B[0][1];
	A[1][1] = B[1][1];

	return A;
}
