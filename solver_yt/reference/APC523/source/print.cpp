#include "print.h"

void print_vector(vec &x)
{
	cout << endl;
	cout << "Start of printing a vector " << endl;
	for (int i = 0; i < x.size(); ++i)
	{
		cout << x[i] << endl;
	}
	cout << "End of printing a vector " << endl;
	cout << endl;
}

void print_matrix(matrix &X)
{
	cout << endl;
	cout << "Start of printing a matrix " << endl;	
	for (int i = 0; i < X.size(); ++i){
		for (int j = 0; j < X[0].size(); ++j)
		{	
			cout << X[i][j] << "  ";
			// cout << fixed << setprecision(3) << X[i][j] << "  ";
		}
		cout << endl;
	}
	cout << "End of printing a matrix " << endl;
	cout << endl;	
}

void print( string title, vec &v )
{
   cout << title << '\n';

   int n = v.size();           
   for ( int i = 0; i < n; i++ )
   {
      double x = v[i];   if ( abs( x ) < NEARZERO ) x = 0.0;
      cout << x << '\t';
   }
   cout << '\n';
}


void print( string title, matrix &A )
{
   cout << title << '\n';
   
   // A is an m x n matrix
   int m = A.size(), n = A[0].size();                     
   for ( int i = 0; i < m; i++ )
   {
      for ( int j = 0; j < n; j++ )
      {
         double x = A[i][j];   if ( abs( x ) < NEARZERO ) x = 0.0;
         cout << x << '\t';
      }
      cout << '\n';
   }
}
