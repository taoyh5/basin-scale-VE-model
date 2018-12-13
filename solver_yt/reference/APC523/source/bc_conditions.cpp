#include "bc_conditions.h"

double Dirichlet_bc_value(vec &x)
{
	double r = sqrt(pow(x[0], 2) + pow(x[1], 2));
	double value = cos(M_PI*r/2.);
	return value;
}

double Neumann_bc_value(vec &x)
{	
	// to be implemented if needed
	return 0;

}


double f_value(vec &x)
{
	double r = sqrt(pow(x[0], 2) + pow(x[1], 2));	
	double value = M_PI/2.*(1./r*sin(M_PI*r/2.) + M_PI/2.*cos(M_PI*r/2.));
	return value;
}


double u_exact(vec &x)
{
	double r = sqrt(pow(x[0], 2) + pow(x[1], 2));
	double value = cos(M_PI*r/2.);
	return value;
}