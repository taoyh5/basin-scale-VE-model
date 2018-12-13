#ifndef FE_VALUE_H
#define FE_VALUE_H

#include "common_definitions.h"
#include "element.h"
#include "operations.h"

/* The class {FeValue} holds many important values that we might need when we sit on some certain
   finite element cell. Each element will have one corresponding {FeValue} object as its support. 

   The way we enuerate degree of freedoms can be shown in the following diagram ({ORDER} = 2)

	6--7--8
	|     |
	3  4  5
	|     |
	0--1--2

*/
class FeValue
{
	public:
		/* {xi} is the quadrature point that this {FeValue} object is focusing on
		   Every computations in this {FeValue} object is based on this quadrature point */
    	vec xi;	

    	/* Constructor. 
    	   @Input: {ele_ptr_} as a pointer to the element because we will need information from the element*/
    	FeValue(Element *ele_ptr_);

    	/* @Input: the quadrature point {xi} and the physical coordinates list {x_list}
    	   @Output: the physical coordinate corresponding to {xi} */
		vec get_physical_crds(vec &xi, matrix &x_list);

		/* Evaluate the {i}th shape function value at the point {xi}. */
	   	double shape_value(int i);

	   	/* Evaluate the {i}th shape function gradient at the point {xi}. */
    	vec shape_gradient(int i);

    	/* Initialize the value of {xi}. */
    	virtual void get_quadrature(int q);

    	/* {JxW} represents the jacobian times the weight as the name indicates.
    		An interesting fact is that the jacobian always comes along with the
    		weight of the quadrature point, so we can always treat them as a whole. */
  		virtual double JxW();
	
	protected:

		/* {ele_ptr_} is the pointer to this {Element} object.
			It is declared as {protected} because its derived class also needs this. */
		Element *ele_ptr;

		/* The deformation gradient (often denoted as F) is useful for computing the jacobian */
		matrix get_deformation_gradient(vec &xi, matrix &x_list);

	private:

		/* indices for quadrature point */
		int q_index_xi_1;
		int q_index_xi_2;

		/* We need to know which quadrature point to use.  */
    	void get_q_index(int q);

    	/* The following two functions {lagrange_basis} and are {d_lagrange_basis} are two
    	   fundamental functions. It uses Lagrange polynomials and returns the value evaluated
    	   at {x} using the {i}th base function (in 1D) */
		double lagrange_basis(vec &nodes, int index, double x);
		double d_lagrange_basis(vec &nodes, int index, double x);

		/* Generate {length} uniformly spaced points starting from {lower_limit} and ends at 
		   {upper_limit}*/
		vec linspace(int length, double lower_limit, double upper_limit);


		/* The following two functions {base} and {grad_base} takes {xi} but returns a value
		   that is actually in the physical domain.
		   Notice that in {grad_base}, the gradient is with respect to the physical coordinates,
		   not {xi}! */
		double base(int base_index, vec &xi);
		vec grad_base(int base_index, vec &xi);

		/* Get the jacobian. */
		double get_jacobian(matrix &F);
};

/* The class {FeFaceValue} is derived from {FeValue}. This class holds the information that is necessary 
   to know for each face. */
class FeFaceValue: public FeValue 
{
	public:
		/* Constructor. */
		FeFaceValue(Element *ele_ptr_, ElementFace *eleface_ptr_);

		/* Rewrite the virtual function {get_quadrature} defined in its base class. */
		void get_quadrature(int q);

		/* Use Nanson's formula to compute the JxW term. This also rewrites the virtual function {JxW}
		   defined in its base class. */
		double JxW();	
				
	private:
		/* The corresponding {eleface_ptr} pointing to the face. */
		ElementFace *eleface_ptr;

		/* We need to know which quadrature point to use.  */
		int q_index_xi;
};

#endif