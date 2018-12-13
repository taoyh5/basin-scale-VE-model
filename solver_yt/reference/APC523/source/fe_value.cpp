#include "fe_value.h"

//**************************************************************************************
// {public} functions of {FeValue} start from here

FeValue::FeValue(Element *ele_ptr_)
{
    ele_ptr = ele_ptr_;
}

vec FeValue::get_physical_crds(vec &xi, matrix &x_list)
{
	vec x(DIM, 0);
	for (int index = 0; index < DOF_PER_ELE; ++index)
	{
		for (int i = 0; i < DIM; ++i)
		{
			x[i] += x_list[index][i]*base(index, xi);
		}
	}
	return x;	
}

double FeValue::shape_value(int i)
{
	return base(i, xi);
}

// Pls pay attention to how I use chain rule to get the gradient w.r.t. physical
// coordinates rather that parametric coordinates.
vec FeValue::shape_gradient(int i)
{
	matrix F(DIM, vec(DIM, 0));
	vec gradient(DIM, 0);
	F = get_deformation_gradient(xi, ele_ptr->x_list);

	matrix matrix_temp1(DIM, vec(DIM, 0));
	matrix matrix_temp2(DIM, vec(DIM, 0));
	vec vec_temp1(DIM, 0);

	matrix_temp1 = inverse(F);
	matrix_temp2 = transpose(matrix_temp1);

	vec_temp1 = grad_base(i, xi);

	gradient = matrix_mul_vector(matrix_temp2 , vec_temp1);
	return gradient;
}

void FeValue::get_quadrature(int q)
{
	xi.clear();
	get_q_index(q);
	xi.push_back(QUADRATURE_TABLE[QUAD_NUM - 1][q_index_xi_1]);
	xi.push_back(QUADRATURE_TABLE[QUAD_NUM - 1][q_index_xi_2]);
}

double FeValue::JxW()
{
	matrix F(DIM, vec(DIM, 0));
	F = get_deformation_gradient(xi, ele_ptr->x_list);
	double J = determinant(F);
	double W1 = WEIGHT_TABLE[QUAD_NUM - 1][q_index_xi_1];
	double W2 = WEIGHT_TABLE[QUAD_NUM - 1][q_index_xi_2];
	return J*W1*W2;
}

// {public} functions of {FeValue} end here
//**************************************************************************************




//**************************************************************************************
// {protected} functions of {FeValue} start from here

matrix FeValue::get_deformation_gradient(vec &xi, matrix &x_list)
{
	matrix F(DIM, vec(DIM, 0));
	vec grad_phi(DIM, 0); 
	vec x(DIM, 0); 

	for (int index = 0; index < DOF_PER_ELE; ++index)
	{
		for (int i = 0; i < DIM; ++i)
		{
			x[i] = x_list[index][i];
		}

		grad_phi = grad_base(index, xi);
				
		for (int i = 0; i < DIM; ++i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				F[i][j] += x[i]*grad_phi[j];
			}
		}
	}

	return F;
}

// {protected} functions of {FeValue} end here
//**************************************************************************************




//**************************************************************************************
// {private} functions of {FeValue} start from here

void FeValue::get_q_index(int q)
{	
	q_index_xi_1 = q%(QUAD_NUM);
	q_index_xi_2 = q/(QUAD_NUM);
}

double FeValue::lagrange_basis(vec &nodes, int index, double x)
{
	double value = 1.;
	for (int i = 0; i < ORDER + 1; ++i)
	{
		if (i == index) continue;
		value *= (x - nodes[i])/(nodes[index] - nodes[i]);
	}

	return value;
}

double FeValue::d_lagrange_basis(vec &nodes, int index, double x)
{
	double d_value = 0;

	for (int i = 0; i < ORDER + 1; ++i)
	{
		double d_value_temp = 1;
		if (i == index) continue;
		for (int j = 0; j < ORDER + 1; ++j)
		{
			if (j == index) continue;
			if (j == i){
				d_value_temp*=1./(nodes[index] - nodes[j]);
			}
			else
			{
				d_value_temp*=(x - nodes[j])/(nodes[index] - nodes[j]);
			} 
		}
		d_value+=d_value_temp;
	}
	return d_value;
}

vec FeValue::linspace(int length, double lower_limit, double upper_limit)
{

	vec array(length); 
	double step = (upper_limit - lower_limit)/(length - 1);
	for (int i = 0; i < length; ++i)
	{
		array[i] = lower_limit + step*i;
	}
	return array; 
}

// Input: base index, xi, eta
// Output: base value
double FeValue::base(int base_index, vec &xi)
{
	vec param_one_dim(ORDER + 1);

	param_one_dim = linspace(param_one_dim.size(), XI_1, XI_2);
	int base_index_xi_1 = base_index%(ORDER + 1);
	int base_index_xi_2 = base_index/(ORDER + 1);

	return lagrange_basis(param_one_dim, base_index_xi_1, xi[0])*
		   lagrange_basis(param_one_dim, base_index_xi_2, xi[1]);

}

vec FeValue::grad_base(int base_index, vec &xi)
{

	vec grad_base_value(DIM); 

	vec param_one_dim(ORDER + 1);
	param_one_dim = linspace(param_one_dim.size(), XI_1, XI_2);
	int base_index_xi_1 = base_index%(ORDER + 1);
	int base_index_xi_2 = base_index/(ORDER + 1);
	grad_base_value[0] = d_lagrange_basis(param_one_dim, base_index_xi_1, xi[0])*
						 lagrange_basis(param_one_dim, base_index_xi_2, xi[1]);
	grad_base_value[1] = lagrange_basis(param_one_dim, base_index_xi_1, xi[0])*
						 d_lagrange_basis(param_one_dim, base_index_xi_2, xi[1]);

	return grad_base_value;

}

double FeValue::get_jacobian(matrix &F)
{
	return determinant(F);
}

// {private} functions of {FeValue} end here
//**************************************************************************************




//**************************************************************************************
// {public} functions of {FeFaceValue} start from here

FeFaceValue::FeFaceValue(Element *ele_ptr_, ElementFace *eleface_ptr_):FeValue(ele_ptr_)
{
    ele_ptr = ele_ptr_;
    eleface_ptr = eleface_ptr_;

}

void FeFaceValue::get_quadrature(int q)
{
	// Quadrature points on a face is essentially a 1-D treatment
	xi.clear();
	q_index_xi = q;
	switch(eleface_ptr->face_id)
    {
        case 0:
        {
			xi.push_back(-1);
			xi.push_back(QUADRATURE_TABLE[QUAD_NUM - 1][q]);
            break;              
        }
        case 1:
        {
			xi.push_back(1);
			xi.push_back(QUADRATURE_TABLE[QUAD_NUM - 1][q]);
            break;              
        }
        case 2:
        {
			xi.push_back(QUADRATURE_TABLE[QUAD_NUM - 1][q]);
			xi.push_back(-1);
            break;              
        }
        case 3:
        {
			xi.push_back(QUADRATURE_TABLE[QUAD_NUM - 1][q]);
			xi.push_back(1);
            break;              
        }      
    }
}

double FeFaceValue::JxW()
{
	matrix F(DIM, vec(DIM, 0));
	F = get_deformation_gradient(xi, ele_ptr->x_list);

	double J = determinant(F);

	matrix matrix_temp1(DIM, vec(DIM, 0));
	matrix matrix_temp2(DIM, vec(DIM, 0));
	vec vec_temp1(DIM, 0);
	vec vec_temp2(DIM, 0);

	matrix_temp1 = inverse(F);
	matrix_temp2 = transpose(matrix_temp1);

	vec_temp1 = eleface_ptr->normal_vector;
	vec_temp2 = matrix_mul_vector(matrix_temp2 , vec_temp1);

	double temp = norm(vec_temp2);
	double W = WEIGHT_TABLE[QUAD_NUM - 1][q_index_xi];
	return temp*J*W;
}

// {public} functions of {FeFaceValue} end here
//**************************************************************************************
