#include <iostream>
#include <fstream>
using namespace std;


#include "getfem/getfem_assembling.h" /* import assembly methods      */
#include "getfem/getfem_export.h"     /* export functions             */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_import.h"
#include "gmm/gmm.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_vector;
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::short_type;  /* = short */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

typedef gmm::row_matrix<sparse_vector> sparse_row_matrix;


//#include "abc_problem.h"


double a_init = 0.0;
double b_init = 0.0;
double a = 0.0;
double b = 0.0;

size_type L = 0;

double alphax = 0.0;
double betax = 0.0;
double nux = 0.0;

double eps = 0.0;

int TYPE_REG = 1;
int TYPE_ACT_FUNC = 1;


size_type nb_dof = 0.0;





/***************************************************************************************************************************/
/******************************************************* THE NOISE *********************************************************/
/***************************************************************************************************************************/

double noise1(const base_node p) {

	double espmu = 0.5;
	double sigmax = 0.2;

	return 0.05*exp(-1.0*pow(p[0]-espmu, 2.0)/(2.0*sigmax*sigmax))/sqrt(2.0*M_PI*sigmax*sigmax);

}


double noise2(const base_node p) {

	double espmu = 0.5;
	double sigmax = 0.2;

	return 0.05*exp(-1.0*pow(p[0]-espmu, 2.0)/(2.0*sigmax*sigmax))/sqrt(2.0*M_PI*sigmax*sigmax);

}



/***************************************************************************************************************************/
/*********************************************** ROUTINES FOR THE SUBPROBLEM ***********************************************/
/***************************************************************************************************************************/


void matlab_export_misfit(const double misfit) {

  std::ofstream expmisfit;
  expmisfit.open("./MATLAB/Misfit.txt", ios::out|ios::app);
  expmisfit << misfit << ",   " ;
  expmisfit.close();

}


void matlab_export_coeff1(const double coeff) {

  std::ofstream evolcoeff;
  evolcoeff.open("./MATLAB/Coeffdata1.txt", ios::out|ios::app);
  evolcoeff << coeff << ",   " ;
  evolcoeff.close();

}

void matlab_export_coeff2(const double coeff) {

  std::ofstream evolcoeff;
  evolcoeff.open("./MATLAB/Coeffdata2.txt", ios::out|ios::app);
  evolcoeff << coeff << ",   " ;
  evolcoeff.close();

}

void matlab_export_coeffbis1(const double coeff) {

  std::ofstream evolcoeff;
  evolcoeff.open("./MATLAB/Coeffdata_validation1.txt", ios::out|ios::app);
  evolcoeff << coeff << ",   " ;
  evolcoeff.close();

}

void matlab_export_coeffbis2(const double coeff) {

  std::ofstream evolcoeff;
  evolcoeff.open("./MATLAB/Coeffdata_validation2.txt", ios::out|ios::app);
  evolcoeff << coeff << ",   " ;
  evolcoeff.close();

}


void matlab_export_coefflearnedNN1(const double coeff) {

  std::ofstream evolcoeff;
  evolcoeff.open("./MATLAB/CoefflearnedNN1.txt", ios::out|ios::app);
  evolcoeff << coeff << ",   " ;
  evolcoeff.close();

}

void matlab_export_coefflearnedNN2(const double coeff) {

  std::ofstream evolcoeff;
  evolcoeff.open("./MATLAB/CoefflearnedNN2.txt", ios::out|ios::app);
  evolcoeff << coeff << ",   " ;
  evolcoeff.close();

}


void matlab_export_supertol(const double Normgrad) {

  std::ofstream supertol;
  supertol.open("./MATLAB/Supertol.txt", ios::out|ios::app);

  supertol << Normgrad << ", " ;

  supertol.close();

}


void matlab_export_weights1(const std::vector<double> weights) {

  std::ofstream finalweights;
  finalweights.open("./MATLAB/Weights1.txt", ios::out|ios::app);

  for (size_type i=0; i<weights.size(); ++i) {

  	finalweights << weights[i] << ", " ;

  }

  finalweights.close();

}


void matlab_export_weights2(const std::vector<double> weights) {

  std::ofstream finalweights;
  finalweights.open("./MATLAB/Weights2.txt", ios::out|ios::app);

  for (size_type i=0; i<weights.size(); ++i) {

  	finalweights << weights[i] << ", " ;

  }

  finalweights.close();

}



void matlab_export_weights_trans1(const std::vector<double> weights) {

  std::ofstream finalweights;
  finalweights.open("./MATLAB/Weights_trans1.txt", ios::out|ios::app);

  for (size_type i=0; i<weights.size(); ++i) {

  	finalweights << weights[i] << ", " ;

  }

  finalweights.close();

}



void matlab_export_weights_trans2(const std::vector<double> weights) {

  std::ofstream finalweights;
  finalweights.open("./MATLAB/Weights_trans2.txt", ios::out|ios::app);

  for (size_type i=0; i<weights.size(); ++i) {

  	finalweights << weights[i] << ", " ;

  }

  finalweights.close();

}



void matlab_export_truth1(const std::vector<double> u) {

  std::ofstream blabla0;
  blabla0.open("./MATLAB/Coeff_truth1.txt", ios::out|ios::trunc);
  blabla0.close();

  std::ofstream blabla;
  blabla.open("./MATLAB/Coeff_truth1.txt", ios::out|ios::app);

  for (size_type i=0; i<u.size(); ++i) {

  	blabla << u[i] << ", " ;

  }

  blabla.close();

}



void matlab_export_truth2(const std::vector<double> u) {

  std::ofstream blabla0;
  blabla0.open("./MATLAB/Coeff_truth2.txt", ios::out|ios::trunc);
  blabla0.close();

  std::ofstream blabla;
  blabla.open("./MATLAB/Coeff_truth2.txt", ios::out|ios::app);

  for (size_type i=0; i<u.size(); ++i) {

  	blabla << u[i] << ", " ;

  }

  blabla.close();

}




double exact_sol(const base_node p) {

	double x = p[0];

	return x*x*sin(M_PI*x);


}




double exact_rhs_a(const base_node p) {

	double x = p[0];

	return b*x*x*sin(M_PI*x) - 1.0*a*( (2.0-M_PI*M_PI*x*x)*sin(M_PI*x) + 4.0*x*M_PI*cos(M_PI*x) );

}


void direct_solve(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &y, const double u1, const double u2) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, u2));
  MS.add_initialized_fixed_size_data("agiven", plain_vector(1, 1.0)); // Diffusion term: Laplace

  MS.add_fem_variable("y", mf);
  plain_vector F(mf.nb_dof());

  cout << "Assembling for problem-a ..." << endl;
  sparse_matrix M(mf.nb_dof(),mf.nb_dof());
  getfem::asm_mass_matrix(M, mim, mf, mf);
  gmm::scale(M, u1);
  getfem::add_explicit_matrix(MS, "y", "y", M);
  getfem::add_generic_elliptic_brick(MS, mim, "y", "a");
  getfem::interpolation_function(mf, F, exact_rhs_a);
  add_Dirichlet_condition_with_simplification(MS, "y", BOUNDARY);

  MS.add_initialized_fem_data("VolumicData", mf, F);
  getfem::add_source_term_brick(MS, mim, "y", "VolumicData");

  cout << "Solving..." << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  gmm::resize(y, mf.nb_dof());
  gmm::copy(MS.real_variable("y"), y);

}



void adjoint_solve(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &p,
                   const double u1, const double u2, const plain_vector y, const plain_vector ydata) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, u2));

  MS.add_fem_variable("p", mf);
  plain_vector F(nb_dof);

  cout << "Assembling adjoint system for problem-a ..." << endl;
  sparse_matrix M(mf.nb_dof(),mf.nb_dof());
  getfem::asm_mass_matrix(M, mim, mf, mf);
  gmm::scale(M, u1);
  getfem::add_explicit_matrix(MS, "p", "p", M);
  getfem::add_generic_elliptic_brick(MS, mim, "p", "a");
  add_Dirichlet_condition_with_simplification(MS, "p", BOUNDARY);

  gmm::copy(ydata, F);
  gmm::scale(F, -1.0);
  gmm::add(y, F);
  gmm::scale(F, -1.0);

  MS.add_initialized_fem_data("VolumicData", mf, F);
  getfem::add_source_term_brick(MS, mim, "p", "VolumicData");

  cout << "Solving the adjoint system..." << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  gmm::resize(p, nb_dof);
  gmm::copy(MS.real_variable("p"), p);

}



void compute_grad_control_L2(double &Gu1, double &Gu2, getfem::mesh &mesh, getfem::mesh_fem &mf, getfem::mesh_im &mim, const double u1, const double u2, const plain_vector ydata) {

  Gu1 = alphax*u1;
  Gu2 = alphax*u2;

  plain_vector y(nb_dof);
  direct_solve(mesh, mim, mf, y, u1, u2);
  plain_vector p(nb_dof);
  adjoint_solve(mesh, mim, mf, p, u1, u2, y, ydata);

  sparse_matrix M(nb_dof, nb_dof);
  gmm::clear(M);

  sparse_matrix Mass(nb_dof, nb_dof);
  gmm::clear(Mass);

  plain_vector prod(nb_dof);
  gmm::clear(prod);

  plain_vector Agiven(nb_dof);
  getfem::model MS;

  cout << "Assembling the gradient for problem-a ..." << endl;
  for (size_type i=0; i<nb_dof; ++i) {
	Agiven[i] = u2;
  }
  getfem::asm_stiffness_matrix_for_laplacian(M, mim, mf, mf, Agiven);
  gmm::mult(M, y, prod);
  Gu2 += gmm::vect_sp(prod, p);

  getfem::asm_mass_matrix(Mass, mim, mf, mf);
  gmm::scale(Mass, u1);
  gmm::clear(prod);
  gmm::mult(Mass, y , prod);
  Gu1 += gmm::vect_sp(prod, p);

  cout << "Gradient with the L2 reg is computed." << endl;

}







/***************************************************************************************************************************/
/********************************************* NEURAL NETWORKS AND SENSITIVITY *********************************************/
/***************************************************************************************************************************/


double func_max(const double x, const double y) {


  return 0.5*(x+y+gmm::abs(x-y));

}


double activation_func(const double p, const int ACT_FUNC) {


  double res = 0.0;
  double epsi = 1.0;

  switch (ACT_FUNC) {

	case 0:

		res = (exp(p) - exp(-1.0*p))/(exp(p)+exp(-1.0*p));
		break;

	case 1:

		res = func_max(0,p);
		break;

	case 2:

		res = log(1.0+exp(p)) - log(2.0);
		break;

	case 3:

		res = p*p*p;//func_max(0,p*p*p);
		break;

	case 4: 

		res = atan(p);
		break;

	case 5:

		if ((p>-1.0*epsi) && (p<epsi)) {
			res = ((p+epsi)*(p+epsi))/(4.0*epsi);
		}
		else if (p>epsi) {
			res = p;
		}
		else {
			res = 0.0;
		}
		break;

	default:

		cout << "Wrong choice of activation function!" << endl;

  }

  
  return res;

}



double activation_func_prime(const double p, const int ACT_FUNC) {


  double res = 0.0;
  double thx = 0.0;

  double epsi = 1.0;

  switch (ACT_FUNC) {

	case 0:

		thx = (exp(p) - exp(-1.0*p))/(exp(p)+exp(-1.0*p));
		res = 1.0 - thx*thx;
		//res = 1.0-activation_func(p, ACT_FUNC)*activation_func(p, ACT_FUNC);
		break;

	case 1:

		res = 1.0;
  		if (p<0.0) {
			res = 0.0;
  		}  
		break;

	case 2:

		res = 1.0/(1.0+exp(-p)) - 0.5;
		break;

	case 3:

		res = 3.0*p*p;
  		/*if (p<0.0) {
			res = 0.0;
  		} */ 
		break;

	case 4:

		res = 1.0/(1.0+p*p);
		break;

	case 5:

		if ((p>-1.0*epsi) && (p<epsi)) {
			res = (p+epsi)/(2.0*epsi);
		}
		else if (p>epsi) {
			res = 1.0;
		}
		else {
			res = 0.0;
		}
		break;

	default:

		cout << "Wrong choice of activation function!" << endl;

  }



  return res;

}


double activation_func_second(const double p, const int ACT_FUNC) {


  double res = 0.0;
  double thx = 0.0;

  double epsi = 1.0;

  switch (ACT_FUNC) {

	case 0:

		thx = (exp(p) - exp(-1.0*p))/(exp(p)+exp(-1.0*p));
		res = -2.0*thx*(1-thx*thx);
		//res = -2.0*activation_func(p, ACT_FUNC)*activation_func_prime(p, ACT_FUNC);
		break;

	case 1:

		res = 0.0;
  		if (p<0.0) {
			res = 0.0;
  		}  
		break;

	case 2:

		res = -1.0*(exp(-p)) / ( (1.0+exp(-p))*(1.0+exp(-p)) );
		break;

	case 3:

		res = 6.0*p;
  		/*if (p<0.0) {
			res = 0.0;
  		}  */
		break;

	case 4:

		res = -2.0*p/((1.0+p*p)*(1.0+p*p));
		break;

	case 5:

		if ((p>-1.0*epsi) && (p<epsi)) {
			res = 1.0/(2.0*epsi);
		}
		else if (p>epsi) {
			res = 0.0;
		}
		else {
			res = 0.0;
		}
		break;


	default:

		cout << "Wrong choice of activation function!" << endl;

  }



  return res;

}



double eval_NN(const size_type LL, const std::vector<double> weights, const std::vector<double> weights_trans, const double control, const int ACT_FUNC) {

  double res = weights[0]*control + weights_trans[0];

  for (size_type l=1; l<LL; ++l) {
	
	res = weights[l]*activation_func(res, ACT_FUNC) + weights_trans[l];

  }


  return res;

}




double reg_prime_u(const double u, const std::vector<double> weights, const std::vector<double> weights_trans, const size_type LL) {

	// !!!! Achtung !!! Everything works only in 1D


	// We assume that L <= weights.size();

	double prod = 1.0;
	if (LL > 1) {
	for (size_type l=0; l<LL-1; ++l) {

		//cout << "ll = " << l << endl;
		prod *= activation_func_prime(eval_NN(l, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC)*weights[l];

	}	
	}

	prod = weights[LL-1]*prod;

	return prod;

}



double reg_prime_w(const double u, const std::vector<double> weights, const std::vector<double> weights_trans, const size_type LL, const size_type s) {


	double prod = 1.0;

	if (LL > s) {

		for (size_type k=s; k<LL-1; ++k) {

			prod *= weights[k+1]*activation_func_prime(eval_NN(k, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);

		}

	}
		//res[l] *= weights[weights.size()-1];

	if (s>0) {

		prod *= activation_func(eval_NN(s-1, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);

	}


	return prod;

}



double reg_prime_w_trans(const double u, const std::vector<double> weights, const std::vector<double> weights_trans, const size_type LL, const size_type s) {


	double prod = 0.0;

	if (LL > s) {

		prod = 1.0;
		for (size_type k=s; k<LL-1; ++k) {

			prod *= weights[k+1]*activation_func_prime(eval_NN(k, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);

		}
	}

	return prod;

}




double reg_second_uu(const double u, const std::vector<double> weights, const std::vector<double> weights_trans) {

	// !!!! Achtung !!! Everything works only in 1D

	double sum = 0.0;
	double prod = 1.0;

	for (size_type l=0; l<weights.size()-1; ++l) {

		prod = 1.0;

		for (size_type k=0; k<l; ++k) {

			prod *= activation_func_prime(eval_NN(k, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
			prod *= weights[k];

		}

		prod *= activation_func_second(eval_NN(l, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
		prod *= reg_prime_u(u, weights, weights_trans, l);
		prod *= weights[l];

		for (size_type k=l+1; k<weights.size()-1; ++k) {

			prod *= activation_func_prime(eval_NN(k, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
			prod *= weights[k];

		}

		prod = weights[L-1]*prod;

		sum += prod;

	} 

	return sum;

}




double reg_second_uw(const double u, const std::vector<double> weights, const std::vector<double> weights_trans, const size_type s) {

	// !!!! Achtung !!! Everything works only in 1D

	double sum = 0.0;
	double prod = 1.0;

	for (size_type k=0; k<weights.size()-1; ++k) {

		prod = 1.0;

		if (k>0) {
		for (size_type l=0; l<k; ++l) {

			prod *= activation_func_prime(eval_NN(l, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
			prod *= weights[l];

		}
		}

		if (s > k) {

			prod = 0.0;

		}

		if (s <= k) {

			prod *= activation_func_second(eval_NN(k, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
			prod *= reg_prime_w(u, weights, weights_trans, k, s);
			prod *= weights[k];

		}

		for (size_type l=k+1; l<weights.size()-1; ++l) {

			prod *= activation_func_prime(eval_NN(l, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
			prod *= weights[l];

		}

		prod = weights[L-1]*prod;

		sum += prod;

	} 

	return sum;

}



double reg_second_uw_trans(const double u, const std::vector<double> weights, const std::vector<double> weights_trans, const size_type s) {

	// !!!! Achtung !!! Everything works only in 1D

	double sum = 0.0;
	double prod = 1.0;

	for (size_type k=0; k<weights.size()-1; ++k) {

		prod = 1.0;

		if (k>0) {
		for (size_type l=0; l<k; ++l) {

			prod *= activation_func_prime(eval_NN(l, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
			prod *= weights[l];

		}
		}

		if (s > k) {

			prod = 0.0;

		}

		if (s <= k) {

			prod *= activation_func_second(eval_NN(k, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
			prod *= reg_prime_w_trans(u, weights, weights_trans, k, s);
			prod *= weights[k];

		}

		for (size_type l=k+1; l<weights.size()-1; ++l) {

			prod *= activation_func_prime(eval_NN(l, weights, weights_trans, u, TYPE_ACT_FUNC), TYPE_ACT_FUNC);
			prod *= weights[l];

		}

		prod = weights[L-1]*prod;

		sum += prod;

	} 

	return sum;

}



double grad_penalization_term_constant_NN(const std::vector<double> weights, const std::vector<double> weights_trans, const size_type l) {

  double sum = 0.0;

  for(double u=1.0; u<1.6; u+=0.01) {

	sum += log(abs(reg_prime_u(u, weights, weights_trans, weights.size())))*reg_second_uw(u, weights, weights_trans, l)/reg_prime_u(u, weights, weights_trans, weights.size());

  }

  sum *= betax*0.01;

  return sum;

} 




/***************************************************************************************************************************/
/********************************************** ROUTINES FOR THE MAIN PROBLEM **********************************************/
/***************************************************************************************************************************/


double mynorm(const std::vector<double> Supergrad) {

  double sum = 0.0;

  for (size_type k=0; k<Supergrad.size(); ++k) {

	sum += std::abs(Supergrad[k]);

  }

  return sum/double(Supergrad.size());

}



void compute_misfit(double &misfit, const std::vector<double> u1, const std::vector<double> u2, const std::vector<double> udata1, const std::vector<double> udata2) {

  std::vector<double> diff1(udata1.size());
  for (size_type k=0; k<udata1.size(); ++k) {

	diff1[k] = u1[k] - udata1[k];

  }
  std::vector<double> diff2(udata2.size());
  for (size_type k=0; k<udata2.size(); ++k) {

	diff2[k] = u2[k] - udata2[k];

  }
  double ratio = (gmm::vect_norm2(diff1)+gmm::vect_norm2(diff2))/(gmm::vect_norm2(udata1)+gmm::vect_norm2(udata2));

  misfit = 100.0*ratio;

}



void direct_solve_prime1(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &y, const double u1, const double u2) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  plain_vector y0(mf.nb_dof());
  direct_solve(mesh, mim, mf, y0, u1, u2);

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, u2));

  MS.add_fem_variable("y", mf);
  plain_vector F(mf.nb_dof());

  cout << "Assembling for problem-a ..." << endl;
  
  sparse_matrix M(mf.nb_dof(), mf.nb_dof());
  getfem::asm_mass_matrix(M, mim, mf, mf);
  gmm::scale(M, u1);
  getfem::add_explicit_matrix(MS, "y", "y", M);

  getfem::add_generic_elliptic_brick(MS, mim, "y", "a");
  gmm::clear(F);
  gmm::add(gmm::scaled(y0, -1.0), F);
  add_Dirichlet_condition_with_simplification(MS, "y", BOUNDARY);

  MS.add_initialized_fem_data("VolumicData", mf, F);
  getfem::add_source_term_brick(MS, mim, "y", "VolumicData");

  cout << "Solving system of direct prime..." << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  gmm::resize(y, mf.nb_dof());
  gmm::copy(MS.real_variable("y"), y);

}



void direct_solve_prime2(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &y, const double u1, const double u2) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  plain_vector y0(mf.nb_dof());
  direct_solve(mesh, mim, mf, y0, u1, u2);

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, -1.0*u2*u2));

  MS.add_fem_variable("y", mf);
  plain_vector F(mf.nb_dof());

  cout << "Assembling for problem-a ..." << endl;
  
  sparse_matrix M(mf.nb_dof(), mf.nb_dof());
  getfem::asm_mass_matrix(M, mim, mf, mf);
  gmm::scale(M, -1.0*u1*u2);
  getfem::add_explicit_matrix(MS, "y", "y", M);

  getfem::add_generic_elliptic_brick(MS, mim, "y", "a");
  getfem::interpolation_function(mf, F, exact_rhs_a);
  gmm::add(gmm::scaled(y0, -1.0*u1), F);
  add_Dirichlet_condition_with_simplification(MS, "y", BOUNDARY);

  MS.add_initialized_fem_data("VolumicData", mf, F);
  getfem::add_source_term_brick(MS, mim, "y", "VolumicData");

  cout << "Solving system of direct prime..." << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  gmm::resize(y, mf.nb_dof());
  gmm::copy(MS.real_variable("y"), y);

}


void direct_solve_second1(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &y, const double u1, const double u2) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

        //plain_vector y0(mf.nb_dof());
        //direct_solve(mesh, mim, mf, y0, u1, u2);
  plain_vector y1(mf.nb_dof());
  direct_solve_prime1(mesh, mim, mf, y1, u1, u2);

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, u2));

  MS.add_fem_variable("y", mf);
  plain_vector F(mf.nb_dof());

  cout << "Assembling for problem-a ..." << endl;

  sparse_matrix M(mf.nb_dof(), mf.nb_dof());
  getfem::asm_mass_matrix(M, mim, mf, mf);
  gmm::scale(M, u1);
  getfem::add_explicit_matrix(MS, "y", "y", M);

  getfem::add_generic_elliptic_brick(MS, mim, "y", "a");
  gmm::clear(F);
  gmm::add(gmm::scaled(y1,-2.0), F);
  add_Dirichlet_condition_with_simplification(MS, "y", BOUNDARY);

  MS.add_initialized_fem_data("VolumicData", mf, F);
  getfem::add_source_term_brick(MS, mim, "y", "VolumicData");

  cout << "Solving system of direct second..." << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  gmm::resize(y, mf.nb_dof());
  gmm::copy(MS.real_variable("y"), y);

}



void direct_solve_second2(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &y, const double u1, const double u2) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  plain_vector y0(mf.nb_dof());
  direct_solve(mesh, mim, mf, y0, u1, u2);
  plain_vector y2(mf.nb_dof());
  direct_solve_prime2(mesh, mim, mf, y2, u1, u2);

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, 0.5*u2*u2*u2));

  MS.add_fem_variable("y", mf);
  plain_vector F(mf.nb_dof());

  cout << "Assembling for problem-a ..." << endl;

  sparse_matrix M(mf.nb_dof(), mf.nb_dof());
  getfem::asm_mass_matrix(M, mim, mf, mf);
  gmm::scale(M, 0.5*u1*u2*u2);
  getfem::add_explicit_matrix(MS, "y", "y", M);

  getfem::add_generic_elliptic_brick(MS, mim, "y", "a");
  getfem::interpolation_function(mf, F, exact_rhs_a);
  gmm::add(gmm::scaled(y0,-1.0*u1), F);
  gmm::add(gmm::scaled(y2,u1*u2), F);
  add_Dirichlet_condition_with_simplification(MS, "y", BOUNDARY);

  MS.add_initialized_fem_data("VolumicData", mf, F);
  getfem::add_source_term_brick(MS, mim, "y", "VolumicData");

  cout << "Solving system of direct second..." << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  gmm::resize(y, mf.nb_dof());
  gmm::copy(MS.real_variable("y"), y);

}




void compute_grad_control_NN(double &Gu1, double &Gu2, getfem::mesh &mesh, getfem::mesh_fem &mf, getfem::mesh_im &mim, 
                             const double u1, const double u2, const plain_vector ydata, 
                             std::vector<double> weights1, const std::vector<double> weights_trans1, std::vector<double> weights2, const std::vector<double> weights_trans2) {


  Gu1 = reg_prime_u(u1, weights1, weights_trans1, L);
  Gu2 = reg_prime_u(u2, weights2, weights_trans2, L);

  plain_vector y(nb_dof);
  direct_solve(mesh, mim, mf, y, u1, u2);
  plain_vector p(nb_dof);
  adjoint_solve(mesh, mim, mf, p, u1, u2, y, ydata);

  sparse_matrix M(nb_dof, nb_dof);
  gmm::clear(M);

  plain_vector prod(nb_dof);
  gmm::clear(prod);

  plain_vector Agiven(nb_dof);
  getfem::model MS;

  cout << "Assembling the gradient for problem-a ..." << endl;
  for (size_type i=0; i<nb_dof; ++i) {
	Agiven[i] = u2;
  }
  getfem::asm_stiffness_matrix_for_laplacian(M, mim, mf, mf, Agiven);
  gmm::mult(M, y, prod);
  Gu2 += gmm::vect_sp(prod, p);

  sparse_matrix Mass(mf.nb_dof(), mf.nb_dof());
  getfem::asm_mass_matrix(Mass, mim, mf, mf);
  gmm::scale(Mass, u1);
  gmm::clear(prod);
  gmm::mult(Mass, y , prod);
  Gu1 += gmm::vect_sp(prod, p);
  cout << "Gradient with the NN reg is computed." << endl;

}


void matlab_export_tol_subproblem(const double absGu) {

  std::ofstream tol;
  tol.open("./MATLAB/Subtol.txt", ios::out|ios::app);

  tol << absGu << ", " ;

  tol.close();

}




void subproblem_Nesterov_L2(double &u1, double &u2, const plain_vector ydata,
			    getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf,
			    const unsigned Nmax, const double small_step, const double tol_nes) {

	// Step n=0
        double su1 = u1;
        double su2 = u2;
        double umem1 = 0.0;
        double umem2 = 0.0;
        double Gu1 = 0.0;
        double Gu2 = 0.0;
	
	unsigned n=1;
	double tol = 1.0;
	double tolmem = 0.0;
	double difftol = 1.0;

	while ((n<=Nmax) && (std::min(tol,difftol)>tol_nes)) {

		// Step n-1
		compute_grad_control_L2(Gu1, Gu2, mesh, mf, mim, su1, su2, ydata);
                umem1 = u1;
                umem2 = u2;
		// Step n
		u1 = su1 - small_step*Gu1;
		su1 = u1 + (double(n-1)/double(n+2))*(u1 - umem1);
		u2 = su2 - small_step*Gu2;
		su2 = u2 + (double(n-1)/double(n+2))*(u2 - umem2);	

		tolmem = tol;
		tol = abs(Gu1) + abs(Gu2);
		cout << "Nesterov for subproblem - step " << n << " ; tol = " << tol << endl;
		matlab_export_tol_subproblem(tol);
		difftol = abs(tol-tolmem);
		matlab_export_tol_subproblem(difftol);	

		n++;

	}

  	std::ofstream tolstream;
  	tolstream.open("./MATLAB/Subtol.txt", ios::out|ios::app);
	tolstream << "; " << endl;// << endl;
	tolstream.close();

	matlab_export_coefflearnedNN1(u1);
	matlab_export_coefflearnedNN2(u2);		

	cout << "Sub-Nesterov terminated after " << Nmax << " iterations; " << "Tol = " << abs(Gu1) << endl;
	cout << "Sub-Nesterov terminated after " << Nmax << " iterations; " << "Tol = " << abs(Gu2) << endl;

}




void subproblem_Nesterov(double &u1, double &u2, const plain_vector ydata, 
                        const std::vector<double> weights1, const std::vector<double> weights_trans1,
                        const std::vector<double> weights2, const std::vector<double> weights_trans2,
			 getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf,
			 const unsigned Nmax, const double small_step, const double tol_nes) {

	// Step n=0
	double su1 = u1;
	double umem1 = 0.0;
	double Gu1 = 0.0;
	double su2 = u2;
	double umem2 = 0.0;
	double Gu2 = 0.0;
	
	unsigned n=1;
	double tol = 1.0;
	double tolmem = 0.0;
	double difftol = 1.0;

	while ((n<=Nmax) && (std::min(tol,difftol)>tol_nes)) {

		// Step n-1
		compute_grad_control_NN(Gu1, Gu2, mesh, mf, mim, su1, su2, ydata, weights1, weights_trans1, weights2, weights_trans2);
		umem1 = u1;
		umem2 = u2;
		// Step n
		u1 = su1 - small_step*Gu1;
		su1 = u1 + (double(n-1)/double(n+2))*(u1 - umem1);
		u2 = su2 - small_step*Gu2;
		su2 = u2 + (double(n-1)/double(n+2))*(u2 - umem2);		

		tolmem = tol;
		tol = abs(Gu1) + abs(Gu2);
		cout << "Nesterov for subproblem - step " << n << " ; tol = " << tol << endl;
		matlab_export_tol_subproblem(tol);
		difftol = abs(tol-tolmem);
		matlab_export_tol_subproblem(difftol);	

		n++;

	}

  	std::ofstream tolstream;
  	tolstream.open("./MATLAB/Subtol.txt", ios::out|ios::app);
	tolstream << "; " << endl;// << endl;
	tolstream.close();

	matlab_export_coefflearnedNN1(u1);	
	matlab_export_coefflearnedNN2(u2);	

	cout << "Sub-Nesterov terminated after " << Nmax << " iterations; " << "Tol = " << abs(Gu1) << endl;
	cout << "Sub-Nesterov terminated after " << Nmax << " iterations; " << "Tol = " << abs(Gu2) << endl;

}




// For solving the subproblem with the neural network as regularizer




void superadjoint_solve(std::vector<double> &mu1, std::vector<double> &mu2, getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, 
			const std::vector<double> weights1, const std::vector<double> weights_trans1,
			const std::vector<double> weights2, const std::vector<double> weights_trans2,
			const std::vector<plain_vector> y, const std::vector<plain_vector> ydata, 
			const std::vector<double> u1, const std::vector<double> u2, 
                       const std::vector<double> udata1, const std::vector<double> udata2,
                       const std::vector<double> utruth1, const std::vector<double> utruth2) {

	gmm::resize(mu1, u1.size()); // We have K data u_k, so dim(Gradient) = K, and dim(mu) = K 
	gmm::resize(mu2, u2.size());

        double Gprime1 = 0.0;
        double Gprime2 = 0.0;

		//size_type nb_dof = mf.nb_dof();

	sparse_matrix M(nb_dof, nb_dof);
	asm_mass_matrix(M, mim, mf, mf);

	plain_vector yprime1(nb_dof);
	plain_vector ysecond1(nb_dof);
	plain_vector yprime2(nb_dof);
	plain_vector ysecond2(nb_dof);
	plain_vector prod(nb_dof);

	for (size_type k=0; k<y.size(); ++k) {

		Gprime1 = 0.0;
                Gprime2 = 0.0;

		gmm::clear(yprime1);
		gmm::clear(yprime2);
			//cout << "UUUUUUUUUUUUU ===== " << u[k] << endl;
		//cout << "Solving direct-prime..." << endl;
		a = utruth2[k];
		b = utruth1[k];

		direct_solve_prime1(mesh, mim, mf, yprime1, u1[k], u2[k]); 
		direct_solve_prime2(mesh, mim, mf, yprime2, u1[k], u2[k]); 
		//cout << "Direct-prime done." << endl;
		gmm::clear(prod);
		gmm::mult(M, yprime1, prod);
		Gprime1 += gmm::vect_sp(prod, yprime1);

		gmm::clear(prod);
		gmm::mult(M, yprime2, prod);
		Gprime2 += gmm::vect_sp(prod, yprime2);

		gmm::clear(prod);
		gmm::copy(ydata[k], prod);
		gmm::scale(prod, -1.0);
		gmm::add(y[k], prod);
		gmm::clear(ysecond1);
		gmm::clear(ysecond2);
		//cout << "Solving direct-second..." << endl;
		direct_solve_second1(mesh, mim, mf, ysecond1, u1[k], u2[k]); 
		direct_solve_second2(mesh, mim, mf, ysecond2, u1[k], u2[k]); 
		//cout << "Direct-second done." << endl;
		Gprime1 += gmm::vect_sp(prod, ysecond1);
		Gprime1 += reg_second_uu(u1[k], weights1, weights_trans1);
		Gprime2 += gmm::vect_sp(prod, ysecond2);
		Gprime2 += reg_second_uu(u2[k], weights2, weights_trans2);

		if (gmm::abs(Gprime1) > 1.e-12) {

			mu1[k] = -(u1[k] - utruth1[k])/Gprime1;

		}

		if (gmm::abs(Gprime2) > 1.e-12) {

			mu2[k] = -(u2[k] - utruth2[k])/Gprime2;

		}

	}

}


void compute_supergradient(std::vector<double> &Supergrad1, std::vector<double> &Supergrad2, 
                          const std::vector<double> weights1, const std::vector<double> weights_trans1,
                          const std::vector<double> weights2, const std::vector<double> weights_trans2,
			   getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf,
			   std::vector<plain_vector> &y, const std::vector<plain_vector> ydata,
			   const unsigned Nmax_subnesterov, const double step_size_nesterov, const double tol_subnesterov, 
			   std::vector<double> &mu1, std::vector<double> &mu2, 
                          std::vector<double> &u1, std::vector<double> &u2, 
                          const std::vector<double> udata1, const std::vector<double> udata2, 
                          const std::vector<double> utruth1, const std::vector<double> utruth2) {

	size_type Ndata1 = utruth1.size();
	size_type Ndata2 = utruth2.size();

	// Solve the equivalent of the state equation
	for (size_type k=0; k<Ndata1; ++k) {

		cout << "k = " << k << endl;
		u2[k] = a_init;
		u1[k] = b_init;
		a = utruth2[k];
		b = utruth1[k];
		subproblem_Nesterov(u1[k], u2[k], ydata[k], weights1, weights_trans1, weights2, weights_trans2, mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
		//subproblem_barzilai_borwein(u[k], ydata[k], weights, mesh, mim, mf);
		direct_solve(mesh, mim, mf, y[k], u1[k], u2[k]);

	}

  	std::ofstream coefflearnedstream1;
  	coefflearnedstream1.open("./MATLAB/CoefflearnedNN1.txt", ios::out|ios::app);
	coefflearnedstream1 << "; " << endl << endl;
	coefflearnedstream1.close();

  	std::ofstream coefflearnedstream2;
  	coefflearnedstream2.open("./MATLAB/CoefflearnedNN2.txt", ios::out|ios::app);
	coefflearnedstream2 << "; " << endl << endl;
	coefflearnedstream2.close();

	cout << "SUPER ADJOINT SOLVE..." << endl;
	superadjoint_solve(mu1, mu2, mesh, mim, mf, weights1, weights_trans1, weights2, weights_trans2, y, ydata, u1, u2, udata1, udata2, utruth1, utruth2);


	gmm::resize(Supergrad1, 2*weights1.size());	
	gmm::clear(Supergrad1);
	gmm::resize(Supergrad2, 2*weights2.size());	
	gmm::clear(Supergrad2);
	
	for (size_type l=0; l<weights1.size(); ++l) {

		Supergrad1[l] = 0.0;

		for (size_type k=0; k<mu1.size(); ++k) {

			Supergrad1[l] += reg_second_uw(u1[k], weights1, weights_trans1, l)*mu1[k];
			//Supergrad[l] += grad_penalization_term_constant_NN(weights, weights_trans, l);

		}

		Supergrad1[l] += nux*weights1[l];

	}

	for (size_type l=0; l<weights2.size(); ++l) {

		Supergrad2[l] = 0.0;

		for (size_type k=0; k<mu2.size(); ++k) {

			Supergrad2[l] += reg_second_uw(u2[k], weights2, weights_trans2, l)*mu2[k];
			//Supergrad[l] += grad_penalization_term_constant_NN(weights, weights_trans, l);

		}

		Supergrad2[l] += nux*weights2[l];

	}

	Supergrad1[weights1.size()-1] = 0.0; // For avoiding trivial NN to be critical points
	Supergrad2[weights2.size()-1] = 0.0; // For avoiding trivial NN to be critical points

	/*for (size_type l=weights.size(); l<Supergrad.size(); ++l) {

		Supergrad[l] = 0.0;

		for (size_type k=0; k<mu.size(); ++k) {

			//cout << mu[k] << endl;
			Supergrad[l] += reg_second_uw_trans(u[k], weights, weights_trans, l)*mu[k];

		}

		Supergrad[l] += nux*weights_trans[l];
			//cout << Supergrad[l] << endl;
	}*/

		//gmm::add(gmm::scaled(weights, nux), Supergrad);

	/*double misfit = 0.0;
	compute_misfit(misfit, u, utruth);
	matlab_export_misfit(misfit);*/


}





void Super_BB_algorithm(std::vector<double> &weights1, std::vector<double> &weights_trans1, std::vector<double> &weights2, std::vector<double> &weights_trans2, 
                       std::vector<double> &Supergrad1, std::vector<double> &Supergrad2,
			std::vector<double> &mu1, std::vector<double> &mu2, 
                       getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, 
			std::vector<plain_vector> &y, const std::vector<plain_vector> ydata, 
			std::vector<double> &u1, const std::vector<double> udata1, const std::vector<double> utruth1,
			std::vector<double> &u2, const std::vector<double> udata2, const std::vector<double> utruth2,
			const unsigned Nmax_subnesterov, const double step_size_nesterov, const double tol_subnesterov, 
			const unsigned Nmax, const double super_step, const double super_tol) {


	std::vector<double> Gw1(weights1.size());
	std::vector<double> Gwmem1(weights1.size());
	std::vector<double> weightsmem1(weights1.size());
	std::vector<double> Sw1(weights1.size());
	std::vector<double> Tw1(weights1.size());

	std::vector<double> Gwt1(weights_trans1.size());
	std::vector<double> Gwtmem1(weights_trans1.size());
	std::vector<double> weightstransmem1(weights_trans1.size());
	std::vector<double> Swt1(weights_trans1.size());
	std::vector<double> Twt1(weights_trans1.size());

	std::vector<double> Gw2(weights2.size());
	std::vector<double> Gwmem2(weights2.size());
	std::vector<double> weightsmem2(weights2.size());
	std::vector<double> Sw2(weights2.size());
	std::vector<double> Tw2(weights2.size());

	std::vector<double> Gwt2(weights_trans2.size());
	std::vector<double> Gwtmem2(weights_trans2.size());
	std::vector<double> weightstransmem2(weights_trans2.size());
	std::vector<double> Swt2(weights_trans2.size());
	std::vector<double> Twt2(weights_trans2.size());


	// First gradient
	compute_supergradient(Supergrad1, Supergrad2, weights1, weights_trans1, weights2, weights_trans2, mesh, mim, mf, y, ydata, Nmax_subnesterov, step_size_nesterov, tol_subnesterov, mu1, mu2, u1, u2, udata1, udata2, utruth1, utruth2);

		//gmm::copy(Gw, Gwmem);
		//gmm::copy(Gwt, Gwtmem);

	for (size_type l=0; l<weights1.size(); ++l) {
		Gw1[l] = Supergrad1[l];
		Gwt1[l] = Supergrad1[l+weights1.size()];
	}
	for (size_type l=0; l<weights2.size(); ++l) {
		Gw2[l] = Supergrad2[l];
		Gwt2[l] = Supergrad2[l+weights2.size()];
	}


	// First naive update
	gmm::copy(weights1, weightsmem1);
	gmm::copy(weights_trans1, weightstransmem1);
	gmm::copy(weights2, weightsmem2);
	gmm::copy(weights_trans2, weightstransmem2);	

	gmm::add(gmm::scaled(Gw1, -1.0*super_step), weights1);
	gmm::add(gmm::scaled(Gwt1, -1.0*super_step), weights_trans1);
	gmm::add(gmm::scaled(Gw2, -1.0*super_step), weights2);
	gmm::add(gmm::scaled(Gwt2, -1.0*super_step), weights_trans2);

	// New gradient
	gmm::copy(Gw1, Gwmem1);
	gmm::copy(Gwt1, Gwtmem1);
	gmm::copy(Gw2, Gwmem2);
	gmm::copy(Gwt2, Gwtmem2);

	compute_supergradient(Supergrad1, Supergrad2, weights1, weights_trans1, weights2, weights_trans2, mesh, mim, mf, y, ydata, Nmax_subnesterov, step_size_nesterov, tol_subnesterov, mu1, mu2, u1, u2, udata1, udata2, utruth1, utruth2);


	for (size_type l=0; l<weights1.size(); ++l) {
		Gw1[l] = Supergrad1[l];
		Gwt1[l] = Supergrad1[l+weights1.size()];
	}
	for (size_type l=0; l<weights2.size(); ++l) {
		Gw2[l] = Supergrad2[l];
		Gwt2[l] = Supergrad2[l+weights2.size()];
	}


	// Initialize the steps
	double a1lpha_w_1 = 0.0;
	double a1lpha_w_2 = 0.0;
	double a1lpha_wt_1 = 0.0; 
	double a1lpha_wt_2 = 0.0;
	double aw_1_up1 = 0.0;
	double aw_1_down1 = 0.0;
	double aw_2_up1 = 0.0;
	double aw_2_down1 = 0.0;
	double awt_1_up1 = 0.0;
	double awt_1_down1 = 0.0;
	double awt_2_up1 = 0.0; 
	double awt_2_down1 = 0.0;

	double a2lpha_w_1 = 0.0;
	double a2lpha_w_2 = 0.0;
	double a2lpha_wt_1 = 0.0; 
	double a2lpha_wt_2 = 0.0;
	double aw_1_up2 = 0.0;
	double aw_1_down2 = 0.0;
	double aw_2_up2 = 0.0;
	double aw_2_down2 = 0.0;
	double awt_1_up2 = 0.0;
	double awt_1_down2 = 0.0;
	double awt_2_up2 = 0.0; 
	double awt_2_down2 = 0.0;

	// Let's go!!!
	double misfit = 0.0;

	double tol = 1.0;
	double tolmem = 0.0;
	double difftol = 1.0;

	unsigned n=1;
	while ((n<=Nmax) && (std::min(tol,100.0*difftol)>super_tol)) {

		// BB coeffs
		
		// For the weights
		gmm::clear(Sw1);
		gmm::clear(Tw1);
		gmm::add(weights1, gmm::scaled(weightsmem1, -1.0), Sw1);
		gmm::add(Gw1, gmm::scaled(Gwmem1, -1.0), Tw1);

		gmm::clear(Sw2);
		gmm::clear(Tw2);
		gmm::add(weights2, gmm::scaled(weightsmem2, -1.0), Sw2);
		gmm::add(Gw2, gmm::scaled(Gwmem2, -1.0), Tw2);

		aw_1_up1 = gmm::vect_sp(Sw1, Sw1); 
		aw_1_down1 = gmm::vect_sp(Sw1, Tw1);
		aw_2_up1 = gmm::vect_sp(Sw1, Tw1);
		aw_2_down1 = gmm::vect_sp(Tw1, Tw1);

		aw_1_up2 = gmm::vect_sp(Sw2, Sw2); 
		aw_1_down2 = gmm::vect_sp(Sw2, Tw2);
		aw_2_up2 = gmm::vect_sp(Sw2, Tw2);
		aw_2_down2 = gmm::vect_sp(Tw2, Tw2);

		if (gmm::abs(aw_1_down1) > 1.e-12) a1lpha_w_1 = aw_1_up1/aw_1_down1;
		if (gmm::abs(aw_2_down1) > 1.e-12) a1lpha_w_2 = aw_2_up1/aw_2_down1;
		if (gmm::abs(aw_1_down2) > 1.e-12) a2lpha_w_1 = aw_1_up2/aw_1_down2;
		if (gmm::abs(aw_2_down2) > 1.e-12) a2lpha_w_2 = aw_2_up2/aw_2_down2;

		// For the weights_trans
		gmm::clear(Swt1);
		gmm::clear(Twt1);
		gmm::add(weights_trans1, gmm::scaled(weightstransmem1, -1.0), Swt1);
		gmm::add(Gwt1, gmm::scaled(Gwtmem1, -1.0), Twt1);

		gmm::clear(Swt2);
		gmm::clear(Twt2);
		gmm::add(weights_trans2, gmm::scaled(weightstransmem2, -1.0), Swt2);
		gmm::add(Gwt2, gmm::scaled(Gwtmem2, -1.0), Twt2);

		awt_1_up1 = gmm::vect_sp(Swt1, Swt1); 
		awt_1_down1 = gmm::vect_sp(Swt1, Twt1);
		awt_2_up1 = gmm::vect_sp(Swt1, Twt1);
		awt_2_down1 = gmm::vect_sp(Twt1, Twt1);
		awt_1_up2 = gmm::vect_sp(Swt2, Swt2); 
		awt_1_down2 = gmm::vect_sp(Swt2, Twt2);
		awt_2_up2 = gmm::vect_sp(Swt2, Twt2);
		awt_2_down2 = gmm::vect_sp(Twt2, Twt2);

		if (gmm::abs(awt_1_down1) > 1.e-12) a1lpha_wt_1 = awt_1_up1/awt_1_down1;
		if (gmm::abs(awt_2_down1) > 1.e-12) a1lpha_wt_2 = awt_2_up1/awt_2_down1;
		if (gmm::abs(awt_1_down2) > 1.e-12) a2lpha_wt_1 = awt_1_up2/awt_1_down2;
		if (gmm::abs(awt_2_down2) > 1.e-12) a2lpha_wt_2 = awt_2_up2/awt_2_down2;

		// Update the control
		gmm::copy(weights1, weightsmem1);
		gmm::copy(weights_trans1, weightstransmem1);
		gmm::copy(weights2, weightsmem2);
		gmm::copy(weights_trans2, weightstransmem2);

		if (n & 1) {
			gmm::add(gmm::scaled(Gw1, -1.0*a1lpha_w_1), weights1);
			gmm::add(gmm::scaled(Gwt1, -1.0*a1lpha_wt_1), weights_trans1);
			gmm::add(gmm::scaled(Gw2, -1.0*a2lpha_w_1), weights2);
			gmm::add(gmm::scaled(Gwt2, -1.0*a2lpha_wt_1), weights_trans2);
		}
    		else {
			gmm::add(gmm::scaled(Gw1, -1.0*a1lpha_w_2), weights1);
			gmm::add(gmm::scaled(Gwt1, -1.0*a1lpha_wt_2), weights_trans1);
			gmm::add(gmm::scaled(Gw2, -1.0*a2lpha_w_2), weights2);
			gmm::add(gmm::scaled(Gwt2, -1.0*a2lpha_wt_2), weights_trans2);
		}
			
		// Update the gradient
		compute_supergradient(Supergrad1, Supergrad2, weights1, weights_trans1, weights2, weights_trans2, mesh, mim, mf, y, ydata, Nmax_subnesterov, step_size_nesterov, tol_subnesterov, mu1, mu2, u1, u2, udata1, udata2, utruth1, utruth2);

		gmm::copy(Gw1, Gwmem1);
		gmm::copy(Gwt1, Gwtmem1);
		gmm::copy(Gw2, Gwmem2);
		gmm::copy(Gwt2, Gwtmem2);
		for (size_type l=0; l<weights1.size(); ++l) {

			Gw1[l] = Supergrad1[l];
			Gwt1[l] = Supergrad1[weights1.size()+l];

		}
		for (size_type l=0; l<weights2.size(); ++l) {

			Gw2[l] = Supergrad2[l];
			Gwt2[l] = Supergrad2[weights2.size()+l];

		}

		// Update the tolerance
		tolmem = tol;
		tol = mynorm(Supergrad1) + mynorm(Supergrad2);
		difftol = std::abs(tol-tolmem);

		cout << "Residual = " << tol << endl;
		matlab_export_supertol(gmm::vect_norm2(Supergrad1)+gmm::vect_norm2(Supergrad2));

		compute_misfit(misfit, u1, u2, utruth1, utruth2);
		matlab_export_misfit(misfit);	

		n++;

	}


}




/***************************************************************************************************************************/
/***************************************************************************************************************************/
/*********************************************                              ************************************************/
/*********************************************   MAIN DATA-DRIVEN PROBLEM   ************************************************/
/*********************************************                              ************************************************/
/***************************************************************************************************************************/
/***************************************************************************************************************************/



int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
    
  // Read parameters.
  bgeot::md_param PARAM;
  PARAM.read_command_line(argc, argv);
  a_init = double(PARAM.real_value("A_INIT", "Initial coefficient for problem-a"));
  b_init = double(PARAM.real_value("B_INIT", "Initial coefficient for problem-b"));
  	//a = double(PARAM.real_value("A", "Ground truth coefficient for problem-a"));

  L = int(PARAM.int_value("L", "Number of layers for the neural network"));

  alphax = double(PARAM.real_value("ALPHAX", "Cost parameter for standard regularization")); 
  betax = double(PARAM.real_value("BETAX", "Coefficient penalizing constant neural networks"));
  nux = double(PARAM.real_value("NUX", "Cost parameter for the neural network weights"));

  //TYPE_REG = int(PARAM.int_value("TYPE_REG", "Choice of regularization type"));

  TYPE_ACT_FUNC = int(PARAM.int_value("TYPE_ACT_FUNC", "Choice of activation functions"));

  eps = double(PARAM.real_value("EPS", "Tolerance for gradient algorithms"));


  // Load the mesh
  getfem::mesh mesh;
  std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
  getfem::import_mesh(MESH_FILE, mesh);

  // Integration method
  std::string IM = PARAM.string_value("IM", "Mesh file");
  getfem::mesh_im mim(mesh);
  mim.set_integration_method(mesh.convex_index(), getfem::int_method_descriptor(IM));

  // Finite element method
  getfem::mesh_fem mf(mesh);
  std::string FEM = PARAM.string_value("FEM", "finite element method");
  mf.set_qdim(bgeot::dim_type(1));
  mf.set_finite_element(mesh.convex_index(), getfem::fem_descriptor(FEM));
  nb_dof = mf.nb_dof();

  getfem::mesh_fem mf_data(mesh);
  std::string FEM_DATA = PARAM.string_value("FEM_DATA", "finite element method");
  mf_data.set_qdim(bgeot::dim_type(1));
  mf_data.set_finite_element(mesh.convex_index(), getfem::fem_descriptor(FEM_DATA));
  	//size_type nb_dof_data = mf_data.nb_dof();




  /**************************************************************************************/
  /********************************* GENERATE DATA **************************************/
  /**************************************************************************************/

  size_type Ndata = 1 + size_type((1.0-0.5)/0.05);
  	//size_type K = Ndata;
  std::vector<double> utruth1(Ndata);
  std::vector<double> utruth2(Ndata);

  for (size_type k=0; k<Ndata; ++k) {
	
	utruth1[k] = 0.5 + double(k)*0.05;
	utruth2[k] = 0.5 + double(k)*0.05;

  }

  matlab_export_truth1(utruth1);
  matlab_export_truth2(utruth2);

  std::vector<plain_vector> ydata(Ndata);

  for (size_type k=0; k<Ndata; ++k) {

	a = utruth2[k];
	b = utruth1[k];
	direct_solve(mesh, mim, mf, ydata[k], utruth1[k], utruth2[k]); // 'ydata' is measured, from 'utruth', namely the ground truth

  }

  // Introduce a noise on the measurements of the state
  std::vector<double> vec_noise(nb_dof);

  getfem::interpolation_function(mf, vec_noise, noise1);
  for (size_type k=0; k<Ndata; ++k) {

	gmm::add(vec_noise, ydata[k]);	

  }

  cout << "Data input generated!" << endl;




  unsigned Nmax_subnesterov = unsigned(PARAM.real_value("NMAX_SUB", "Max nb iterations for sub-Nesterov"));
  double step_size_nesterov = double(PARAM.real_value("STEP_NESTEROV", "Step size for sub-Nesterov"));
  double tol_subnesterov = double(PARAM.real_value("TOL_SUB", "Max nb iterations for sub-Nesterov"));



  /**************************************************************************************/
  /***************************** GENERATE ARTIFICIAL DATA *******************************/
  /**************************************************************************************/


  std::ofstream tolstream;
  tolstream.open("./MATLAB/Subtol.txt", ios::out|ios::trunc);
  tolstream.close();

  std::ofstream expmisfit;
  expmisfit.open("./MATLAB/Misfit.txt", ios::out|ios::trunc);
  expmisfit.close();

  std::ofstream coeffstream1;
  coeffstream1.open("./MATLAB/Coeffdata1.txt", ios::out|ios::trunc);
  coeffstream1.close();
  std::ofstream coeffstream2;
  coeffstream2.open("./MATLAB/Coeffdata2.txt", ios::out|ios::trunc);
  coeffstream2.close();

  std::ofstream coefflearnedNN01;
  coefflearnedNN01.open("./MATLAB/CoefflearnedNN1.txt", ios::out|ios::trunc);
  coefflearnedNN01.close();
  std::ofstream coefflearnedNN02;
  coefflearnedNN02.open("./MATLAB/CoefflearnedNN2.txt", ios::out|ios::trunc);
  coefflearnedNN02.close();


  std::ofstream coefflearnedL20;
  coefflearnedL20.open("./MATLAB/CoefflearnedL2.txt", ios::out|ios::trunc);
  coefflearnedL20.close();

  std::vector<double> udata1(Ndata);
  std::vector<double> udata2(Ndata);

  for (size_type k=0; k<Ndata; ++k) {

	a = utruth2[k];
	udata2[k] = a_init;//utruth[k];//a_init;
	b = utruth1[k];
	udata1[k] = b_init;//utruth[k];//a_init;

	subproblem_Nesterov_L2(udata1[k], udata2[k], ydata[k], mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
		//simple_algorithm_barzilai_borwein(udata[k], ydata[k], mesh, mim, mf);
  	matlab_export_coeff1(udata1[k]);
  	matlab_export_coeff2(udata2[k]);

		// Also export in a specific file
  		std::ofstream evolcoeffNN033;
  		evolcoeffNN033.open("./MATLAB/CoefflearnedL2.txt", ios::out|ios::app);
  		evolcoeffNN033 << udata1[k] << ",   " ;
  		evolcoeffNN033.close();

  }

  std::ofstream evolcoeffNN011;
  evolcoeffNN011.open("./MATLAB/CoefflearnedNN1.txt", ios::out|ios::app);
  evolcoeffNN011 << " ; " << endl << endl;
  evolcoeffNN011.close();
  std::ofstream evolcoeffNN012;
  evolcoeffNN012.open("./MATLAB/CoefflearnedNN2.txt", ios::out|ios::app);
  evolcoeffNN012 << " ; " << endl << endl;
  evolcoeffNN012.close();
  
  cout << "Data output generated!" << endl;


  /**************************************************************************************/
  /****************************** SOLVE THE MAIN PROBLEM ********************************/
  /**************************************************************************************/

  std::ofstream supertolstream;
  supertolstream.open("./MATLAB/Supertol.txt", ios::out|ios::trunc);
  supertolstream.close();

  std::vector<double> weights1(L);
  for (size_type l=0; l<L; ++l) {
	weights1[l] = -1.0*(double(l&1)-0.5)*(-2.0);
  }
  std::vector<double> weights_trans1(L);
  for (size_type l=0; l<L; ++l) {
	weights_trans1[l] = 0.0;//(double(l&1)-0.5)*(1.0);
  }

  std::vector<double> weights2(L);
  for (size_type l=0; l<L; ++l) {
	weights2[l] = -1.0*(double(l&1)-0.5)*(-2.0);
  }
  std::vector<double> weights_trans2(L);
  for (size_type l=0; l<L; ++l) {
	weights_trans2[l] = 0.0;//(double(l&1)-0.5)*(1.0);
  }

  std::vector<double> mu1(Ndata);
  std::vector<double> mu2(Ndata);
  std::vector<plain_vector> y(Ndata);
  std::vector<double> u1(Ndata);
  std::vector<double> u2(Ndata);
  std::vector<double> Supergrad1(2*L);
  std::vector<double> Supergrad2(2*L);
  	//gmm::clear(Supergrad);

  cout << "First gradient..." << endl;
  compute_supergradient(Supergrad1, Supergrad2, weights1, weights_trans1, weights2, weights_trans2, mesh, mim, mf, y, ydata, Nmax_subnesterov, step_size_nesterov, tol_subnesterov, mu1, mu2, u1, u2, udata1, udata2, utruth1, utruth2);


  cout << "Super Gradient method is starting..." << endl;

  unsigned Nmax = unsigned(PARAM.real_value("NMAX", "Max nb iterations for Super-Nesterov"));
  double super_step = double(PARAM.real_value("SUPER_STEP", "Step size for Super-Nesterov"));
  double super_tol = double(PARAM.real_value("SUPER_TOL", "Tolerance for the Super-Nesterov"));

  /*double tol = 1.0;
  double tolmem = 0.0;
  double difftol = 1.0;

  unsigned n=1;
  while ((n<=Nmax) && (std::min(tol,100.0*difftol)>super_tol)) {
  //for (unsigned it=0; it<Nmax; ++it) {

	//cout << "ITER = " << it << endl;
		//gmm::add(gmm::scaled(Supergrad, -1.0*super_step), weights); 
	for (size_type l=0; l<weights.size(); ++l) {

		weights[l] -= super_step*Supergrad[l];
		weights_trans[l] -= super_step*Supergrad[weights.size()+l];

	}

	// Update y and u from the subproblem with NN as regularizer
	for (size_type k=0; k<Ndata; ++k) {

		//cout << "k = " << k << endl;
		u[k] = a_init;
		a = utruth[k];
		subproblem_Nesterov(u[k], ydata[k], weights, weights_trans, mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
		//subproblem_barzilai_borwein(u[k], ydata[k], weights, mesh, mim, mf);
		direct_solve(mesh, mim, mf, y[k], u[k]);

	}

  	std::ofstream coefflearnedstream2;
  	coefflearnedstream2.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::app);
	coefflearnedstream2 << "; " << endl << endl;
	coefflearnedstream2.close();

	// Compute the supergradient
	superadjoint_solve(mu, mesh, mim, mf, weights, weights_trans, y, ydata, u, udata, utruth);
	compute_supergradient(Supergrad, weights, weights_trans, mesh, mim, mf, y, ydata, Nmax_subnesterov, step_size_nesterov, tol_subnesterov, mu, u, udata, utruth);


	tolmem = tol;
	tol = mynorm(Supergrad);
	difftol = std::abs(tol-tolmem);

	cout << "Residual = " << tol << endl;
	matlab_export_supertol(gmm::vect_norm2(Supergrad));

	compute_misfit(misfit, u, utruth);
	matlab_export_misfit(misfit);

	++n;

  }*/

  /*Super_Nesterov_algorithm(weights, weights_trans, Supergrad, mu, mesh, mim, mf, 
			   y, ydata, u, udata, utruth,
			   Nmax_subnesterov, step_size_nesterov, tol_subnesterov, 
			   Nmax, super_step, super_tol);*/

  Super_BB_algorithm(weights1, weights_trans1, weights2, weights_trans2, Supergrad1, Supergrad2, mu1, mu2, mesh, mim, mf, 
		     y, ydata, u1, udata1, utruth1, u2, udata2, utruth2,
		     Nmax_subnesterov, step_size_nesterov, tol_subnesterov, 
		     Nmax, super_step, super_tol);

  std::ofstream coeffweights1;
  coeffweights1.open("./MATLAB/Weights1.txt", ios::out|ios::trunc);
  coeffweights1.close();
  std::ofstream coeffweights2;
  coeffweights2.open("./MATLAB/Weights2.txt", ios::out|ios::trunc);
  coeffweights2.close();

  std::ofstream coeffweightstrans1;
  coeffweightstrans1.open("./MATLAB/Weights_trans1.txt", ios::out|ios::trunc);
  coeffweightstrans1.close();
  std::ofstream coeffweightstrans2;
  coeffweightstrans2.open("./MATLAB/Weights_trans12.txt", ios::out|ios::trunc);
  coeffweightstrans2.close();

  matlab_export_weights1(weights1);
  matlab_export_weights_trans1(weights_trans1);
  matlab_export_weights2(weights2);
  matlab_export_weights_trans2(weights_trans2);



  /**************************************************************************************/
  /************************************** VALIDATION ************************************/
  /**************************************************************************************/


  // New exact data
  size_type Ndatabis = Ndata;
  std::vector<double> utruthbis1(Ndatabis);
  std::vector<double> utruthbis2(Ndatabis);
  std::ofstream coefflearnedL204;
  coefflearnedL204.open("./MATLAB/Coeffdata_validation1.txt", ios::out|ios::trunc);
  coefflearnedL204.close();
  std::ofstream coefflearnedL2042;
  coefflearnedL2042.open("./MATLAB/Coeffdata_validation2.txt", ios::out|ios::trunc);
  coefflearnedL2042.close();
  for (size_type k=0; k<Ndata; ++k) {
        utruthbis1[k] = 0.475 + double(k)*0.05;
        utruthbis2[k] = 0.475 + double(k)*0.05;
        matlab_export_coeffbis1(utruth1[k]);
        matlab_export_coeffbis2(utruth2[k]);
  }
  std::vector<plain_vector> ydata2(Ndatabis);
  for (size_type k=0; k<Ndatabis; ++k) {
        a = utruthbis2[k];
        b = utruthbis1[k];
        direct_solve(mesh, mim, mf, ydata2[k], utruthbis1[k], utruthbis2[k]); // 'ydata' is measured, from 'utruth', namely the ground truth
  }
  cout << "Validation Data generated!" << endl;

  // New artificial data
  std::ofstream coefflearnedL202;
  coefflearnedL202.open("./MATLAB/CoefflearnedNoise_validation1.txt", ios::out|ios::trunc);
  coefflearnedL202.close();
  std::ofstream coefflearnedL2022;
  coefflearnedL2022.open("./MATLAB/CoefflearnedNoise_validation2.txt", ios::out|ios::trunc);
  coefflearnedL2022.close();
  std::vector<double> udatabis1(Ndatabis);
  std::vector<double> udatabis2(Ndatabis);

  std::vector<double> vec_noise2(nb_dof);
  getfem::interpolation_function(mf, vec_noise2, noise2);
  for (size_type k=0; k<Ndatabis; ++k) {
	gmm::add(vec_noise2, ydata2[k]);	
  }


  for (size_type k=0; k<Ndatabis; ++k) {

	a = utruthbis2[k];
	udatabis2[k] = a_init;
	b = utruthbis1[k];
	udatabis1[k] = b_init;
	subproblem_Nesterov_L2(udatabis1[k], udatabis2[k], ydata2[k], mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);

  	//matlab_export_coeff2(udata2[k]);

  		std::ofstream evolcoeffNN02;
  		evolcoeffNN02.open("./MATLAB/CoefflearnedNoise_validation1.txt", ios::out|ios::app);
  		evolcoeffNN02 << udatabis1[k] << ",   " ;
  		evolcoeffNN02.close();
  		std::ofstream evolcoeffNN022;
  		evolcoeffNN022.open("./MATLAB/CoefflearnedNoise_validation2.txt", ios::out|ios::app);
  		evolcoeffNN022 << udatabis2[k] << ",   " ;
  		evolcoeffNN022.close();

  }
  cout << "Validation Data output generated!" << endl;


  // Output of the NN from these new (artificial) data
  std::vector<double> ubis1(Ndatabis);
  std::vector<double> ubis2(Ndatabis);
  for (size_type k=0; k<Ndatabis; ++k) {

        cout << "k = " << k << endl;	
        ubis2[k] = a_init;
        a = utruthbis2[k];
        ubis1[k] = b_init;
        b = utruthbis1[k];
        subproblem_Nesterov(ubis1[k], ubis2[k], ydata2[k], weights1, weights_trans1, weights2, weights_trans2, mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);

  }

  // Export the result
  std::ofstream coefflearnedNN03;
  coefflearnedNN03.open("./MATLAB/CoefflearnedNN_validation1.txt", ios::out|ios::trunc);
  coefflearnedNN03.close();
  std::ofstream coefflearnedNN032;
  coefflearnedNN032.open("./MATLAB/CoefflearnedNN_validation2.txt", ios::out|ios::trunc);
  coefflearnedNN032.close();
  for (size_type k=0; k<Ndatabis; ++k) {
        std::ofstream evolcoeffNN035;
        evolcoeffNN035.open("./MATLAB/CoefflearnedNN_validation1.txt", ios::out|ios::app);
        evolcoeffNN035 << ubis1[k] << ",   " ;
        evolcoeffNN035.close();
        std::ofstream evolcoeffNN036;
        evolcoeffNN036.open("./MATLAB/CoefflearnedNN_validation2.txt", ios::out|ios::app);
        evolcoeffNN036 << ubis2[k] << ",   " ;
        evolcoeffNN036.close();
  }

  cout << "Validation terminated!" << endl;



}















































