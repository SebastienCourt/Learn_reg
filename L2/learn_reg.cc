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
double a = 0.0;

size_type L = 0;

double alphax = 0.0;
double nux = 0.0;

double eps = 0.0;

int TYPE_REG = 1;
int TYPE_ACT_FUNC = 1;


size_type nb_dof = 0.0;




/***************************************************************************************************************************/
/*********************************************** ROUTINES FOR THE SUBPROBLEM ***********************************************/
/***************************************************************************************************************************/


void matlab_export_misfit(const double misfit) {

  std::ofstream expmisfit;
  expmisfit.open("./MATLAB/Misfit.txt", ios::out|ios::app);
  expmisfit << misfit << ",   " ;
  expmisfit.close();

}


void matlab_export_coeff(const double coeff) {

  std::ofstream evolcoeff;
  evolcoeff.open("./MATLAB/Coeffdata.txt", ios::out|ios::app);
  evolcoeff << coeff << ",   " ;
  evolcoeff.close();

}


void matlab_export_coefflearnedNN(const double coeff) {

  std::ofstream evolcoeff;
  evolcoeff.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::app);
  evolcoeff << coeff << ",   " ;
  evolcoeff.close();

}


void matlab_export_supertol(const double Normgrad) {

  std::ofstream supertol;
  supertol.open("./MATLAB/Supertol.txt", ios::out|ios::app);

  supertol << Normgrad << ", " ;

  supertol.close();

}


void matlab_export_weights(const std::vector<double> weights) {

  std::ofstream finalweights;
  finalweights.open("./MATLAB/Weights.txt", ios::out|ios::app);

  for (size_type i=0; i<weights.size(); ++i) {

  	finalweights << weights[i] << ", " ;

  }

  finalweights.close();

}



void matlab_export_weights_trans(const std::vector<double> weights) {

  std::ofstream finalweights;
  finalweights.open("./MATLAB/Weights_trans.txt", ios::out|ios::app);

  for (size_type i=0; i<weights.size(); ++i) {

  	finalweights << weights[i] << ", " ;

  }

  finalweights.close();

}









double exact_rhs_a(const base_node p) {

	double x = p[0];

	return -1.0*a*( (2.0-M_PI*M_PI*x*x)*sin(M_PI*x) + 4.0*x*M_PI*cos(M_PI*x) );

}


void direct_solve(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &y, const double u) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, u));
  MS.add_initialized_fixed_size_data("agiven", plain_vector(1, 1.0)); // Diffusion term: Laplace

  MS.add_fem_variable("y", mf);
  plain_vector F(mf.nb_dof());

  cout << "Assembling for problem-a ..." << endl;
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
		   const double u, const plain_vector y, const plain_vector ydata) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, u));

  MS.add_fem_variable("p", mf);
  plain_vector F(nb_dof);

  cout << "Assembling adjoint system for problem-a ..." << endl;
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



void compute_grad_control_L2(double &Gu, getfem::mesh &mesh, getfem::mesh_fem &mf, getfem::mesh_im &mim, const double u, const plain_vector ydata) {

  Gu = alphax*u;

  plain_vector y(nb_dof);
  direct_solve(mesh, mim, mf, y, u);
  plain_vector p(nb_dof);
  adjoint_solve(mesh, mim, mf, p, u, y, ydata);

  sparse_matrix M(nb_dof, nb_dof);
  gmm::clear(M);

  plain_vector prod(nb_dof);
  gmm::clear(prod);

  plain_vector Agiven(nb_dof);
  getfem::model MS;

  cout << "Assembling the gradient for problem-a ..." << endl;
  for (size_type i=0; i<nb_dof; ++i) {
	Agiven[i] = 1.0;
  }
  getfem::asm_stiffness_matrix_for_laplacian(M, mim, mf, mf, Agiven);
  gmm::mult(M, y, prod);
  Gu += gmm::vect_sp(prod, p);
  cout << "Gradient with the L2 reg is computed." << endl;

}



void simple_algorithm_barzilai_borwein(double &u, const plain_vector ydata,
				       getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf) {


  scalar_type alpha_u_1 = 0.0;
  scalar_type alpha_u_2 = 0.0;

  double Gu = 0.0;//base_small_vector Gu(u.size()); 

  double umem = 0.0;//base_small_vector umem(u.size());
  double Gumem = 0.0;//base_small_vector Gumem(u.size());

  double Su = 0.0;//base_small_vector Su(u.size());
  double Tu = 0.0;//base_small_vector Tu(Gu.size());

  scalar_type au_1_up, au_1_down, au_2_up, au_2_down;

  int compteur = 0;
  double tol = 0.0;
  double tol_mem = 0.0;

  // Initializing BB
  compute_grad_control_L2(Gu, mesh, mf, mim, u, ydata);
  umem = u;//gmm::copy(u, umem);
  Gumem = Gu;//gmm::copy(Gu, Gumem);
  tol_mem = abs(Gumem);//gmm::vect_norm2(Gumem);

  u -= 0.0001*Gu;
  compute_grad_control_L2(Gu, mesh, mf, mim, u, ydata);
  tol = abs(Gu);//gmm::vect_norm2(Gu);

  cout << "Simple-Barzilai-Borwein is starting..." << endl;

  while (((tol > eps) && ( compteur < 2000) && (gmm::abs(tol-tol_mem) > eps)) || (compteur < 3)) {

	compteur++;

	// Compute scalar-products
	Su = u - umem;//gmm::add(u, gmm::scaled(umem, -1.0), Su);
	Tu = Gu - Gumem;//gmm::add(Gu, gmm::scaled(Gumem, -1.0), Tu);

	au_1_up = Su*Su;//gmm::vect_sp(Su, Su); 
	au_1_down = Su*Tu;//gmm::vect_sp(Su, Tu);
	au_2_up = Su*Tu;//gmm::vect_sp(Su, Tu);
	au_2_down = Tu*Tu;//gmm::vect_sp(Tu, Tu);

	if (gmm::abs(au_1_down) > 1.e-12) alpha_u_1 = au_1_up/au_1_down;
	if (gmm::abs(au_2_down) > 1.e-12) alpha_u_2 = au_2_up/au_2_down;


	// Save the previous quantities
	umem = u;//gmm::copy(u, umem);
	Gumem = Gu;//gmm::copy(Gu, Gumem);
	

	// Update of the unknowns -- control parameters
	if (compteur & 1) {
		u -= alpha_u_1*Gu;//gmm::add(gmm::scaled(Gu, -1.0*alpha_u_1), u);
	}
    	else {
		u -= alpha_u_2*Gu;//gmm::add(gmm::scaled(Gu, -1.0*alpha_u_2), u);
	}


	compute_grad_control_L2(Gu, mesh, mf, mim, u, ydata);
	tol_mem = tol;
	tol = abs(Gu);//gmm::vect_norm2(Gu);

	cout << "Tol Simple-Barzilai-Borwein-" << compteur << " = " << tol;

  }


  cout << " END " << endl;
  cout << "Simple-Barzilai-Borwein done, after " << compteur << " iterations." << endl;


}



/***************************************************************************************************************************/
/********************************************* NEURAL NETWORKS AND SENSITIVITY *********************************************/
/***************************************************************************************************************************/


double func_max(const double x, const double y) {


  return 0.5*(x+y+gmm::abs(x-y));

}


double activation_func(const double p, const int ACT_FUNC) {


  double res = 0.0;

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

	default:

		cout << "Wrong choice of activation function!" << endl;

  }

  
  return res;

}



double activation_func_prime(const double p, const int ACT_FUNC) {


  double res = 0.0;
  double thx = 0.0;

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

	default:

		cout << "Wrong choice of activation function!" << endl;

  }



  return res;

}


double activation_func_second(const double p, const int ACT_FUNC) {


  double res = 0.0;
  double thx = 0.0;

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



void compute_misfit(double &misfit, const std::vector<double> u, const std::vector<double> udata) {

  std::vector<double> diff(udata.size());
  for (size_type k=0; k<udata.size(); ++k) {

	diff[k] = u[k] - udata[k];

  }
  double ratio = gmm::vect_norm2(diff)/gmm::vect_norm2(udata);

  misfit = 100.0*ratio;

}



void direct_solve_prime(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &y, const double u) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, -1.0*u*u));

  MS.add_fem_variable("y", mf);
  plain_vector F(mf.nb_dof());

  cout << "Assembling for problem-a ..." << endl;
  getfem::add_generic_elliptic_brick(MS, mim, "y", "a");
  getfem::interpolation_function(mf, F, exact_rhs_a);
  add_Dirichlet_condition_with_simplification(MS, "y", BOUNDARY);

  MS.add_initialized_fem_data("VolumicData", mf, F);
  getfem::add_source_term_brick(MS, mim, "y", "VolumicData");

  cout << "Solving system of direct prime..." << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  gmm::resize(y, mf.nb_dof());
  gmm::copy(MS.real_variable("y"), y);

}


void direct_solve_second(getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, plain_vector &y, const double u) {


  int BOUNDARY = 1;
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh, r);
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	mesh.region(BOUNDARY).add(i.cv(),i.f());
  }

  getfem::model MS;

  MS.add_initialized_fixed_size_data("a", plain_vector(1, 0.5*u*u*u));

  MS.add_fem_variable("y", mf);
  plain_vector F(mf.nb_dof());

  cout << "Assembling for problem-a ..." << endl;
  getfem::add_generic_elliptic_brick(MS, mim, "y", "a");
  getfem::interpolation_function(mf, F, exact_rhs_a);
  add_Dirichlet_condition_with_simplification(MS, "y", BOUNDARY);

  MS.add_initialized_fem_data("VolumicData", mf, F);
  getfem::add_source_term_brick(MS, mim, "y", "VolumicData");

  cout << "Solving system of direct second..." << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  gmm::resize(y, mf.nb_dof());
  gmm::copy(MS.real_variable("y"), y);

}




void compute_grad_control_NN(double &Gu, getfem::mesh &mesh, getfem::mesh_fem &mf, getfem::mesh_im &mim, 
			     const double u, const plain_vector ydata, std::vector<double> weights, const std::vector<double> weights_trans) {


  Gu = reg_prime_u(u, weights, weights_trans, L);

  plain_vector y(nb_dof);
  direct_solve(mesh, mim, mf, y, u);
  plain_vector p(nb_dof);
  adjoint_solve(mesh, mim, mf, p, u, y, ydata);

  sparse_matrix M(nb_dof, nb_dof);
  gmm::clear(M);

  plain_vector prod(nb_dof);
  gmm::clear(prod);

  plain_vector Agiven(nb_dof);
  getfem::model MS;

  cout << "Assembling the gradient for problem-a ..." << endl;
  for (size_type i=0; i<nb_dof; ++i) {
	Agiven[i] = 1.0;
  }
  getfem::asm_stiffness_matrix_for_laplacian(M, mim, mf, mf, Agiven);
  gmm::mult(M, y, prod);
  Gu += gmm::vect_sp(prod, p);
  cout << "Gradient with the NN reg is computed." << endl;

}


void matlab_export_tol_subproblem(const double absGu) {

  std::ofstream tol;
  tol.open("./MATLAB/Subtol.txt", ios::out|ios::app);

  tol << absGu << ", " ;

  tol.close();

}




void subproblem_Nesterov_L2(double &u, const plain_vector ydata,
			    getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf,
			    const unsigned Nmax, const double small_step, const double tol_nes) {

	// Step n=0
	double su = u;
	double umem = 0.0;
	double Gu = 0.0;
	
	unsigned n=1;
	double tol = 1.0;
	double tolmem = 0.0;
	double difftol = 1.0;

	while ((n<=Nmax) && (std::min(tol,difftol)>tol_nes)) {

		// Step n-1
		compute_grad_control_L2(Gu, mesh, mf, mim, su, ydata);
		umem = u;
		// Step n
		u = su - small_step*Gu;
		su = u + (double(n-1)/double(n+2))*(u - umem);	

		tolmem = tol;
		tol = abs(Gu);
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

	matlab_export_coefflearnedNN(u);	

	cout << "Sub-Nesterov terminated after " << Nmax << " iterations; " << "Tol = " << abs(Gu) << endl;

}




void subproblem_Nesterov(double &u, const plain_vector ydata, const std::vector<double> weights, const std::vector<double> weights_trans,
			 getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf,
			 const unsigned Nmax, const double small_step, const double tol_nes) {

	// Step n=0
	double su = u;
	double umem = 0.0;
	double Gu = 0.0;
	
	unsigned n=1;
	double tol = 1.0;
	double tolmem = 0.0;
	double difftol = 1.0;

	while ((n<=Nmax) && (std::min(tol,difftol)>tol_nes)) {

		// Step n-1
		compute_grad_control_NN(Gu, mesh, mf, mim, su, ydata, weights, weights_trans);
		umem = u;
		// Step n
		u = su - small_step*Gu;
		su = u + (double(n-1)/double(n+2))*(u - umem);	

		tolmem = tol;
		tol = abs(Gu);
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

	matlab_export_coefflearnedNN(u);	

	cout << "Sub-Nesterov terminated after " << Nmax << " iterations; " << "Tol = " << abs(Gu) << endl;

}




// For solving the subproblem with the neural network as regularizer
void subproblem_barzilai_borwein(double &u, const plain_vector ydata, const std::vector<double> weights, const std::vector<double> weights_trans,
				 getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf) {


  scalar_type alpha_u_1 = 0.0;
  scalar_type alpha_u_2 = 0.0;

  double umem = 0.0;//base_small_vector umem(u.size());
  double Gumem = 0.0;//base_small_vector Gumem(Gu.size());

  double Su = 0.0;//base_small_vector Su(u.size());
  double Tu = 0.0;//base_small_vector Tu(Gu.size());

  scalar_type au_1_up, au_1_down, au_2_up, au_2_down;

  int compteur = 0;
  double tol = 0.0;
  double tol_mem = 0.0;

  // Initializing BB
  double Gu = 0.0;
  compute_grad_control_NN(Gu, mesh, mf, mim, u, ydata, weights, weights_trans);
  umem = u;//gmm::copy(u, umem);
  Gumem = Gu;//gmm::copy(Gu, Gumem);
  tol_mem = abs(Gumem);//gmm::vect_norm2(Gumem);

  u -= 0.0001*Gu;//gmm::add(gmm::scaled(Gu, -1.0*0.01), u);
  compute_grad_control_NN(Gu, mesh, mf, mim, u, ydata, weights, weights_trans);
  tol = abs(Gu);//gmm::vect_norm2(Gu);

  cout << "NN-Barzilai-Borwein is starting..." << endl;

  while (((tol > eps) && ( compteur < 2000) && (gmm::abs(tol-tol_mem) > eps)) || (compteur < 3)) {

	compteur++;

	// Compute scalar-products
	Su = u - umem;//gmm::add(u, gmm::scaled(umem, -1.0), Su);
	Tu = Gu - Gumem;//gmm::add(Gu, gmm::scaled(Gumem, -1.0), Tu);

	au_1_up = Su*Su;//gmm::vect_sp(Su, Su); 
	au_1_down = Su*Tu;//gmm::vect_sp(Su, Tu);
	au_2_up = Su*Tu;//gmm::vect_sp(Su, Tu);
	au_2_down = Tu*Tu;//gmm::vect_sp(Tu, Tu);

	if (gmm::abs(au_1_down) > 1.e-12) alpha_u_1 = au_1_up/au_1_down;
	if (gmm::abs(au_2_down) > 1.e-12) alpha_u_2 = au_2_up/au_2_down;


	// Save the previous quantities
	umem = u;//gmm::copy(u, umem);
	Gumem = Gu;//gmm::copy(Gu, Gumem);
	

	// Update of the unknowns -- control parameters
	if (compteur & 1) {
		u -= alpha_u_1*Gu;//gmm::add(gmm::scaled(Gu, -1.0*alpha_u_1), u);
	}
    	else {
		u -= alpha_u_2*Gu;//gmm::add(gmm::scaled(Gu, -1.0*alpha_u_2), u);
	}


	compute_grad_control_NN(Gu, mesh, mf, mim, u, ydata, weights, weights_trans);
	tol_mem = tol;
	tol = abs(Gu);//gmm::vect_norm2(Gu);

	cout << "Tol NN-Barzilai-Borwein-" << compteur << " = " << tol;

  }


  cout << " END " << endl;
  cout << "NN-Borwein done, after " << compteur << " iterations." << endl;


}




void superadjoint_solve(std::vector<double> &mu, getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, 
			const std::vector<double> weights, const std::vector<double> weights_trans,
			const std::vector<plain_vector> y, const std::vector<plain_vector> ydata, 
			const std::vector<double> u, const std::vector<double> udata, const std::vector<double> utruth) {

	gmm::resize(mu, u.size()); // We have K data u_k, so dim(Gradient) = K, and dim(mu) = K 

	double Gprime = 0.0;

		//size_type nb_dof = mf.nb_dof();

	sparse_matrix M(nb_dof, nb_dof);
	asm_mass_matrix(M, mim, mf, mf);

	plain_vector yprime(nb_dof);
	plain_vector ysecond(nb_dof);
	plain_vector prod(nb_dof);

	for (size_type k=0; k<y.size(); ++k) {

		Gprime = 0.0;

		gmm::clear(yprime);
			//cout << "UUUUUUUUUUUUU ===== " << u[k] << endl;
		//cout << "Solving direct-prime..." << endl;
		a = utruth[k];

		direct_solve_prime(mesh, mim, mf, yprime, u[k]); 
		//cout << "Direct-prime done." << endl;
		gmm::clear(prod);
		gmm::mult(M, yprime, prod);
		Gprime += gmm::vect_sp(prod, yprime);

		gmm::clear(prod);
		gmm::copy(ydata[k], prod);
		gmm::scale(prod, -1.0);
		gmm::add(y[k], prod);
		gmm::clear(ysecond);
		//cout << "Solving direct-second..." << endl;
		direct_solve_second(mesh, mim, mf, ysecond, u[k]); 
		//cout << "Direct-second done." << endl;
		Gprime += gmm::vect_sp(prod, ysecond);

		Gprime += reg_second_uu(u[k], weights, weights_trans);

		if (gmm::abs(Gprime) > 1.e-12) {

			mu[k] = -(u[k] - udata[k])/Gprime;

		}

	}

}


void compute_supergradient(std::vector<double> &Supergrad, const std::vector<double> weights, const std::vector<double> weights_trans,
			   const std::vector<double> mu, const std::vector<double> u) {


	gmm::resize(Supergrad, 2*weights.size());	
	gmm::clear(Supergrad);
	
	for (size_type l=0; l<weights.size(); ++l) {

		Supergrad[l] = 0.0;

		for (size_type k=0; k<mu.size(); ++k) {

			//cout << mu[k] << endl;
			Supergrad[l] += reg_second_uw(u[k], weights, weights_trans, l)*mu[k];

		}

		Supergrad[l] += nux*weights[l];
			//cout << Supergrad[l] << endl;
	}

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


}



void Super_Nesterov_algorithm(std::vector<double> &weights, std::vector<double> &weights_trans, std::vector<double> &Supergrad, 
			      std::vector<double> &mu, getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, 
			      std::vector<plain_vector> &y, const std::vector<plain_vector> ydata, 
			      std::vector<double> &u, const std::vector<double> udata, const std::vector<double> utruth,
			      const unsigned Nmax_subnesterov, const double step_size_nesterov, const double tol_subnesterov, 
			      const unsigned Nmax, const double super_step, const double super_tol) {


	std::vector<double> sw(weights.size());
	std::vector<double> sw_trans(weights_trans.size());
	// Step n=0
	for (size_type l=0; l<weights.size(); ++l) {
		sw[l] = weights[l];
		sw_trans[l] = weights_trans[l];
	}

	std::vector<double> weightsmem(weights.size());
	std::vector<double> weightstransmem(weights_trans.size());

	double misfit = 0.0;

	double tol = 1.0;
	double tolmem = 0.0;
	double difftol = 1.0;

	unsigned n=1;
	while ((n<=Nmax) && (std::min(tol,100.0*difftol)>super_tol)) {
		//for (unsigned n=1; n<=Nmax; n++) {


		// Gradient
		for (size_type k=0; k<utruth.size(); ++k) {

			//cout << "k = " << k << endl;
			u[k] = a_init;
			a = utruth[k];
			subproblem_Nesterov(u[k], ydata[k], sw, sw_trans, mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
			//subproblem_barzilai_borwein(u[k], ydata[k], weights, mesh, mim, mf);
			direct_solve(mesh, mim, mf, y[k], u[k]);

		}
		superadjoint_solve(mu, mesh, mim, mf, sw, sw_trans, y, ydata, u, udata, utruth);
		compute_supergradient(Supergrad, sw, sw_trans, mu, u);


		tolmem = tol;
		tol = mynorm(Supergrad);
		difftol = std::abs(tol-tolmem);

			// Update (quick and cheap)
			for (size_type l=0; l<weights.size(); ++l) {

				// Step n-1
				weightsmem[l] = weights[l];
				weightstransmem[l] = weights_trans[l];
				// Step n
				weights[l] = sw[l] - super_step*Supergrad[l];
				weights_trans[l] = sw_trans[l] - super_step*Supergrad[weights.size()+l];
				sw[l] = weights[l] + (double(n-1)/double(n+2))*(weights[l] - weightsmem[l]);
				sw_trans[l] = weights_trans[l] + (double(n-1)/double(n+2))*(weights_trans[l] - weightstransmem[l]);
			
			}
/*for (size_type l=0; l<Supergrad.size(); ++l) {
  
matlab_export_coefflearnedNN(Supergrad[l]);
cout << Supergrad[l] << endl;
}*/

		cout << "Residual = " << tol << endl;
		matlab_export_supertol(gmm::vect_norm2(Supergrad));

  		std::ofstream coefflearnedstream2;
  		coefflearnedstream2.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::app);
		coefflearnedstream2 << "; " << endl << endl;
		coefflearnedstream2.close();

		compute_misfit(misfit, u, udata);
		matlab_export_misfit(misfit);	

		++n;

	}


}




void Super_BB_algorithm(std::vector<double> &weights, std::vector<double> &weights_trans, std::vector<double> &Supergrad, 
			std::vector<double> &mu, getfem::mesh &mesh, getfem::mesh_im &mim, getfem::mesh_fem &mf, 
			std::vector<plain_vector> &y, const std::vector<plain_vector> ydata, 
			std::vector<double> &u, const std::vector<double> udata, const std::vector<double> utruth,
			const unsigned Nmax_subnesterov, const double step_size_nesterov, const double tol_subnesterov, 
			const unsigned Nmax, const double super_step, const double super_tol) {


	std::vector<double> Gw(weights.size());
	std::vector<double> Gwmem(weights.size());
	std::vector<double> weightsmem(weights.size());
	std::vector<double> Sw(weights.size());
	std::vector<double> Tw(weights.size());

	std::vector<double> Gwt(weights_trans.size());
	std::vector<double> Gwtmem(weights_trans.size());
	std::vector<double> weightstransmem(weights_trans.size());
	std::vector<double> Swt(weights_trans.size());
	std::vector<double> Twt(weights_trans.size());


/*std::ofstream mys;
mys.open("./MATLAB/Weights.txt", ios::out|ios::app);
mys << "Firstgradient BB ;" ;
mys.close();*/
	// First gradient
	for (size_type k=0; k<utruth.size(); ++k) {

		//cout << "k = " << k << endl;
		u[k] = a_init;
		a = utruth[k];
		subproblem_Nesterov(u[k], ydata[k], weights, weights_trans, mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
		//subproblem_barzilai_borwein(u[k], ydata[k], weights, mesh, mim, mf);
		direct_solve(mesh, mim, mf, y[k], u[k]);

	}
	superadjoint_solve(mu, mesh, mim, mf, weights, weights_trans, y, ydata, u, udata, utruth);
	compute_supergradient(Supergrad, weights, weights_trans, mu, u);

		//gmm::copy(Gw, Gwmem);
		//gmm::copy(Gwt, Gwtmem);
/*std::ofstream mys1;
mys1.open("./MATLAB/Weights.txt", ios::out|ios::app);
mys1 << "Maybe here BB ;" ;
mys1.close();*/
	for (size_type l=0; l<weights.size(); ++l) {
		Gw[l] = Supergrad[l];
		Gwt[l] = Supergrad[l+weights.size()];
	}

  		std::ofstream coefflearnedstream2000;
  		coefflearnedstream2000.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::app);
		coefflearnedstream2000 << "; " << endl << endl;
		coefflearnedstream2000.close();


/*std::ofstream mys2;
mys2.open("./MATLAB/Weights.txt", ios::out|ios::app);
mys2 << "First naive update BB ;" ;
mys2.close();*/
	// First naive update
	gmm::copy(weights, weightsmem);
	gmm::copy(weights_trans, weightstransmem);
	
	gmm::add(gmm::scaled(Gw, -1.0*super_step), weights);
	gmm::add(gmm::scaled(Gwt, -1.0*super_step), weights_trans);

/*std::ofstream mys3;
mys3.open("./MATLAB/Weights.txt", ios::out|ios::app);
mys3 << "New gradient BB ;" ;
mys3.close();*/
	// New gradient
	gmm::copy(Gw, Gwmem);
	gmm::copy(Gwt, Gwtmem);

	for (size_type k=0; k<utruth.size(); ++k) {

		//cout << "k = " << k << endl;
		u[k] = a_init;
		a = utruth[k];
		subproblem_Nesterov(u[k], ydata[k], weights, weights_trans, mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
		//subproblem_barzilai_borwein(u[k], ydata[k], weights, mesh, mim, mf);
		direct_solve(mesh, mim, mf, y[k], u[k]);

	}
	superadjoint_solve(mu, mesh, mim, mf, weights, weights_trans, y, ydata, u, udata, utruth);
	compute_supergradient(Supergrad, weights, weights_trans, mu, u);

  		std::ofstream coefflearnedstream200;
  		coefflearnedstream200.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::app);
		coefflearnedstream200 << "; " << endl << endl;
		coefflearnedstream200.close();
/*std::ofstream mys35;
mys35.open("./MATLAB/Weights.txt", ios::out|ios::app);
mys35 << "New gradient BB2 ;" ;
mys35.close();*/
	for (size_type l=0; l<weights.size(); ++l) {
		Gw[l] = Supergrad[l];
		Gwt[l] = Supergrad[l+weights.size()];
	}


/*std::ofstream mys4;
mys4.open("./MATLAB/Weights.txt", ios::out|ios::app);
mys4 << "Init BB ;" ;
mys4.close();*/
	// Initialize the steps
	double alpha_w_1 = 0.0;
	double alpha_w_2 = 0.0;
	double alpha_wt_1 = 0.0; 
	double alpha_wt_2 = 0.0;
	double aw_1_up = 0.0;
	double aw_1_down = 0.0;
	double aw_2_up = 0.0;
	double aw_2_down = 0.0;
	double awt_1_up = 0.0;
	double awt_1_down = 0.0;
	double awt_2_up = 0.0; 
	double awt_2_down = 0.0;

	// Let's go!!!
	double misfit = 0.0;

	double tol = 1.0;
	double tolmem = 0.0;
	double difftol = 1.0;

	unsigned n=1;
	while ((n<=Nmax) && (std::min(tol,100.0*difftol)>super_tol)) {

		// BB coeffs
		
		// For the weights
		gmm::clear(Sw);
		gmm::clear(Tw);
		gmm::add(weights, gmm::scaled(weightsmem, -1.0), Sw);
		gmm::add(Gw, gmm::scaled(Gwmem, -1.0), Tw);

		aw_1_up = gmm::vect_sp(Sw, Sw); 
		aw_1_down = gmm::vect_sp(Sw, Tw);
		aw_2_up = gmm::vect_sp(Sw, Tw);
		aw_2_down = gmm::vect_sp(Tw, Tw);

		if (gmm::abs(aw_1_down) > 1.e-12) alpha_w_1 = aw_1_up/aw_1_down;
		if (gmm::abs(aw_2_down) > 1.e-12) alpha_w_2 = aw_2_up/aw_2_down;

		// For the weights_trans
		gmm::clear(Swt);
		gmm::clear(Twt);
		gmm::add(weights_trans, gmm::scaled(weightstransmem, -1.0), Swt);
		gmm::add(Gwt, gmm::scaled(Gwtmem, -1.0), Twt);

		awt_1_up = gmm::vect_sp(Swt, Swt); 
		awt_1_down = gmm::vect_sp(Swt, Twt);
		awt_2_up = gmm::vect_sp(Swt, Twt);
		awt_2_down = gmm::vect_sp(Twt, Twt);

		if (gmm::abs(awt_1_down) > 1.e-12) alpha_wt_1 = awt_1_up/awt_1_down;
		if (gmm::abs(awt_2_down) > 1.e-12) alpha_wt_2 = awt_2_up/awt_2_down;

		// Update the control
		gmm::copy(weights, weightsmem);
		gmm::copy(weights_trans, weightstransmem);

		if (n & 1) {
			gmm::add(gmm::scaled(Gw, -1.0*alpha_w_1), weights);
			gmm::add(gmm::scaled(Gwt, -1.0*alpha_wt_1), weights_trans);
		}
    		else {
			gmm::add(gmm::scaled(Gw, -1.0*alpha_w_2), weights);
			gmm::add(gmm::scaled(Gwt, -1.0*alpha_wt_2), weights_trans);
		}
			
		// Update the gradient
		for (size_type k=0; k<utruth.size(); ++k) {

			//cout << "k = " << k << endl;
			u[k] = a_init;
			a = utruth[k];
			subproblem_Nesterov(u[k], ydata[k], weights, weights_trans, mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
			//subproblem_barzilai_borwein(u[k], ydata[k], weights, mesh, mim, mf);
			direct_solve(mesh, mim, mf, y[k], u[k]);

		}
		superadjoint_solve(mu, mesh, mim, mf, weights, weights_trans, y, ydata, u, udata, utruth);
		compute_supergradient(Supergrad, weights, weights_trans, mu, u);

		gmm::copy(Gw, Gwmem);
		gmm::copy(Gwt, Gwtmem);
		for (size_type l=0; l<weights.size(); ++l) {

			Gw[l] = Supergrad[l];
			Gwt[l] = Supergrad[weights.size()+l];

		}

		// Update the tolerance
		tolmem = tol;
		tol = mynorm(Supergrad);
		difftol = std::abs(tol-tolmem);

		cout << "Residual = " << tol << endl;
		matlab_export_supertol(gmm::vect_norm2(Supergrad));

  		std::ofstream coefflearnedstream2;
  		coefflearnedstream2.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::app);
		coefflearnedstream2 << "; " << endl << endl;
		coefflearnedstream2.close();

		compute_misfit(misfit, u, udata);
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
  	//a = double(PARAM.real_value("A", "Ground truth coefficient for problem-a"));

  L = int(PARAM.int_value("L", "Number of layers for the neural network"));

  alphax = double(PARAM.real_value("ALPHAX", "Cost parameter for standard regularization"));
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

  size_type Ndata = size_type((1.55-1.05)/0.05);
  	//size_type K = Ndata;
  std::vector<double> utruth(Ndata);

  for (size_type k=0; k<Ndata; ++k) {
	
	utruth[k] = 0.5 + double(k)*0.05;

  }

  std::vector<plain_vector> ydata(Ndata);

  for (size_type k=0; k<Ndata; ++k) {

	a = utruth[k];
	direct_solve(mesh, mim, mf, ydata[k], utruth[k]); // 'ydata' is measured, from 'utruth', namely the ground truth

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

  std::ofstream coeffstream;
  coeffstream.open("./MATLAB/Coeffdata.txt", ios::out|ios::trunc);
  coeffstream.close();

  std::ofstream coefflearnedNN0;
  coefflearnedNN0.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::trunc);
  coefflearnedNN0.close();

  std::ofstream coefflearnedL20;
  coefflearnedL20.open("./MATLAB/CoefflearnedL2.txt", ios::out|ios::trunc);
  coefflearnedL20.close();

  std::vector<double> udata(Ndata);

  for (size_type k=0; k<Ndata; ++k) {

	a = utruth[k];
	udata[k] = a_init;//utruth[k];//a_init;
	subproblem_Nesterov_L2(udata[k], ydata[k], mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
		//simple_algorithm_barzilai_borwein(udata[k], ydata[k], mesh, mim, mf);
  	matlab_export_coeff(udata[k]);

  		std::ofstream evolcoeffNN0;
  		evolcoeffNN0.open("./MATLAB/CoefflearnedL2.txt", ios::out|ios::app);
  		evolcoeffNN0 << udata[k] << ",   " ;
  		evolcoeffNN0.close();

  }

  std::ofstream evolcoeffNN01;
  evolcoeffNN01.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::app);
  evolcoeffNN01 << " ; " << endl << endl;
  evolcoeffNN01.close();
  
  cout << "Data output generated!" << endl;


  /**************************************************************************************/
  /****************************** SOLVE THE MAIN PROBLEM ********************************/
  /**************************************************************************************/

  std::ofstream supertolstream;
  supertolstream.open("./MATLAB/Supertol.txt", ios::out|ios::trunc);
  supertolstream.close();

  std::vector<double> weights(L);
  for (size_type l=0; l<L; ++l) {
	weights[l] = (double(l&1)-0.5)*(-2.0);
  }
  std::vector<double> weights_trans(L);
  for (size_type l=0; l<L; ++l) {
	weights_trans[l] = (double(l&1)-0.5)*(1.0);
  }

  std::vector<double> mu(Ndata);
  std::vector<plain_vector> y(Ndata);
  std::vector<double> u(Ndata);
  std::vector<double> Supergrad(2*L);
  	//gmm::clear(Supergrad);

  double misfit = 0.0;

  cout << "First gradient..." << endl;

	for (size_type k=0; k<Ndata; ++k) {

		cout << "k = " << k << endl;
		u[k] = a_init;
		a = utruth[k];
		subproblem_Nesterov(u[k], ydata[k], weights, weights_trans, mesh, mim, mf, Nmax_subnesterov, step_size_nesterov, tol_subnesterov);
		//subproblem_barzilai_borwein(u[k], ydata[k], weights, mesh, mim, mf);
		direct_solve(mesh, mim, mf, y[k], u[k]);

	}

  	std::ofstream coefflearnedstream;
  	coefflearnedstream.open("./MATLAB/CoefflearnedNN.txt", ios::out|ios::app);
	coefflearnedstream << "; " << endl << endl;
	coefflearnedstream.close();

	cout << "SUPER ADJOINT SOLVE..." << endl;
	superadjoint_solve(mu, mesh, mim, mf, weights, weights_trans, y, ydata, u, udata, utruth);
	cout << "COMPUTE SUPER GRADIENT..." << endl;
  	compute_supergradient(Supergrad, weights, weights_trans, mu, u);

	compute_misfit(misfit, u, udata);
	matlab_export_misfit(misfit);


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
	compute_supergradient(Supergrad, weights, weights_trans, mu, u);


	tolmem = tol;
	tol = mynorm(Supergrad);
	difftol = std::abs(tol-tolmem);

	cout << "Residual = " << tol << endl;
	matlab_export_supertol(gmm::vect_norm2(Supergrad));

	compute_misfit(misfit, u, udata);
	matlab_export_misfit(misfit);

	++n;

  }*/

  /*Super_Nesterov_algorithm(weights, weights_trans, Supergrad, mu, mesh, mim, mf, 
			   y, ydata, u, udata, utruth,
			   Nmax_subnesterov, step_size_nesterov, tol_subnesterov, 
			   Nmax, super_step, super_tol);*/

  Super_BB_algorithm(weights, weights_trans, Supergrad, mu, mesh, mim, mf, 
		     y, ydata, u, udata, utruth,
		     Nmax_subnesterov, step_size_nesterov, tol_subnesterov, 
		     Nmax, super_step, super_tol);

  std::ofstream coeffweights;
  coeffweights.open("./MATLAB/Weights.txt", ios::out|ios::trunc);
  coeffweights.close();

  std::ofstream coeffweightstrans;
  coeffweightstrans.open("./MATLAB/Weights_trans.txt", ios::out|ios::trunc);
  coeffweightstrans.close();

  matlab_export_weights(weights);
  matlab_export_weights_trans(weights_trans);






}















































