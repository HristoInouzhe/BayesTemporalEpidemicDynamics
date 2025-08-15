// HaiCS (Hamiltonians in Computational Statistics)
//            Developed in the group of
//         Professor Elena Akhmatskaya by:
//     Elena Akhmatskaya, Mario Fernández-Pendás,
//    Hristo Inouzhe, Felix Müller, Lorenzo Nagar, 
//       Jorge Pérez Heredia, Tijana Radivojević, 
//           María Xosé Rodríguez-Álvarez.
//   For references to this package, please cite:
//     * Radivojevic, T. Enhancing Sampling in
//     Computational Statistics Using Modified
//     Hamiltonians. PhD Thesis, UPV/EHU (2016)
//  * Radivojevic, T., Akhmatskaya, E. Modified
// Hamiltonian Monte Carlo for Bayesian inference.
//          Stat Comput 30, 377–404 (2020).
//
// This version provides the code for running SIR/SEIR-like
// probabilistic models (for Incidence data) with a time-dependent
// transmission rate under a Bayesian approach to infer the
// epidemiological parameters using HaiCS HMC-based samplers.

// Definitions of general constants values and symbols

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <stdbool.h>

#define TRUE 1
#define FALSE 0

// Relational operators
#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define DIST(x1,y1,x2,y2) (sqrt(((x1)-(x2))*((x1)-(x2))-((y1)-(y2))*((y1)-(y2))))

#define EPSILON 1e-10  // tolerance
#define DOUBLE_EQ(x,v) (x > ((v - EPSILON)) && (x < (v + EPSILON)))
#define DOUBLE_GE(x,v) (x > (v - EPSILON))
#define DOUBLE_LE(x,v) (x < (v + EPSILON))

#define FPRINT(token) (fprintf(lfp, #token "\n")) //obsolete


// Store/add a value to/read a matrix element in COLUMN major order
// band matrix
#define CBM_STORE(matrix_band_upper_sym_packed, lda, i, j, vij)  (*(&(matrix_band_upper_sym_packed[0][0]) + (j)*lda+(i)) = (vij))

#define CBM_ADD(matrix_band_upper_sym_packed, lda, i, j, vij)  (*(&(matrix_band_upper_sym_packed[0][0]) + (j)*lda+(i)) += (vij))

#define CBM_READ(matrix_band_upper_sym_packed, lda, i, j)  (*(&(matrix_band_upper_sym_packed[0][0]) + (j)*lda+(i)))

// symmetric packed
#define CSPM_STORE(matrix_upper_sym_packed, i, j, vij)  (*(&(matrix_upper_sym_packed[0]) + (j)*(j+1)/2+(i)) = (vij))

#define CSPM_ADD(matrix_upper_sym_packed, i, j, vij)  (*(&(matrix_upper_sym_packed[0]) + (j)*(j+1)/2+(i)) += (vij))

#define CSPM_READ(matrix_upper_sym_packed, i, j)  (*(&(matrix_upper_sym_packed[0]) + (j)*(j+1)/2+(i)))

// symmetric
#define CSM_STORE(matrix_upper_sym, dim, i, j, vij)  (*(&(matrix_upper_sym[0][0]) + (j)*dim+(i)) = (vij))

#define CSM_ADD(matrix_upper_sym, dim, i, j, vij)  (*(&(matrix_upper_sym[0][0]) + (j)*dim+(i)) += (vij))

#define CSM_READ(matrix_upper_sym, dim, i, j)  (*(&(matrix_upper_sym[0][0]) + (j)*dim+(i)))

// Data structures
// store all simulation parameters
typedef struct {
  char model[50];
  char data[50];
  char method[50];
  int seed;
  int iter_sampling;
  int iter_burn_in;
  char integrator[50];
  int t_L;
  int L;
  int t_stepsize;
  double stepsize;
  int t_varphi;
  double varphi;
  int thinning;
  // ATune features
  double scaling_value;
  double stepsize_delta;
  int iter_tune;
  double AR_target;
  double delta_AR_target;
  char fitting_factor_approach[50];
  int opt_param_setting;
  // entries for SIKR/SEIKR model
  int num_basis_spline;
  int spline_pol_degree;
  double tfin_spl_bas;
  int num_out_cvodes;
  int initialCountsFixed;  // Flag indicating if initial counts are fixed
  double S0_fixed_value;   // Fixed value for S0
  double E0_fixed_value;   // Fixed value for E0
  double I0_fixed_value;   // Fixed value for I0
  double dirich_alpha_0;
  double dirich_alpha_1;
  double dirich_alpha_2;
  double dirich_alpha_3;
  char alpha_prior[50];
  double alpha_param_1;
  double alpha_param_2;
  char gamma_prior[50];
  double gamma_param_1;
  double gamma_param_2;
  char phi_inv_prior[50];
  double phi_inv_param_1;
  double phi_inv_param_2;
  char tau_prior[50];
  double tau_param_1;
  double tau_param_2;
  int num_comp;
  int num_comp_E;
  int do_under_report;
  int T0_under;
  double U0_under;
  int T1_under;
  double U1_under;
  int do_grad_desc;
  double learning_rate;
  int max_iter;
  int gammaFixed;
  double gamma_fixed_value;
  int gammaBounded;
  double gammaUpper;
  double gammaLower;
  char beta_prior[50];
  double beta_param_1;
  double beta_param_2;
} TypeInputSettings;

//Parts of the Hamiltonian function
typedef struct {
  double ham; //hamiltonian
  double pot; //potential
  double kin; //kinetic
  double grad_term; //term with gradient
  double order2; //h^2/24*(...)
  //double hes_term; //term with hessian
  double order4; //h^4/720*(...)
} TypeHamiltonian;

//states
typedef struct {
  double *pos; //position
  double *mom; //momenta
  double *grad; //gradient
  double **hes; //hessian matrix
  //double beta;
} TypeState;

//miscellaneous
typedef struct {
  int dim1_hes; //dimension 1 of the hessian
  int dim2_hes; //dimension 2 of the hessian
} TypeMisc;


typedef struct {
  double **solution; //ODESol
} TypeODESol;

typedef struct {
  double p[50];           /* problem parameters */
  int spline_degree;
  int number_extended_knots;
  int number_basis_elements;
  double *knots;
  double *basis;
} *UserData;

typedef struct {
  double p[50];           /* problem parameters */
  int spline_degree;
  int number_extended_knots;
  int number_basis_elements;
  double *knots;
  double *basis;
} *UserDataSEIR;

#endif
