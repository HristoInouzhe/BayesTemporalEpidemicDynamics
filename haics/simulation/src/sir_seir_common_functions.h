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

#include <stdio.h>
#include "Definitions.h"

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */

#ifndef SIR_SEIR_COMMON_FUNCTIONS_H
#define SIR_SEIR_COMMON_FUNCTIONS_H

// Enumeration for different types of prior distributions
typedef enum {
    PRIOR_NORMAL,
    PRIOR_GAMMA,
    PRIOR_BETA,
    PRIOR_LOGNORMAL,
    PRIOR_WEIBULL,
    PRIOR_KUMARASWAMY,
    PRIOR_LORENTZ,
    PRIOR_INVGAMMA,
    PRIOR_EXPONENTIAL,
    PRIOR_CHISQUARED,
    PRIOR_UNIFORM,
    PRIOR_CONSTANT
} PriorDistributionType;

// Structure to hold prior settings for each parameter
typedef struct {
    PriorDistributionType dist_type;  // Type of prior distribution
    double param1;                    // First parameter (e.g., mean, shape, rate)
    double param2;                    // Second parameter (e.g., sd, scale); unused for one-parameter distributions
} PriorSetting;

// Define macros for initializing PriorSetting structures for all prior distributions

// Existing macros
#define PRIOR_SETTING_NORMAL(mu, sigma)          ((PriorSetting){PRIOR_NORMAL, mu, sigma})
#define PRIOR_SETTING_GAMMA(k, theta)            ((PriorSetting){PRIOR_GAMMA, k, theta})
#define PRIOR_SETTING_EXPONENTIAL(lambda)        ((PriorSetting){PRIOR_EXPONENTIAL, lambda, 0})
#define PRIOR_SETTING_INVGAMMA(alpha, beta)      ((PriorSetting){PRIOR_INVGAMMA, alpha, beta})
#define PRIOR_SETTING_UNIFORM(a, b)              ((PriorSetting){PRIOR_UNIFORM, a, b})
#define PRIOR_SETTING_BETA(a, b)                 ((PriorSetting){PRIOR_BETA, a, b})
#define PRIOR_SETTING_LOGNORMAL(mu, sigma)       ((PriorSetting){PRIOR_LOGNORMAL, mu, sigma})
#define PRIOR_SETTING_WEIBULL(k, lambda)         ((PriorSetting){PRIOR_WEIBULL, k, lambda})
#define PRIOR_SETTING_KUMARASWAMY(a, b)          ((PriorSetting){PRIOR_KUMARASWAMY, a, b})
#define PRIOR_SETTING_LORENTZ(x0, gamma)         ((PriorSetting){PRIOR_LORENTZ, x0, gamma})
#define PRIOR_SETTING_CHISQUARED(k)              ((PriorSetting){PRIOR_CHISQUARED, k, 0})
#define PRIOR_SETTING_CONSTANT(c)              ((PriorSetting){PRIOR_CONSTANT, c, 0})

/* Prototypes of private functions */

void ProcessArgs(char *argv[], booleantype *sensi, int *sensi_meth, booleantype *err_con);
void WrongArgs(char *name);
//static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
//static void PrintOutputS(N_Vector *uS);
//static void PrintFinalStats(void *cvode_mem, booleantype sensi);
 
int check_retval(void *returnvalue, const char *funcname, int opt);

int mindex(int i, int j, int nrow);

double bs(int nknots, int nspline, int updegree, double x, double * knots);

void splinebasis(int *d, int *n, int *m, double * x, double * knots, double * basis);

double splinefunction(double * beta, int *d, int *n, int *m, double * x, double * knots, double * basis, int *nf);

double vectors_dot_prod(double *x, double *y, int n);

double underReport(int t, double u0, double u1 , int t0, int t1);

void matrix_vector_mult(double **mat, double *vec, double *result, int rows, int cols);

void matrix_matrix_mult(double **matrix_left, double **matrix_right, double **result_matrix, int rows, int cols, int prod_dim);

double lbeta(double a, double b);

double binomial_coefficient(int n, int k);

PriorDistributionType get_prior_type_from_string(const char* prior_type_str);

PriorSetting initialize_prior(const char* prior_type_str, double param1, double param2);

double log_prior_density(double x, PriorSetting *prior);

double dirichlet_logdensity_4_param(double x0, double x1, double x2, const double alpha[4]);

double dirichlet_logdensity_3_param(double x0, double x1, const double alpha[3]);

double dirichlet_logdensity_4_param_derivative(double x0, double x1, double x2, const double alpha[4], int i);

double dirichlet_logdensity_3_param_derivative(double x0, double x1, const double alpha[3], int i);

double log_prior_density_derivative(double x, PriorSetting *prior);

double sample_from_prior(PriorSetting *prior);

void sample_dirichlet_4_param(const double alpha[4], double x[4]);

void sample_dirichlet_3_param(const double alpha[3], double x[3]);

#endif
