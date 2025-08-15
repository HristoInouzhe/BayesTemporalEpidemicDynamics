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

//#define NOUT  100     /* number of output times */ //substituted by g_inpset->num_out_cvodes

//#define NSpB  12 
//#define NP_SEIR    3 + NSpB           /* number of problem parameters */
//#define NS_SEIR    3 + NSpB            /* number of sensitivities computed */

#ifndef SIR_MODEL_FUNCTIONS_H
#define SIR_MODEL_FUNCTIONS_H

/* Type : UserDataSEIR */

/* typedef struct {
  realtype p[50]; 
  int spline_degree;
  int number_extended_knots;
  int number_basis_elements;
  double *knots;
  double *basis;
} *UserDataSEIR; */


/* Prototypes of functions by CVODES */

int f_PREV(realtype t, N_Vector y, N_Vector ydot, void *user_data);

int Jac_PREV(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int fS_PREV(int Ns, realtype t, N_Vector y, N_Vector ydot,
              int iS, N_Vector yS, N_Vector ySdot,
              void *user_data, N_Vector tmp1, N_Vector tmp2);

void SIR_CVODES(double results_y[][num_out_cvodes], double results_s[][num_out_cvodes][3 + g_inpset->num_basis_spline], double state_params[3 + g_inpset->num_basis_spline],
                 double degree, double m, double nf, double *knots, double *basis, char *sensitivity_parameters[]);

void printODE(FILE *odefp);

#endif
