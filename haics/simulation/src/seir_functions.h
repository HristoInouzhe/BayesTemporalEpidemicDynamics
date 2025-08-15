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

#ifndef SEIR_FUNCTIONS_H
#define SEIR_FUNCTIONS_H

/* Prototypes of functions by CVODES */

int f_SEIR(realtype t, N_Vector y, N_Vector ydot, void *user_data);

int Jac_SEIR(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int fS_SEIR(int Ns, realtype t, N_Vector y, N_Vector ydot,
              int iS, N_Vector yS, N_Vector ySdot,
              void *user_data, N_Vector tmp1, N_Vector tmp2);

void SEIR_CVODES(double results_y[][num_out_cvodes], double results_s[][num_out_cvodes][kNumParam - 2], double state_params[kNumParam - 2],
                 double degree, double m, double nf, double *knots, double *basis, char *sensitivity_parameters[]);

void printODE_SEIR(FILE *odefp);

#endif
