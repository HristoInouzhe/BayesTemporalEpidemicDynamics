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

#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include <stdio.h>
#include "Definitions.h"

void MetropolisMove(void (*gradient_function)(TypeState *st), int dim, TypeState **st, double ss, double *coeff, int L); //added by florian

void VV(void (*gradient_function)(TypeState *st), int dim, TypeState **st, double ss, double *coeff, int L);

void V2S(void (*gradient_function)(TypeState *st), int dim, TypeState **st, double ss, double *coeff, int L);

void V3S(void (*gradient_function)(TypeState *st), int dim, TypeState **st, double ss, double *coeff, int L);

void fnAssignIntegratorCoeffs(double **mat_g_coeff, int start_assign, int end_assign, int flag_burn_in); //added by Felix

void fnAssignStepsizes(double *stepsize, double stepsize_value, int type_random_stepsize, int start_iter_assign, int end_iter_assign, double random_step); //added by Felix

void fnAssignTrajectoryLengths(int *traj_length, int *L_value, int type_random_L, int start_iter_assign, int end_iter_assign, int L_lower); //added by Felix

void fnAssignPhi(double *varphi, double varphi_value, int type_random_varphi, int start_iter_assign, int end_iter_assign, int kNumParam); //added by Felix

#endif
