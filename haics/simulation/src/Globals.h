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

// Global variables shared among all files

#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#include "Definitions.h"

extern TypeInputSettings *g_inpset;

extern gsl_rng *g_rng;

//files
extern FILE *red_lfp; //reduced logfile
extern FILE *tfp; //trajectories
extern FILE *lpdfp; //logPD
extern FILE *hfp; //hamiltonians
extern FILE *tsfp; //time steps

extern FILE *tlfp; //trajectory lengths
extern FILE *artfp; //art.txt for metrics during production
extern FILE *artburnin_fp;//art_burnin.txt for metrics during burn-in

extern FILE *fp_startingpoint; //starting point
extern FILE *final_point_tune_fp; //final point from tuning
extern FILE *final_point_burnin_fp; //final point from burn-in

extern FILE *scfp; //store the scaling factor calculated during the burn-in stage

extern FILE *odefp; //ODE solution file

//variables
extern char *id_global;
extern int hbar_to_b_num_disc_steps;
extern int hbar_to_b_num_disc_steps_3s;
extern double delta_t;

extern TypeMisc *global_misc;

extern int kNumParam;
extern int kNumObs; //dimension of observation set
extern double num_pop; //population size of the dataset

extern int num_out_cvodes;
extern int tfin_spl_bas;

extern double *stepsize;//declared in read_input.c
extern int *traj_length;
extern int *L_value; //added for normalization time
extern int L_lower;
extern double *varphi;
extern double fitting_factor;
extern double h_lower;
extern double h_upper;
extern double hsl_low;
extern double hsl_upp;
extern double hsl_simulation;
extern int dt_random_scheme;
extern double step_randomization;
extern int L_random_scheme;
extern int phi_random_scheme;

extern double **g_coeff;//use in hmc.c, read_input, integrators

extern int *accept_arr; //store the acceptances (1) and the rejections (0)
extern double *dH_array; //store the energy errors
extern int *num_proposals_dHneg_arr; //store the number of proposals with negative energy errors

enum step {back2, back1, curr, forw1, forw2};

//extern TypeState *state;
extern TypeState *state[5];
extern TypeState *temporalStatePlusH[5], *temporalStateMinusH[5];
extern TypeODESol *ODESol[5];

//extern TypeState *state_prev;
//extern TypeState *state_next;
extern TypeHamiltonian *Ham;

extern double **g_obs;//declared in _model.c
extern double **rho_basis;
//extern int *t; //response variable, declared in blr_model.c

extern double (*Prior)();//declared in _model.c
extern double (*Loglikelihood)();
extern double (*LogJacobian)();
extern void (*JacobianMat)();
extern void (*Gradient)(TypeState *state);
extern void (*Hessian)();
extern void (*Transf)();
extern void (*InvTransf)();

extern void (*MDMove)(void (*grad_func)(TypeState *state), int dim, TypeState **st, double ss, double *coeff, int L);

extern const double PI;

#endif
