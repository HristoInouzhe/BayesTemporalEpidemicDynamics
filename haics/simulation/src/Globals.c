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

#include "Globals.h"

TypeInputSettings *g_inpset;

gsl_rng *g_rng;

//files
FILE *red_lfp; //reduced logfile
FILE *tfp; //trajectories
FILE *lpdfp; //logPD
FILE *hfp; //hamiltonians
FILE *tsfp; //time steps

FILE *tlfp; //trajectory lengths
FILE *artfp; //art.txt for metrics during production
FILE *artburnin_fp;//art_burnin.txt for metrics during burn-in

FILE *fp_startingpoint; //starting point
FILE *final_point_tune_fp; //final point from tuning
FILE *final_point_burnin_fp; //final point from burn-in

FILE *scfp; //store the scaling factor calculated during the burn-in stage

FILE *odefp; //ODE solution file

//variables
char *id_global; //instance id, needed globally for assigning integrator coeffs
int hbar_to_b_num_disc_steps = 40000; //this is the number of discretization steps in the "h_bar to b" 2-stage table
int hbar_to_b_num_disc_steps_3s = 60000; //this is the number of discretization steps in the "h_bar to b" 3-stage table
double delta_t; //this is the stepsize used for advancing the integrator. g_inpset->stepsize can also refer to hbar

TypeMisc *global_misc;

int kNumParam;
int kNumObs;
double num_pop;

double *stepsize;
int *traj_length;
int *L_value;  //added for normalization time
int L_lower;
double *varphi;
double fitting_factor;
double h_lower;
double h_upper;
double hsl_low;
double hsl_upp;
double hsl_simulation;
int dt_random_scheme;
double step_randomization;
int L_random_scheme;
int phi_random_scheme;

double **g_coeff;

int *accept_arr; //store the acceptances (1) and the rejections (0)
double *dH_array; //store the energy errors
int *num_proposals_dHneg_arr; //store the number of proposals with negative energy errors

//TypeState *state;
TypeState *state[5] = {NULL,NULL,NULL,NULL,NULL};
TypeODESol *ODESol[5] = {NULL,NULL,NULL,NULL,NULL};
TypeState *temporalStatePlusH[5] = {NULL,NULL,NULL,NULL,NULL}, *temporalStateMinusH[5] = {NULL,NULL,NULL,NULL,NULL};
//TypeState *state_prev;
//TypeState *state_next;
TypeHamiltonian *Ham;

double **g_obs;
double **rho_basis;

int num_out_cvodes;
int tfin_spl_bas;

double (*Prior)();
double (*Loglikelihood)();
double (*LogJacobian)();
void (*JacobianMat)();
void (*Gradient)(TypeState *state);
void (*Hessian)();
void (*Transf)();
void (*InvTransf)();

void (*MDMove)(void (*grad_func)(TypeState *state), int dim, TypeState **st, double ss, double *coeff, int L);//prev and next used only if S&H shadow Hamiltonian

const double PI = 3.1415926535897932384626433832795;
