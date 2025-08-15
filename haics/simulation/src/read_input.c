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
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctype.h>

#include "Globals.h"
#include "Definitions.h"
#include "utils.h"
#include "read_input.h"
#include "integrators.h"
#include "adapt_integrator.h"
#include "sir_model_synthetic.h"
#include "seir_model_synthetic.h"
#include "sir_stand_model_synthetic.h"
#include "sir_stand_model_prevalence_synthetic.h"

//GetValue()
static int GetValue(FILE* fp, char *str, char *fmt,  void *val) {
  char c[128];

  if (1 != fscanf(fp, "%s", c)) {
    fprintf(stderr, "read_input.c - Error reading input file at %s\n", str);
    exit(1);
  }

  if (strcmp(str, c) != 0) {
    fprintf(stderr, "read_input.c - Error reading input file - expected %s found %s.\n Check the input file.\n", str, c);
    exit(1);
  }

  if (1 != fscanf(fp, fmt, val)) {
    fprintf(stderr, "read_input.c - Error reading input file at %s\n", str);
    fprintf(stderr, "read_input.c - Cannot read value format %s\n", fmt);
    exit(1);
  }

  return fscanf(fp, "%*[^\n]");
}

//ReadInput()
int ReadInput(char *id) {

  FILE *ifp;
  unsigned long int seed;
  const gsl_rng_type *T;
  int iter_total;
  char filename[80];

  //global structure for misc variables (allocation + initialisation)
  global_misc = (TypeMisc *)malloc(sizeof(TypeMisc));
  global_misc->dim1_hes = 1;
  global_misc->dim2_hes = 1;

  /////////////////////////
  //inputfile: opening, reading in, closing

  sprintf(filename, "./output/%s/inputfile_%s.txt", id, id);
  ifp = fopen(filename, "r");
  if (ifp == NULL) {
    fprintf(stderr, "read_input.c - Error: Cannot open input file %s\n", filename);
    exit(EXIT_FAILURE);
  }

  //allocate memory for inputfile structure
  g_inpset = (TypeInputSettings *)malloc(sizeof(TypeInputSettings));

  fprintf(red_lfp, "------ Inputfile parameters\n");

  //model: compartmental model
  GetValue(ifp, "model", "%s", &g_inpset->model);
  fprintf(red_lfp, "model: %s\n", g_inpset->model);
  printf("model: %s\n", g_inpset->model);

  //data: set of data to used for Bayesian inference
  GetValue(ifp, "data", "%s", &g_inpset->data);
  fprintf(red_lfp, "data: %s\n",  g_inpset->data);
  printf("data: %s\n", g_inpset->data);

  //method: sampler to use, HMC/GHMC
  GetValue(ifp, "method", "%s", &g_inpset->method);
  fprintf(red_lfp, "method: %s\n", g_inpset->method);

  //seed
  GetValue(ifp, "seed", "%i", &g_inpset->seed);
  fprintf(red_lfp, "seed: %i\n", g_inpset->seed);

  //iter_sampling: number of iterations for the production stage
  GetValue(ifp, "iter_sampling", "%i", &g_inpset->iter_sampling);
  fprintf(red_lfp, "iter_sampling: %i\n", g_inpset->iter_sampling);

  //iter_burn_in: number of iterations for the burn-in stage
  GetValue(ifp, "iter_burn_in", "%i", &g_inpset->iter_burn_in);
  fprintf(red_lfp, "iter_burn_in: %i\n", g_inpset->iter_burn_in);

  //integrator: numerical scheme for integrating Hamiltonian dynamics
  GetValue(ifp, "integrator", "%s", &g_inpset->integrator);
  fprintf(red_lfp, "integrator: %s\n", g_inpset->integrator);

  //t_L: randomization scheme for the number of integration steps per iteration L
  GetValue(ifp, "t_L", "%i", &g_inpset->t_L);
  fprintf(red_lfp, "t_L: %i\n", g_inpset->t_L);

  //L: number of integration steps per iteration
  GetValue(ifp, "L", "%i", &g_inpset->L);
  fprintf(red_lfp, "L: %i\n", g_inpset->L);

  //t_stepsize: randomization scheme for the integration step size
  GetValue(ifp, "t_stepsize", "%i", &g_inpset->t_stepsize);
  fprintf(red_lfp, "t_stepsize: %i\n", g_inpset->t_stepsize);

  //stepsize: integration step size
  GetValue(ifp, "stepsize", "%lf", &g_inpset->stepsize);
  fprintf(red_lfp, "stepsize: %1.4f\n", g_inpset->stepsize);
  printf("stepsize: %1.4f\n", g_inpset->stepsize);

  //thinning: how much the (G)HMC chains should be thinned out before storing them
  GetValue(ifp, "thinning", "%i", &g_inpset->thinning);
  fprintf(red_lfp, "thinning: %i\n", g_inpset->thinning);

  //t_varphi: randomization scheme for the noise in the Partial Momentum Update for GHMC
  GetValue(ifp, "t_Phi", "%i", &g_inpset->t_varphi);
  fprintf(red_lfp, "t_Phi: %i\n", g_inpset->t_varphi);

  //Phi: noise in the Partial Momentum Update for GHMC
  GetValue(ifp, "Phi", "%lf", &g_inpset->varphi);
  fprintf(red_lfp, "Phi: %1.4lf\n", g_inpset->varphi);

  //scaling_value: for mapping dt->h for s-AIA integrator
  GetValue(ifp, "scaling_value", "%lf", &g_inpset->scaling_value);
  fprintf(red_lfp, "scaling_value: %1.4lf\n", g_inpset->scaling_value);

  //stepsize_delta: extra parameter needed for stepsize randomization schemes
  GetValue(ifp, "stepsize_delta", "%lf", &g_inpset->stepsize_delta);
  fprintf(red_lfp, "stepsize_delta: %1.4lf\n", g_inpset->stepsize_delta);

  //iter_tune: number of iterations for tuning a stepsize for AR_target during the burn-in stage
  GetValue(ifp, "iter_tune", "%i", &g_inpset->iter_tune);
  fprintf(red_lfp, "iter_tune: %i\n", g_inpset->iter_tune);

  //AR_target: target acceptance rate for tuning a stepsize
  GetValue(ifp, "AR_target", "%lf", &g_inpset->AR_target);
  fprintf(red_lfp, "AR_target: %1.4lf\n", g_inpset->AR_target);

  //delta_AR_target: sensibility of AR_target for tuning a stepsize
  GetValue(ifp, "delta_AR_target", "%lf", &g_inpset->delta_AR_target);
  fprintf(red_lfp, "delta_AR_target: %1.4lf\n", g_inpset->delta_AR_target);

  //fitting factor approach
  GetValue(ifp, "fitting_factor_approach", "%s", &g_inpset->fitting_factor_approach);
  fprintf(red_lfp, "fitting_factor_approach: %s\n", g_inpset->fitting_factor_approach);

  //ATune: use (1) or not (0)
  GetValue(ifp, "opt_param_setting", "%i", &g_inpset->opt_param_setting);
  fprintf(red_lfp, "opt_param_setting: %1.4i\n", g_inpset->opt_param_setting);

  //num_basis_spline: the number of basis in the beta spline approximation
  GetValue(ifp, "num_basis_spline", "%i", &g_inpset->num_basis_spline);
  fprintf(red_lfp, "num_basis_spline: %i\n", g_inpset->num_basis_spline);

  //spline_pol_degree: the degree of the polynomials used in the splines
  GetValue(ifp, "spline_pol_degree", "%i", &g_inpset->spline_pol_degree);
  fprintf(red_lfp, "spline_pol_degree: %i\n", g_inpset->spline_pol_degree);

  //tfin_spl_bas: the final time of the spline basis
  GetValue(ifp, "tfin_spl_bas", "%lf", &g_inpset->tfin_spl_bas);
  fprintf(red_lfp, "tfin_spl_bas: %1.4lf\n", g_inpset->tfin_spl_bas);

  //num_out_cvodes: number of time points returned by the ODE's solver
  GetValue(ifp, "num_out_cvodes", "%i", &g_inpset->num_out_cvodes);
  fprintf(red_lfp, "num_out_cvodes: %i\n", g_inpset->num_out_cvodes);

  //Fixed Initial Counts:
  GetValue(ifp, "initialCountsFixed", "%d", &g_inpset->initialCountsFixed);
  fprintf(red_lfp, "initialCountsFixed: %d\n", g_inpset->initialCountsFixed);

  //S0_fixed_value:
  GetValue(ifp, "S0_fixed_value", "%lf", &g_inpset->S0_fixed_value);
  fprintf(red_lfp, "S0_fixed_value: %1.4lf\n", g_inpset->S0_fixed_value);

  //E0_fixed_value:
  GetValue(ifp, "E0_fixed_value", "%lf", &g_inpset->E0_fixed_value);
  fprintf(red_lfp, "E0_fixed_value: %1.4lf\n", g_inpset->E0_fixed_value);

  //I0_fixed_value:
  GetValue(ifp, "I0_fixed_value", "%lf", &g_inpset->I0_fixed_value);
  fprintf(red_lfp, "I0_fixed_value: %1.4lf\n", g_inpset->I0_fixed_value);

  //dirich_alpha_0:
  GetValue(ifp, "dirich_alpha_0", "%lf", &g_inpset->dirich_alpha_0);
  fprintf(red_lfp, "dirich_alpha_0: %1.4lf\n", g_inpset->dirich_alpha_0);

  //dirich_alpha_1:
  GetValue(ifp, "dirich_alpha_1", "%lf", &g_inpset->dirich_alpha_1);
  fprintf(red_lfp, "dirich_alpha_1: %1.4lf\n", g_inpset->dirich_alpha_1);

  //dirich_alpha_2:
  GetValue(ifp, "dirich_alpha_2", "%lf", &g_inpset->dirich_alpha_2);
  fprintf(red_lfp, "dirich_alpha_2: %1.4lf\n", g_inpset->dirich_alpha_2);

  //dirich_alpha_3:
  GetValue(ifp, "dirich_alpha_3", "%lf", &g_inpset->dirich_alpha_3);
  fprintf(red_lfp, "dirich_alpha_3: %1.4lf\n", g_inpset->dirich_alpha_3);

  //alpha_prior: 
  GetValue(ifp, "alpha_prior", "%s", &g_inpset->alpha_prior);
  fprintf(red_lfp, "alpha_prior: %s\n", g_inpset->alpha_prior);

  //alpha_param_1: mean of the gaussian prior for the rate of transition between the E and I compartments
  GetValue(ifp, "alpha_param_1", "%lf", &g_inpset->alpha_param_1);
  fprintf(red_lfp, "alpha_param_1: %1.4lf\n", g_inpset->alpha_param_1);

  //alpha_param_2: mean of the gaussian prior for the rate of transition between the E and I compartments
  GetValue(ifp, "alpha_param_2", "%lf", &g_inpset->alpha_param_2);
  fprintf(red_lfp, "alpha_param_2: %1.4lf\n", g_inpset->alpha_param_2);

  //gamma_prior: 
  GetValue(ifp, "gamma_prior", "%s", &g_inpset->gamma_prior);
  fprintf(red_lfp, "gamma_prior: %s\n", g_inpset->gamma_prior);

  //gamma_param_1: 
  GetValue(ifp, "gamma_param_1", "%lf", &g_inpset->gamma_param_1);
  fprintf(red_lfp, "gamma_param_1: %1.4lf\n", g_inpset->gamma_param_1);

  //gamma_param_2: 
  GetValue(ifp, "gamma_param_2", "%lf", &g_inpset->gamma_param_2);
  fprintf(red_lfp, "gamma_param_2: %1.4lf\n", g_inpset->gamma_param_2);

  //phi_inv_prior: 
  GetValue(ifp, "phi_inv_prior", "%s", &g_inpset->phi_inv_prior);
  fprintf(red_lfp, "phi_inv_prior: %s\n", g_inpset->phi_inv_prior);

  //phi_inv_param_1:
  GetValue(ifp, "phi_inv_param_1", "%lf", &g_inpset->phi_inv_param_1);
  fprintf(red_lfp, "phi_inv_param_1: %1.4lf\n", g_inpset->phi_inv_param_1);

  //phi_inv_param_2:
  GetValue(ifp, "phi_inv_param_2", "%lf", &g_inpset->phi_inv_param_2);
  fprintf(red_lfp, "phi_inv_param_2: %1.4lf\n", g_inpset->phi_inv_param_2);

  //tau_prior: 
  GetValue(ifp, "tau_prior", "%s", &g_inpset->tau_prior);
  fprintf(red_lfp, "tau_prior: %s\n", g_inpset->tau_prior);

  //tau_param_1:
  GetValue(ifp, "tau_param_1", "%lf", &g_inpset->tau_param_1);
  fprintf(red_lfp, "tau_param_1: %1.4lf\n", g_inpset->tau_param_1);

  //tau_param_2: 
  GetValue(ifp, "tau_param_2", "%lf", &g_inpset->tau_param_2);
  fprintf(red_lfp, "tau_param_2: %1.4lf\n", g_inpset->tau_param_2);

  //num_comp: number of infectious compartments (K) for the SIKR/SEMIKR models
  GetValue(ifp, "num_comp", "%i", &g_inpset->num_comp);
  fprintf(red_lfp, "num_comp: %i\n", g_inpset->num_comp);

  //num_comp_E: number of latent compartments (M) for SEMIKR model
  GetValue(ifp, "num_comp_E", "%i", &g_inpset->num_comp_E);
  fprintf(red_lfp, "num_comp_E: %i\n", g_inpset->num_comp_E);

  //do_under_report: option to do (1) or not (0) the under report  
  GetValue(ifp, "do_under_report", "%i", &g_inpset->do_under_report);
  fprintf(red_lfp, "do_under_report: %i\n", g_inpset->do_under_report);

  //T0_under: first time point where underreport was measured
  GetValue(ifp, "T0_under", "%i", &g_inpset->T0_under);
  fprintf(red_lfp, "T0_under: %i\n", g_inpset->T0_under);

  //U0_under: fraction of total cases that was reported at the first time point (between 0 and 1)
  GetValue(ifp, "U0_under", "%lf", &g_inpset->U0_under);
  fprintf(red_lfp, "U0_under: %1.4lf\n", g_inpset->U0_under);

  //T1_under: second time point where underreport was measured
  GetValue(ifp, "T1_under", "%i", &g_inpset->T1_under);
  fprintf(red_lfp, "T1_under: %i\n", g_inpset->T1_under);

  //U1_under: fraction of total cases that was reported at the second time point (between 0 and 1)
  GetValue(ifp, "U1_under", "%lf", &g_inpset->U1_under);
  fprintf(red_lfp, "U1_under: %1.4lf\n", g_inpset->U1_under);

  //do_grad_desc: option to do (1) or not (0) the gradient descent
  GetValue(ifp, "do_grad_desc", "%i", &g_inpset->do_grad_desc);
  fprintf(red_lfp, "do_grad_desc: %i\n", g_inpset->do_grad_desc);

  //learning_rate: for the gradient descent algorithm
  GetValue(ifp, "learning_rate", "%lf", &g_inpset->learning_rate);
  fprintf(red_lfp, "learning_rate: %1.4lf\n", g_inpset->learning_rate);

  //max_iter: amount of steps that the algorithm follows the gradient
  GetValue(ifp, "max_iter", "%i", &g_inpset->max_iter);
  fprintf(red_lfp, "max_iter: %i\n", g_inpset->max_iter);

  //gammaFixed: 
  GetValue(ifp, "gammaFixed", "%d", &g_inpset->gammaFixed);
  fprintf(red_lfp, "gammaFixed: %d\n", g_inpset->gammaFixed);

  //gamma_fixed_value: 
  GetValue(ifp, "gamma_fixed_value", "%lf", &g_inpset->gamma_fixed_value);
  fprintf(red_lfp, "gamma_fixed_value: %1.4lf\n", g_inpset->gamma_fixed_value);

  //gammaBounded: 
  GetValue(ifp, "gammaBounded", "%d", &g_inpset->gammaBounded);
  fprintf(red_lfp, "gammaBounded: %d\n", g_inpset->gammaBounded);

  //gammaUpper: upper bound for the gamma (inverse of the average infectious period) parameter
  GetValue(ifp, "gammaUpper", "%lf", &g_inpset->gammaUpper);
  fprintf(red_lfp, "gammaUpper: %1.4lf\n", g_inpset->gammaUpper);

  //gammaLower: lower bound for the gamma (inverse of the average infectious period) parameter
  GetValue(ifp, "gammaLower", "%lf", &g_inpset->gammaLower);
  fprintf(red_lfp, "gammaLower: %1.4lf\n", g_inpset->gammaLower);

  //beta_prior: 
  GetValue(ifp, "beta_prior", "%s", &g_inpset->beta_prior);
  fprintf(red_lfp, "beta_prior: %s\n", g_inpset->beta_prior);

  //beta_param_1: param_1 of the prior for the transmission rate distribution in the standard SIR model
  GetValue(ifp, "beta_param_1", "%lf", &g_inpset->beta_param_1);
  fprintf(red_lfp, "beta_param_1: %1.4lf\n", g_inpset->beta_param_1);

  //beta_param_2: param 2 of the prior for the transmission rate distribution in the standard SIR model
  GetValue(ifp, "beta_param_2", "%lf", &g_inpset->beta_param_2);
  fprintf(red_lfp, "beta_param_2: %1.4lf\n", g_inpset->beta_param_2);

  fclose(ifp); // close inputfile stream after reading 

  if((g_inpset->iter_burn_in > 0) && (g_inpset->iter_tune == 0)){
    printf("ERROR: read_input.c - iter_tune cannot be 0 when running burn-in.\n");
    exit(EXIT_FAILURE);
  }

  //////////////////////////
  // opening of files

  //files that are always written
  sprintf(filename, "./output/%s%s", id, "/trajectories.txt");
  tfp = fopen(filename, "w+");

  sprintf(filename, "./output/%s%s", id, "/logPD.txt");
  lpdfp = fopen(filename, "w+");

  sprintf(filename, "./output/%s%s", id, "/hamiltonian.txt");
  hfp = fopen(filename, "w+");

  sprintf(filename, "./output/%s%s", id, "/timestep.txt");
  tsfp = fopen(filename, "w+");

  sprintf(filename, "./output/%s%s", id, "/traj_lengths.txt");
  tlfp = fopen(filename, "w+");

  sprintf(filename, "./output/%s%s", id,"/art.txt");
  artfp = fopen(filename, "w+");

  sprintf(filename, "./output/%s%s", id,"/art_burnin.txt");
  artburnin_fp = fopen(filename, "w+");

  if (g_inpset->do_grad_desc == 1){
    sprintf(filename, "./output/%s%s", id, "/starting_point.txt");
    fp_startingpoint = fopen(filename, "w+");
  }

  sprintf(filename, "./output/%s%s", id, "/scaling_all.txt");//store the scaling factor calculated during the burn-in stage
  scfp = fopen(filename, "w+");

  sprintf(filename, "./output/%s%s", id, "/ode_sol.txt");
  odefp = fopen(filename, "w+");

  if(g_inpset->iter_burn_in > 0){
    sprintf(filename, "./output/%s%s", id, "/finalpoint_tune.txt");
    final_point_tune_fp = fopen(filename, "w+");

    sprintf(filename, "./output/%s%s", id, "/finalpoint_burnin.txt");
    final_point_burnin_fp = fopen(filename, "w+");
  }

  // end of opening of files
  //////////////////////////

  ////////////////
  //stuff depending on inputfile parameters

  //calculate total number of iterations
  iter_total = g_inpset->iter_burn_in + g_inpset->iter_sampling;

  //find number of stages as first chracter of g_inpset->integrator for assignment later
  char strNumStages[256];
  strncpy(strNumStages, g_inpset->integrator, 1);
  strNumStages[1] = 0; // null terminate destination

  //get the number of this run from the ID. these lines find the last number in the ID (which has the format *_<number>).
  int run_number = -1;
  char* copy_id_global = id_global;
  while (*copy_id_global) {
    if (isdigit(*copy_id_global)) {
      run_number = strtol(copy_id_global, &copy_id_global, 10);
    } else {
      copy_id_global++;
    }
  }

  if(run_number < 0){
    fprintf(stderr, "read_input.c - Problem finding out run_number from id_global in read_input.c, run_number is %i, id_global is %s\nPlease put a number in the ID. This is important for properly assigning the seed.", run_number, id_global);
    exit(EXIT_FAILURE);
  }

  ////////////////////////////////////
  //initialize random number generator

  //read the environment variables - GSL_RNG_TYPE and GSL_RNG_SEED uses these values to set the corresponding library variables gsl_rng_default and gsl_rng_default_seed
  gsl_rng_env_setup();

  //depending on input generate a random seed (if seed==0), otherwise use the given seed (and add the number of the run)
  if (g_inpset->seed==0){
    seed = time (NULL) * getpid();
  } else{
    seed = g_inpset->seed + run_number - 1;
    //-1 so it stays the same for the first run (i.e. run_number = 1)
  }

  fprintf(red_lfp, "seed for calculation: %li\n", seed);

  gsl_rng_default_seed = seed;
  T = gsl_rng_default;
  g_rng = gsl_rng_alloc(T);
  gsl_rng_set(g_rng, seed);

  //reading model: define number of parameters, functions for prior, likelihood, grad, hessian, transformations, read observations and initialize state
  if (strcmp(g_inpset->model, "SIKR_Incidence") == 0) {
    ReadModelSIKR_Incidence();
  }else if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0) {
    ReadModelSEMIKR_Incidence();
  }else if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0) {
    ReadModelSIKR_stand_Incidence();
  } else if (strcmp(g_inpset->model, "SIR_standard_Prevalence") == 0) {
    ReadModelSIKR_stand_Prevalence();
  }else{
    fprintf(stderr, "ERROR in read_input.c - no match found in if-elseif for model. model is %s\n", g_inpset->model);;
    exit(EXIT_FAILURE);
  }

  return 0;
}
