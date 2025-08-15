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

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>

#include "hmc.h"
#include "Globals.h"
#include "integrators.h"
#include "adapt_integrator.h" 
#include "Definitions.h"
#include "utils.h"
#include "sir_seir_common_functions.h"
#include "sir_model_functions.h" // Cote: check whether this can be done in a different way (need function printODE() from sir_model.c)
#include "sir_functions.h"
#include "seir_functions.h"
#include "sir_stand_functions.h"
#include "sir_stand_model_prevalence_synthetic.h"
#define S0_constant_1 1
#define I0_constant_1 1

void Hamiltonian(double logpd, TypeState **st, int dimension, double stepsize, TypeHamiltonian *H, int iter);

//florian--wayne: metropolis random walk
void HamiltonianMetropolis(double logpd, TypeState **st, int dimension, double stepsize, TypeHamiltonian *H, int iter);

// pointer to the Hamiltonian function (shadow or true)
void (*HamiltonianF)(double /*logpd*/, TypeState **/*state*/, int /*dimension*/, double /*stepsize*/, TypeHamiltonian *H, int iter);

//florian--wayne: pointer to the Hamiltonian function to be called in the PMU step. This is different for Metropolis as the kinetic energy in the Metropolis RW is set to 0
void (*HamiltonianPMU)(double /*logpd*/, TypeState **/*state*/, int /*dimension*/, double /*stepsize*/, TypeHamiltonian *H, int iter);

static void MomFlip(double *mom, int dimension);
static void MomFlipNo(double *mom, int dimension);
/* momentum flip function */
static void (*MomentumFlipF)(double */*mom*/, int /*dimension*/);

void MetropolisTest(double dH, bool *acc, bool *acc_ed, bool *prev_acc, double *acc_prob_prev_inv, int *numRF);

void AlwaysAcceptTest(double dH, bool *acc, bool *acc_ed, bool *PREV_ACC, double *acc_prob_prev_inv, int *numRF); //added by Felix, 12.07.2018

static void (*AcceptanceTest)(double /*dH*/, bool */*acc*/, bool */*acc_ed*/, bool */*prev_acc*/, double */*accProbPrevInv*/, int */*numRF*/);

static void PMU(double phi, TypeState **st, double *mom_old, int dimension, double logpd, double ss, bool *Accepted, int *numAM, TypeHamiltonian *H, int iter);

static void (*MomentumUpdateF)(double /*phi*/, TypeState **st, double */*mom_old*/, int /*dimension*/, double /*logpd*/, double /*ss*/, bool */*mom_acc*/, int */*numAM*/, TypeHamiltonian *H, int iter);

static void NoWeights(TypeHamiltonian Ham);
static void (*CalculateWeights)(TypeHamiltonian Ham);

/* book keeping */
static int num_proposed, num_accepted, num_accepted_ed, num_prop_mom, num_acc_mom, num_red_flips;
static int num_acc_tune, iter_tuning, iter_tune_extra = 0; // for tuning the time step
static double ar_tuning; //for tuning the time step
static double dt_propos_tuning; //for tuning the time step
static int flag_dt_tuned = 0; // flag for checking if a stepsize has been tuned
static double dt_tuned; //tuned time step for the scaling burn-in
static bool ACCEPTED, ACCEPTED_ed, MOM_ACCEPTED, PREV_ACCEPTED;
static double acc_prob_prev_inv;
static int iter;

//variables for normalization time
double t0 = 0.;
double t1 = 0.;
double time_iter = 0.;
double NormalizationTime = 0.;

//variables for TUNE method
double ar_last = 0.;

TypeState *state_old[5] = {NULL,NULL,NULL,NULL,NULL};

// Hamiltonian()
// mom[]: momenta
// dimension: number of elements in the momenta vector
void Hamiltonian(double logpd, TypeState **st, int dimension, double stepsize, TypeHamiltonian *H, int iter) {

  H->pot = -logpd;
  
  H->kin = 0.5*cblas_ddot(dimension, st[curr]->mom, 1, st[curr]->mom, 1);

  H->ham = H->pot + H->kin;

}

// florian--wayne:
/*************************************************************
 *     hamiltonian for metropolis RW
 *     kinetic energy is 0
 *     mom[]: momenta
 *     dimension: number of elements in the momenta vector
 ************************************************************/

void HamiltonianMetropolis(double logpd, TypeState **st, int dimension, double stepsize, TypeHamiltonian *H, int iter) {

  H->pot = -logpd;

  H->kin = 0;

  H->ham = H->pot;

}

//momentum flip - automatic
static void MomFlip(double *mom, int dimension) {

  double alpha = -1.0;
  cblas_dscal(dimension, alpha, mom, 1);

}

//momentum flip - no flip
static void MomFlipNo(double *mom, int dimension) {
}

//metropolis test
void MetropolisTest(double dH, bool *acc, bool *acc_ed, bool *prev_acc, double *acc_prob_prev_inv, int *numRF) {

  double exp_dH,
  r = -1;

  if (dH <= 0) {
    *acc = TRUE;
    *acc_ed = TRUE;
  } else {
    *acc_ed = FALSE;
    r = gsl_rng_uniform(g_rng);
    exp_dH = exp(-dH);
    if (r < exp_dH) {
      *acc = TRUE;
    } else {
      *acc = FALSE;
    }
  }
}

//AlwaysAcceptTest
void AlwaysAcceptTest(double dH, bool *acc, bool *acc_ed, bool *prev_acc, double *acc_prob_prev_inv, int *numRF) {

  *acc = TRUE;

} //end of function AlwaysAcceptTest

//partial momentum update step
static void PMU(double phi, TypeState **st, double *mom_old, int dimension, double logpd, double ss, bool *Accepted, int *numAM, TypeHamiltonian *H, int iter) {

  int i;
  double aa = 0;
  double r = -1;
  double rphi;
  double rphicomp;
  double u[dimension];
  int INCX = 1,
  INCY = 1;

  for (i = 0; i < dimension; i++){
    u[i] = gsl_ran_gaussian_ziggurat(g_rng, 1.);
  }

  //check if phi is adaptive
  if ((g_inpset->t_varphi == 4) || (g_inpset->t_varphi == 5) || (g_inpset->t_varphi == 6)) {
    aa = cblas_ddot(dimension, mom_old, INCX, mom_old, INCY) - cblas_ddot(dimension, u, INCX, u, INCY);
  }
  
  //asign value if phi is adaptive
  if (g_inpset->t_varphi == 4) { //adaptive maximal MD AR
    if (aa < 0) {
      phi = 1 - phi;
    } else if (aa == 0) {
      r = gsl_rng_uniform(g_rng);
      if (r < 0.5) {
        phi = 1 - phi;
      }
    }
  } else if (g_inpset->t_varphi == 5) { //adaptive minimal MD AR
    if (aa > 0) {
      phi = 1 - phi;
    } else if (aa == 0) {
      r = gsl_rng_uniform(g_rng);
      if (r < 0.5) {
        phi = 1 - phi;
      }
    }
  } else if (g_inpset->t_varphi == 6) {}//adaptive optimal MD AR
    
  /// perform the momentum update
  if ( DOUBLE_EQ(phi,1) ) {
    cblas_dcopy(dimension, u, INCX, st[curr]->mom, INCY);//mom=u
  } else {
    rphi = sqrt(phi);
    rphicomp = sqrt(1-phi);
    for (i = 0; i < dimension; i++){ st[curr]->mom[i] =  rphicomp * mom_old[i] + rphi * u[i]; }
  }
  
  *numAM += 1;
  *Accepted = TRUE;
  
  //florian--wayne: 
  HamiltonianPMU(logpd, st, dimension, ss, H, iter);

}

//no weights needed
static void NoWeights(TypeHamiltonian Ham) {

}

//hmc_update()
//Update a state of the Markov chain
static void HMCUpdate(int iter) {

  TypeHamiltonian Ham_curr;

  fflush(stdout);

  double logpd = 0.,
  logpd_old = 0.,
  deltaH = 0.;
  int i;

  ////////////////////////////
  //MODEL PARAMETERS

  logpd_old = logpd = deltaH = 0.0;
  
  // keep the old position, momenta, gradient
  CopyState(state_old[curr], state[curr], kNumParam);

  logpd_old = Prior() + Loglikelihood() + LogJacobian();

  num_prop_mom++;


  // printf("Position and Momentum for gamma prior momentum update: %f, %f\n", state[curr]->pos[1], state[curr]->mom[1]);
  //update momenta and calculate current Hamiltonian; Hamiltonian here does not inlcude Jacobian, as it is constant for the PMMC step
  MomentumUpdateF(varphi[iter], state, state_old[curr]->mom, kNumParam, logpd_old, stepsize[iter], &MOM_ACCEPTED, &num_acc_mom, Ham, iter);

  //store the current Ham into the global Ham; it used some values from the previous state of the Ham
  Ham_curr = *Ham;

  // keep the current momenta, needed in case of rejection
  memcpy(state_old[curr]->mom, state[curr]->mom, kNumParam*sizeof(double));

  fprintf(tsfp, "%f\n", stepsize[iter]);
  fprintf(tlfp, "%d\n", traj_length[iter]);

  // For fixed gamma in SIR-like models we need to supress the momentum from the Hamiltonian calculation.
  // For the logposterior it is handled internaly.
  if(g_inpset->gammaFixed == 1){
    if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0 || strcmp(g_inpset->model, "SIKR_Incidence") == 0){
      state[curr]->mom[0] = 0;
    }
    if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0){
      state[curr]->mom[1] = 0;
    }
  }
  if(S0_constant_1 == 1){
    if(strcmp(g_inpset->model, "SIR_standard_Prevalence") == 0){
      state[curr]->mom[1] = 0;
    }
  }
  if(I0_constant_1 == 1){
    if(strcmp(g_inpset->model, "SIR_standard_Prevalence") == 0){
      state[curr]->mom[2] = 0;
    }
  }
  if(g_inpset->initialCountsFixed == 1){
    if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0 || strcmp(g_inpset->model, "SIKR_Incidence") == 0){
      // state[curr]->mom[0] = 0;
    }
    if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0){
      state[curr]->mom[2] = 0;
      state[curr]->mom[3] = 0;
      state[curr]->mom[4] = 0;
    }
  }

  // run the trajectory / make the proposal
  MDMove(Gradient, kNumParam, state, stepsize[iter], &(g_coeff[iter][0]), traj_length[iter]);
  
  // SIKR/SEMIKR models: Check if the proposed point fulfils that all states sum to N (the size of the poulation.)
  if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0){
    InvTransf();
    double R_compartment = num_pop - state[curr]-> pos[1] - state[curr]-> pos[2];
    if (R_compartment < 0){
      state[curr]-> pos[1] = num_pop - state[curr]-> pos[2];
    }
    if(R_compartment > num_pop){
      if(state[curr]-> pos[1] < 0){
        state[curr]-> pos[1] = 0;
      } else if (state[curr]-> pos[2] < 0){
        state[curr]-> pos[2] = 0;
      }    
    }
    Transf();
  } else if (strcmp(g_inpset->model, "SIKR_Incidence") == 0){
    InvTransf();
    double R_compartment = num_pop - state[curr]-> pos[1] - state[curr]-> pos[2];
    if (R_compartment < 0){
      state[curr]-> pos[1] = num_pop - state[curr]-> pos[2] ;
    }
    if(R_compartment > num_pop){
      if(state[curr]-> pos[1] < 0){
        state[curr]-> pos[1] = 0;
      } else if (state[curr]-> pos[2] < 0){
        state[curr]-> pos[2] = 0;
      }    
    }
    Transf();
  } else if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0){
    InvTransf();
    double R_compartment = num_pop - state[curr]-> pos[2] - state[curr]-> pos[3] - state[curr]-> pos[4];
    if (R_compartment < 0){
      state[curr]-> pos[2] = num_pop - state[curr]-> pos[3] - state[curr]-> pos[4];
    }
    if(R_compartment > num_pop){
      if(state[curr]-> pos[2] < 0){
          state[curr]-> pos[2] = 0;
      } else if (state[curr]-> pos[3] < 0){
          state[curr]-> pos[3] = 0;
      } else {
          state[curr]-> pos[4] = 0;
      }        
    }
    Transf();
  }

  num_proposed++;

  logpd = Prior() + Loglikelihood() + LogJacobian();

  HamiltonianF(logpd, state, kNumParam, stepsize[iter], Ham, iter);//proposed Ham is stored in Ham global

  deltaH = (Ham->ham) - (Ham_curr.ham);
  if(deltaH < 0){
    num_proposals_dHneg_arr[iter] = 1; //store the number of proposals with negative energy error
  }else{
    num_proposals_dHneg_arr[iter] = 0;
  }

  // acceptance test
  AcceptanceTest(deltaH, &ACCEPTED, &ACCEPTED_ed, &PREV_ACCEPTED, &acc_prob_prev_inv, &num_red_flips);

  //HessianSIKR_Incidence();

  if (ACCEPTED) {
    num_accepted++;
    accept_arr[iter] = 1;
    if((flag_dt_tuned == 0) && (iter < g_inpset->iter_burn_in)){ //during the burn-in, AR is used both for tuning a time step and for calculating fitting factors
      num_acc_tune++;
    }
    if (ACCEPTED_ed){
      num_accepted_ed++;
    }
    CalculateWeights(*Ham);
    fprintf(hfp, "%f\n", Ham->ham);
  } else {
    accept_arr[iter] = 0;
    CopyState(state[curr], state_old[curr], kNumParam);
    logpd = logpd_old;
    *Ham = Ham_curr;
    MomentumFlipF(state[curr]->mom, kNumParam);
    CalculateWeights(*Ham);
    fprintf(hfp, "%f\n", Ham->ham);
  }

  // tuning the time step during the burn-in stage
  if(flag_dt_tuned == 0){  
    if((iter+1)%(g_inpset->iter_tune) == 0){
      ar_tuning = (double)(num_acc_tune)/(g_inpset->iter_tune); // corrected acceptance rate
      num_acc_tune = 0;
      iter_tuning = iter+1+iter_tune_extra;
      if(ar_tuning < g_inpset->AR_target - g_inpset->delta_AR_target){ //if the AR is smaller than the target, lower a stepsize
        dt_propos_tuning = stepsize[iter] - stepsize[iter]/10 * (1 + 1./iter_tuning);
        for (i = iter + 1; i < g_inpset->iter_burn_in + g_inpset->iter_sampling; i = i + 1){
          stepsize[i] = stepsize[iter] - stepsize[iter]/10 * (1 + 1./iter_tuning);
        }
      }else if(ar_tuning > g_inpset->AR_target + g_inpset->delta_AR_target){ //if the AR is bigger than the target, increase a stepsize
        dt_propos_tuning = stepsize[iter] + stepsize[iter]/10 * (1 + 1./iter_tuning);
        for (i = iter + 1; i < g_inpset->iter_burn_in + g_inpset->iter_sampling; i = i + 1){
          stepsize[i] = stepsize[iter] + stepsize[iter]/10 * (1 + 1./iter_tuning);
        }
      }else{
        iter_tuning = iter+1+iter_tune_extra;
        flag_dt_tuned = 1;
        for (i = 0; i < kNumParam; i = i + 1){
          InvTransf();
          fprintf(final_point_tune_fp, "%f ", state[curr]->pos[i]);
          Transf();
        }
        //fprintf(final_point_tune_fp, "\n");
        dt_tuned = stepsize[iter];
        for(i = 0; i < g_inpset->iter_burn_in + g_inpset->iter_sampling; i = i + 1){
          stepsize[i] = dt_tuned;
        }
        //print the last AR tuned
        printf("Tuned AR: %lf\n", ar_tuning);
        fprintf(red_lfp, "Tuned AR: %lf\n", ar_tuning);
        fprintf(artfp, "AR_tuned %lf\n", ar_tuning);

        //print the tuned dt
        printf("Tuned dt: %lf\n", dt_tuned);
        fprintf(red_lfp, "Tuned dt: %lf\n", dt_tuned);
        fprintf(artfp, "dt_tuned %lf\n", dt_tuned);

        //number of iterations needed for tuning
        printf("Iterations needed for tuning: %i\n", iter_tuning);
        fprintf(red_lfp, "Iterations needed for tuning: %i\n", iter_tuning);
        fprintf(artfp, "iter_tuning %i\n", iter_tuning);
      }
    }
  }

  fprintf(lpdfp, "%f\n", logpd);

} //end of function "HMCUpdate"

//hmc()
double HMC()
{

  int i;
  int ii;
  int tot_iter = g_inpset->iter_burn_in + g_inpset->iter_sampling;
  clock_t start0, end0, startT, endT;
  double elapsed0, elapsedT;

  int accept_burnin = 0;
  int accept_dHpos_burnin = 0;
  int accept_dHneg_burnin = 0;
  int num_proposal_burnin = 0;
  int num_proposal_dHpos_burnin = 0;
  double dH_avg_burnin = 0.;
  double dH_pos_avg_burnin = 0.;
  double dH2_avg_burnin = 0.;
  int L_target = 1;
  double phi_target = 1;

  ACCEPTED = ACCEPTED_ed = MOM_ACCEPTED = FALSE;
  PREV_ACCEPTED = FALSE;

  acc_prob_prev_inv = 0;

  num_proposed = 0;
  num_accepted = 0;
  num_accepted_ed = 0;
  num_prop_mom = 0;
  num_acc_mom = 0;
  num_red_flips = 0;

  printf("The dimension of the system is %i\n", kNumParam);
  fprintf(red_lfp, "The dimension of the system is %i\n", kNumParam);
  fprintf(artfp, "D %i\n", kNumParam);

  // allocate memory for the old state keeping
  for (i = back2; i <= forw2; i++) {
    SAlloc(&(state_old[i]), kNumParam);
  }

  // setup of functions
  if (strcmp(g_inpset->method, "HMC") == 0){

    AcceptanceTest = MetropolisTest;
    MomentumFlipF = MomFlipNo;

  } else if (strcmp(g_inpset->method, "GHMC") == 0){
    
    AcceptanceTest = MetropolisTest;
    MomentumFlipF = MomFlip;

  }else{

    fprintf(stderr, "hmc.c - no match found in if-elseif for method. Exiting.");
    exit(EXIT_FAILURE);

  } //end of if-elseif for method

  CalculateWeights = NoWeights;

  //allocating memory for stepsize, traj_length, g_coeff, g_SH_coeff, varphi
  //stepsize = (double *)malloc(tot_iter*sizeof(double)); //allocate and init with zeroes
  //traj_length = (int *)malloc(tot_iter*sizeof(int)); //allocate and init with zeroes
  //g_coeff = MatrixAlloc(tot_iter, 8); //allocate and init with zeroes
  //varphi = (double *)malloc(tot_iter*sizeof(double)); //allocate and init with zeroes

  //allocating memory for stepsize, traj_length, g_coeff, g_SH_coeff, varphi
  stepsize = (double *)malloc(tot_iter*sizeof(double)); //allocate and init with zeroes
  traj_length = (int *)malloc(tot_iter*sizeof(int)); //allocate and init with zeroes
  g_coeff = MatrixAlloc(tot_iter, 8); //allocate and init with zeroes
  varphi = (double *)malloc(tot_iter*sizeof(double)); //allocate and init with zeroes

  //allocation memory for acceptance/rejection arrays
  accept_arr = (int *)malloc(tot_iter*sizeof(int)); //allocate and init with zeroes
  dH_array = (double *)malloc(tot_iter*sizeof(double)); //allocate and init with zeroes
  num_proposals_dHneg_arr = (int *)malloc(tot_iter*sizeof(int)); //allocate and init with zeroes

  Transf();//transformation of the model parameters at the beginning
  start0 = clock();

  ////////////////////////////////////////////////////////////
  // burn-in phase

  HamiltonianF = Hamiltonian;
  HamiltonianPMU = Hamiltonian;
  MomentumUpdateF = PMU;
  
  ///////////////////////////////////////////////////////
  //assign quantities for burn-in phase
  //1. stepsize = 1/D, 2. number of integration steps L = 1, 3. integrator coefficients of 1-stage VV, 4. Shadow Hamiltonian coefficients, 5. phi chosen from the input

double delta_t_start = 0.;

if(g_inpset->iter_burn_in == 0){
  flag_dt_tuned = 1;
  printf("The number of burn-in iterations chosen is 0, only production stage.\n");
}else{
  if(g_inpset->opt_param_setting == 1){ //option to calculate optimal parameters setting
    
    //perform tuning
    flag_dt_tuned = 0;

    //assign step size for tuning
    delta_t_start = 1. / kNumParam; // dt for starting the burn-in
    if ((strcmp(g_inpset->model, "SIR_standard_Incidence") == 0) || (strcmp(g_inpset->model, "SIKR_Incidence") == 0) || (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0)){
      if(strcmp(g_inpset->method, "GHMC") == 0){
        delta_t_start = 0.002;
      }else{
        delta_t_start = 0.001;
      }
    }

    //assign optimal phi for PMU
    phi_random_scheme = 2;

  }else if(g_inpset->opt_param_setting == 0){

    //do not perform tuning
    flag_dt_tuned = 1;

    //assign step size from input
    delta_t_start = g_inpset->stepsize;

    //assign phi for PMU from input
    phi_random_scheme = g_inpset->t_varphi;
    phi_target = g_inpset->varphi;

  }else{
    fprintf(stderr, "hmc.c - Error in phi assignation for the burn-in, opt_param_setting does not fit the if-else, it has to be 0 or 1 and it is %i.\n", g_inpset->opt_param_setting);
    exit(EXIT_FAILURE);
  }

  //assignation for burnin
  for (i = 0; i < g_inpset->iter_burn_in; i++) {
    stepsize[i] = delta_t_start;
    traj_length[i] = 1;
    g_coeff[i][0] = 0.5; //1-stage VV
    g_coeff[i][1] = 1.0; //1-stage VV
  }

  fprintf(artburnin_fp, "dt_start %lf\n", delta_t_start);
  fprintf(artburnin_fp, "L_burnin %i\n", 1);


  fnAssignPhi(varphi, phi_target, phi_random_scheme, 0, g_inpset->iter_burn_in - 1, kNumParam); //choice of phi given in the inputfile, if HMC is set to be constant equal to 1 in the fnAssignPhi function in integrators.c

  if(flag_dt_tuned == 0){
    
    printf("The first time step for the tuning is %lf\n", stepsize[0]);
    fprintf(red_lfp, "The first time step for the tuning is %lf\n", stepsize[0]);

    printf("The trajectory used for the tuning is %d\n", traj_length[0]);
    fprintf(red_lfp, "The trajectory used for the tuning is %d\n", traj_length[0]);
  }

  //assign MDMove for the burn-in, i.e. 1-stage Velocity Verlet
  printf("Assigning MDMove for the burn-in stage\n");
  fprintf(red_lfp, "Assigning MDMove for the burn-in stage\n");

  MDMove = VV;
  printf("MDMove assigned: VV\n");

  iter = 0;

  printf("------------------------------\n");
  printf("---- Starting tuning phase ---\n");
  printf("------------------------------\n");
  fprintf(red_lfp, "------------------------------\n");
  fprintf(red_lfp, "---- Starting tuning phase ---\n");
  fprintf(red_lfp, "------------------------------\n");

  //tuning the time step
  while (flag_dt_tuned == 0) {

    HMCUpdate(iter);

    iter++;

    if((iter == g_inpset->iter_burn_in) && (flag_dt_tuned == 0)){
      iter = iter - g_inpset->iter_tune;
      for (i = iter; i < g_inpset->iter_burn_in + g_inpset->iter_sampling; i = i + 1){
          stepsize[i] = dt_propos_tuning;
      }
      iter_tune_extra += g_inpset->iter_tune;
    }
  }

  printf("Tuning finished.\n");
  printf("------------------------------\n");
  fprintf(red_lfp, "Tuning finished\n");
  fprintf(red_lfp, "------------------------------\n");

  printf("------------------------------\n");
  printf("--- Starting burn-in phase ---\n");
  printf("------------------------------\n");
  fprintf(red_lfp, "------------------------------\n");
  fprintf(red_lfp, "--- Starting burn-in phase ---\n");
  fprintf(red_lfp, "------------------------------\n");

  //burn-in
  for (iter = 0; iter < g_inpset->iter_burn_in; iter++) {
    
    if (iter % ((g_inpset->iter_burn_in) / 10) == 0 && iter > 0){
      printf("%d%% of burn-in completed\n", 100 * iter / g_inpset->iter_burn_in);
      printf("------------------------------\n");
    }

    HMCUpdate(iter);

    //write trajectory file
    if ((iter+1)%(g_inpset->thinning) == 0) {
      //main trajectory
      for (i = 0; i < kNumParam; i++) { fprintf(tfp, "%f ", state[curr]->pos[i]); }
      fprintf(tfp, "\n");

      //SIKR Incidence
      if (strcmp(g_inpset->model, "SIKR_Incidence") == 0) {
        printODE_SIR(odefp);
      }
      //SEMIKR synthetic
     if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0) {
        printODE_SEIR(odefp);
      }
      //SIR standard Incidence
      if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0) {
        printODE_SIR_stand(odefp);
      }
      //SIR standard Prevalence
      if (strcmp(g_inpset->model, "SIR_standard_Prevalence") == 0) {
        printODE_SIR_stand(odefp);
      }
    } //end of "if"-clause for iter
  }

  printf("Burn-in finished.\n");

  end0 = clock();
  elapsed0 = ((double)(end0-start0))/CLOCKS_PER_SEC;

  //store the final point of the burn-in stage
  InvTransf();
  for (i = 0; i < kNumParam; i++) { fprintf(final_point_burnin_fp, "%f ", state[curr]->pos[i]); }
  //fprintf(final_point_burnin_fp, "\n");
  Transf();

  //SIKR_Prevalence
  if (strcmp(g_inpset->model, "SIKR_Prevalence") == 0) {
    printODE(odefp);
  }
  //SIKR Incidence
  if (strcmp(g_inpset->model, "SIKR_Incidence") == 0) {
    printODE_SIR(odefp);
  }
  //SEMIKR synthetic
  if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0) {
    printODE_SEIR(odefp);
  }
  //SIR standard Incidence
  if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0) {
    printODE_SIR_stand(odefp);
  }
  //SIR standard Prevalence
  if (strcmp(g_inpset->model, "SIR_standard_Prevalence") == 0) {
    printODE_SIR_stand(odefp);
  }

  //time elapsed in burn-in
  printf("Elapsed time burnin: %lf sec\n", elapsed0);
  fprintf(red_lfp, "Elapsed time burnin: %lf\n", elapsed0);
  fprintf(artfp, "TimeBurnin %lf\n", elapsed0);

  printf("--- Printing statistics for burn-in. --\n");
  fprintf(red_lfp, "--- Printing statistics to art_burnin file. --\n");

  if(flag_dt_tuned == 0){
    fprintf(stderr, "hmc.c - Stepsize cannot be tuned with the burn-in iterations chosen in input. Exiting.\n");
    fprintf(stderr, "The last AR obtained is: %lf\n", ar_tuning);
    fprintf(stderr, "The last dt used is: %lf\n", stepsize[iter-1]);
    exit(EXIT_FAILURE);
  }

  //number of iterations in the burn-in with dt_tuned - N
  num_proposal_burnin = g_inpset->iter_burn_in;
  fprintf(red_lfp, "Proposals burn-in: %i\n", num_proposal_burnin);
  fprintf(artburnin_fp, "N_burnin %i\n", num_proposal_burnin);

  //number of accepted proposals during the burn-in - N_acc
  for(ii=0; ii < g_inpset->iter_burn_in; ii++){
    accept_burnin += accept_arr[ii];
  }
  //printf("Accepted proposals burn-in: %i\n", accept_burnin);
  fprintf(red_lfp, "Accepted proposals burn-in: %i\n", accept_burnin);
  fprintf(artburnin_fp, "N_acc_burnin %i\n", accept_burnin);

  //number of accepted proposals with negative energy error during the burn-in - N-
  for(ii=0; ii < g_inpset->iter_burn_in; ii++){
    accept_dHneg_burnin += num_proposals_dHneg_arr[ii];
  }
  //printf("Accepted proposals with dH<0 burn-in: %i\n", accept_dHneg_burnin);
  fprintf(red_lfp, "Accepted proposals with dH<0 burn-in: %i\n", accept_dHneg_burnin);
  fprintf(artburnin_fp, "N_dHneg_burnin %i\n", accept_dHneg_burnin);

  //number of proposals with positive energy error during the burn-in - N+
  num_proposal_dHpos_burnin = num_proposal_burnin - accept_dHneg_burnin;
  fprintf(red_lfp, "Proposals with dH>0 burn-in: %i\n", num_proposal_dHpos_burnin);
  fprintf(artburnin_fp, "N_dHpos_burnin %i\n", num_proposal_dHpos_burnin);

  //number of accepted proposals with positive energy error during the burn-in - N+_acc
  accept_dHpos_burnin = accept_burnin-accept_dHneg_burnin;
  fprintf(red_lfp, "Accepted proposals with dH>0 burn-in: %i\n", accept_dHpos_burnin);
  fprintf(artburnin_fp, "N_acc_dHpos_burnin %i\n", accept_dHpos_burnin);

  //average acceptance rate over the burn-in stage, excluding the initial tuning stage - AR = N_acc / N
  double AR_burnin = (double)accept_burnin/num_proposal_burnin;
  printf("Acceptance rate burn-in: %lf\n", AR_burnin);
  fprintf(red_lfp, "Acceptance rate burn-in: %lf\n", AR_burnin);
  fprintf(artburnin_fp, "AR_burnin %lf\n", AR_burnin);

  //average filtered acceptance rate over the burn-in stage (only cases accepted using the energy difference over proposals with positive energy error) - filtAR = N+_acc / N+
  double AR_filt_burnin = (double)accept_dHpos_burnin/num_proposal_dHpos_burnin;
  fprintf(red_lfp, "Filtered AR burn-in: %lf\n", AR_filt_burnin);
  fprintf(artburnin_fp, "AR_filt_burnin %lf\n", AR_filt_burnin);

  //average new filtered acceptance rate over the burn-in stage (only cases accepted using the energy difference over the total number of proposals) - filtAR_new = N+acc / N
  double AR_filtNew_burnin = (double)accept_dHpos_burnin/num_proposal_burnin;
  fprintf(red_lfp, "New Filtered AR burn-in: %lf\n", AR_filtNew_burnin);
  fprintf(artburnin_fp, "AR_filtNew_burnin %lf\n", AR_filtNew_burnin);

  //average energy error over the burn-in stage (after tuning) - dH = sum_i^N dH_i / N
  for(ii=0; ii < g_inpset->iter_burn_in; ii++){
    dH_avg_burnin += dH_array[ii];
  }
  dH_avg_burnin = dH_avg_burnin/num_proposal_burnin;
  fprintf(red_lfp, "Average Energy Error burn-in: %lf\n", dH_avg_burnin);
  fprintf(artburnin_fp, "dH_avg_burnin %lf\n", dH_avg_burnin);

  //average positive energy error over the burn-in stage (after tuning, averaged over the total number of proposals) - dH+ = sum_i^N dH+_i / N
  for(ii=0; ii < g_inpset->iter_burn_in; ii++){
    if(dH_array[ii] > 0){
      dH_pos_avg_burnin += dH_array[ii];
    }
  }
  dH_pos_avg_burnin = dH_pos_avg_burnin/num_proposal_burnin;
  fprintf(red_lfp, "Average positive Energy Error burnin: %lf\n", dH_pos_avg_burnin);
  fprintf(artburnin_fp, "dH_pos_avg_burnin %lf\n", dH_pos_avg_burnin);

  //average positive energy error over the burnin stage (averaged over the proposals with positive energy error) - H+_av = sum_1^N+ dH+_i / N+
  double H_pos_av_burnin = dH_pos_avg_burnin/num_proposal_dHpos_burnin;
  fprintf(red_lfp, "H+_av burnin: %lf\n", H_pos_av_burnin);
  fprintf(artburnin_fp, "H_pos_av_burnin %lf\n", H_pos_av_burnin);

  //average positive energy error over the burnin stage (averaged over the accepted proposals with positive energy error) - H+_acc_av = sum_1^N+_acc dH+_i / N+_acc
  double H_pos_acc_av_burnin = dH_pos_avg_burnin/accept_dHpos_burnin;
  fprintf(red_lfp, "H+_acc_av burnin: %lf\n", H_pos_acc_av_burnin);
  fprintf(artburnin_fp, "H_pos_acc_av_burnin %lf\n", H_pos_acc_av_burnin);

  //average squared energy error over the burn-in stage (after tuning) - dH^2 = sum_i^N dH^2_i / N
  for(ii=0; ii < g_inpset->iter_burn_in; ii++){
    dH2_avg_burnin += dH_array[ii]*dH_array[ii];
  }
  dH2_avg_burnin = dH2_avg_burnin/num_proposal_burnin;
  fprintf(red_lfp, "Average Squared Energy Error burn-in: %lf\n", dH2_avg_burnin);
  fprintf(artburnin_fp, "dH2_avg_burnin %lf\n", dH2_avg_burnin);

  //standard dev of energy error over the burn-in stage
  double dH_sd_burnin = sqrt(dH2_avg_burnin - dH_avg_burnin*dH_avg_burnin);
  fprintf(red_lfp, "StdDevEnergyError burn-in: %lf\n", dH_sd_burnin);
  fprintf(artfp, "dH_sd_burnin %lf\n", dH_sd_burnin);

  //mean of angle phi over burn-in
  double mean_phi_burnin = 0;
  for (ii=0; ii < g_inpset->iter_burn_in; ii++){ mean_phi_burnin += varphi[ii]; }
  mean_phi_burnin = mean_phi_burnin / num_proposal_burnin;
  fprintf(red_lfp, "Average angle phi during burn-in: %lf\n", mean_phi_burnin);
  fprintf(artburnin_fp, "mean_phi_burnin %lf\n", mean_phi_burnin);

  //standard deviation of angle phi during burn-in
  double sd_phi_burnin = 0;
  for (ii=0; ii < g_inpset->iter_burn_in; ii++){ sd_phi_burnin += pow(varphi[ii] - mean_phi_burnin, 2.0); }
  sd_phi_burnin = sqrt(sd_phi_burnin / num_proposal_burnin);
  fprintf(red_lfp, "standard deviation phi during burn-in: %lf\n", sd_phi_burnin);
  fprintf(artburnin_fp, "sd_phi_burnin %lf\n", sd_phi_burnin);

  //calculate and store the fitting factor with all the approaches
  fprintf(scfp, "FiltAR %lf\n", fnCalculateFittingFactor("filtAR", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "FiltAR_new %lf\n", fnCalculateFittingFactor("filtAR_new", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "Gupta %lf\n", fnCalculateFittingFactor("Gupta", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "TightUpper %lf\n", fnCalculateFittingFactor("TightUpper", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "dH_average %lf\n", fnCalculateFittingFactor("dH_average", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "dHpos_average %lf\n", fnCalculateFittingFactor("dHpos_average", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "H_pos_avg %lf\n", fnCalculateFittingFactor("H_pos_avg", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "H_pos_acc_avg %lf\n", fnCalculateFittingFactor("H_pos_acc_avg", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "Gupta_pos %lf\n", fnCalculateFittingFactor("Gupta_pos", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "dH2_average %lf\n", fnCalculateFittingFactor("dH2_average", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "mean %lf\n", fnCalculateFittingFactor("mean", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
  fprintf(scfp, "mean_new %lf\n", fnCalculateFittingFactor("mean_new", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin));
} //end of the loop for iter burn-in > 0
// end of the loop for the real burn-in phase

///////////////////////////
//start production phase loop

//florian: determine number of stages (normal and latent). Needed to assign coefficients of shadow Hamiltonian
char strNumStages[256];
strncpy(strNumStages, g_inpset->integrator, 1);
strNumStages[1] = 0; // null terminate destination

if(g_inpset->iter_sampling == 0){

  printf("The number of production iterations chosen is 0, only burnin stage.\n");

}else{

  //skip the tuning part in HMCUpdate
  flag_dt_tuned = 1;
  
  //if-else for the optimal parameters setting procedure
  if(g_inpset->opt_param_setting == 1){

    printf("ATune: Optimal parameters setting procedure.\n");

    //optimal parameters setting procedure cannot applied without burn-in
    if(g_inpset->iter_burn_in == 0){
      fprintf(stderr, "hmc.c - Error using the optimal parameter setting approach. Burn-in is mandatory but has %i iterations.", g_inpset->iter_burn_in);
      exit(EXIT_FAILURE);
      }

    fitting_factor = fnCalculateFittingFactor("Gupta", 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin);

    printf("The fitting factor after the optimal parameters procedure is: %lf\n", fitting_factor);
    fprintf(artfp, "fitting_factor %lf\n", fitting_factor);
    
    //assigning h_lower
    h_lower = 2.077236689777599565332/3;

    //assigning h_upper
    h_upper = 1;

    //calculate upper and lower bounds in dimensional units
    hsl_low = 1 / fitting_factor * h_lower;
    hsl_upp = 1 / fitting_factor * h_upper;

    //transform into k-stages units (for the moment k = 2, 3)
    if(strcmp(strNumStages, "2") == 0){
      hsl_low = hsl_low * 2;
      hsl_upp = hsl_upp * 2;
      hsl_simulation = hsl_simulation * 2;
    }else if(strcmp(strNumStages, "3") == 0){
      hsl_low = hsl_low * 3;
      hsl_upp = hsl_upp * 3;
      hsl_simulation = hsl_simulation * 3;
    }

    printf("The lower HSL is: %lf\n", hsl_low);
    fprintf(artfp, "hsl_low %lf\n", hsl_low);

    printf("The upper HSL is: %lf\n", hsl_upp);
    fprintf(artfp, "hsl_upp %lf\n", hsl_upp);

    printf("The randomization interval for the stepsize is (%lf, %lf).\n", hsl_low, hsl_upp);

    //set randomization scheme for dt
    delta_t = hsl_upp;
    dt_random_scheme = 20;
    step_randomization = hsl_upp - hsl_low;

    //set randomization scheme for L
    L_random_scheme = 2;
    L_target = 6;
    L_lower = 2;

    //set randomization scheme for phi
    phi_random_scheme = 2;

    if(step_randomization < 0){
      fprintf(stderr, "hmc.c - Error in step_randomization calculation for the optimal parameters setting procedure, it should be positive and it is %lf\n", step_randomization);
      exit(EXIT_FAILURE);
    }

  }else if (g_inpset->opt_param_setting == 0){

    //No optimal parameters settings, everything is given in input

    //stepsize
    delta_t = g_inpset->stepsize;

    //assign randomization scheme for the stepsize and stepsize_delta
    dt_random_scheme = g_inpset->t_stepsize;
    step_randomization = g_inpset->stepsize_delta;

    //assign randomization schemes for L
    L_random_scheme = g_inpset->t_L;
    L_target = g_inpset->L;
    L_lower = g_inpset->L;

    //assign randomization schemes for phi
    phi_random_scheme = g_inpset->t_varphi;
    phi_target = g_inpset->varphi;

    //check for s-AIA integrator, the fitting factor must be provided to enable the h->b_opt map
    char integrator_compare[50];
    if ((strcmp(integrator_compare, "2sAIA-HMC") == 0) || (strcmp(integrator_compare, "3sAIA-HMC") == 0)){
      printf("You chose Adaptive Integration Approach.\n");

      //compute fitting factor from burn-in data
      if(g_inpset->iter_burn_in > 0){
        printf("Compute fitting factor from the burn-in data with the approach chosen in input.\n");
        fitting_factor = fnCalculateFittingFactor(g_inpset->fitting_factor_approach, 0, dt_tuned, 1, kNumParam, kNumParam, accept_burnin, accept_dHpos_burnin, num_proposal_burnin, num_proposal_dHpos_burnin, dH_avg_burnin, dH_pos_avg_burnin, dH2_avg_burnin);
        printf("The fitting factor after the burn-in is: %lf\n", fitting_factor);
      }else{//give it in input
        printf("You did not run burn-in. s-AIA will use the fitting factor value indicated in the inputfile.\n");
        fitting_factor = g_inpset->scaling_value;
        printf("The chosen fitting factor is: %lf\n", fitting_factor);
      }
      fprintf(artfp, "fitting_factor %lf\n", fitting_factor);
    }

  }else{
    fprintf(stderr, "hmc.c - Error in optimal parameters setting procedure, opt_param_setting does not fit the if-else, it has to be 0 or 1 and it is %i.\n", g_inpset->opt_param_setting);
    exit(EXIT_FAILURE);
  }//end if-else for the optimal parameters setting procedure

  //assign (G)HMC hyperparameters for sampling
    
  fnAssignStepsizes(stepsize, delta_t, dt_random_scheme, g_inpset->iter_burn_in, tot_iter - 1, step_randomization);
  fnAssignTrajectoryLengths(traj_length, &L_target, L_random_scheme, g_inpset->iter_burn_in, tot_iter - 1, L_lower);
  fnAssignIntegratorCoeffs(g_coeff, g_inpset->iter_burn_in, tot_iter - 1, 0);
  fnAssignPhi(varphi, phi_target, phi_random_scheme, g_inpset->iter_burn_in, tot_iter - 1, kNumParam);

  //assign MDMove for the production stage depending on the number of stages of the integrator chosen in input
  printf("Assigning MDMove for the production stage\n");
  fprintf(red_lfp, "Assigning MDMove for the production stage\n");

  //find number of stages as first chracter of g_inpset->integrator for assignment later
  char strNumStages[256];
  strncpy(strNumStages, g_inpset->integrator, 1);
  strNumStages[1] = 0; // null terminate destination

  if(strcmp(strNumStages, "1") == 0){ //1-stage integrators
    MDMove = VV;
  }else if(strcmp(strNumStages, "2") == 0){ //2-stage integrators
    MDMove = V2S;
  }else if(strcmp(strNumStages, "3") == 0){ //3-stage integrators
    MDMove = V3S;
  }else{
    fprintf(stderr, "HMC.c - No match found in cases for integrator for assigning MDMove. Exiting.");
    exit(EXIT_FAILURE);
  }

  ///////////////////////////
  //start production phase loop

  printf("--- Starting Production phase (collecting posterior samples...)---\n");
  printf("------------------------------------------------------------------\n");
  fprintf(red_lfp, "--- Starting Production phase (collecting posterior samples...)---\n");
  fprintf(red_lfp, "------------------------------------------------------------------\n");

  startT = clock();

  num_proposed  = 0;
  num_accepted  = 0;
  num_prop_mom = 0;
  num_acc_mom = 0;
  num_red_flips = 0;

  for (iter = g_inpset->iter_burn_in; iter < tot_iter; iter++) {
    if ((iter - g_inpset->iter_burn_in) % ((tot_iter - g_inpset->iter_burn_in) / 10) == 0 && (iter - g_inpset->iter_burn_in) > 0){
      printf("%d%% of production completed\n", 100 * (iter - g_inpset->iter_burn_in) / (tot_iter - g_inpset->iter_burn_in));
      printf("------------------------------\n");
    }

    HMCUpdate(iter);

    //write trajectory file
    if ((iter+1)%(g_inpset->thinning) == 0) {
      //main trajectory
      for (i = 0; i < kNumParam; i++) { fprintf(tfp, "%f ", state[curr]->pos[i]); }
      fprintf(tfp, "\n");

      //SIKR Incidence
      if (strcmp(g_inpset->model, "SIKR_Incidence") == 0) {
        printODE_SIR(odefp);
      }
      //SEMIKR synthetic
      if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0) {
        printODE_SEIR(odefp);
      }
      //SIR standard Incidence
      if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0) {
        printODE_SIR_stand(odefp);
      }
      //SIR standard Prevalence
      if (strcmp(g_inpset->model, "SIR_standard_Prevalence") == 0) {
        printODE_SIR_stand(odefp);
      }
    }
  }
    
  endT = clock();
  elapsedT = ((double)(endT-startT))/CLOCKS_PER_SEC;

  /////////////////////////
  //print results to files

  printf("--- Printing statistics.---\n");

  fprintf(red_lfp, "--- Printing statistics to art file.---\n");

  //number of proposals in the production stage: N
  int num_proposal_prod = tot_iter - g_inpset->iter_burn_in;
  fprintf(red_lfp, "Proposals production: %i\n", num_proposal_prod);
  fprintf(artfp, "N_prod %i\n", num_proposal_prod);

  //number of accepted proposals during the production - N_acc
  int accept_prod = 0;
  for(ii=g_inpset->iter_burn_in; ii < tot_iter; ii++){
    accept_prod += accept_arr[ii];
  }
  fprintf(red_lfp, "Accepted proposals production: %i\n", accept_prod);
  fprintf(artfp, "N_acc_prod %i\n", accept_prod);

  //number of accepted proposals with negative energy error during the production - N-
  int accept_dHneg_prod = 0;
  for(ii=g_inpset->iter_burn_in; ii < tot_iter; ii++){
    accept_dHneg_prod += num_proposals_dHneg_arr[ii];
  }
  fprintf(red_lfp, "Accepted proposals with dH<0 production: %i\n", accept_dHneg_prod);
  fprintf(artfp, "N_dHneg_prod %i\n", accept_dHneg_prod);

  //number of proposals with positive energy error during the production - N+
  int num_proposal_dHpos_prod = num_proposal_prod - accept_dHneg_prod;
  fprintf(red_lfp, "Proposals with dH>0 production: %i\n", num_proposal_dHpos_prod);
  fprintf(artfp, "N_dHpos_prod %i\n", num_proposal_dHpos_prod);

  //number of accepted proposals with positive energy error during the production - N+_acc
  int accept_dHpos_prod = accept_prod-accept_dHneg_prod;
  fprintf(red_lfp, "Accepted proposals with dH>0 production: %i\n", accept_dHpos_prod);
  fprintf(artfp, "N_acc_dHpos_prod %i\n", accept_dHpos_prod);

  //average acceptance rate over the production stage - AR = N_acc / N
  double AR = (double)accept_prod/num_proposal_prod;
  printf("Acceptance rate production: %lf\n", AR);
  fprintf(red_lfp, "Acceptance rate production: %lf\n", AR);
  fprintf(artfp, "AR_prod %lf\n", AR);

  //average acceptance rate
  AR = (double)num_accepted/(double)num_proposed;
  printf("Acceptance rate: %lf\n", AR);
  fprintf(red_lfp, "Acceptance rate: %lf\n", AR);
  fprintf(artfp, "AR %lf\n", AR);

  //average filtered acceptance rate over the production stage (only cases accepted using the energy difference over proposals with positive energy error) - filtAR = N+acc / N+
  double AR_filt_prod = (double)accept_dHpos_prod/num_proposal_dHpos_prod;
  fprintf(red_lfp, "Filtered AR production: %lf\n", AR_filt_prod);
  fprintf(artfp, "AR_filt_prod %lf\n", AR_filt_prod);

  //average new filtered acceptance rate over the burn-in stage (only cases accepted using the energy difference over the total number of proposals) - filtAR_new = N+_acc / N
  double AR_filtNew_prod = (double)accept_dHpos_prod/num_proposal_prod;
  fprintf(red_lfp, "New Filtered AR production: %lf\n", AR_filtNew_prod);
  fprintf(artfp, "AR_filtNew_prod %lf\n", AR_filtNew_prod);

  //average energy error over the production stage - dH = sum_1^N dH_i /N
  double dH_avg_prod = 0.;
  for(ii=g_inpset->iter_burn_in; ii < tot_iter; ii++){
    dH_avg_prod += dH_array[ii];
  }
  dH_avg_prod = dH_avg_prod/num_proposal_prod;
  fprintf(red_lfp, "Average Energy Error production: %lf\n", dH_avg_prod);
  fprintf(artfp, "dH_avg_prod %lf\n", dH_avg_prod);

  //average positive energy error over the production stage (averaged over the total number of proposals) - dH+ = sum_1^N dH+_i / N
  double dH_pos_avg_prod = 0.;
  for(ii=g_inpset->iter_burn_in; ii < tot_iter; ii++){
    if(dH_array[ii] > 0){
      dH_pos_avg_prod += dH_array[ii];
    }
  }
  dH_pos_avg_prod = dH_pos_avg_prod/num_proposal_prod;
  fprintf(red_lfp, "Average positive Energy Error production: %lf\n", dH_pos_avg_prod);
  fprintf(artfp, "dH_pos_avg_prod %lf\n", dH_pos_avg_prod);

  //average positive energy error over the production stage (averaged over the proposals with positive energy error) - H+_av = sum_1^N+ dH+_i / N+
  double H_pos_av_prod = dH_pos_avg_prod/num_proposal_dHpos_prod;
  fprintf(red_lfp, "H+_av production: %lf\n", H_pos_av_prod);
  fprintf(artfp, "H_pos_av_prod %lf\n", H_pos_av_prod);

  //average positive energy error over the production stage (averaged over the accepted proposals with positive energy error) - H+_acc_av = sum_1^N+_acc dH+_i / N+_acc
  double H_pos_acc_av_prod = dH_pos_avg_prod/accept_dHpos_prod;
  fprintf(red_lfp, "H+_acc_av production: %lf\n", H_pos_acc_av_prod);
  fprintf(artfp, "H_pos_acc_av_prod %lf\n", H_pos_acc_av_prod);

  //average squared energy error over the production stage - dH^2 = sum_i^N dH^2_i / N
  double dH2_avg_prod = 0.;
  for(ii=g_inpset->iter_burn_in; ii < tot_iter; ii++){
    dH2_avg_prod += dH_array[ii]*dH_array[ii];
  }
  dH2_avg_prod = dH2_avg_prod/num_proposal_prod;
  fprintf(red_lfp, "Average Squared Energy Error production: %lf\n", dH2_avg_prod);
  fprintf(artfp, "dH2_avg_prod %lf\n", dH2_avg_prod);

  //standard dev of energy error over the production stage
  double dH_sd_prod = sqrt(dH2_avg_prod - dH_avg_prod*dH_avg_prod);
  fprintf(red_lfp, "StdDevEnergyError production: %lf\n", dH_sd_prod);
  fprintf(artfp, "dH_sd_prod %lf\n", dH_sd_prod);

  //time elapsed total
  printf("Elapsed time: %lf sec\n", elapsedT);
  fprintf(red_lfp, "Elapsed time: %lf\n", elapsedT);
  fprintf(artfp, "TimeTotal %lf\n", elapsedT);

  //mean of integrator coefficient b1
  double mean_g_coeff_i_0 = 0;
  for (ii = 0; ii < tot_iter; ii++){ mean_g_coeff_i_0 += g_coeff[ii][0]; }
  mean_g_coeff_i_0 = mean_g_coeff_i_0 / tot_iter;
  printf("Average b: %lf\n", mean_g_coeff_i_0);
  fprintf(red_lfp, "Average b: %lf\n", mean_g_coeff_i_0);
  fprintf(artfp, "mean_b %lf\n", mean_g_coeff_i_0);

  //standard deviation of integrator coefficient b1
  double sd_g_coeff_i_0 = 0;
  for (ii = 0; ii < tot_iter; ii++){ sd_g_coeff_i_0 += pow(g_coeff[ii][0] - mean_g_coeff_i_0, 2.0); }
  sd_g_coeff_i_0 = sqrt(sd_g_coeff_i_0 / tot_iter);
  fprintf(red_lfp, "standard deviation b: %lf\n", sd_g_coeff_i_0);
  fprintf(artfp, "sd_b %lf\n", sd_g_coeff_i_0);

  //mean of integrator coefficient a1
  double mean_g_coeff_i_1 = 0;
  for (ii = g_inpset->iter_burn_in; ii < tot_iter; ii++){ mean_g_coeff_i_1 += g_coeff[ii][1]; }
  mean_g_coeff_i_1 = mean_g_coeff_i_1 / (tot_iter - g_inpset->iter_burn_in);
  printf("Average a: %lf\n", mean_g_coeff_i_1);
  fprintf(red_lfp, "Average a: %lf\n", mean_g_coeff_i_1);
  fprintf(artfp, "mean_a %lf\n", mean_g_coeff_i_1);

  //standard deviation of integrator coefficient a1
  double sd_g_coeff_i_1 = 0;
  for (ii = g_inpset->iter_burn_in; ii < tot_iter; ii++){ sd_g_coeff_i_1 += pow(g_coeff[ii][1] - mean_g_coeff_i_1, 2.0); }
  sd_g_coeff_i_1 = sqrt(sd_g_coeff_i_1 / (tot_iter - g_inpset->iter_burn_in));
  fprintf(red_lfp, "standard deviation a: %lf\n", sd_g_coeff_i_1);
  fprintf(artfp, "sd_a %lf\n", sd_g_coeff_i_1);

  //mean of angle phi
  double mean_phi = 0;
  for (ii = 0; ii < tot_iter; ii++){ mean_phi += varphi[ii]; }
  mean_phi = mean_phi / tot_iter;
  printf("Average angle phi: %lf\n", mean_phi);
  fprintf(red_lfp, "Average angle phi: %lf\n", mean_phi);
  fprintf(artfp, "mean_phi %lf\n", mean_phi);

  //standard deviation of angle phi
  double sd_phi = 0;
  for (ii = 0; ii < tot_iter; ii++){ sd_phi += pow(varphi[ii] - mean_phi, 2.0); }
  sd_phi = sqrt(sd_phi / tot_iter);
  fprintf(red_lfp, "standard deviation phi: %lf\n", sd_phi);
  fprintf(artfp, "sd_phi %lf\n", sd_phi);

  //mean of time step
  double mean_dt = 0;
  for (ii = g_inpset->iter_burn_in; ii < tot_iter; ii++){ mean_dt += stepsize[ii]; }
  mean_dt = mean_dt / (tot_iter - g_inpset->iter_burn_in);
  printf("Average time step: %lf\n", mean_dt);
  fprintf(red_lfp, "Average time step: %lf\n", mean_dt);
  fprintf(artfp, "mean_timestep %lf\n", mean_dt);

  //mean of L
  double mean_L = 0;
  for (ii = g_inpset->iter_burn_in; ii < tot_iter; ii++){ mean_L += traj_length[ii]; }
  mean_L = mean_L / (tot_iter - g_inpset->iter_burn_in);
  printf("Average L: %lf\n", mean_L);
  fprintf(red_lfp, "Average L: %lf\n", mean_L);
  fprintf(artfp, "mean_L %lf\n", mean_L);

  }
  //end of main sampling loop

  ///////////////////////////////////////////////
  //close opened file pointers
  fclose(tfp);
  fclose(lpdfp);
  fclose(hfp);
  fclose(tsfp);
  fclose(artfp);
  fclose(artburnin_fp);
  fclose(scfp);
  if(g_inpset->iter_burn_in > 0){
    fclose(final_point_tune_fp);
    fclose(final_point_burnin_fp);
  }
  fclose(odefp);

  ///////////////////////////////////////////////
  //free allocated memory
  for (i = back2; i <= forw2; i++) {
    SFree(&state[i]);
    SFree(&state_old[i]);
  }

  MatrixFree(g_obs);

  free(stepsize);
  free(traj_length);
  free(varphi);

  gsl_rng_free(g_rng);

  free(g_inpset);

  return 0;
}
