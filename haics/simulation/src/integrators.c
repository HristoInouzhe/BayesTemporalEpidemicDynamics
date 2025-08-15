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
#include <string.h>
//#include "mkl.h"
#include <gsl/gsl_cblas.h>
//#include <cblas.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "integrators.h"
#include "adapt_integrator.h"
#include "Globals.h"
#include "utils.h"
#include "hmc.h"

/*************************************************************
 *     MovePos()
 *     drift
 *     position = position + stepsize*mom
 *     pos: position
 *     mom: momenta
 ***********************************************************/
static void MovePos(int dim, TypeState *st, double coeff, double ss) {
  double sss = coeff*ss;
  
  cblas_daxpy(dim, sss, st->mom, 1, st->pos, 1);
  
}

/*************************************************************
 *   MoveMom()
 *   kick
 *   momenta = momenta + stepsize*grad
 *   gradient_function (force): gradient_model
 *************************************************************/

static void MoveMom(int dim, TypeState *st, double coeff, double ss) {
  double sss = coeff*ss;

  cblas_daxpy(dim, sss, st->grad, 1, st->mom, 1);

}

// florian-wayne:
/************************************************************
 *   MetropolisMove()
 *   To be used for Metropolis-Hastings MCMC
 ************************************************************/
void MetropolisMove(void (*GradientFunction)(TypeState *st), int dim, TypeState **st, double ss, double *coeff, int L) {

    /* The stepsize ss here will not be used in the same way as with the other integrators.
     * For Metropolis MCMC, the stepsize array stepsize[] as initialized in read_input has all elements set to 1.
     * The actual stepsize that will be used here will be taken from g_inpset->stepsize which is specified in
     * the input file. */

    // Note that g_coeff in read_input has elements of value 0 for Metropolis MCMC

    cblas_daxpy(dim,g_inpset->stepsize,st[curr]->mom,1,st[curr]->pos,1);
}

/************************************************************
 *   VV()
 *   Verlet integrator (velocity) (used for HMC)
 *   Integration of a trajectory of length L coefficients for kick and drift moves as:
 *   coeff = { 0.5,  1, 0, 0,  0,  0,  0,  0}

 ************************************************************/
void VV(void (*GradientFunction)(TypeState *st), int dim, TypeState **st, double ss, double *coeff, int L) {

  int istep;

  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[0], ss); //0.5

  for (istep = 0; istep < L-1; istep++) {
    MovePos(dim, st[curr], coeff[1], ss); //1
    GradientFunction(st[curr]);
    MoveMom(dim, st[curr], 2.*coeff[0], ss); //1
  }

  MovePos(dim, st[curr], coeff[1], ss); //1
  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[0], ss); //0.5
}

/************************************************************
 *   V2S()
 *   Two stages integrator (velocity) for HMC
 *   Integration of a trajectory of length L coefficients for kick and drift moves as:
 *   coeff = { b1, a1, b2, a1,  0,  0,  0,  0}

 ************************************************************/
void V2S(void (*GradientFunction)(TypeState *st), int dim, TypeState **st, double ss, double *coeff, int L) {

  int istep;

  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[0], ss); //b1

  for (istep = 0; istep < L-1; istep++) {

    MovePos(dim, st[curr], coeff[1], ss); //a1
    GradientFunction(st[curr]);
    MoveMom(dim, st[curr], coeff[2], ss); //b2
    MovePos(dim, st[curr], coeff[3], ss); //a1
    GradientFunction(st[curr]);
    MoveMom(dim, st[curr], 2.*coeff[0], ss);//2*b1
  }

  MovePos(dim, st[curr], coeff[3], ss); //a1
  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[2], ss); //b2
  MovePos(dim, st[curr], coeff[1], ss); //a1
  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[0], ss); //b1
}

/************************************************************
 *   V3S()
 *   Three stages integrator (velocity) for HMC
 *   Integration of a trajectory of length L with coefficients for kick and drift moves as:
 *   coeff = { b1, a1, b2, a2, b2, a1,  0,  0}
 ************************************************************/
void V3S(void (*GradientFunction)(TypeState *st), int dim, TypeState **st, double ss, double *coeff, int L) {

  int istep;

  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[0], ss); //b1

  for (istep = 0; istep < L-1; istep++) {
    MovePos(dim, st[curr], coeff[1], ss); //a1
    GradientFunction(st[curr]);
    MoveMom(dim, st[curr], coeff[2], ss); //b2
    MovePos(dim, st[curr], coeff[3], ss); //a2
    GradientFunction(st[curr]);
    MoveMom(dim, st[curr], coeff[4], ss); //b2
    MovePos(dim, st[curr], coeff[5], ss); //a1
    GradientFunction(st[curr]);
    MoveMom(dim, st[curr], 2*coeff[0], ss);//2*b1
  }

  MovePos(dim, st[curr], coeff[5], ss); //a1
  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[4], ss); //b2
  MovePos(dim, st[curr], coeff[3], ss); //a2
  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[2], ss); //b2
  MovePos(dim, st[curr], coeff[1], ss); //a1
  GradientFunction(st[curr]);
  MoveMom(dim, st[curr], coeff[0], ss); //b1
}

/************************************************************
 *  fnAssignIntegratorCoeffs()
 *  assigns coefficients into mat_g_coeff depending on chosen integrator
 *  Arguments:
 ** mat_g_coeff -  structure where coefficients should be written into
 ** start and stop of iteration number to assign (in array mat_g_coeff)
 ** flag_burn_in -  0: assignment for main sampling phase, 1: for burn-in phase
 ***********************************************************/

void fnAssignIntegratorCoeffs(double **mat_g_coeff, int start_iter_assign, int end_iter_assign, int flag_burn_in){

  //printing to logfiles
  fprintf(red_lfp, "-- fnAssignIntegratorCoeffs: from iteration %i to %i, flag_burn_in is %i\n", start_iter_assign+1, end_iter_assign+1, flag_burn_in);

  //declaring
  FILE *bfp;
  FILE *b3sfp;
  char filename[80];
  char filename_b[80];
  double b_aia[hbar_to_b_num_disc_steps];
  double b_aia_3s[hbar_to_b_num_disc_steps_3s];
  double num;
  int i;
  char integrator_compare[50];

  //assign integrator string 
  strcpy(integrator_compare, g_inpset->integrator);

  //find number of stages as first chracter of g_inpset->integrator
  char strNumStages_compare[256];
  strncpy(strNumStages_compare, integrator_compare, 1);
  strNumStages_compare[1] = 0; // null terminate destination

  //assign data string
  char data_compare[50];
  strcpy(data_compare, g_inpset->data);

  //assign integrator's coefficients
  if (strcmp(integrator_compare, "1sVV") == 0) {//1-stage integrator Velocity Verlet
    double b1 = 0.5;

    double a1 = 1.0;

    for (i = start_iter_assign; i <= end_iter_assign; i++) {
      mat_g_coeff[i][0] = b1;
      mat_g_coeff[i][1] = a1;
    }

  } else if (strcmp(integrator_compare, "2sVV") == 0) {//2-stage integrator Velocity Verlet
    double b1 = 0.25;
    double b2 = 1 - 2 * b1;

    double a1 = 0.5;

    for (i = start_iter_assign; i <= end_iter_assign; i++) {
      mat_g_coeff[i][0] = b1;
      mat_g_coeff[i][1] = a1;
      mat_g_coeff[i][2] = b2;
      mat_g_coeff[i][3] = a1;
    }

  } else if (strcmp(integrator_compare, "2sBCSS-HMC") == 0) {//2-stage integrator (min rho) for HMC - see Blanes, Casas, Sanz-Serna (2014)
    double b1 = 0.211781;
    double b2 = 1 - 2 * b1;

    double a1 = 0.5;

    for (i = start_iter_assign; i <= end_iter_assign; i++) {
      mat_g_coeff[i][0] = b1;
      mat_g_coeff[i][1] = a1;
      mat_g_coeff[i][2] = b2;
      mat_g_coeff[i][3] = a1;
    }

  } else if (strcmp(integrator_compare, "2sMinError-HMC") == 0) {//2-stage minimum error integrator for HMC - see McLachlan (1995)
    double b1 = 0.1932;
    double b2 = 1 - 2 * b1;

    double a1 = 0.5;

    for (i = start_iter_assign; i <= end_iter_assign; i++) {
      mat_g_coeff[i][0] = b1;
      mat_g_coeff[i][1] = a1;
      mat_g_coeff[i][2] = b2;
      mat_g_coeff[i][3] = a1;
    }


  } else if (strcmp(integrator_compare, "2sAIA-HMC") == 0) {//2-stage s-AIA - see Nagar et al. (2023)
    //this if clause catches integrators which are a) 2 stage and b) use AIA

    //assign b_aia input file
    sprintf(filename, "./input/b_2sAIA-HMC.txt");

    //check for errors
    bfp = fopen(filename, "r");
    if (bfp == NULL) {
      fprintf(stderr, "integrators.c - Cannot open file %s \n", filename);
      exit(1);
    }

    //read coefficients from file into variable
    for (i = 0; i < hbar_to_b_num_disc_steps; i++) {
      if (fscanf(bfp, "%lf", &num) != EOF) {
        b_aia[i] = num;
      }
    }

    //close file
    fclose(bfp);

    //b1 and b2 are calculated within the for loop
    double a1 = 0.5;

    for (i = start_iter_assign; i <= end_iter_assign; i++){
      double b1 = FindIntegrator2stage(b_aia, stepsize[i]);
      double b2 = 1 - 2 * b1;

      mat_g_coeff[i][0] = b1;
      mat_g_coeff[i][1] = a1;
      mat_g_coeff[i][2] = b2;
      mat_g_coeff[i][3] = a1;

    } //end of for loop
  } else if (strcmp(integrator_compare, "3sBCSS-HMC") == 0) {//3-stage integrator (min rho) for HMC - see Blanes, Casas, Sanz-Serna (2014)
    double b1 = 0.11888010966548;
    double b2 = 0.5 - b1;

    double a1 = 0.29619504261126;
    double a2 = 1 - 2 * a1;

    for (i = start_iter_assign; i <= end_iter_assign; i++) {
      mat_g_coeff[i][0] = b1;
      mat_g_coeff[i][1] = a1;
      mat_g_coeff[i][2] = b2;
      mat_g_coeff[i][3] = a2;
      mat_g_coeff[i][4] = b2;
      mat_g_coeff[i][5] = a1;
    }

  } else if (strcmp(integrator_compare, "3sAIA-HMC") == 0) {//3-stage s-AIA - see Nagar et al. (2023)

    //assign b_aia input file
    sprintf(filename_b, "./input/b_3sAIA-HMC.txt");

    //open and check for errors
    b3sfp = fopen(filename_b, "r");
    if (b3sfp == NULL) {
      fprintf(stderr, "integrators.c - Cannot open file %s \n", filename_b);
      exit(EXIT_FAILURE);
    }

    //read coefficients from file into variable
    for (i = 0; i < hbar_to_b_num_disc_steps_3s; i++) {
      if (fscanf(b3sfp, "%lf", &num) != EOF) {
        b_aia_3s[i] = num;
      }
    }

    //close file
    fclose(b3sfp);
    
    for (i = start_iter_assign; i <= end_iter_assign; i++) {
      double b1 = FindIntegrator3stage(b_aia_3s, stepsize[i]);
      double b2 = 0.5 - b1;

      double a1 = (1-2*b1)/(4*(1-3*b1)); //optimal stability hyperbola for (a, b), from Eq. (29) in Radivojevic et al. (2018)
      double a2 = 1 - 2 * a1;
        
      mat_g_coeff[i][0] = b1;
      mat_g_coeff[i][1] = a1;
      mat_g_coeff[i][2] = b2;
      mat_g_coeff[i][3] = a2;
      mat_g_coeff[i][4] = b2;
      mat_g_coeff[i][5] = a1;
    }
  } else if (strcmp(integrator_compare, "3sMinError-HMC") == 0) {//3-stage minimum error integrator for HMC - see Predescu et al. (2012)
    double b1 = 0.108991425403425322; // a in eq. (26)
    double b2 = 0.5 - b1;

    double a1 = 0.290485609075128726; // b in eq. (26)
    double a2 = 1 - 2 * a1;

    for (i = start_iter_assign; i <= end_iter_assign; i++) {
      mat_g_coeff[i][0] = b1;
      mat_g_coeff[i][1] = a1;
      mat_g_coeff[i][2] = b2;
      mat_g_coeff[i][3] = a2;
      mat_g_coeff[i][4] = b2;
      mat_g_coeff[i][5] = a1;
    }
  
  } else {
    fprintf(stderr, "integrators.c - Integrator option not found. Choose one of 1sVV, 2sVV, 2sBCSS-HMC, 3sBCSS-HMC,, 2sAIA-HMC, 3sAIA-HMC, 2sMinError-HMC, 3sMinError-HMC.\n");
    exit(1);
  }
}//end of function fnAssignIntegratorCoeffs


/************************************************************
 *  fnAssignStepsizes()
 *  assigns stepsizes into array stepsize depending on t_stepsize
 *  Arguments:
 ** stepsize - structure where stepsizes should be written into
 ** stepsize_value - the value that should be taken as the baseline for assigning stepsizes
 ** type_random_stepsize - the type of randomization that should be used
 ** start and stop of iteration number to assign (in array stepsize)
 ***********************************************************/

void fnAssignStepsizes(double *stepsize, double stepsize_value, int type_random_stepsize, int start_iter_assign, int end_iter_assign, double random_step){

  //printing to logfiles
  fprintf(red_lfp, "-- fnAssignStepsizes: from iteration %i to %i, stepsize_value is %lf, t_stepsize is %i\n", start_iter_assign+1, end_iter_assign+1, stepsize_value, type_random_stepsize);

  //declaring
  int ii;
  double r;
  double std;

  //find number of stages as first chracter of g_inpset->integrator (for safety check later)
  char strNumStages[256];
  strncpy(strNumStages, g_inpset->integrator, 1);
  strNumStages[1] = 0; //null terminate destination

  switch(type_random_stepsize){

  case 0: //constant

    for(ii=start_iter_assign; ii<=end_iter_assign; ii++){
      stepsize[ii] = stepsize_value;
    }

    break;
  case 1: //random from N(stepsize, 0.0025*stepsize^2)

    for(ii=start_iter_assign; ii<=end_iter_assign; ii++){
      std = 0.05 * stepsize_value;
      r = gsl_ran_gaussian_ziggurat(g_rng, std);
      stepsize[ii] = (0.8 + r * 0.4) * stepsize_value;
    }

    break;
  case 20: //random from U(stepsize_value - stepsize_delta, stepsize_value)

    for(ii=start_iter_assign; ii<=end_iter_assign; ii++){
      r = gsl_rng_uniform(g_rng);
      stepsize[ii] = stepsize_value - random_step*r;
    }

    break;
  case 21: //random from U(stepsize_value, stepsize_value + stepsize_delta)

    for(ii=start_iter_assign; ii<=end_iter_assign; ii++){
      r = gsl_rng_uniform(g_rng);
      stepsize[ii] = stepsize_value + g_inpset->stepsize_delta*r;
    }

    break;        
  case 3: //random from U(stepsize_value - stepsize_delta, stepsize_value + stepsize_delta)    

    for(ii=start_iter_assign; ii<=end_iter_assign; ii++){
      r = gsl_rng_uniform(g_rng);
      stepsize[ii] = stepsize_value - g_inpset->stepsize_delta + 2*r*g_inpset->stepsize_delta;
    }

    break;
  default:
    fprintf(stderr, "integrators.c - type_random_stepsize has to be one of 0, 1, 20, 21, 3.\n");
    exit(EXIT_FAILURE);
  break;
  } //end switch type_random_stepsize

} //end of function fnAssignStepsizes


/************************************************************
 *  fnAssignTrajectoryLengths()
 *  assigns trajectory lengths into array traj_length depending on t_L
 *  Arguments:
 ** traj_length - structure where trajectory lengths should be written into
 ** L_value - the value that should be taken as the baseline for assigning trajectory lengths
 ** type_random_L - the type of randomization that should be used
 ** start and stop of iteration number to assign (in array stepsize)
 ***********************************************************/

void fnAssignTrajectoryLengths(int *traj_length, int *L_value, int type_random_L, int start_iter_assign, int end_iter_assign, int L_lower){

  //printing to logfiles
  fprintf(red_lfp, "-- fnAssignTrajectoryLengths: from iteration %i to %i, L_value is %i, t_L is %i\n", start_iter_assign+1, end_iter_assign+1, *L_value, type_random_L);

  //declaring
  int ii;
  double r;

  switch(type_random_L){
  case 0: //constant

    for(ii=start_iter_assign; ii<=end_iter_assign; ii++){
      traj_length[ii] = *L_value;
    }

    break;
  case 1: //random from U{1,L}

    for(ii=start_iter_assign; ii<=end_iter_assign; ii++){
      r = gsl_rng_uniform_int(g_rng, (*L_value));
      traj_length[ii] = 1 + r;
    }

    break;
  case 2: //random from U{L_lower,L_value}

    for(ii=start_iter_assign; ii<=end_iter_assign; ii++){
      r = gsl_rng_uniform(g_rng);
      traj_length[ii] = round(L_lower + r * (*L_value - L_lower));
    }

    printf("Trajectory length for the production L_prod in U{%d, %d}.\n", L_lower, *L_value);

    fprintf(artfp, "L_lower %d\n", L_lower);
    fprintf(artfp, "L_upper %d\n", *L_value);

    printf("The theoretical average of L is %f\n", (L_lower + *L_value)/2.);
    fprintf(artfp, "avL %f\n", (L_lower + *L_value)/2.);

    break;
  default:
    fprintf(stderr, "integrators.c - type_random_L has to be one of 0, 1.\n");
    exit(EXIT_FAILURE);
  break;  
  } //end of switch type_random_L

} //end of function fnAssignTrajectoryLengths

//////////////////////////////////////////////////////
// fnAssignPhi()
// assigns phi into the array for phi depending on t_varphi
// Arguments:
// varphi - structure where phi should be written into
// varphi_value - the value that should be taken as the baseline for assigning phi
// type_random_varphi - the type of randomization that should be used
// start and stop of iteration number to assign (in array varphi)

void fnAssignPhi(double *varphi, double varphi_value, int type_random_varphi, int start_iter_assign, int end_iter_assign, int kNumParam){

  //printing to logfiles
  fprintf(red_lfp, "-- fnAssignPhi: from iteration %i to %i, varphi_value is %lf, t_varphi is %i\n", start_iter_assign+1, end_iter_assign+1, varphi_value, type_random_varphi);

  //declaring
  int ii;
  double r;
  double h_phi, h_low, lambda_phi, lambda_low;

  char strNumStages[256];
  strncpy(strNumStages, g_inpset->integrator, 1);
  strNumStages[1] = 0; //null terminate destination

  if ((strcmp(g_inpset->method, "HMC") == 0) || (strcmp(g_inpset->method, "TUNE") == 0)) {
    varphi_value = 1;
    type_random_varphi = 0;
  }

  switch(type_random_varphi){
  case 0: //constant

    for (ii=start_iter_assign; ii<=end_iter_assign; ii++){
      varphi[ii] = varphi_value;
    }

    break;
  case 1: //random from U(0,Phi)

    for (ii=start_iter_assign; ii<=end_iter_assign; ii++){
      r = gsl_rng_uniform(g_rng);
      varphi[ii] = varphi_value * r;
    }

    break;
  case 2: //optimal phi for GHMC, randomize in (phi(h_CoLSI), phi(h_lower)), see Akhmatskaya et al. (2025)
    
    //calculate h and lambda depending on a number of stages
    if(strcmp(strNumStages, "1") == 0){
      h_phi = 1;
      h_low = 0.5193;
      double b_phi = 0.5;
      double b_low = 0.5;
      lambda_phi = (6 * b_phi - 1) / 24;
      lambda_low = (6 * b_low - 1) / 24;
    }else if(strcmp(strNumStages, "2") == 0){
      h_phi = 2;
      h_low = 1.0386;
      double b_phi = 0.211781;
      double b_low = 0.1958;
      lambda_phi = (6 * b_phi - 1) / 24;
      lambda_low = (6 * b_low - 1) / 24;
    }else if(strcmp(strNumStages, "3") == 0){
      h_phi = 3;
      h_low = 2.0772;
      double b_phi = 0.11888010966548;
      double a_phi = (1-2*b_phi)/(4*(1-3*b_phi));
      double b_low = 0.113252;
      double a_low = (1-2*b_low)/(4*(1-3*b_low));
      lambda_phi = (1 - 6 * a_phi * (1 - a_phi) * (1 - 2 * b_phi)) / 12;
      lambda_low = (1 - 6 * a_low * (1 - a_low) * (1 - 2 * b_low)) / 12;
    }

    //calculate boundaries for phi
    double varphi_1 = fmin(1, - log(0.999) * (1 + 2 * h_phi * h_phi * lambda_phi) / (2 * kNumParam * h_phi * h_phi * h_phi * h_phi * lambda_phi * lambda_phi));
    varphi_1 = fmax(0.001, varphi_1);

    double varphi_2 = fmin(1, - log(0.999) * (1 + 2 * h_low * h_low * lambda_low) / (2 * kNumParam * h_low * h_low * h_low * h_low * lambda_low * lambda_low));

    //assign phi at each iteration
    for (ii=start_iter_assign; ii<=end_iter_assign; ii++){
      r = gsl_rng_uniform_pos(g_rng);
      varphi[ii] = varphi_1 + r * (varphi_2 - varphi_1);
    }//for loop end

    printf("The randomization scheme for Phi is (%lf, %lf).\n", varphi_1, varphi_2);
    fprintf(artfp, "phi_low %lf\n", varphi_1);
    fprintf(artfp, "phi_upp %lf\n", varphi_2);

    break;
  default:
    fprintf(stderr, "integrators.c - type_random_varphi has to be one of 0, 1.\n");
    exit(EXIT_FAILURE);
  break;
  } //end switch type_random_varphi

} //end of function fnAssignPhi
