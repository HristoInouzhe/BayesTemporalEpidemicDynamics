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
//
// SUNDIALS functions for the SEMIKR models
/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is an example for the simple SIR model, with the coding
 * needed for its solution by CVODES for Forward Sensitivity
 * Analysis. The problem is defined as (S, I, R) = (y1, y2, y3) and
 * (beta, gamma, I(t0)) = (p1, p2, p3)
 * with the following three ordinary diferential equations:
 *    dy1/dt = -p1*y1*y2/NPOP
 *    dy2/dt =  p1*y1*y2/NPOP - p2*y2
 *    dy3/dt =  p2*y2
 * on the interval from t = 0.0 to t = 20, with initial
 * conditions y1 = NPOP - p3, y2 = p3 and y3 = 0. NPOP = y1 + y2 + y3 = S + I + R is a constant.
 * The initial parameters are: p1 = 2, p2 = 1, and p3 = 1. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the DENSE linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 * Output is printed in days from t = 1 to t = 20.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters p1, p2 and p3.
 * The sensitivity right hand side is given analytically through the
 * user routine fS (of type SensRhs1Fn).
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on SUNTRUE or SUNFALSE,
 * respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % cvsRoberts_FSA_dns -nosensi
 * If sensitivities are to be computed:
 *    % cvsRoberts_FSA_dns -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>


#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */
//#include <legendre/legendre_polynomial.h>
/* Accessor macros */

#include "Globals.h"
#include "utils.h"
#include "seir_functions.h"
#include "read_input.h"
#include "sir_seir_common_functions.h"


#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

/* Problem Constants */

#define Y3    RCONST(0.0)
#define RTOL  RCONST(1e-10)  /* scalar relative tolerance */
#define ATOL1 RCONST(1e-10)  /* vector absolute tolerance components */
#define MAXSTEPS RCONST(1e+10)
#define T0    RCONST(0.0)   /* initial time */
#define T1    RCONST(1.0)   /* first output time */
#define TAD RCONST(1.0)  /* output time factor */

#define X0 0.0
#define gamma0 1e-10
#define phiInv0 1e-10
#define alpha0 1e-10

//#define gammaL 0.03333333333333 
//#define gammaU 1.0

#define ZERO  RCONST(0.0)

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

void SEIR_CVODES(double results_y[][num_out_cvodes], double results_s[][num_out_cvodes][kNumParam - 2], double state_params[kNumParam - 2],
                 double degree, double m, double nf, double *knots, double *basis, char *sensitivity_parameters[])
//void SIR_CVODES(double **results_y, double ***results_s, double *state_params,
//                 double degree, double m, double nf, double *knots, double *basis, char *sensitivity_parameters[])
{
  //choose the number of components with respect to the model
  int NICOMP = g_inpset->num_comp;
  int NECOMP = g_inpset->num_comp_E;
  int NEQ = 3 + NECOMP + NICOMP;


  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  UserDataSEIR data;
  realtype t, tout;
  N_Vector y;
  int iout, retval, j, jj;
  //static double results_y[NEQ][num_out_cvodes];
  realtype *y_data;
  //double *y_data;

  //realtype pbar[3 + g_inpset->num_basis_spline];
  int is;
  N_Vector *yS;
  booleantype sensi, err_con;
  int sensi_meth;
  //double results_s[NEQ][num_out_cvodes][3 + g_inpset->num_basis_spline];
  realtype *s_data;

  cvode_mem = NULL;
  data      = NULL;
  y         = NULL;
  yS        = NULL;
  A         = NULL;
  LS        = NULL;

  /* Process arguments */
  ProcessArgs(sensitivity_parameters, &sensi, &sensi_meth, &err_con);

  /* User data structure */
  data = (UserDataSEIR) malloc(sizeof *data);
  if (check_retval((void *)data, "malloc", 2)) {
    printf("Error: UserDataSEIR\n");
    exit(EXIT_FAILURE);
  }
  for(j = 0; j < kNumParam - 2; j++){
    data->p[j] = state_params[j];
  }
  data->spline_degree = degree;
  data->number_extended_knots = m;
  data->number_basis_elements = nf;
  data->knots = knots;
  data->basis = basis;

  /* Initial conditions */
  y = N_VNew_Serial(NEQ);
  if (check_retval((void *)y, "N_VNew_Serial", 0)){
    printf("Error: N_VNew_Serial\n");
    exit(EXIT_FAILURE);
  }

  Ith(y,1) = data->p[2];
  Ith(y,2) = data->p[3];
  for(j = 1; j <= (NECOMP - 1); j++){
    Ith(y, 2 + j) = RCONST(0.0);
  }
  Ith(y, NECOMP + 2) = data->p[4];
  for(j = 1; j <= (NICOMP - 1); j++){
    Ith(y, NECOMP + 2 + j) = RCONST(0.0);
  }
  Ith(y, NEQ - 1) = num_pop - data->p[2] - data->p[3] - data->p[4];
  Ith(y, NEQ) = num_pop - data->p[2];


  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) {
      printf("Error: check_retval\n");
      exit(EXIT_FAILURE);
  }

  /* Allocate space for CVODES */
  retval = CVodeInit(cvode_mem, f_SEIR, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) {
      printf("Error: CVodeInit\n");
      exit(EXIT_FAILURE);
  }

  /* Use private function to compute error weights */
  /*retval = CVodeWFtolerances(cvode_mem, ewt);*/

  realtype rtol, atol;

  rtol = RTOL;
  atol = ATOL1;
//  atol[1] = ATOL2;
//  atol[2] = ATOL3;

  retval = CVodeSStolerances(cvode_mem, rtol, atol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) {
    printf("Error: CVodeSStolerances\n");
    exit(EXIT_FAILURE);
  }

  /* Attach user data */
  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval(&retval, "CVodeSetUserData", 1)) {
    printf("Error: CVodeSetUserData\n");
    exit(EXIT_FAILURE);
  }

  /* Create dense SUNMatrix */
  A = SUNDenseMatrix(NEQ, NEQ);
  if (check_retval((void *)A, "SUNDenseMatrix", 0)) {
    printf("Error: SUNDenseMatrix\n");
    exit(EXIT_FAILURE);
  }

  /* Create dense SUNLinearSolver */
  LS = SUNLinSol_Dense(y, A);
  if (check_retval((void *)LS, "SUNDenseLinearSolver", 0)) {
    printf("Error: SUNDenseLinearSolver\n");
    exit(EXIT_FAILURE);
  }

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) {
    printf("Error: CVodeSetLinearSolver\n");
    exit(EXIT_FAILURE);
  }

  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, Jac_SEIR);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) {
    printf("Error: CVodeSetJacFn\n");
    exit(EXIT_FAILURE);
  }

  retval = CVodeSetMaxStep(cvode_mem, MAXSTEPS);

  //printf("\nSIR with variable initial conditions\n");

  /* Sensitivity-related settings */
  if (sensi) {

    /* Set parameter scaling factor */
    /*for(j = 0; j < 3 + g_inpset->num_basis_spline; j++){
      pbar[j] = data->p[j];
    }*/

    /*pbar[0] = 0.1;
    pbar[1] = 0.1;
    pbar[2] = 0.1;*/

    /* Set sensitivity initial conditions */
    yS = N_VCloneVectorArray(kNumParam - 2, y);
    if (check_retval((void *)yS, "N_VCloneVectorArray", 0)){
      printf("Error: N_VCloneVectorArray\n");
      exit(EXIT_FAILURE);
    }

    N_VConst(ZERO, yS[0]);
    N_VConst(ZERO, yS[1]);

    Ith(yS[2], 1) = RCONST(1.0);
    for(j = 0; j <= (NECOMP - 1 + NICOMP); j++){
      Ith(yS[2], 2 + j) = RCONST(0.0);
    }
    Ith(yS[2], NEQ - 1) = -RCONST(1.0);
    Ith(yS[2], NEQ) = -RCONST(1.0);

    Ith(yS[3], 1) = RCONST(0.0);
    Ith(yS[3], 2) = RCONST(1.0);
    for(j = 1; j <= (NECOMP - 1 + NICOMP); j++){
      Ith(yS[3], 2 + j) = RCONST(0.0);
    }
    Ith(yS[3], NEQ - 1) = -RCONST(1.0);
    Ith(yS[3], NEQ) = RCONST(0.0);

    Ith(yS[4], 1) = RCONST(0.0);
    Ith(yS[4], 2) = RCONST(0.0);
    for(j = 1; j <= (NECOMP - 1); j++){
      Ith(yS[4], 2 + j) = RCONST(0.0);
    }
    Ith(yS[4], NECOMP + 2) = RCONST(1.0);
    for(j = 1; j <= (NICOMP - 1); j++){
      Ith(yS[4], NECOMP + 2 + j) = RCONST(0.0);
    }
    Ith(yS[4], NEQ - 1) = -RCONST(1.0);
    Ith(yS[4], NEQ) = RCONST(0.0);

    for (is = 5; is < kNumParam - 2; is++) N_VConst(ZERO, yS[is]);

    /* Call CVodeSensInit1 to activate forward sensitivity computations
       and allocate internal memory for COVEDS related to sensitivity
       calculations. Computes the right-hand sides of the sensitivity
       ODE, one at a time */
    retval = CVodeSensInit1(cvode_mem, kNumParam - 2, sensi_meth, fS_SEIR, yS);
    if(check_retval(&retval, "CVodeSensInit", 1)) {
      printf("Error: CVodeSensInit\n");
      exit(EXIT_FAILURE);
    }

    /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
       variables based on the rolerances supplied for states variables and
       the scaling factor pbar */
    /*
    realtype rtol, atol[3];

    rtol    = RTOL;
    atol[0] = ATOL1;
    atol[1] = ATOL2;
    atol[2] = ATOL3;

    retval = CVodeSensSStolerances(cvode_mem, rtol, atol);
    */
    retval = CVodeSensEEtolerances(cvode_mem);
    if(check_retval(&retval, "CVodeSensEEtolerances", 1)){
      printf("Error: CVodeSensEEtolerances\n");
      exit(EXIT_FAILURE);
    }

    /* Set sensitivity analysis optional inputs */
    /* Call CVodeSetSensErrCon to specify the error control strategy for
       sensitivity variables */
    retval = CVodeSetSensErrCon(cvode_mem, err_con);
    if (check_retval(&retval, "CVodeSetSensErrCon", 1)){
      printf("Error: CVodeSetSensErrCon\n");
      exit(EXIT_FAILURE);
    }

    /* Call CVodeSetSensParams to specify problem parameter information for
       sensitivity calculations */
    /*retval = CVodeSetSensParams(cvode_mem, NULL, pbar, NULL);
    if (check_retval(&retval, "CVodeSetSensParams", 1)) return(1);*/

    /*printf("Sensitivity: YES ");
    if(sensi_meth == CV_SIMULTANEOUS)
      printf("( SIMULTANEOUS +");
    else
      if(sensi_meth == CV_STAGGERED) printf("( STAGGERED +");
      else                           printf("( STAGGERED1 +");
    if(err_con) printf(" FULL ERROR CONTROL )");
    else        printf(" PARTIAL ERROR CONTROL )");*/

  } else {

    //printf("Sensitivity: NO ");

  }

  /* In loop over output points, call CVode, print results, test for error */

  /*printf("\n\n");
  printf("===========================================");
  printf("============================\n");
  printf("     T     Q       H      NST           y1");
  printf("           y2           y3    \n");
  printf("===========================================");
  printf("============================\n");*/

  for (iout=1, tout=T1; iout <= num_out_cvodes; iout++, tout += TAD) {

    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) break;
    y_data = N_VGetArrayPointer(y);
    for(j = 0; j < NEQ; j++){
    	results_y[j][iout - 1] = y_data[j];
        //results_y[j][iout - 1] = Ith(y,j);
        //results_y[j][iout - 1] = *(y_data + j);
    }
//    printf("results_y: %12.9e %12.9e %12.9e \n", results_y[0][iout-1], results_y[1][iout-1],results_y[2][iout-1]);

    //PrintOutput(cvode_mem, t, y);

    /* Call CVodeGetSens to get the sensitivity solution vector after a
       successful return from CVode */
    if (sensi) {
      retval = CVodeGetSens(cvode_mem, &t, yS);
      if (check_retval(&retval, "CVodeGetSens", 1)) break;
      for (jj = 0; jj < kNumParam - 2; jj++) {
        s_data = N_VGetArrayPointer(yS[jj]);
        for(j = 0; j < NEQ; j++){
          results_s[j][iout-1][jj] = s_data[j];
        }
//        printf("results_s_%d: %12.9e %12.9e %12.9e \n", jj, results_s[0][iout-1][jj], results_s[1][iout-1][jj],results_s[2][iout-1][jj]);
      }

//      printf("results_s_2: %12.9e %12.9e %12.9e \n", results_s[0][iout-1][1], results_s[1][iout-1][1],results_s[2][iout-1][1]);
//      printf("results_s_3: %12.9e %12.9e %12.9e \n", results_s[0][iout-1][2], results_s[1][iout-1][2],results_s[2][iout-1][2]);
      //PrintOutputS(yS);
    }
//    printf("-----------------------------------------");
//    printf("------------------------------\n");

  }

  //return results_y;

  /* Print final statistics */
//  PrintFinalStats(cvode_mem, sensi);

  /* Free memory */

  N_VDestroy(y);                    /* Free y vector */
  if (sensi) {
    N_VDestroyVectorArray(yS, kNumParam - 2);  /* Free yS vector */
  }
  free(data);                              /* Free user data */
  CVodeFree(&cvode_mem);                   /* Free CVODES memory */
  SUNLinSolFree(LS);                       /* Free the linear solver memory */
  SUNMatDestroy(A);                        /* Free the matrix memory */

//  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
 */

int f_SEIR(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  //choose the number of components with respect to the model
  int NICOMP = g_inpset->num_comp;
  int NECOMP = g_inpset->num_comp_E;
  int NEQ = 3 + NECOMP + NICOMP;


  realtype y1, y2, y_tot = 0.0;
  UserDataSEIR data;
  realtype p1, p2;
  data = (UserDataSEIR) user_data;
  int iii;
  int n = 1, *np;
  np = &n;
  int *degreep;
  degreep = &data->spline_degree;
  int *mp;
  mp = &data->number_extended_knots;
  int *nfp;
  nfp = &data->number_basis_elements;
  double beta[data->number_basis_elements], beta_t, *tp = &t;

  y1 = Ith(y, 1); y2 = Ith(y, 2);
  p1 = data->p[0];
  p2 = data->p[1];
  for(iii = 5; iii < kNumParam - 2; iii++){
    beta[iii - 5] = data->p[iii];
  }

  beta_t = exp(splinefunction(beta, degreep, np, mp, tp, data->knots, data->basis, nfp));

  for (iii = (1 + NECOMP); iii <= (NECOMP + NICOMP); iii++){
    y_tot += Ith(y, 1 + iii);
  }
  Ith(ydot, 1) = -beta_t * y1 * y_tot / num_pop;
  Ith(ydot, 2) = beta_t * y1 * y_tot / num_pop - NECOMP * p1 * y2;
  for (iii = 2; iii <= NECOMP; iii++){
    Ith(ydot, 1 + iii) = NECOMP * p1 * Ith(y, iii) - NECOMP * p1 * Ith(y, 1 + iii);
  }
  Ith(ydot, 2 + NECOMP) = NECOMP * p1 * Ith(y, 1 + NECOMP) - NICOMP * p2 * Ith(y, 2 + NECOMP);
  for (iii = (2 + NECOMP); iii <= (NECOMP + NICOMP); iii++){
    Ith(ydot, 1 + iii) = NICOMP * p2 * Ith(y, iii) - NICOMP * p2 * Ith(y, 1 + iii);
  }
  Ith(ydot, NEQ - 1) = NICOMP * p2 * Ith(y, 1 + NECOMP + NICOMP);
  Ith(ydot, NEQ) = beta_t * y1 * y_tot / num_pop;

  return(0);
}


/*
 * Jacobian routine. Compute J(t,y).
 */

int Jac_SEIR(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  //choose the number of components with respect to the model
  int NICOMP = g_inpset->num_comp;
  int NECOMP = g_inpset->num_comp_E;
  int NEQ = 3 + NECOMP + NICOMP;

  realtype y1, y_tot = 0.0;
  UserDataSEIR data;
  realtype p1, p2;
  data = (UserDataSEIR) user_data;
  int iii;
  int n = 1, *np;
  np = &n;
  int *degreep;
  degreep = &data->spline_degree;
  int *mp;
  mp = &data->number_extended_knots;
  int *nfp;
  nfp = &data->number_basis_elements;
  double beta[data->number_basis_elements], beta_t, *tp = &t;

  y1 = Ith(y,1);
  p1 = data->p[0];
  p2 = data->p[1];
  for(iii = 5; iii < kNumParam - 2; iii++){
    beta[iii - 5] = data->p[iii];
  }

  beta_t = exp(splinefunction(beta, degreep, np, mp, tp, data->knots, data->basis, nfp));
  for (iii = (1 + NECOMP); iii <= (NECOMP + NICOMP); iii++){
    y_tot += Ith(y, 1 + iii);
  }

  IJth(J, 1, 1) = -beta_t * y_tot / num_pop;
  IJth(J, 2, 1) =  beta_t * y_tot / num_pop;
  IJth(J, NEQ, 1) =  beta_t * y_tot / num_pop;

  for (iii = 2; iii <= (1 + NECOMP); iii++){
    IJth(J, iii, iii) = - NECOMP * p1;
    IJth(J, 1 + iii, iii) = NECOMP * p1;
  }

  for (iii = (1 + NECOMP); iii <= (NECOMP + NICOMP); iii++){
    IJth(J, 1, 1 + iii) = -beta_t * y1 / num_pop;
    IJth(J, 2, 1 + iii) = beta_t * y1 / num_pop;
    IJth(J, 1 + iii, 1 + iii) = -NICOMP * p2;
    IJth(J, 2 + iii, 1 + iii) = NICOMP * p2;
    IJth(J, NEQ, 1 + iii) = beta_t * y1 / num_pop;
  }

  return(0);
}

/*
 * fS routine. Compute sensitivity r.h.s.
 */

int fS_SEIR(int Ns, realtype t, N_Vector y, N_Vector ydot,
              int iS, N_Vector yS, N_Vector ySdot,
              void *user_data, N_Vector tmp1, N_Vector tmp2)
{
  //choose the number of components with respect to the model
  int NICOMP = g_inpset->num_comp;
  int NECOMP = g_inpset->num_comp_E;
  int NEQ = 3 + NECOMP + NICOMP;


  UserDataSEIR data;
  realtype p1, p2;
  realtype y1, y2, y_tot = 0.0;
  realtype s1, s2;
  realtype sd[NEQ];
  data = (UserDataSEIR) user_data;
  int iii;
  int n = 1, *np;
  np = &n;
  int *degreep;
  degreep = &data->spline_degree;
  int *mp;
  mp = &data->number_extended_knots;
  int *nfp;
  nfp = &data->number_basis_elements;
  double beta[data->number_basis_elements], beta_t, *tp = &t;

  p1 = data->p[0];
  p2 = data->p[1];
  for(iii = 5; iii < kNumParam - 2; iii++){
    beta[iii - 5] = data->p[iii];
  }

  y1 = Ith(y, 1);  y2 = Ith(y, 2);
  s1 = Ith(yS, 1); s2 = Ith(yS, 2);

  beta_t = exp(splinefunction(beta, degreep, np, mp, tp, data->knots, data->basis, nfp));
  splinebasis(degreep, np, mp, tp, data->knots, data->basis);
  for (iii = (1 + NECOMP); iii <= (NECOMP + NICOMP); iii++){
    y_tot += Ith(y, 1 + iii);
  }

  sd[0] = -beta_t * y_tot * s1 / num_pop;
  sd[1] = beta_t * y_tot * s1 / num_pop  -  NECOMP * p1 * s2;
  sd[NECOMP + 1] = NECOMP * p1 * Ith(yS, 1 + NECOMP) - NICOMP * p2 * Ith(yS, 2 + NECOMP);
  sd[NEQ - 2] = NICOMP * p2 * Ith(yS, 1 + NECOMP + NICOMP);
  sd[NEQ - 1] = beta_t * y_tot * s1/ num_pop;
  for (iii = (2 + NECOMP); iii <= (1 + NECOMP + NICOMP); iii++){
    sd[0] += - (beta_t * y1 / num_pop) * Ith(yS, iii);
    sd[1] += (beta_t * y1 / num_pop) * Ith(yS, iii);
    sd[NEQ - 1] += (beta_t * y1 / num_pop) * Ith(yS, iii);
  }

  for(iii = 2; iii <= NECOMP; iii++){
    sd[iii] = NECOMP * p1 * Ith(yS, iii) - NECOMP * p1 * Ith(yS, 1 + iii);
  }

  for(iii = (NECOMP + 2); iii <= (NECOMP + NICOMP); iii++){
    sd[iii] = NICOMP * p2 * Ith(yS, iii) - NICOMP * p2 * Ith(yS, 1 + iii);
  }

  if (iS == 0){
    sd[1] += - NECOMP * y2;
    sd[NECOMP + 1] += NECOMP * Ith(y, 1 + NECOMP);
    for (iii = 2; iii <= NECOMP; iii++){
      sd[iii] += NECOMP * Ith(y, iii) - NECOMP * Ith(y, iii + 1);
    }
  } else{
    if(iS == 1){
      sd[NECOMP + 1] += - NICOMP * Ith(y, 2 + NECOMP);
      for(iii = (1 + NECOMP); iii < (NICOMP + NECOMP); iii++){
        sd[1 + iii] += NICOMP * Ith(y, 1 + iii) - NICOMP * Ith(y, 2 + iii);
      }
      sd[NEQ - 2] += NICOMP * Ith(y, 1 + NECOMP + NICOMP);
    } else{
      for (iii = 5; iii < kNumParam - 2; iii++){
        if(iS == iii){
          sd[0] += -beta_t * data->basis[iii - 5] * y1 * y_tot / num_pop;
          sd[1] += beta_t * data->basis[iii - 5] * y1 * y_tot / num_pop;
          sd[NEQ - 1] += beta_t * data->basis[iii - 5] * y1 * y_tot / num_pop;
        }
      }
    }
  }

  for(iii = 0; iii < NEQ; iii++){
    Ith(ySdot, iii + 1) = sd[iii];
  }

  return(0);
}


void printODE_SEIR(FILE *odefp) {
    //choose the number of components with respect to the model
    int NICOMP = g_inpset->num_comp;
    int NECOMP = g_inpset->num_comp_E;
    int NEQ = 3 + NECOMP + NICOMP;

    int i, j;

    double t_alpha, t_gamma, t_s0, t_e0, t_i0;

    t_alpha = exp(state[curr]->pos[0]) + alpha0;
    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma from the state vector to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[1]));
        }
        else{
            // Transform the phi inverse parameter using an exponential function and a shift
            t_gamma =  exp(state[curr]->pos[1]) + gamma0;
        }
        
    }
    t_s0 = 0.0 + (num_pop - 0)/(1 + exp(-state[curr]->pos[2]));
    t_e0 = 0.0 + (num_pop - 0)/(1 + exp(-state[curr]->pos[3]));
    t_i0 = 0.0 + (num_pop - 0)/(1 + exp(-state[curr]->pos[4]));

    char *sensitivities[] = {"-nosensi"};
    double parameters_odes[kNumParam - 2];
    double results_y[NEQ][kNumObs];
    double results_s[NEQ][kNumObs][kNumParam - 2];
    parameters_odes[0] = t_alpha;
    parameters_odes[1] = t_gamma;
    parameters_odes[2] = t_s0;
    parameters_odes[3] = t_e0;
    parameters_odes[4] = t_i0;
    for(i = 5; i < kNumParam - 2; i++){
        parameters_odes[i] = state[curr]->pos[i + 2];
    }

    int m = (g_inpset->num_basis_spline - g_inpset->spline_pol_degree - 1) + 2*(g_inpset->spline_pol_degree + 1);
    int nf = g_inpset->num_basis_spline;
    double dx, knots[m], basis[g_inpset->num_basis_spline];
    dx = (tfin_spl_bas - X0)/(g_inpset->num_basis_spline - g_inpset->spline_pol_degree);
    for (i = 0; i < m; i++){
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }

    SEIR_CVODES(results_y, results_s, parameters_odes, g_inpset->spline_pol_degree, m, nf, knots, basis, sensitivities);

    for (i = 0; i < kNumObs; i++) {
      for(j = 0; j < NEQ; j++) {
        fprintf(odefp, "%f ", results_y[j][i]);
      }
    }
    fprintf(odefp, "\n");
}
