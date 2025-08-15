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

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

/* Problem Constants */

#define NEQ   3                /* number of equations  */
//#define Y1    RCONST(99.0)   /* initial y components */
//#define Y2    RCONST(1.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1e-10)   /* scalar relative tolerance */
#define ATOL1 RCONST(1e-10)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1e-10)
#define ATOL3 RCONST(1e-10)
#define T0    RCONST(0.0)      /* initial time */
#define T1    RCONST(1.0)      /* first output time */
//#define TMULT RCONST(10.0)  /* output time factor */
#define TAD RCONST(1.0)       /* output time factor */
#define NOUT  14              /* number of output times */

#define NP    3               /* number of problem parameters */
#define NS    3               /* number of sensitivities computed */
//#define NPOP    RCONST(763.0)             /* number of TOTAL POPULATION */

#define ZERO  RCONST(0.0)


/* Type : UserData */

typedef struct {
  realtype p[3];           /* problem parameters */
} *UserData0;

/* Prototypes of functions by CVODES */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot,
              int iS, N_Vector yS, N_Vector ySdot,
              void *user_data, N_Vector tmp1, N_Vector tmp2);

/*static int ewt(N_Vector y, N_Vector w, void *user_data);*/

/* Prototypes of private functions */
static void ProcessArgs(char *argv[],
                        booleantype *sensi, int *sensi_meth,
                        booleantype *err_con);
static void WrongArgs(char *name);
/*static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
static void PrintOutputS(N_Vector *uS);
static void PrintFinalStats(void *cvode_mem, booleantype sensi);*/
static int check_retval(void *returnvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

void SIR_CVODES(double results_y[][NOUT], double results_s[][NOUT][NS], double state_params[NP], char *sensitivity_parameters[])
{

  //printf("Entraaaaaaaaaaaaaaaaaaaaaaaaaaamos\n");

  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  UserData0 data;
  realtype t, tout;
  N_Vector y;
  int iout, retval, j, jj;
  //static double results_y[NEQ][NOUT];
  realtype *y_data;
  //double *y_data;

  //realtype pbar[NS];
  int is;
  N_Vector *yS;
  booleantype sensi, err_con;
  int sensi_meth;
  //double results_s[NEQ][NOUT][NS];
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
  data = (UserData0) malloc(sizeof *data);
  if (check_retval((void *)data, "malloc", 2)) { 
    printf("Error: UserData\n");
    exit(EXIT_FAILURE);
  }


  //printf("pbar1: %f \n", state_params[0]);
  //printf("pbar2: %f \n", state_params[1]);
  //printf("pbar3: %f \n", state_params[2]);

  data->p[0] = state_params[0];
  data->p[1] = state_params[1];
  data->p[2] = state_params[2];

  /* Initial conditions */
  y = N_VNew_Serial(NEQ);
  if (check_retval((void *)y, "N_VNew_Serial", 0)){ 
    printf("Error: N_VNew_Serial\n");
    exit(EXIT_FAILURE);
  }

  Ith(y,1) = num_pop - data->p[2];
  Ith(y,2) = data->p[2];
  Ith(y,3) = Y3;

  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) {
      printf("Error: check_retval\n");
      exit(EXIT_FAILURE);
  }

  /* Allocate space for CVODES */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { 
      printf("Error: CVodeInit\n");
      exit(EXIT_FAILURE);
  }

  /* Use private function to compute error weights */
  /*retval = CVodeWFtolerances(cvode_mem, ewt);
  if (check_retval(&retval, "CVodeSetEwtFn", 1)) {
    printf("Error: CVodeWFtolerances\n");
    exit(EXIT_FAILURE);
  }*/
  realtype rtol, atol[3];

  rtol    = RTOL;
  atol[0] = ATOL1;
  atol[1] = ATOL2;
  atol[2] = ATOL3;

  retval = CVodeSStolerances(cvode_mem, rtol, atol[0]);
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
  LS = SUNDenseLinearSolver(y, A);
  //LS = SUNLinSol_Dense(y, A);
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
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) { 
    printf("Error: CVodeSetJacFn\n");
    exit(EXIT_FAILURE);
  }


  /* To test: Stan*/
  double init_step = 0;
  retval = CVodeSetInitStep(cvode_mem, init_step);

  long int max_err_test_fails = 20;  // NOLINT(runtime/int)
  retval = CVodeSetMaxErrTestFails(cvode_mem, max_err_test_fails);

  long int max_conv_fails = 50;  // NOLINT(runtime/int)
  retval = CVodeSetMaxConvFails(cvode_mem, max_conv_fails);

  retval = CVodeSetMaxStep(cvode_mem, 1e+3);

  //printf("\nSIR with variable initial conditions\n");

  /* Sensitivity-related settings */
  if (sensi) {

    /* Set parameter scaling factor */
    /* pbar[0] = data->p[0];
    pbar[1] = data->p[1];
    pbar[2] = data->p[2];*/

    //printf("pbar1: %f \n", data->p[0]);
    //printf("pbar2: %f \n", data->p[1]);
    //printf("pbar3: %f \n", data->p[2]);

    /* Set sensitivity initial conditions */
    yS = N_VCloneVectorArray(NS, y);
    if (check_retval((void *)yS, "N_VCloneVectorArray", 0)){ 
      printf("Error: N_VCloneVectorArray\n");
      exit(EXIT_FAILURE);
    }

    for (is = 0; is<NS;is++) N_VConst(ZERO, yS[is]);
    Ith(yS[NS - 1],1) = RCONST(-1.0);
    Ith(yS[NS - 1],2) = RCONST(1.0);
    Ith(yS[NS - 1],3) = RCONST(0.0);

    /* Call CVodeSensInit1 to activate forward sensitivity computations
       and allocate internal memory for COVEDS related to sensitivity
       calculations. Computes the right-hand sides of the sensitivity
       ODE, one at a time */
    retval = CVodeSensInit1(cvode_mem, NS, sensi_meth, fS, yS);
    if(check_retval(&retval, "CVodeSensInit", 1)) { 
      printf("Error: CVodeSensInit\n");
      exit(EXIT_FAILURE);
    }

    /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
       variables based on the rolerances supplied for states variables and
       the scaling factor pbar */
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
    retval = CVodeSetSensParams(cvode_mem, NULL, pbar, NULL);
    if (check_retval(&retval, "CVodeSetSensParams", 1)){ 
      printf("Error: CVodeSetSensParams\n");
      exit(EXIT_FAILURE);
    }

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

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += TAD) {

    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) break;
    y_data = N_VGetArrayPointer(y);
    for(j = 0; j < NEQ; j++){
      results_y[j][iout - 1] = y_data[j];
        //results_y[j][iout - 1] = Ith(y,j);
        //results_y[j][iout - 1] = *(y_data + j);
    }
    //printf("results_y: %12.4e %12.4e %12.4e \n", results_y[0][iout-1], results_y[1][iout-1],results_y[2][iout-1]);

    //PrintOutput(cvode_mem, t, y);

    /* Call CVodeGetSens to get the sensitivity solution vector after a
       successful return from CVode */
    if (sensi) {
      retval = CVodeGetSens(cvode_mem, &t, yS);
      if (check_retval(&retval, "CVodeGetSens", 1)) break;
      for (jj = 0; jj < NS; jj++) {
        s_data = N_VGetArrayPointer(yS[jj]);
        for(j = 0; j < NEQ; j++){
          results_s[j][iout-1][jj] = s_data[j];
        }
      }
      /*printf("results_s_1: %12.4e %12.4e %12.4e \n", results_s[0][iout-1][0], results_s[1][iout-1][0],results_s[2][iout-1][0]);
      printf("results_s_2: %12.4e %12.4e %12.4e \n", results_s[0][iout-1][1], results_s[1][iout-1][1],results_s[2][iout-1][1]);
      printf("results_s_3: %12.4e %12.4e %12.4e \n", results_s[0][iout-1][2], results_s[1][iout-1][2],results_s[2][iout-1][2]);
      //PrintOutputS(yS);*/
    }
    /*printf("-----------------------------------------");
    printf("------------------------------\n");*/

  }

  //return results_y;

  /* Print final statistics */
  /*PrintFinalStats(cvode_mem, sensi);*/

  /* Free memory */

  N_VDestroy(y);                    /* Free y vector */
  if (sensi) {
    N_VDestroyVectorArray(yS, NS);  /* Free yS vector */
  }
  free(data);                              /* Free user data */
  CVodeFree(&cvode_mem);                   /* Free CVODES memory */
  SUNLinSolFree(LS);                       /* Free the linear solver memory */
  SUNMatDestroy(A);                        /* Free the matrix memory */

  /*return(0); */
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, y3, yd1, yd2;
  UserData0 data;
  realtype p1, p2;

  y1 = Ith(y,1);
  y2 = Ith(y,2);
  y3 = Ith(y,3);

  data = (UserData0) user_data;
  p1 = data->p[0];
  p2 = data->p[1];

  yd1 = Ith(ydot,1) = -p1*y1*y2/num_pop;
  yd2 = Ith(ydot,2) = p1*y1*y2/num_pop - p2*y2;
        Ith(ydot,3) = -yd1 - yd2;

  return(0);
}


/*
 * Jacobian routine. Compute J(t,y).
 */

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2;
  UserData0 data;
  realtype p1, p2;

  y1 = Ith(y,1); 
  y2 = Ith(y,2);
  data = (UserData0) user_data;
  p1 = data->p[0]; 
  p2 = data->p[1];

  IJth(J,1,1) = -p1*y2/num_pop;
  IJth(J,1,2) = -p1*y1/num_pop;
  IJth(J,1,3) = 0;
  IJth(J,2,1) =  p1*y2/num_pop;
  IJth(J,2,2) = p1*y1/num_pop-p2;
  IJth(J,2,3) = 0;
  IJth(J,3,1) = 0;
  IJth(J,3,2) = p2;
  IJth(J,3,3) = 0;

  return(0);
}

/*
 * fS routine. Compute sensitivity r.h.s.
 */

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot,
              int iS, N_Vector yS, N_Vector ySdot,
              void *user_data, N_Vector tmp1, N_Vector tmp2)
{
  UserData0 data;
  realtype p1, p2;
  realtype y1, y2;
  realtype s1, s2;
  realtype sd1, sd2, sd3;

  data = (UserData0) user_data;
  p1 = data->p[0];
  p2 = data->p[1];

  y1 = Ith(y,1);
  y2 = Ith(y,2);
  s1 = Ith(yS,1); 
  s2 = Ith(yS,2);

  sd1 = -p1*y2*s1/num_pop - p1*y1*s2/num_pop;
  sd2 = p1*y2*s1/num_pop + (p1*y1/num_pop - p2)*s2;
  sd3 = p2*s2;

  switch (iS) {
  case 0:
    sd1 += -y1*y2/num_pop;
    sd2 +=  y1*y2/num_pop;
    break;
  case 1:
    sd2 += -y2;
    sd3 += y2;
    break;
  case 2:
    break;
  }

  Ith(ySdot,1) = sd1;
  Ith(ySdot,2) = sd2;
  Ith(ySdot,3) = sd3;

  return(0);
}

/*
 * EwtSet function. Computes the error weights at the current solution.
 */
/*
static int ewt(N_Vector y, N_Vector w, void *user_data)
{
  int i;
  realtype yy, ww, rtol, atol[3];

  rtol    = RTOL;
  atol[0] = ATOL1;
  atol[1] = ATOL2;
  atol[2] = ATOL3;

  for (i=1; i<=3; i++) {
    yy = Ith(y,i);
    ww = rtol * SUNRabs(yy) + atol[i-1];
    if (ww <= 0.0) return (-1);
    Ith(w,i) = 1.0/ww;
  }

  return(0);
}
*/
/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Process and verify arguments to cvsfwddenx.
 */

static void ProcessArgs(char *argv[], booleantype *sensi, int *sensi_meth, booleantype *err_con)
{
  *sensi = SUNFALSE;
  *sensi_meth = -1;
  *err_con = SUNFALSE;

//  if (argc < 2) WrongArgs(argv[0]);

  if (strcmp(argv[0],"-nosensi") == 0)
    *sensi = SUNFALSE;
  else if (strcmp(argv[0],"-sensi") == 0)
    *sensi = SUNTRUE;
  else
    WrongArgs("SIR_CVODES");

  if (*sensi) {

    if (strcmp(argv[1],"sim") == 0)
      *sensi_meth = CV_SIMULTANEOUS;
    else if (strcmp(argv[1],"stg") == 0)
      *sensi_meth = CV_STAGGERED;
    else if (strcmp(argv[1],"stg1") == 0)
      *sensi_meth = CV_STAGGERED1;
    else
      WrongArgs("SIR_CVODES");

    if (strcmp(argv[2],"t") == 0)
      *err_con = SUNTRUE;
    else if (strcmp(argv[2],"f") == 0)
      *err_con = SUNFALSE;
    else
      WrongArgs("SIR_CVODES");
  }

}

static void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = t or f\n");

    exit(0);
}

/*
 * Print current t, step count, order, stepsize, and solution.
 */

/*
static void PrintOutput(void *cvode_mem, realtype t, N_Vector u)
{
  long int nst;
  int qu, retval;
  realtype hu, *udata;

  udata = N_VGetArrayPointer(u);

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetLastOrder(cvode_mem, &qu);
  check_retval(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetLastStep(cvode_mem, &hu);
  check_retval(&retval, "CVodeGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#else
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#endif

  printf("                  Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", udata[0], udata[1], udata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#endif

}
*/

/*
 * Print sensitivities.
*/
/*
static void PrintOutputS(N_Vector *uS)
{
  realtype *sdata;

  sdata = N_VGetArrayPointer(uS[0]);
  printf("                  Sensitivity 1  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e  \n", sdata[0], sdata[1][2]);
#endif

  sdata = N_VGetArrayPointer(uS[1]);
  printf("                  Sensitivity 2  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

  sdata = N_VGetArrayPointer(uS[2]);
  printf("                  Sensitivity 3  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
}
*/

/*
 * Print some final statistics from the CVODES memory.
 */
/*
static void PrintFinalStats(void *cvode_mem, booleantype sensi)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nje, nfeLS;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  if (sensi) {
    retval = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    check_retval(&retval, "CVodeGetSensNumRhsEvals", 1);
    retval = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    check_retval(&retval, "CVodeGetNumRhsEvalsSens", 1);
    retval = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    check_retval(&retval, "CVodeGetSensNumLinSolvSetups", 1);
    retval = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
    check_retval(&retval, "CVodeGetSensNumErrTestFails", 1);
    retval = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
    check_retval(&retval, "CVodeGetSensNumNonlinSolvIters", 1);
    retval = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
    check_retval(&retval, "CVodeGetSensNumNonlinSolvConvFails", 1);
  }

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  //printf("\nFinal Statistics\n\n");
  //printf("nst     = %5ld\n\n", nst);
  //printf("nfe     = %5ld\n",   nfe);
  //printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  //printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  //if(sensi) {
  //  printf("\n");
  //  printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
  //  printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
  //  printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  //}

  //printf("\n");
  //printf("nje    = %5ld    nfeLS     = %5ld\n", nje, nfeLS);
  

}
*/

/*
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns an integer value so check if
 *             retval < 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
        funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
    return(1); }

  return(0);
}
