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

#include <stdlib.h>
#include <time.h>

#include "utils.h"
#include "Definitions.h"

/*************************************************
 Copy state
 ************************************************/
void CopyState(TypeState *st_copy, TypeState *st, int dim) {
  memcpy(st_copy->mom, st->mom, dim*sizeof(double));
  memcpy(st_copy->pos, st->pos, dim*sizeof(double));
  memcpy(st_copy->grad, st->grad, dim*sizeof(double));
  
}

/*************************************************
 Allocate state
 ************************************************/
void SAlloc(TypeState **st, int dim) {
  *st = (TypeState *)malloc(sizeof(TypeState));
  (*st)->pos  = (double *)malloc(dim*sizeof(double));
  (*st)->mom  = (double *)malloc(dim*sizeof(double));
  (*st)->grad = (double *)malloc(dim*sizeof(double));
}

/*************************************************
 Free state
 ************************************************/
void SFree(TypeState **st) {
  free((*st)->pos);
  free((*st)->mom);
  free((*st)->grad);
  free(*st);
}

/*************************************************
 MatrixAlloc
 routine to allocate memory for a matrix with r rows and
 c columns, with consecutive elements in memory (row major)
 ************************************************/

double ** MatrixAlloc(int r, int c) {
  
  double **matrix;
  int i, j;
  
  matrix = (double**)malloc(r*sizeof(double *));//allocates a pointer to r double pointers
  matrix[0] = (double*)malloc(r*c*sizeof(double));//allocates r x c doubles
  for (i = 1; i < r; i++){
    matrix[i] = matrix[i-1] + c;
  }

  //init with zeroes
  for(i = 0; i<r; i++){
    for(j = 0; j<c; j++){
      matrix[i][j] = 0.;
    }
  }
  
  return matrix;
}

/*************************************************
 MatrixFree
 routine to free memory for a matrix with consecutive elements in memory
 ************************************************/

void MatrixFree(double **matrix) {
  
  free(matrix[0]);
  free(matrix);
}