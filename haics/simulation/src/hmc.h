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

#ifndef HMC_H
#define HMC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "Globals.h"
#include "integrators.h"

double HMC();

//function for calculating eigenvalues of a general matrix, link: http://www.netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen_ga66e19253344358f5dee1e60502b9e96f.html
//extern void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA, double* WR, double* WI, double* VL, int* LDVL, double* VR, int* LDVR, double* WORK, int* LWORK, int* INFO);

//function for calculating eigenvalues of a symmetric matrix, link: http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html
extern void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W, double* WORK, int* LWORK, int* INFO);

#endif
