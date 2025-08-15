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

#include <math.h>
#include <time.h>
//#include <gsl/gsl_rng.h>

#include "read_input.h"
#include "hmc.h"
#include "Globals.h"


//main()

int main(int argc, char *argv[])
{
  char timestr[100];
  time_t now;
  char filename[80];

  if (argc != 2) {
    fprintf(stderr, "main.c - ID needed\n");
    fprintf(stderr, "main.c - Usage: %s <ID> \n", argv[0]);
    exit(1);
  }

  now = time(0); //take time for logging
  strftime(timestr, 100, "%Y-%m-%d %H:%M:%S.", localtime(&now));

  //reduced logfile
  sprintf(filename, "./output/%s%s%s%s", argv[1], "/reduced_logfile_", argv[1],".txt");
  red_lfp = fopen(filename, "w+");
  fprintf(red_lfp, "Created on %s\n", timestr);

  // get the parameters from the input file and set up the whole system
  id_global = argv[1];
  ReadInput(argv[1]);

  // perform HMC sampling
  HMC();

  now = time(0);
  strftime(timestr, 100, "%Y-%m-%d %H:%M:%S.", localtime(&now));

  fprintf(red_lfp, "Finished on %s\n", timestr);
  fclose(red_lfp);

  return 0;
}
