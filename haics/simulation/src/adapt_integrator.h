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

#ifndef HaMonCaSa_adapt_integrator_h
#define HaMonCaSa_adapt_integrator_h

double scale_timestep(double ss);

double FindIntegrator2stage(double *b_coeff, double ss);

double FindIntegrator3stage(double *b_coeff, double ss);

extern double fnCalculateFittingFactor(char *fitting_factor_approach, int use_freqs, double dt_vv, double max_freqs, double sum_freqs, int D, int num_accepted, int num_accepted_pos, int num_proposed, int num_proposed_pos, double dH_avg, double dHpos_avg, double dH2_avg);

#endif
