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

#include <string.h>
#include <math.h>

#include "Globals.h"
#include "Definitions.h"

// scale user given time step
// scale user given time step
double scale_timestep(double ss) {
  double ss_scaled = 0.;

  ss_scaled = fitting_factor * ss;

  return ss_scaled;
}

double FindIntegrator2stage(double *b_coeff, double ss) {
  double ss_scaled, b;
  int ind;

  ss_scaled = scale_timestep(ss);
  
  if (ss_scaled >= 4){
    b = 0.25;
    fprintf(red_lfp, "WARNING: Scaled timestep (%f) exceeds stability limit. b set to 0.25 (VV).\n", ss_scaled);
  } else if (ss_scaled > 0){

    ind = round(ss_scaled * 10000);
    
    //safety checks
    if(ind > hbar_to_b_num_disc_steps){
      ind = hbar_to_b_num_disc_steps;
    }else if(ind == 0){
      ind = 1;
    }

    //assign b
    b = b_coeff[ind-1];
    
  } else{
    perror("adapt_integrator.c - There was a problem in the else-if clauses while assigning the coefficient b in the FindIntegrator2stage function.\n");
    exit(EXIT_FAILURE);
  }

  return b;
}

double FindIntegrator3stage(double *b_coeff, double ss) {
  double ss_scaled, b;
  int ind;

  ss_scaled = scale_timestep(ss);
  
  if (ss_scaled >= 6){
    b = 1.0f/6.0f;
    fprintf(red_lfp, "WARNING: Scaled timestep (%f) exceeds stability limit. b set to 1/6 (Strang).\n", ss_scaled);
  } else if (ss_scaled > 0){
    
    ind = round(ss_scaled * 10000);

    //safety checks
    if(ind > hbar_to_b_num_disc_steps_3s){
      ind = hbar_to_b_num_disc_steps_3s;
    }else if(ind == 0){
      ind = 1;
    }

    //assign b
    b = b_coeff[ind-1];
    
  } else{
    perror("adapt_integrator.c - There was a problem in the else-if clauses while assigning the coefficient b in the FindIntegrator3stage function.\n");
    exit(EXIT_FAILURE);
  }

  return b;
}

double fnCalculateFittingFactor(char *fitting_factor_approach, int use_freqs, double dt_vv, double max_freqs, double sum_freqs, int D, int num_accepted, int num_accepted_pos, int num_proposed, int num_proposed_pos, double dH_avg, double dHpos_avg, double dH2_avg){

  double exp_dH_scaling = 0.;
  double fitting_filtAR, fitting_Gupta, exp_dH_filtAR, exp_dH_Gupta; //for mean and mean_new cases
  double AR = (double)num_accepted/num_proposed;
  double AR_filt = (double)num_accepted_pos/num_proposed_pos;
  double AR_filt_new = (double)num_accepted_pos/num_proposed;
  double fitting_factor = 1.;

  if(strcmp(fitting_factor_approach, "filtAR") == 0){
    exp_dH_scaling = - log(AR_filt);
  }else if(strcmp(fitting_factor_approach, "filtAR_new") == 0){
    exp_dH_scaling = - log(AR_filt_new);
  }else if(strcmp(fitting_factor_approach, "Gupta") == 0){
    exp_dH_scaling = 4 * PI * (1 - AR) * (1 - AR);
  }else if(strcmp(fitting_factor_approach, "TightUpper") == 0){
    exp_dH_scaling = - log(AR_filt)*num_proposed_pos/num_proposed;
  }else if(strcmp(fitting_factor_approach, "dH_average") == 0){
    exp_dH_scaling = dH_avg;
  }else if(strcmp(fitting_factor_approach, "dHpos_average") == 0){
    exp_dH_scaling = dHpos_avg;
  }else if(strcmp(fitting_factor_approach, "H_pos_avg") == 0){
    exp_dH_scaling = dHpos_avg*num_proposed/num_proposed_pos;
  }else if(strcmp(fitting_factor_approach, "H_pos_acc_avg") == 0){
    exp_dH_scaling = dHpos_avg*num_proposed/num_accepted_pos;
  }else if(strcmp(fitting_factor_approach, "Gupta_pos") == 0){
    exp_dH_scaling = 4 * PI * (1 - 2*num_proposed_pos/num_proposed + 2*dHpos_avg*num_proposed/num_proposed_pos) * (1 - 2*num_proposed_pos/num_proposed + 2*dHpos_avg*num_proposed/num_proposed_pos);
  }else if(strcmp(fitting_factor_approach, "dH2_average") == 0){
    exp_dH_scaling = dH2_avg;
  }else if(strcmp(fitting_factor_approach, "mean") == 0){
    exp_dH_filtAR = - log(AR_filt);
    exp_dH_Gupta = 4 * PI * (1 - AR) * (1 - AR);
    if(use_freqs == 0){
      fitting_filtAR = fmax(1, 1 / (dt_vv * max_freqs) * pow(32*exp_dH_filtAR/D, 1.0/6.0));
      fitting_Gupta = fmax(1, 1 / (dt_vv * max_freqs) * pow(32*exp_dH_Gupta/D, 1.0/6.0));
    }else if(use_freqs == 1){
      fitting_filtAR = fmax(1, 1 / dt_vv * pow(32*exp_dH_filtAR/sum_freqs, 1.0/6.0));
      fitting_Gupta = fmax(1, 1 / dt_vv * pow(32*exp_dH_Gupta/sum_freqs, 1.0/6.0));
    }else{
      fprintf(stderr, "adapt_integrator.c - fnCalculateFittingFactor, use_freqs_fitting has to be 0 or 1, it is %i\n", use_freqs);
      exit(EXIT_FAILURE);
    }
    fitting_factor = 2 * fitting_Gupta * fitting_filtAR / (fitting_Gupta + fitting_filtAR);
  }else if(strcmp(fitting_factor_approach, "mean_new") == 0){
    exp_dH_filtAR = - log(AR_filt_new);
    exp_dH_Gupta = 4 * PI * (1 - AR) * (1 - AR);
    if(use_freqs == 0){
      fitting_filtAR = fmax(1, 1 / (dt_vv * max_freqs) * pow(32*exp_dH_filtAR/D, 1.0/6.0));
      fitting_Gupta = fmax(1, 1 / (dt_vv * max_freqs) * pow(32*exp_dH_Gupta/D, 1.0/6.0));
    }else if(use_freqs == 1){
      fitting_filtAR = fmax(1, 1 / dt_vv * pow(32*exp_dH_filtAR/sum_freqs, 1.0/6.0));
      fitting_Gupta = fmax(1, 1 / dt_vv * pow(32*exp_dH_Gupta/sum_freqs, 1.0/6.0));
    }else{
      fprintf(stderr, "adapt_integrator.c - fnCalculateFittingFactor, use_freqs_fitting has to be 0 or 1, it is %i\n", use_freqs);
      exit(EXIT_FAILURE);
    }
    fitting_factor = 2 * fitting_Gupta * fitting_filtAR / (fitting_Gupta + fitting_filtAR);
  }else{
    fprintf(stderr, "adapt_integrator.c - fnCalculateFittingFactor, fitting_factor_approach has to be filtAR, filtAR_new, Gupta, TightUpper, dH_average, dHpos_average, H_pos_avg, H_pos_acc_avg, Gupta_pos, dH2_average, mean, mean_new; it is %s\n", fitting_factor_approach);
    exit(EXIT_FAILURE);
  }

  if((strcmp(fitting_factor_approach, "mean") != 0) && (strcmp(fitting_factor_approach, "mean_new") != 0)){
    if(use_freqs == 0){
      fitting_factor = fmax(1, 1 / (dt_vv * max_freqs) * pow(32*exp_dH_scaling/D, 1.0/6.0));
    }else if(use_freqs == 1){
      fitting_factor = fmax(1, 1 / dt_vv * pow(32*exp_dH_scaling/sum_freqs, 1.0/6.0));
    }else{
      fprintf(stderr, "adapt_integrator.c - fnCalculateFittingFactor, use_freqs_fitting has to be 0 or 1, it is %i\n", use_freqs);
      exit(EXIT_FAILURE);
    }
  }

  return fitting_factor;
}
