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
#include <gsl/gsl_sf_psi.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "sir_stand_model_synthetic.h"
#include "sir_seir_common_functions.h"
#include "Globals.h"
#include "utils.h"
#include "sir_stand_functions.h"
#include "read_input.h"

#define M_PI 3.14159265358979323846
#define X0 0.0

#define gamma0 1e-2
#define phiInv0 1e-2
#define beta0 1e-2
#define n_chains 10
#define I0_mean 1.0
#define I0_sd 1.0
#define S0_constant 1
#define I0_low 1.0
#define I0_constant 1
#define I0_value 1.0

double PriorSIKR_stand_Prevalence(){

    double t_gamma, t_s0, t_i0, t_phi_inv, t_beta;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma from the state vector to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[0]));
        }
        else{
            // Transform the phi inverse parameter using an exponential function and a shift
            t_gamma =  exp(state[curr]->pos[0]) + gamma0;
        }
        
    }

    // Transform the initial susceptible count parameter to stay within [0, num_pop]
    t_s0 = 0.0 + (num_pop - 0)/(1 + exp(-state[curr]->pos[1]));
    // printf("S_0: %f\n", t_s0);


    // Transform the initial infected count parameter to stay within [0, num_pop]
    t_i0 = I0_low + (num_pop - I0_low)/(1 + exp(-state[curr]->pos[2]));

    // Transform the phi inverse parameter using an exponential function and a shift
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;

    // Transform the beta parameter using an exponential function
    t_beta = exp(state[curr]->pos[4]) + beta0;

    // printf("------------------------------\n");
    // printf("HMC proposal transformed\n");
    // printf("%f %f %f %f ", t_gamma, t_i0, t_phi_inv, t_beta);
    // printf("\n");
    // printf("------------------------------\n");

    double logp = 0.0;
    double params[5];

    // Assign transformed parameters to the params array
    params[0] = t_gamma;
    // params[1] = t_s0 / num_pop;
    // params[2] = t_i0 / num_pop;
    params[1] = t_s0;
    params[2] = t_i0;
    params[3] = t_phi_inv;
    params[4] = t_beta;
    
    // Initialize prior_settings array dynamically
    PriorSetting prior_settings[5];
    if(I0_constant == 0){
        prior_settings[0] = initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
        prior_settings[1] = initialize_prior("PRIOR_CONSTANT", num_pop - t_i0, 0);
        prior_settings[2] = initialize_prior("PRIOR_NORMAL", I0_mean, I0_sd);
        prior_settings[3] = initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
        prior_settings[4] = initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
    } else{
        prior_settings[0] = initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
        prior_settings[1] = initialize_prior("PRIOR_CONSTANT", num_pop - t_i0, 0);
        prior_settings[2] = initialize_prior("PRIOR_CONSTANT", I0_value, 0);
        prior_settings[3] = initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
        prior_settings[4] = initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
    }

    // Calculate the log-prior for gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, we can skip its prior contribution
    }
    else{
        logp += log_prior_density(params[0], &prior_settings[0]);
    }

    // printf("R_0: %f\n", num_pop - params[1] - params[2]);
    // Calculate the log-prior for initial susceptible count parameter (t_s0)
    // logp += log_prior_density(num_pop - params[1] - params[2], &prior_settings[1]);

    // // Calculate the log-prior for initial infected count parameter (t_i0)
    // // logp += log_prior_density(params[2], &prior_settings[2]);
    // double dirich_alpha[3];
    // dirich_alpha[0] = g_inpset->dirich_alpha_0;
    // dirich_alpha[1] = g_inpset->dirich_alpha_1;
    // dirich_alpha[2] = g_inpset->dirich_alpha_2;
    // logp += dirichlet_logdensity_3_param(params[1], params[2], dirich_alpha);
    
    // Calculate the log-prior for phi inverse parameter (t_s0)
    logp += log_prior_density(params[1], &prior_settings[1]);

    // Calculate the log-prior for phi inverse parameter (t_i0)
    logp += log_prior_density(params[2], &prior_settings[2]);

    // Calculate the log-prior for phi inverse parameter (t_phi_inv)
    logp += log_prior_density(params[3], &prior_settings[3]);
    
    // Calculate the log-prior for beta parameter (t_beta)
    logp += log_prior_density(params[4], &prior_settings[4]);
    
    return logp;
}

double PriorSIKR_stand_Prevalence_2(double *H){

    double t_gamma, t_s0, t_i0, t_phi_inv, t_beta;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma from the state vector to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-(state[curr]->pos[0] + H[0])));
        }
        else{
            // Transform the phi inverse parameter using an exponential function and a shift
            t_gamma =  exp(state[curr]->pos[0]+ H[0]) + gamma0;
        }
        
    }

    // Transform the initial susceptible count parameter to stay within [0, num_pop]
    t_s0 = 0.0 + (num_pop - 0)/(1 + exp(-(state[curr]->pos[1] + H[1])));
    // printf("S_0: %f\n", t_s0);


    // Transform the initial infected count parameter to stay within [0, num_pop]
    t_i0 = I0_low + (num_pop - I0_low)/(1 + exp(-(state[curr]->pos[2] + H[2])));

    // Transform the phi inverse parameter using an exponential function and a shift
    t_phi_inv = exp(state[curr]->pos[3] + H[3]) + phiInv0;

    // Transform the beta parameter using an exponential function
    t_beta = exp(state[curr]->pos[4] + H[4]) + beta0;

    // printf("------------------------------\n");
    // printf("HMC proposal transformed\n");
    // printf("%f %f %f %f ", t_gamma, t_i0, t_phi_inv, t_beta);
    // printf("\n");
    // printf("------------------------------\n");

    double logp = 0.0;
    double params[5];

    // Assign transformed parameters to the params array
    params[0] = t_gamma;
    // params[1] = t_s0 / num_pop;
    // params[2] = t_i0 / num_pop;
    params[1] = t_s0;
    params[2] = t_i0;
    params[3] = t_phi_inv;
    params[4] = t_beta;
    
    // Initialize prior_settings array dynamically
    PriorSetting prior_settings[5];
    if(I0_constant == 0){
        prior_settings[0] = initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
        prior_settings[1] = initialize_prior("PRIOR_CONSTANT", num_pop - t_i0, 0);
        prior_settings[2] = initialize_prior("PRIOR_NORMAL", I0_mean, I0_sd);
        prior_settings[3] = initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
        prior_settings[4] = initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
    } else{
        prior_settings[0] = initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
        prior_settings[1] = initialize_prior("PRIOR_CONSTANT", num_pop - t_i0, 0);
        prior_settings[2] = initialize_prior("PRIOR_CONSTANT", I0_value, 0);
        prior_settings[3] = initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
        prior_settings[4] = initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
    }

    // Calculate the log-prior for gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, we can skip its prior contribution
    }
    else{
        logp += log_prior_density(params[0], &prior_settings[0]);
    }

    // printf("R_0: %f\n", num_pop - params[1] - params[2]);
    // Calculate the log-prior for initial susceptible count parameter (t_s0)
    // logp += log_prior_density(num_pop - params[1] - params[2], &prior_settings[1]);

    // // Calculate the log-prior for initial infected count parameter (t_i0)
    // // logp += log_prior_density(params[2], &prior_settings[2]);
    // double dirich_alpha[3];
    // dirich_alpha[0] = g_inpset->dirich_alpha_0;
    // dirich_alpha[1] = g_inpset->dirich_alpha_1;
    // dirich_alpha[2] = g_inpset->dirich_alpha_2;
    // logp += dirichlet_logdensity_3_param(params[1], params[2], dirich_alpha);
    
    // Calculate the log-prior for phi inverse parameter (t_s0)
    logp += log_prior_density(params[1], &prior_settings[1]);

    // Calculate the log-prior for phi inverse parameter (t_i0)
    logp += log_prior_density(params[2], &prior_settings[2]);

    // Calculate the log-prior for phi inverse parameter (t_phi_inv)
    logp += log_prior_density(params[3], &prior_settings[3]);
    
    // Calculate the log-prior for beta parameter (t_beta)
    logp += log_prior_density(params[4], &prior_settings[4]);
    
    return logp;
}

double LoglikelihoodSIKR_stand_Prevalence(){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp;

    double t_gamma, t_s0, t_i0, t_phi_inv, t_beta;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma from the state vector to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[0]));
        }
        else{
            // Transform gamma using an exponential function and a shift
            t_gamma = exp(state[curr]->pos[0]) + gamma0;
        }
    }

    // Transform the initial susceptible count parameter to stay within [0, num_pop]
    t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-state[curr]->pos[1]));

    // Transform the initial infected count parameter to stay within [0, num_pop]
    t_i0 = I0_low + (num_pop - I0_low) / (1 + exp(-state[curr]->pos[2]));

    // Calculate the dispersion parameter phi as the inverse of t_phi_inv
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;

    // Calculate the dispersion parameter phi as the inverse of t_phi_inv
    double t_phi = 1 / t_phi_inv;

    // Transform beta using an exponential function and a shift
    t_beta = exp(state[curr]->pos[4]) + beta0;

    int i;

    // Set up sensitivities (sensitivities are required here)
    char *sensitivities[] = {"-nosensi"};

    // Prepare the parameters for the ODE solver
    double parameters_odes[kNumParam - 1];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][kNumParam - 1];

    // Assign gamma and initial infected count to the parameters array
    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_s0;
    parameters_odes[2] = t_i0;
    parameters_odes[3] = t_beta;

    SIR_CVODES_stand_Incidence(
        results_y,       // Output array for state variables
        results_s,       // Output array for sensitivities
        parameters_odes, // Input parameters for the ODEs
        sensitivities    // Sensitivity settings
    );

    double I[kNumObs];
    for (i = 0; i < kNumObs; i++){
        I[i] = results_y[1][i];
    }

    double u_rep;
    int T0_under_rep, T1_under_rep;
    double U0_under_rep, U1_under_rep;
    if (g_inpset->do_under_report == 0){
        T0_under_rep = 0;
        T1_under_rep = tfin_spl_bas;
        U0_under_rep = U1_under_rep = 1.0;
    } else {
        T0_under_rep = g_inpset-> T0_under;
        T1_under_rep = g_inpset-> T1_under;
        U0_under_rep = g_inpset-> U0_under;
        U1_under_rep = g_inpset-> U1_under;
    }
    double logL = 0.0;
    for (i = 0; i < kNumObs; i++){
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        logL += lgamma(g_obs[0][i] / u_rep + t_phi) - lgamma(t_phi) - lgamma(g_obs[0][i] / u_rep + 1.0);
        logL += (g_obs[0][i] / u_rep) * (log(I[i]) - log(I[i] + t_phi));
        logL += t_phi*(log(t_phi) - log(I[i] + t_phi));
    }

    return logL;
}

double LoglikelihoodSIKR_stand_Prevalence_2(double *H){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp;

    double t_gamma, t_s0, t_i0, t_phi_inv, t_beta;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma from the state vector to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-(state[curr]->pos[0] + H[0])));
        }
        else{
            // Transform gamma using an exponential function and a shift
            t_gamma = exp(state[curr]->pos[0] + H[0]) + gamma0;
        }
    }

    // Transform the initial susceptible count parameter to stay within [0, num_pop]
    t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-(state[curr]->pos[1] + H[1])));

    // Transform the initial infected count parameter to stay within [0, num_pop]
    t_i0 = I0_low + (num_pop - I0_low) / (1 + exp(-(state[curr]->pos[2] + H[2])));

    // Calculate the dispersion parameter phi as the inverse of t_phi_inv
    t_phi_inv = exp(state[curr]->pos[3] + H[3]) + phiInv0;

    // Calculate the dispersion parameter phi as the inverse of t_phi_inv
    double t_phi = 1 / t_phi_inv;

    // Transform beta using an exponential function and a shift
    t_beta = exp(state[curr]->pos[4] + H[4]) + beta0;

    int i;

    // Set up sensitivities (sensitivities are required here)
    char *sensitivities[] = {"-nosensi"};

    // Prepare the parameters for the ODE solver
    double parameters_odes[kNumParam - 1];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][kNumParam - 1];

    // Assign gamma and initial infected count to the parameters array
    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_s0;
    parameters_odes[2] = t_i0;
    parameters_odes[3] = t_beta;

    SIR_CVODES_stand_Incidence(
        results_y,       // Output array for state variables
        results_s,       // Output array for sensitivities
        parameters_odes, // Input parameters for the ODEs
        sensitivities    // Sensitivity settings
    );

    double I[kNumObs];
    for (i = 0; i < kNumObs; i++){
        I[i] = results_y[1][i];
    }

    double u_rep;
    int T0_under_rep, T1_under_rep;
    double U0_under_rep, U1_under_rep;
    if (g_inpset->do_under_report == 0){
        T0_under_rep = 0;
        T1_under_rep = tfin_spl_bas;
        U0_under_rep = U1_under_rep = 1.0;
    } else {
        T0_under_rep = g_inpset-> T0_under;
        T1_under_rep = g_inpset-> T1_under;
        U0_under_rep = g_inpset-> U0_under;
        U1_under_rep = g_inpset-> U1_under;
    }
    double logL = 0.0;
    for (i = 0; i < kNumObs; i++){
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        logL += lgamma(g_obs[0][i] / u_rep + t_phi) - lgamma(t_phi) - lgamma(g_obs[0][i] / u_rep + 1.0);
        logL += (g_obs[0][i] / u_rep) * (log(I[i]) - log(I[i] + t_phi));
        logL += t_phi*(log(t_phi) - log(I[i] + t_phi));
    }

    return logL;
}

void GradientSIKR_stand_Prevalence(TypeState *st){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp;

    int i,j;

    double t_gamma, t_s0, t_i0, t_phi_inv, t_beta;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma from the state vector to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-st->pos[0]));
        }
        else{
            // Transform gamma using an exponential function and a shift, [gamma0, Infty] for stability
            t_gamma = exp(st->pos[0]) + gamma0;
        }
    }

    // Transform the initial susceptible count parameter to stay within [0, num_pop]
    t_s0 = 0.0 + (num_pop - 0)/(1 + exp(-st->pos[1]));

    // Transform the initial infected count parameter to stay within [0, num_pop]
    t_i0 = I0_low + (num_pop - I0_low)/(1 + exp(-st->pos[2]));

    // Transform the phi inverse parameter using an exponential function and a shift, [phiInv0, Infty] for stability
    t_phi_inv = exp(st->pos[3]) + phiInv0;

    // Calculate the dispersion parameter phi as the inverse of t_phi_inv
    double t_phi = 1 / t_phi_inv;

    // Transform beta using an exponential function and a shift, [beta0, Infty] for stability
    t_beta = exp(st->pos[4]) + beta0;

    // Jacobian and its derivative for parameter transformations
    double Jacobian[kNumParam];
    double DJacobian[kNumParam];

    // Calculate Jacobian and its derivative for t_gamma
    if (g_inpset->gammaFixed == 1){
        Jacobian[0] = 0.0;
        DJacobian[0] = 0.0;
    }
    else{
        if (g_inpset->gammaBounded == 1){
            double exp_neg_pos0 = exp(-st->pos[0]);
            double denom = (1 + exp_neg_pos0) * (1 + exp_neg_pos0);
            Jacobian[0] = (g_inpset->gammaUpper - g_inpset->gammaLower) * exp_neg_pos0 / denom;
            DJacobian[0] = (1 - exp(st->pos[0])) / (1 + exp(st->pos[0]));
        }
        else{
            Jacobian[0] = t_gamma - gamma0;  // Since t_gamma = exp(st->pos[0]) + gamma0
            DJacobian[0] = 1.0;
        }
    }

    // Calculate Jacobian and its derivative for t_s0
    if (S0_constant == 1){
        Jacobian[1] = 0.0;
        DJacobian[1] = 0.0;
    } else{
        double exp_neg_pos1 = exp(-st->pos[1]);
        double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
        Jacobian[1] = (num_pop - 0) * exp_neg_pos1 / denom1;
        DJacobian[1] = (1 - exp(st->pos[1])) / (1 + exp(st->pos[1]));
    }

    // Calculate Jacobian and its derivative for t_i0
    if(I0_constant == 1){
        Jacobian[2] = 0.0;
        DJacobian[2] = 0.0;
    } else{
        double exp_neg_pos2 = exp(-st->pos[2]);
        double denom2 = (1 + exp_neg_pos2) * (1 + exp_neg_pos2);
        Jacobian[2] = (num_pop - I0_low) * exp_neg_pos2 / denom2;
        DJacobian[2] = (1 - exp(st->pos[2])) / (1 + exp(st->pos[2]));
    }

    // Calculate Jacobian and its derivative for t_phi_inv
    Jacobian[3] = t_phi_inv - phiInv0;  // Since t_phi_inv = exp(st->pos[2]) + phiInv0
    DJacobian[3] = 1.0;
    
    // Calculate Jacobian and its derivative for t_beta
    Jacobian[4] = t_beta - beta0;  // Since t_beta = exp(st->pos[3]) + beta0
    DJacobian[4] = 1.0;

    // Set up sensitivities (sensitivities are required here)
    char *sensitivities[] = {"-sensi", "stg", "t"};

     // Prepare the parameters for the ODE solver
    double parameters_odes[kNumParam - 1];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][kNumParam - 1];


    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_s0;
    parameters_odes[2] = t_i0;
    parameters_odes[3] = t_beta;

    // Solve the SIR model ODEs using the CVODES solver
    SIR_CVODES_stand_Incidence(
        results_y,       // Output array for state variables
        results_s,       // Output array for sensitivities
        parameters_odes, // Input parameters for the ODEs
        sensitivities    // Sensitivity settings
    );

    // Calculate the model-predicted incidence (C_I) from the cumulative counts
    double I[kNumObs];
    for (i = 0; i < kNumObs; i++){
        I[i] = results_y[1][i];
    }

    // Under-reporting adjustment variables
    double u_rep;
    int T0_under_rep, T1_under_rep;
    double U0_under_rep, U1_under_rep;

    // Check if under-reporting is considered
    if (g_inpset->do_under_report == 0){
        T0_under_rep = 0;
        T1_under_rep = tfin_spl_bas;
        U0_under_rep = U1_under_rep = 1.0;
    } else {
        T0_under_rep = g_inpset-> T0_under;
        T1_under_rep = g_inpset-> T1_under;
        U0_under_rep = g_inpset-> U0_under;
        U1_under_rep = g_inpset-> U1_under;
    }

    // Initialize prior_settings array dynamically
    PriorSetting prior_settings[5];
    if(I0_constant == 0){
        prior_settings[0] = initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
        prior_settings[1] = initialize_prior("PRIOR_CONSTANT", num_pop - t_i0, 0);
        prior_settings[2] = initialize_prior("PRIOR_NORMAL", I0_mean, I0_sd);
        prior_settings[3] = initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
        prior_settings[4] = initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
    } else{
        prior_settings[0] = initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
        prior_settings[1] = initialize_prior("PRIOR_CONSTANT", num_pop - t_i0, 0);
        prior_settings[2] = initialize_prior("PRIOR_CONSTANT", I0_value, 0);
        prior_settings[3] = initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
        prior_settings[4] = initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
    }

    // Compute the gradient for t_gamma
    if (g_inpset->gammaFixed != 1) {
        j = 0;  // Index for t_gamma
        st->grad[j] = 0.0;

        for (i = 0; i < kNumObs; i++) {
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            double obs_adj = g_obs[0][i] / u_rep;
            double I_i = I[i];
            double denom = I_i + t_phi;
            double term = (obs_adj / I_i) - (obs_adj + t_phi) / denom;
            st->grad[j] += results_s[1][i][j] * term;
        }

        // Add the derivative of the log-prior density
        double prior_derivative = log_prior_density_derivative(parameters_odes[j], &prior_settings[j]);

        // Adjust the gradient with Jacobian and derivative
        st->grad[j] += prior_derivative;
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }

    if (S0_constant != 1){
        // Compute the gradient for t_s0
        j = 1;  // Index for t_s0
        st->grad[j] = 0.0;

        for (i = 0; i < kNumObs; i++) {
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            double obs_adj = g_obs[0][i] / u_rep;
            double I_i = I[i];
            double denom = I_i + t_phi;

            double term = (obs_adj / I_i) - (obs_adj + t_phi) / denom;
            st->grad[j] += results_s[1][i][j] * term;
        }
        
        // Add the derivative of the log-prior density
        // double prior_derivative = log_prior_density_derivative(num_pop - parameters_odes[j] - t_i0, &prior_settings[j]);
        double dirich_alpha[3];
        dirich_alpha[0] = g_inpset->dirich_alpha_0;
        dirich_alpha[1] = g_inpset->dirich_alpha_1;
        dirich_alpha[2] = g_inpset->dirich_alpha_2;
        double prior_derivative = dirichlet_logdensity_3_param_derivative(t_s0 / num_pop, t_i0 / num_pop, dirich_alpha, 0) / num_pop;


        // Adjust the gradient with Jacobian and derivative
        st->grad[j] += prior_derivative;
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }

    // Compute the gradient for t_i0
    if(I0_constant != 1){
        j = 2;  // Index for t_i0
        st->grad[j] = 0.0;

        for (i = 0; i < kNumObs; i++) {
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            double obs_adj = g_obs[0][i] / u_rep;
            double I_i = I[i];
            double denom = I_i + t_phi;

            double term = (obs_adj / I_i) - (obs_adj + t_phi) / denom;
            st->grad[j] += results_s[1][i][j] * term;
        }
        
        // Add the derivative of the log-prior density
        double prior_derivative = log_prior_density_derivative(parameters_odes[j], &prior_settings[j]);
        // double prior_derivative = dirichlet_logdensity_3_param_derivative(t_s0 / num_pop, t_i0 / num_pop, dirich_alpha, 1) / num_pop;

        // Adjust the gradient with Jacobian and derivative
        st->grad[j] += prior_derivative;
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }

    // Compute the gradient for t_phi_inv
    j = 3;  // Index for t_phi_inv
    st->grad[j] = 0.0;

    for (i = 0; i < kNumObs; i++) {
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        double obs_adj = g_obs[0][i] / u_rep;
        double I_i = I[i];
        double denom = I_i + t_phi;

        st->grad[j] += gsl_sf_psi(obs_adj + t_phi) - gsl_sf_psi(t_phi);
        st->grad[j] += (I_i - obs_adj) / denom;
        st->grad[j] += log(t_phi) - log(denom);
    }

    // Multiply by derivative of t_phi with respect to t_phi_inv
    st->grad[j] *= -1 / (t_phi_inv * t_phi_inv);

    // Add the derivative of the log-prior density
    double prior_derivative = log_prior_density_derivative(t_phi_inv, &prior_settings[j]);

    // Adjust the gradient with Jacobian and derivative
    st->grad[j] += prior_derivative;
    st->grad[j] *= Jacobian[j];
    st->grad[j] += DJacobian[j];

    // Compute the gradient for t_beta
    j = 4;  // Index for t_beta
    st->grad[j] = 0.0;

    for (i = 0; i < kNumObs; i++) {
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        double obs_adj = g_obs[0][i] / u_rep;
        double I_i = I[i];
        double denom = I_i + t_phi;
        double term = (obs_adj / I_i) - (obs_adj + t_phi) / denom;
        st->grad[j] += results_s[1][i][j - 1] * term;
    }

    // Add the derivative of the log-prior density
    prior_derivative = log_prior_density_derivative(t_beta, &prior_settings[j]);

    // Adjust the gradient with Jacobian and derivative
    st->grad[j] += prior_derivative;
    st->grad[j] *= Jacobian[j];
    st->grad[j] += DJacobian[j];

    // // (Optional) Print the gradient for debugging
    
    // printf("------------------------------\n");
    // printf("Gradient: \n");
    // for (i = 0; i < kNumParam; i++){
    //     printf("%f ", st->grad[i]);
    // }
    // printf("\n");
    // printf("------------------------------\n");
}

double LogJacobianSIKR_stand_Prevalence(){
    
    double Jacobian[kNumParam];
    double t_gamma, t_beta, t_phi_inv;
    
    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // Since t_gamma is fixed, its derivative with respect to st->pos[0] is zero
        // To avoid affecting the Jacobian determinant, we set Jacobian[0] = 1.0
        Jacobian[0] = 1.0;  // No change in volume for fixed parameter
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Derivative of t_gamma with respect to st->pos[0]
            double exp_neg_pos0 = exp(-state[curr]->pos[0]);
            double denom = (1 + exp_neg_pos0) * (1 + exp_neg_pos0);
            double J0 = (g_inpset->gammaUpper - g_inpset->gammaLower) * exp_neg_pos0 / denom;
            // Check for NaN or Inf
            if (isnan(J0) || isinf(J0)) {
                J0 = 0.0;
            }
            Jacobian[0] = J0;
        }
        else{
            // Derivative of t_gamma with respect to st->pos[0]
            Jacobian[0] = exp(state[curr]->pos[0]);
        }
    }

    // Transform the initial infected count parameter t_s0
    if (S0_constant == 1){
        Jacobian[1] = 1.0;  // No change in volume for fixed parameter
    } else {
        double t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-state[curr]->pos[1]));
        // Derivative of t_i0 with respect to st->pos[1]
        double exp_neg_pos1 = exp(-state[curr]->pos[1]);
        double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
        double J1 = (num_pop - 0) * exp_neg_pos1 / denom1;
        // Check for NaN or Inf
        if (isnan(J1) || isinf(J1)) {
            J1 = 0.0;
        }
        Jacobian[1] = J1;
    }

    // Transform the initial infected count parameter t_i0
    if(I0_constant == 1){
        Jacobian[2] = 1.0;  // No change in volume for fixed parameter
    } else{
        double t_i0 = I0_low + (num_pop - I0_low) / (1 + exp(-state[curr]->pos[2]));
        // Derivative of t_i0 with respect to st->pos[1]
        double exp_neg_pos2 = exp(-state[curr]->pos[2]);
        double denom1 = (1 + exp_neg_pos2) * (1 + exp_neg_pos2);
        double J2 = (num_pop - I0_low) * exp_neg_pos2 / denom1;
        // Check for NaN or Inf
        if (isnan(J2) || isinf(J2)) {
            J2 = 0.0;
        }
        Jacobian[2] = J2;
    }

    // Transform the phi inverse parameter t_phi_inv
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;
    // Derivative of t_phi_inv with respect to st->pos[2]
    Jacobian[3] = exp(state[curr]->pos[3]);

    // Transform the beta inverse parameter t_beta
    t_beta = exp(state[curr]->pos[4]) + beta0;
    // Derivative of t_beta with respect to st->pos[3]
    Jacobian[4] = exp(state[curr]->pos[4]);

    int i;
    for (i = 5; i < (kNumParam); i++){
        Jacobian[i] = 1;
    }

    double logJ = 1.0;
    for(i = 0; i < 5; i++){
        logJ *= Jacobian[i];
    }
    logJ = log(logJ);

    return logJ;
}

double LogJacobianSIKR_stand_Prevalence_2(double *H){
    
    double Jacobian[kNumParam];
    double t_gamma, t_beta, t_phi_inv;
    
    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // Since t_gamma is fixed, its derivative with respect to st->pos[0] is zero
        // To avoid affecting the Jacobian determinant, we set Jacobian[0] = 1.0
        Jacobian[0] = 1.0;  // No change in volume for fixed parameter
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Derivative of t_gamma with respect to st->pos[0]
            double exp_neg_pos0 = exp(-(state[curr]->pos[0] + H[0]));
            double denom = (1 + exp_neg_pos0) * (1 + exp_neg_pos0);
            double J0 = (g_inpset->gammaUpper - g_inpset->gammaLower) * exp_neg_pos0 / denom;
            // Check for NaN or Inf
            if (isnan(J0) || isinf(J0)) {
                J0 = 0.0;
            }
            Jacobian[0] = J0;
        }
        else{
            // Derivative of t_gamma with respect to st->pos[0]
            Jacobian[0] = exp(state[curr]->pos[0] + H[0]);
        }
    }

    // Transform the initial infected count parameter t_s0
    if (S0_constant == 1){
        Jacobian[1] = 1.0;  // No change in volume for fixed parameter
    } else {
        double t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-(state[curr]->pos[1] + H[1])));
        // Derivative of t_i0 with respect to st->pos[1]
        double exp_neg_pos1 = exp(-(state[curr]->pos[1] + H[1]));
        double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
        double J1 = (num_pop - 0) * exp_neg_pos1 / denom1;
        // Check for NaN or Inf
        if (isnan(J1) || isinf(J1)) {
            J1 = 0.0;
        }
        Jacobian[1] = J1;
    }

    // Transform the initial infected count parameter t_i0
    if(I0_constant == 1){
        Jacobian[2] = 1.0;  // No change in volume for fixed parameter
    } else {
        double t_i0 = I0_low + (num_pop - I0_low) / (1 + exp(-(state[curr]->pos[2] + H[2])));
        // Derivative of t_i0 with respect to st->pos[1]
        double exp_neg_pos1 = exp(-(state[curr]->pos[2] + H[2]));
        double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
        double J2 = (num_pop - I0_low) * exp_neg_pos1 / denom1;
        // Check for NaN or Inf
        if (isnan(J2) || isinf(J2)) {
            J2 = 0.0;
        }
        Jacobian[2] = J2;
    }

    // Transform the phi inverse parameter t_phi_inv
    t_phi_inv = exp(state[curr]->pos[3] + H[3]) + phiInv0;
    // Derivative of t_phi_inv with respect to st->pos[2]
    Jacobian[3] = exp(state[curr]->pos[3] + H[3]);

    // Transform the beta inverse parameter t_beta
    t_beta = exp(state[curr]->pos[4] + H[4]) + beta0;
    // Derivative of t_beta with respect to st->pos[3]
    Jacobian[4] = exp(state[curr]->pos[4] + H[4]);

    int i;
    for (i = 5; i < (kNumParam); i++){
        Jacobian[i] = 1;
    }

    double logJ = 1.0;
    for(i = 0; i < 5; i++){
        logJ *= Jacobian[i];
    }
    logJ = log(logJ);

    return logJ;
}

void JacobianMatSIKR_stand_Prevalence(double *jac_mat) {

    double t_gamma, t_phi_inv, t_beta;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, the Jacobian entry is 1
        jac_mat[0] = 1.0;  // No change in volume for fixed parameter
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Derivative of t_gamma with respect to st->pos[0]
            double exp_neg_pos0 = exp(-state[curr]->pos[0]);
            double denom = (1 + exp_neg_pos0) * (1 + exp_neg_pos0);
            double J0 = (g_inpset->gammaUpper - g_inpset->gammaLower) * exp_neg_pos0 / denom;
            // Check for NaN or Inf
            if (isnan(J0) || isinf(J0) || J0 == 0.0) {
                jac_mat[0] = 0.0;  // Avoid division by zero
            } else {
                jac_mat[0] = 1.0 / J0;
            }
        }
        else{
            // Derivative of t_gamma with respect to st->pos[0]
            double J0 = exp(state[curr]->pos[0]);
            if (isnan(J0) || isinf(J0) || J0 == 0.0) {
                jac_mat[0] = 0.0;
            } else {
                jac_mat[0] = 1.0 / J0;
            }
        }
    }

    // Transform the initial susceptible count parameter t_s0
    if (S0_constant == 1){
        jac_mat[1] = 1.0;  // No change in volume for fixed parameter
    } else {
        double t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-state[curr]->pos[1]));
        // Derivative of t_i0 with respect to st->pos[1]
        double exp_neg_pos1 = exp(-state[curr]->pos[1]);
        double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
        double J1 = (num_pop - 0) * exp_neg_pos1 / denom1;
        // Check for NaN or Inf
        if (isnan(J1) || isinf(J1) || J1 == 0.0) {
            jac_mat[1] = 0.0;  // Avoid division by zero
        } else {
            jac_mat[1] = 1.0 / J1;
        }
    }

    // Transform the initial infected count parameter t_i0
    if (I0_constant == 1){
        jac_mat[2] = 1.0;  // No change in volume for fixed parameter
    } else {
        double t_i0 = I0_low + (num_pop - I0_low) / (1 + exp(-state[curr]->pos[2]));
        // Derivative of t_i0 with respect to st->pos[1]
        double exp_neg_pos1 = exp(-state[curr]->pos[2]);
        double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
        double J2 = (num_pop - I0_low) * exp_neg_pos1 / denom1;
        // Check for NaN or Inf
        if (isnan(J2) || isinf(J2) || J2 == 0.0) {
            jac_mat[2] = 0.0;  // Avoid division by zero
        } else {
            jac_mat[2] = 1.0 / J2;
        }
    }

    // Transform the phi inverse parameter t_phi_inv
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;
    // Derivative of t_phi_inv with respect to st->pos[2]
    double J3 = exp(state[curr]->pos[3]);
    if (isnan(J3) || isinf(J3) || J3 == 0.0) {
        jac_mat[3] = 0.0;
    } else {
        jac_mat[3] = 1.0 / J3;
    }

   // Transform the beta parameter t_beta
    t_beta = exp(state[curr]->pos[4]) + beta0;
    // Derivative of t_beta with respect to st->pos[3]
    double J4 = exp(state[curr]->pos[4]);
    if (isnan(J4) || isinf(J4) || J4 == 0.0) {
        jac_mat[4] = 0.0;
    } else {
        jac_mat[4] = 1.0 / J4;
    }

    int i;
    for (i = 5; i < (kNumParam); i++){
        jac_mat[i] = 1;
    }
}

void TransfSIKR_stand_Prevalence() {
    // Transform gamma parameter (st->pos[0])
    if (g_inpset->gammaFixed == 1) {
        // Gamma is fixed; no transformation needed
        // st->pos[0] remains unchanged
    } else {
        if (g_inpset->gammaBounded == 1) {
            // Bounded gamma; apply inverse logistic transformation
            // Original transformation: t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-st->pos[0]))
            // Inverse transformation: st->pos[0] = log((t_gamma - g_inpset->gammaLower) / (g_inpset->gammaUpper - t_gamma))
            double t_gamma = state[curr]->pos[0];  // Current value of t_gamma
            state[curr]->pos[0] = log((t_gamma - g_inpset->gammaLower) / (g_inpset->gammaUpper - t_gamma));
        } else {
            // Unbounded gamma; apply logarithm transformation
            // Original transformation: t_gamma = exp(st->pos[0]) + gamma0
            // Inverse transformation: st->pos[0] = log(t_gamma - gamma0)
            double t_gamma = state[curr]->pos[0];  // Current value of t_gamma
            state[curr]->pos[0] = log(t_gamma - gamma0);
        }
    }

    // Transform initial susceptible count parameter t_s0 (st->pos[1])
    // Original transformation: t_s0 = 0.0 + (num_pop - 0.0) / (1 + exp(-st->pos[1]))
    // Inverse transformation: st->pos[1] = log((t_s0 - 0.0) / (num_pop - t_s0))
    double t_s0 = state[curr]->pos[1];  // Current value of t_s0
    state[curr]->pos[1] = log((t_s0 - 0.0) / (num_pop - t_s0));

    // Transform initial infected count parameter t_i0 (st->pos[2])
    // Original transformation: t_i0 = 0.0 + (num_pop - 0.0) / (1 + exp(-st->pos[2]))
    // Inverse transformation: st->pos[2] = log((t_i0 - 0.0) / (num_pop - t_i0))
    double t_i0 = state[curr]->pos[2];  // Current value of t_i0
    state[curr]->pos[2] = log((t_i0 - I0_low) / (num_pop - t_i0));

    // Transform phi inverse parameter t_phi_inv (st->pos[3])
    // Original transformation: t_phi_inv = exp(st->pos[3]) + phiInv0
    // Inverse transformation: st->pos[3] = log(t_phi_inv - phiInv0)
    double t_phi_inv = state[curr]->pos[3];  // Current value of t_phi_inv
    state[curr]->pos[3] = log(t_phi_inv - phiInv0);

    // Transform beta inverse parameter t_beta (st->pos[4])
    // Original transformation: t_beta = exp(st->pos[4]) + beta0
    // Inverse transformation: st->pos[4] = log(t_beta - beta0)
    double t_beta = state[curr]->pos[4];  // Current value of t_beta
    state[curr]->pos[4] = log(t_beta - beta0);
}

void InvTransfSIKR_stand_Prevalence() {
    // Transform gamma parameter (st->pos[0])
    if (g_inpset->gammaFixed == 1) {
        // Gamma is fixed; no transformation needed
        // st->pos[0] remains unchanged
    } else {
        if (g_inpset->gammaBounded == 1) {
            // Bounded gamma; apply logistic transformation
            // Transformation: t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-st->pos[0]))
            state[curr]->pos[0] = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1.0 + exp(-state[curr]->pos[0]));
        } else {
            // Unbounded gamma; apply exponential transformation
            // Transformation: t_gamma = exp(st->pos[0]) + gamma0
            state[curr]->pos[0] = exp(state[curr]->pos[0]) + gamma0;
        }
    }

    // Transform initial infected count parameter t_s0 (st->pos[1])
    // Transformation: t_s0 = 0.0 + (num_pop - 0.0) / (1.0 + exp(-st->pos[1]))
    state[curr]->pos[1] = 0.0 + (num_pop - 0.0) / (1.0 + exp(-state[curr]->pos[1]));

    // Transform initial infected count parameter t_i0 (st->pos[2])
    // Transformation: t_i0 = 0.0 + (num_pop - 0.0) / (1.0 + exp(-st->pos[2]))
    state[curr]->pos[2] = I0_low + (num_pop - I0_low) / (1.0 + exp(-state[curr]->pos[2]));

    // Transform phi inverse parameter t_phi_inv (st->pos[3])
    // Transformation: t_phi_inv = exp(st->pos[3]) + phiInv0
    state[curr]->pos[3] = exp(state[curr]->pos[3]) + phiInv0;

    // Transform beta parameter t_beta (st->pos[4])
    // Transformation: t_beta = exp(st->pos[4]) + beta0
    state[curr]->pos[4] = exp(state[curr]->pos[4]) + beta0;
}

void HessianSIKR_stand_Prevalence(){
    
}

void GradientDescentSIKR_stand_Prevalence(){
    int i, j, ii;
    int init_attempts=0, init_max_attempts=100;
    int l_rate_attempts=0, l_rate_attempts_max = 20;
    int save_counter = 0;  // Counter for saved positions
    char filename[100];
    double logPosterior;
    double st_temp[kNumParam];
    PriorSetting prior_settings[kNumParam - 1];
    printf("Starting Gradient Descent\n");
    printf("------------------------------\n");
    printf("LogPrior: %f\n", PriorSIKR_stand_Prevalence());
    printf("LogJacobian: %f\n", LogJacobianSIKR_stand_Prevalence());
    printf("LogLIkelihood: %f\n", LoglikelihoodSIKR_stand_Prevalence());

    logPosterior = PriorSIKR_stand_Prevalence() + LoglikelihoodSIKR_stand_Prevalence() + LogJacobianSIKR_stand_Prevalence();
    // printf("%d%%\n", 100 * 0 / g_inpset->max_iter);
    // printf("------------------------------\n");
    // printf("Log-Posterior: %.8f \n", logPosterior);
    // printf("------------------------------\n");
    double h = 1e-7;
    double H_values[1 + kNumParam][kNumParam];
    double LogPosterior_values[1 + kNumParam];
    double Num_Derivative[kNumParam];

    for (i = 0; i < (1 + kNumParam); i++){
        if (i == 0){
            for (j = 0; j < (kNumParam); j++){
                H_values[i][j] = 0;
            }
        } else{
            for (j = 0; j < (i - 1); j++){
                H_values[i][j] = 0;
            }
            H_values[i][i - 1] = h;
            for (j = i; j < (kNumParam); j++){
                H_values[i][j] = 0;
            }
        }
    }
    
    while (!isfinite(logPosterior)){
        if(init_attempts < init_max_attempts){
            init_attempts += 1;
            printf("Initial point is not computable (We will try another initialization)\n");
            
            // Initialize t_gamma prior setting
            prior_settings[0] =initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
            // // Initialize t_s0 prior setting
            // prior_settings[1] =initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2);
            // Initialize t_i0 prior setting
            if (I0_constant != 1){
                prior_settings[1] =initialize_prior("PRIOR_NORMAL", I0_mean, I0_sd);
            }
            // Initialize t_phi_inv prior setting
            prior_settings[2] =initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
            // Initialize t_beta prior setting
            prior_settings[3] =initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
            

            // // Initialize the global random number generator (if not already initialized)
            // if (g_rng == NULL) {
            //     const gsl_rng_type * T;
            //     gsl_rng_env_setup();
            //     T = gsl_rng_default;
            //     g_rng = gsl_rng_alloc(T);
            //     gsl_rng_set(g_rng, time(NULL)); // Seed with current time or a fixed seed
            // }

            // Sample initial values for gamma, S_0, I_0, phi_inv, and tau from their priors
            // For gamma (state[curr]->pos[0])
            if (g_inpset->gammaFixed == 1) {
                // If gamma is fixed, use the fixed value
                state[curr]->pos[0] = g_inpset->gamma_fixed_value;
            } else {
                state[curr]->pos[0] = sample_from_prior(&prior_settings[0]);  // gamma
            }

            // For I_0 (state[curr]->pos[2])
            if(I0_constant == 0){
                state[curr]->pos[2] = sample_from_prior(&prior_settings[1]);  // I_0
            } else{
                state[curr]->pos[2] = I0_value;
            }

            // For S_0 (state[curr]->pos[1])
            state[curr]->pos[1] = num_pop - state[curr]->pos[2];  // S_0
            // double dirich_alpha[3], dirich_x[3];
            // dirich_alpha[0] = g_inpset->dirich_alpha_0;
            // dirich_alpha[1] = g_inpset->dirich_alpha_1;
            // dirich_alpha[2] = g_inpset->dirich_alpha_2;

            // sample_dirichlet_3_param(dirich_alpha, dirich_x);
            // state[curr]->pos[1] = num_pop * dirich_x[0];
            // state[curr]->pos[2] = num_pop * dirich_x[1];

            // For phi_inv (state[curr]->pos[3])
            state[curr]->pos[3] = sample_from_prior(&prior_settings[2]);  // phi_inv

            // For beta (state[curr]->pos[4])
            state[curr]->pos[4] = sample_from_prior(&prior_settings[3]);  // beta

            printf("New Random Initial Starting Point (%d/%d)\n", init_attempts, init_max_attempts);
            for(i=0;i<kNumParam;i++){
                printf("%f ", state[curr]->pos[i]);
            }
            printf("\n");
            printf("S0, I0, R0: %f, %f, %f\n", state[curr]->pos[1], state[curr]->pos[2], num_pop - state[curr]->pos[1] - state[curr]->pos[2]);

            TransfSIKR_stand_Prevalence();

            // printf("LogPrior: %f\n", PriorSIKR_stand_Incidence());
            // printf("LogJacobian: %f\n", LogJacobianSIKR_stand_Incidence());
            // printf("LogLIkelihood: %f\n", LoglikelihoodSIKR_stand_Incidence());
            logPosterior = PriorSIKR_stand_Prevalence() + LoglikelihoodSIKR_stand_Prevalence() + LogJacobianSIKR_stand_Prevalence();
            } else{
                printf("Initial point is not computable (All initialization attempts failed). Try changing the priors. Exiting...\n");
                exit(EXIT_FAILURE);
            }  
    }
    printf("Initial point is computable.\n"); 

    // Compute the indices at which to save positions during gradient descent
    // We aim to save positions evenly spaced along the iterations
    int *save_points;
    if (n_chains < g_inpset->max_iter){
        save_points = malloc(n_chains * sizeof(int)); // Allocate for n_chains
        int k;
        for (k = 0; k < n_chains; k++) {
            // Calculate the index (iteration number) at which to save the k-th position
            // The '-1' adjusts for zero-based indexing in C (iterations start from 0)
            save_points[k] = (int)(((double)(k + 1) * g_inpset->max_iter) / n_chains) - 1;
            // Ensure the index is within bounds (non-negative)
            if (save_points[k] < 0) {
                save_points[k] = 0;
            }
        }
    } else{
        save_points = malloc(g_inpset->max_iter * sizeof(int)); // Allocate for max_iter
        int k;
        for (k = 0; k < g_inpset->max_iter; k++) {
            // Calculate the index (iteration number) at which to save the k-th position
            // The '-1' adjusts for zero-based indexing in C (iterations start from 0)
            save_points[k] = k;
        }
    }
    
    i=0;
    // Perform gradient descent for the specified number of iterations
    while (i < g_inpset->max_iter) {
        i += 1;
        // Compute the gradient at the current position and store it in state[curr]->grad[]
        GradientSIKR_stand_Prevalence(state[curr]);

        for(j = 0; j < (1 + kNumParam); j++){
            LogPosterior_values[j] = PriorSIKR_stand_Prevalence_2(H_values[j]) + LoglikelihoodSIKR_stand_Prevalence_2(H_values[j]) + LogJacobianSIKR_stand_Prevalence_2(H_values[j]);
        }


        for(j = 1; j < (1 + kNumParam); j++){
            Num_Derivative[j - 1] = (LogPosterior_values[j] - LogPosterior_values[0]) / h;
        }
        printf("Numerical Derivative\n");
        for (j = 0; j < (kNumParam); j++){
            printf("%f ", Num_Derivative[j]);
        }
        printf("\n");
        printf("------------------------------\n");
        printf("gradient: \n");
        for (j = 0; j < (kNumParam); j++){
            printf("%f ", state[curr]->grad[j]);
        }
        printf("\n");
        printf("------------------------------\n");

        // Update each parameter by moving in the direction of the gradient
        for (j = 0; j < kNumParam; j++) {
            st_temp[j] = state[curr]->pos[j];
            state[curr]->pos[j] += g_inpset->learning_rate * state[curr]->grad[j];
        }

        // Limit the gradient step to fulfill S + I + R = N. If R = N - S  - I < 0, we truncate S so R = 0. 
        InvTransfSIKR_stand_Prevalence();
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

        // printf("New point after gradient descent step %d\n", i);
        // for(jj=0;jj<kNumParam;jj++){
        //     printf("%f ", state[curr]->pos[jj]);
        // }
        // printf("\n");

        TransfSIKR_stand_Prevalence();

        logPosterior = PriorSIKR_stand_Prevalence() + LoglikelihoodSIKR_stand_Prevalence() + LogJacobianSIKR_stand_Prevalence();
        if (!isfinite(logPosterior)){
            l_rate_attempts+=1;
            // i = 0;
            g_inpset->learning_rate = g_inpset->learning_rate / 3.16228;
            for (j = 0; j < kNumParam; j++) {
                state[curr]->pos[j] = st_temp[j];
            }
            i += -1;
            printf("Learing rate may be too big for current state. Reducing by a factor of sqrt(10) (current laerning rate / 3.16228).\n");
            printf("New lerning rate is %e\n", g_inpset->learning_rate);
        }
        if(l_rate_attempts>l_rate_attempts_max){
            printf("Max attempts of refining the learning rate reached. Try changing the priors. Exiting...");
            exit(EXIT_FAILURE);
        }

        // Check if the current iteration is one of the save points
        if (save_counter < n_chains && i == save_points[save_counter]) {
            // Save the current position into positions[save_counter]
            // Save the final position from gradient descent to the starting point file
            sprintf(filename, "./output/%s/startingpoint_%s_%d.txt", id_global, id_global, save_counter);
            FILE *fp_startingpoint = fopen(filename, "w");
            if (fp_startingpoint == NULL) {
                fprintf(stderr, "ReadModelSIKR_stand_Prevalence - Failed to open starting point file '%s'\n", filename);
                exit(EXIT_FAILURE);
            }
            InvTransfSIKR_stand_Prevalence();
            for (ii = 0; ii < kNumParam; ii++) {
                fprintf(fp_startingpoint, "%lf\n", state[curr]->pos[ii]);
            }
            fclose(fp_startingpoint);
            save_counter++;
            TransfSIKR_stand_Prevalence();
        }
        if (g_inpset->max_iter <= 10){
            printf("LogPrior: %f\n", PriorSIKR_stand_Prevalence());
            printf("LogJacobian: %f\n", LogJacobianSIKR_stand_Prevalence());
            printf("LogLIkelihood: %f\n", LoglikelihoodSIKR_stand_Prevalence());

            // logPosterior = PriorSIKR_stand_Incidence() + LoglikelihoodSIKR_stand_Incidence() + LogJacobianSIKR_stand_Incidence();
            printf("%d\n", i);
            printf("------------------------------\n");
            printf("Log-Posterior: %.8f \n", logPosterior);
            printf("------------------------------\n");
        } else{
            if (i % (g_inpset->max_iter / 10) == 0 && i > 0){
                // logPosterior = PriorSIKR_stand_Incidence() + LoglikelihoodSIKR_stand_Incidence() + LogJacobianSIKR_stand_Incidence();
                printf("%d%%\n", 100 * i / g_inpset->max_iter);
                printf("------------------------------\n");
                printf("Log-Posterior: %.8f \n", logPosterior);
                printf("------------------------------\n");
            }
        }  
    }
}

int ReadModelSIKR_stand_Prevalence() {

    char filename[100];
    FILE *fptr_obs;
    FILE *fptr_initialpoint;
    int i = 0, ii = 0;
    int symbol;
    double num;

    kNumObs = 0;

    //open file with dataset
    sprintf(filename, "./benchmarks/epidemics/%s_data.txt", g_inpset->data);

    fptr_obs = fopen(filename, "r");
    if((fptr_obs) == NULL){
        fprintf(stderr, "sir_stand_model_synthetic.c - Failed to open the dataset file '%s'", filename);
        exit(EXIT_FAILURE);
    }

    //count newlines in the whole file (that is the number of observations kNumObs)
    while(!feof(fptr_obs)){
        symbol = fgetc(fptr_obs);
        if(symbol == '\n'){ kNumObs++; }
    }

    rewind(fptr_obs); //go to the beginning of the file

    //allocate memory for the data, only one column (new infectious)
    g_obs = MatrixAlloc(1, kNumObs);

    //read in the file and write data in g_obs
    for(ii = 0; ii < kNumObs; ii++){
        if(fscanf(fptr_obs, "%lf", &num) != EOF){
            g_obs[0][ii] = num;
        }
    }

    num_pop = g_obs[0][kNumObs - 1];

    kNumObs = kNumObs - 1;

    num_out_cvodes = kNumObs;
    tfin_spl_bas = kNumObs;

    fprintf(red_lfp, "Number of observations %d\n", kNumObs);

    fprintf(red_lfp, "Population size %f\n", num_pop);

    rewind(fptr_obs); //go to the beginning of the file

    /* assign number of model parameters and functions */
    kNumParam = 4 + 1;
    Prior = PriorSIKR_stand_Prevalence;
    Loglikelihood = LoglikelihoodSIKR_stand_Prevalence;
    Gradient = GradientSIKR_stand_Prevalence;
    Transf = TransfSIKR_stand_Prevalence;
    InvTransf = InvTransfSIKR_stand_Prevalence;
    LogJacobian = LogJacobianSIKR_stand_Prevalence;
    JacobianMat = JacobianMatSIKR_stand_Prevalence;
    Hessian = HessianSIKR_stand_Prevalence;

    /** initialize state of model parameters ***/
    for (i = back2; i <= forw2; i++) {
        SAlloc(&(state[i]), kNumParam);
        SAlloc(&(temporalStatePlusH[i]), kNumParam);
        SAlloc(&(temporalStateMinusH[i]), kNumParam);
    }

    ///////////////////////////////
    //initialize position variables

    //check if initialpoint file exists
    int bool_initialpointfile_exists = -1;
    sprintf(filename, "./output/%s/initialpoint_%s.txt", id_global, id_global);
    fptr_initialpoint = fopen(filename, "r");
    if ((fptr_initialpoint) == NULL) {
        bool_initialpointfile_exists = 0;
    } else{
        bool_initialpointfile_exists = 1;
    }

    fprintf(red_lfp, "bool_initialpointfile_exists is %i\n", bool_initialpointfile_exists);

    //depending on inputfile generate a random starting point or read in from file
    if(bool_initialpointfile_exists == 0){
        // fprintf(stderr, "Error in SIR model input: initialpoints file not found. \n");
        // exit(EXIT_FAILURE);
        // exit(EXIT_FAILURE);
        // Define prior settings for each parameter
        int total_params = 4;
        PriorSetting prior_settings[total_params];

        // Initialize prior_settings array dynamically
        
        // Initialize t_gamma prior setting
        prior_settings[0] =initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
        // // Initialize t_s0 prior setting
        // prior_settings[1] =initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2);
        // Initialize t_i0 prior setting
        if (I0_constant == 0){
            prior_settings[1] =initialize_prior("PRIOR_NORMAL", I0_mean, I0_sd);
        } else {
            prior_settings[1] =initialize_prior("PRIOR_CONSTANT", I0_value, 0);
        }
        // Initialize t_phi_inv prior setting
        prior_settings[2] =initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
        // Initialize t_beta prior setting
        prior_settings[3] =initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
        

        // // Initialize the global random number generator (if not already initialized)
        // if (g_rng == NULL) {
        //     const gsl_rng_type * T;
        //     gsl_rng_env_setup();
        //     T = gsl_rng_default;
        //     g_rng = gsl_rng_alloc(T);
        //     gsl_rng_set(g_rng, time(NULL)); // Seed with current time or a fixed seed
        // }

        // Sample initial values for gamma, S_0, I_0, phi_inv, and tau from their priors
        // For gamma (state[curr]->pos[0])
        if (g_inpset->gammaFixed == 1) {
            // If gamma is fixed, use the fixed value
            state[curr]->pos[0] = g_inpset->gamma_fixed_value;
        } else {
            state[curr]->pos[0] = sample_from_prior(&prior_settings[0]);  // gamma
        }

        // For I_0 (state[curr]->pos[2])
        if(I0_constant == 0){
            state[curr]->pos[2] = sample_from_prior(&prior_settings[1]);  // I_0
        } else{
            state[curr]->pos[2] = I0_value;
        }
        // For S_0 (state[curr]->pos[1])
        state[curr]->pos[1] = num_pop - state[curr]->pos[2];  // S_0
        // double dirich_alpha[3], dirich_x[3];
        // dirich_alpha[0] = g_inpset->dirich_alpha_0;
        // dirich_alpha[1] = g_inpset->dirich_alpha_1;
        // dirich_alpha[2] = g_inpset->dirich_alpha_2;

        // sample_dirichlet_3_param(dirich_alpha, dirich_x);
        // state[curr]->pos[1] = num_pop * dirich_x[0];
        // state[curr]->pos[2] = num_pop * dirich_x[1];

        // For phi_inv (state[curr]->pos[3])
        state[curr]->pos[3] = sample_from_prior(&prior_settings[2]);  // phi_inv

        // For beta (state[curr]->pos[4])
        state[curr]->pos[4] = sample_from_prior(&prior_settings[3]);  // beta
        printf("Random Initial Starting Point\n");
        for(i=0;i<kNumParam;i++){
            printf("%f ", state[curr]->pos[i]);
        }
        printf("\n");
        printf("S0, I0, R0: %f, %f, %f\n", state[curr]->pos[1], state[curr]->pos[2], num_pop - state[curr]->pos[1] - state[curr]->pos[2]);
    } else if(bool_initialpointfile_exists == 1){ //read in from file

        // count number of spaces in the first line
        int numSizeInitialPoint = 0;
        symbol = 1; //initialise symbol
        while(!feof(fptr_initialpoint) && symbol != '\n'){
            symbol = fgetc(fptr_initialpoint);
            if(isspace(symbol)){ numSizeInitialPoint++; }
        }

        fprintf(red_lfp, "numSizeInitialPoint is %i\n", numSizeInitialPoint);

        if(kNumParam != numSizeInitialPoint){
            fprintf(stderr, "sir_stand_model_synthetic.c - kNumParam (%d) is not equal to numSizeInitialPoint (%i) in sir_stand_model_synthetic.c.\n", kNumParam, numSizeInitialPoint);
            fprintf(stderr, "sir_stand_model_synthetic.c - a random initialization point was drawn from the priors.\n");

            // exit(EXIT_FAILURE);
            // Define prior settings for each parameter
            int total_params = 4;
            PriorSetting prior_settings[total_params];

            // Initialize prior_settings array dynamically
            
            // Initialize t_gamma prior setting
            prior_settings[0] =initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
            // // Initialize t_s0 prior setting
            // prior_settings[1] =initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2);
            // Initialize t_i0 prior setting
            if(I0_constant == 0){
                prior_settings[1] = initialize_prior("PRIOR_NORMAL", I0_mean, I0_sd);
            } else {
                prior_settings[1] = initialize_prior("PRIOR_CONSTANT", I0_value, 0);
            }
            // Initialize t_phi_inv prior setting
            prior_settings[2] =initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
            // Initialize t_beta prior setting
            prior_settings[3] =initialize_prior(g_inpset->beta_prior, g_inpset->beta_param_1, g_inpset->beta_param_2);
            

            // // Initialize the global random number generator (if not already initialized)
            // if (g_rng == NULL) {
            //     const gsl_rng_type * T;
            //     gsl_rng_env_setup();
            //     T = gsl_rng_default;
            //     g_rng = gsl_rng_alloc(T);
            //     gsl_rng_set(g_rng, time(NULL)); // Seed with current time or a fixed seed
            // }

            // Sample initial values for gamma, S_0, I_0, phi_inv, and tau from their priors
            // For gamma (state[curr]->pos[0])
            if (g_inpset->gammaFixed == 1) {
                // If gamma is fixed, use the fixed value
                state[curr]->pos[0] = g_inpset->gamma_fixed_value;
            } else {
                state[curr]->pos[0] = sample_from_prior(&prior_settings[0]);  // gamma
            }

            // For I_0 (state[curr]->pos[2])
            state[curr]->pos[2] = sample_from_prior(&prior_settings[1]);  // I_0

            // For S_0 (state[curr]->pos[1])
            state[curr]->pos[1] = num_pop - state[curr]->pos[2];  // S_0
            // double dirich_alpha[3], dirich_x[3];
            // dirich_alpha[0] = g_inpset->dirich_alpha_0;
            // dirich_alpha[1] = g_inpset->dirich_alpha_1;
            // dirich_alpha[2] = g_inpset->dirich_alpha_2;

            // sample_dirichlet_3_param(dirich_alpha, dirich_x);
            // state[curr]->pos[1] = num_pop * dirich_x[0];
            // state[curr]->pos[2] = num_pop * dirich_x[1];

            // For phi_inv (state[curr]->pos[3])
            state[curr]->pos[3] = sample_from_prior(&prior_settings[2]);  // phi_inv

            // For beta (state[curr]->pos[4])
            state[curr]->pos[4] = sample_from_prior(&prior_settings[3]);  // beta
            printf("Random Initial Starting Point\n");
            for(i=0;i<kNumParam;i++){
                printf("%f ", state[curr]->pos[i]);
            }
            printf("\n");
            printf("S0, I0, R0: %f, %f, %f\n", state[curr]->pos[1], state[curr]->pos[2], num_pop - state[curr]->pos[1] - state[curr]->pos[2]);
        } else {
            //go to beginning of file
            rewind(fptr_initialpoint);

            //read initial points in one by one
            for (ii = 0; ii < numSizeInitialPoint; ii++) {
                if (fscanf(fptr_initialpoint, "%lf", &num) != EOF) {
                    state[curr]->pos[ii] = num;
                }
            }

            fclose(fptr_initialpoint);

            fprintf(red_lfp, "initialpoints are ");
            for (i = 0; i < numSizeInitialPoint; i++) {
                fprintf(red_lfp, "%lf ", state[curr]->pos[i]);
            }
            fprintf(red_lfp, "\n");
        }

    } else{
        fprintf(stderr, "sir_stand_model_synthetic.c - If-else not matched for initial position\n");
        exit(EXIT_FAILURE);
    }
    // end initialize position variables

     // Perform gradient descent if enabled
    if (g_inpset->do_grad_desc == 1) {
        // Transform parameters to the real line before gradient descent
        TransfSIKR_stand_Prevalence();

        // Perform gradient descent and get positions
        GradientDescentSIKR_stand_Prevalence(); // positions

        // Transform parameters back from the real line
        InvTransfSIKR_stand_Prevalence();

        // Save the final position from gradient descent to the starting point file
        sprintf(filename, "./output/%s/startingpoint_%s.txt", id_global, id_global);
        FILE *fp_startingpoint = fopen(filename, "w");
        if (fp_startingpoint == NULL) {
            fprintf(stderr, "ReadModelSIKR_stand_Incidence - Failed to open starting point file '%s'\n", filename);
            exit(EXIT_FAILURE);
        }
        for (ii = 0; ii < kNumParam; ii++) {
            fprintf(fp_startingpoint, "%lf\n", state[curr]->pos[ii]);
        }
        fclose(fp_startingpoint);
    }

    /*** initialize momenta vectors from a gaussian distr ****/
    if (state[curr]->mom == NULL) {
        printf("sir_stand_model_synthetic.c - Error in allocating memory (state mom)");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < kNumParam; i++) {
        state[curr]->mom[i] = gsl_ran_gaussian_ziggurat(g_rng, 1.);
    }

    /***** END initialization ******/
    Ham = (TypeHamiltonian *)malloc(sizeof(TypeHamiltonian));

    return 0;
}
