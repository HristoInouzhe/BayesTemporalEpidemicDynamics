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
//#include <lapacke.h>

#include "sir_model_synthetic.h"
#include "sir_seir_common_functions.h"
#include "Globals.h"
#include "utils.h"
#include "sir_functions.h"
#include "read_input.h"

#define M_PI 3.14159265358979323846
#define X0 0.0

#define gamma0 1e-10
#define phiInv0 1e-10
#define n_chains 10

// #define dirich_alpha_0 99997      
// #define dirich_alpha_1 2
// #define dirich_alpha_2 1

double PriorSIKR_Incidence() {

    double t_gamma, t_s0, t_i0, t_phi_inv;

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
    t_i0 = 0.0 + (num_pop - 0)/(1 + exp(-state[curr]->pos[2]));

    // Transform the phi inverse parameter using an exponential function and a shift
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;

    // Transform the tau parameter using an exponential function
    double t_tau = exp(state[curr]->pos[4]);
    int i, j;

    // // Print transformed parameters for debugging or logging purposes
    // printf("------------------------------\n");
    // printf("HMC proposal transformed\n");
    // printf("%f %f %f %f %f ", t_gamma, t_s0, t_i0, t_phi_inv, t_tau);
    // for (i = 5; i < kNumParam; i++){
    //     printf("%f ", state[curr]->pos[i]);
    // }
    // printf("\n");
    // printf("------------------------------\n");

    // Extract spline coefficients from the state vector
    int nf = g_inpset->num_basis_spline;  // Number of spline basis functions
    double spline_coordinates[nf];
    for (i = 0; i < nf; i++){
        spline_coordinates[i] = state[curr]->pos[i + 5];
    }

    double logp = 0.0;  // Initialize log-prior
    double params[4];   // Array of transformed parameters

    // Assign transformed parameters to the params array
    params[0] = t_gamma;
    params[1] = t_s0 / num_pop;
    params[2] = t_i0 / num_pop;
    params[3] = t_phi_inv;

    // Allocate memory for the penalty matrices KK and KKT (transposed KK)
    // The size of KK depends on the difference order

    int penalty_order = g_inpset->spline_pol_degree - 1;  
    int penalty_size = nf - penalty_order;
    double **KK = malloc(penalty_size * sizeof(double *));
    for (i = 0; i < penalty_size; i++) {
        KK[i] = calloc(nf, sizeof(double));
    }

    // Construct the finite difference penalty matrix KK of specified order
    for (i = 0; i < penalty_size; i++) {
        int k;
        double sign = 1.0;
        for (k = 0; k <= penalty_order; k++) {
            // Compute binomial coefficients with alternating signs
            KK[i][i + k] = sign * binomial_coefficient(penalty_order, k);
            sign *= -1.0;
        }
    }

    // Compute KKT (transpose of KK) and result_matrix (KKT * KK)
    double **KKT = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++) {
        KKT[i] = calloc(penalty_size, sizeof(double));
    }
    double **result_matrix = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++) {
        result_matrix[i] = calloc(nf, sizeof(double));
    }

    // Transpose KK to get KKT
    for (i = 0; i < nf; i++) {
        int k;
        for (k = 0; k < penalty_size; k++) {
            KKT[i][k] = KK[k][i];
        }
    }

    // Compute the penalty matrix: result_matrix = KKT * KK
    matrix_matrix_mult(KKT, KK, result_matrix, nf, nf, penalty_size);

    // Compute the product of result_matrix and spline_coordinates
    double result_mv[nf];
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);

    // Initialize prior_settings array dynamically
    PriorSetting prior_settings[3] = {
        // Initialize t_gamma prior setting
        initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2),

        // // Initialize t_s0 prior setting
        // initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2),

        // // Initialize t_i0 prior setting
        // initialize_prior(g_inpset->I0_prior, g_inpset->I0_param_1, g_inpset->I0_param_2),

        // Initialize t_phi_inv prior setting
        initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2),

        // Initialize t_tau prior setting
        initialize_prior(g_inpset->tau_prior, g_inpset->tau_param_1, g_inpset->tau_param_2)
    };

    // Calculate the log-prior for gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, we can skip its prior contribution
    }
    else{
        logp += log_prior_density(params[0], &prior_settings[0]);
    }
    // printf("LogPrior after gamma: %e\n", logp);
    // printf("S_0, I_0, R_0: %f, %f, %f\n", num_pop * params[1], num_pop * params[2], num_pop * (1 - params[1] - params[2]));
    // // Calculate the log-prior for initial susceptible count parameter (t_s0)
    // logp += log_prior_density(num_pop - params[1] - params[2], &prior_settings[1]);

    // // printf("LogPrior after R0: %e\n", logp);
    // // Calculate the log-prior for initial infected count parameter (t_i0)
    // logp += log_prior_density(params[2], &prior_settings[2]);
    // // printf("LogPrior after I0: %e\n", logp);
    double dirich_alpha[3];
    dirich_alpha[0] = g_inpset->dirich_alpha_0;
    dirich_alpha[1] = g_inpset->dirich_alpha_1;
    dirich_alpha[2] = g_inpset->dirich_alpha_2;
    logp += dirichlet_logdensity_3_param(params[1], params[2], dirich_alpha);

    // Calculate the log-prior for phi inverse parameter (t_phi_inv)
    logp += log_prior_density(params[3], &prior_settings[1]);
    // printf("LogPrior after phi_inv: %e\n", logp);

    // Calculate the log-prior for tau parameter
    logp += log_prior_density(t_tau, &prior_settings[2]);
    // printf("LogPrior after tau: %e\n", logp); 

    // Add the penalty term for spline smoothness
    logp += -(1 / (2 * t_tau * t_tau)) * vectors_dot_prod(spline_coordinates, result_mv, nf);
    // printf("LogPrior after spline-smoothness: %e\n", logp);

    #ifdef DEBUG
    // If debugging, print the log-prior value to the log file
    fprintf(lfp, "Prior is %f\n", logp);
    #endif

    /*
    // Optionally print the log-prior to the console
    printf("------------------------------\n");
    printf("Log-prior: %.8f \n", logp);
    printf("------------------------------\n");
    */

    // Free allocated memory
    for (i = 0; i < penalty_size; i++){
        free(KK[i]);
    }
    free(KK);
    for (i = 0; i < nf; i++){
        free(KKT[i]);
    }
    free(KKT);
    for (i = 0; i < nf; i++){
        free(result_matrix[i]);
    }
    free(result_matrix);

    // Return the calculated log-prior value
    return logp;
}

double PriorSIKR_Incidence_2(double *H) {

    double t_gamma, t_s0, t_i0, t_phi_inv;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma from the state vector to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-(state[curr]->pos[0]+ H[0])));
        }
        else{
            // Transform the phi inverse parameter using an exponential function and a shift
            t_gamma =  exp(state[curr]->pos[0] + H[0]) + gamma0;
        }
        
    }

    // Transform the initial susceptible count parameter to stay within [0, num_pop]
    t_s0 = 0.0 + (num_pop - 0)/(1 + exp(-(state[curr]->pos[1] + H[1])));
    // printf("S_0: %f\n", t_s0);

    // Transform the initial infected count parameter to stay within [0, num_pop]
    t_i0 = 0.0 + (num_pop - 0)/(1 + exp(-(state[curr]->pos[2] + H[2])));

    // Transform the phi inverse parameter using an exponential function and a shift
    t_phi_inv = exp((state[curr]->pos[3] + H[3])) + phiInv0;

    // Transform the tau parameter using an exponential function
    double t_tau = exp((state[curr]->pos[4] + H[4]));
    int i, j;

    // Print transformed parameters for debugging or logging purposes
    // printf("------------------------------\n");
    // printf("HMC proposal transformed\n");
    // printf("%f %f %f %f ", t_gamma, t_i0, t_phi_inv, t_tau);
    // for (i = 4; i < kNumParam; i++){
    //     printf("%f ", state[curr]->pos[i]);
    // }
    // printf("\n");
    // printf("------------------------------\n");

    // Extract spline coefficients from the state vector
    int nf = g_inpset->num_basis_spline;  // Number of spline basis functions
    double spline_coordinates[nf];
    for (i = 0; i < nf; i++){
        spline_coordinates[i] = (state[curr]->pos[i + 5] + H[i + 5]);
    }

    double logp = 0.0;  // Initialize log-prior
    double params[4];   // Array of transformed parameters

    // Assign transformed parameters to the params array
    params[0] = t_gamma;
    params[1] = t_s0 / num_pop;
    params[2] = t_i0 / num_pop;
    params[3] = t_phi_inv;

    // Allocate memory for the penalty matrices KK and KKT (transposed KK)
    // The size of KK depends on the difference order

    int penalty_order = g_inpset->spline_pol_degree - 1;  
    int penalty_size = nf - penalty_order;
    double **KK = malloc(penalty_size * sizeof(double *));
    for (i = 0; i < penalty_size; i++) {
        KK[i] = calloc(nf, sizeof(double));
    }

    // Construct the finite difference penalty matrix KK of specified order
    for (i = 0; i < penalty_size; i++) {
        int k;
        double sign = 1.0;
        for (k = 0; k <= penalty_order; k++) {
            // Compute binomial coefficients with alternating signs
            KK[i][i + k] = sign * binomial_coefficient(penalty_order, k);
            sign *= -1.0;
        }
    }

    // Compute KKT (transpose of KK) and result_matrix (KKT * KK)
    double **KKT = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++) {
        KKT[i] = calloc(penalty_size, sizeof(double));
    }
    double **result_matrix = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++) {
        result_matrix[i] = calloc(nf, sizeof(double));
    }

    // Transpose KK to get KKT
    for (i = 0; i < nf; i++) {
        int k;
        for (k = 0; k < penalty_size; k++) {
            KKT[i][k] = KK[k][i];
        }
    }

    // Compute the penalty matrix: result_matrix = KKT * KK
    matrix_matrix_mult(KKT, KK, result_matrix, nf, nf, penalty_size);

    // Compute the product of result_matrix and spline_coordinates
    double result_mv[nf];
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);

    // Initialize prior_settings array dynamically
    PriorSetting prior_settings[3] = {
        // Initialize t_gamma prior setting
        initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2),

        // // Initialize t_s0 prior setting
        // initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2),

        // // Initialize t_i0 prior setting
        // initialize_prior(g_inpset->I0_prior, g_inpset->I0_param_1, g_inpset->I0_param_2),

        // Initialize t_phi_inv prior setting
        initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2),

        // Initialize t_tau prior setting
        initialize_prior(g_inpset->tau_prior, g_inpset->tau_param_1, g_inpset->tau_param_2)
    };

    // Calculate the log-prior for gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, we can skip its prior contribution
    }
    else{
        logp += log_prior_density(params[0], &prior_settings[0]);
    }
    // printf("LogPrior after gamma: %e\n", logp);
    // printf("R_0: %f\n", num_pop - params[1] - params[2]);
    // Calculate the log-prior for initial susceptible count parameter (t_s0)
    // logp += log_prior_density(num_pop - params[1] - params[2], &prior_settings[1]);

    // printf("LogPrior after R0: %e\n", logp);
    // Calculate the log-prior for initial infected count parameter (t_i0)
    // logp += log_prior_density(params[2], &prior_settings[2]);
    // printf("LogPrior after I0: %e\n", logp);
    double dirich_alpha[3];
    dirich_alpha[0] = g_inpset->dirich_alpha_0;
    dirich_alpha[1] = g_inpset->dirich_alpha_1;
    dirich_alpha[2] = g_inpset->dirich_alpha_2;
    logp += dirichlet_logdensity_3_param(params[1], params[2], dirich_alpha);

    // Calculate the log-prior for phi inverse parameter (t_phi_inv)
    logp += log_prior_density(params[3], &prior_settings[1]);
    // printf("LogPrior after phi_inv: %e\n", logp);

    // Calculate the log-prior for tau parameter
    logp += log_prior_density(t_tau, &prior_settings[2]);
    // printf("LogPrior after tau: %e\n", logp); 

    // Add the penalty term for spline smoothness
    logp += -(1 / (2 * t_tau * t_tau)) * vectors_dot_prod(spline_coordinates, result_mv, nf);
    // printf("LogPrior after spline-smoothness: %e\n", logp);

    #ifdef DEBUG
    // If debugging, print the log-prior value to the log file
    fprintf(lfp, "Prior is %f\n", logp);
    #endif

    /*
    // Optionally print the log-prior to the console
    printf("------------------------------\n");
    printf("Log-prior: %.8f \n", logp);
    printf("------------------------------\n");
    */

    // Free allocated memory
    for (i = 0; i < penalty_size; i++){
        free(KK[i]);
    }
    free(KK);
    for (i = 0; i < nf; i++){
        free(KKT[i]);
    }
    free(KKT);
    for (i = 0; i < nf; i++){
        free(result_matrix[i]);
    }
    free(result_matrix);

    // Return the calculated log-prior value
    return logp;
}
double LoglikelihoodSIKR_Incidence() {
    // Choose the number of components based on the model
    int NNEQ = 3 + g_inpset->num_comp;

    double t_gamma, t_s0, t_i0, t_phi_inv;

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
    t_i0 = 0.0 + (num_pop - 0) / (1 + exp(-state[curr]->pos[2]));

    // Transform the phi inverse parameter using an exponential function and a shift
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;

    // Calculate the dispersion parameter phi as the inverse of t_phi_inv
    double t_phi = 1 / t_phi_inv;

    int i;

    // Set up sensitivities (sensitivities are required here)
    char *sensitivities[] = {"-nosensi"};

    // Prepare the parameters for the ODE solver
    int nf = g_inpset->num_basis_spline; // Number of spline basis functions
    double parameters_odes[kNumParam - 2];

    // Assign gamma and initial infected count to the parameters array
    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_s0;
    parameters_odes[2] = t_i0;
    // Extract spline coefficients from the state vector and assign to parameters array
    for (i = 3; i < kNumParam - 2; i++) {
        parameters_odes[i] = state[curr]->pos[i + 2];
    }

    // Calculate the number of knots for the spline basis functions
    int m = (nf - g_inpset->spline_pol_degree - 1) + 2 * (g_inpset->spline_pol_degree + 1);

    // Declare variables for spline basis and knots
    double dx = (tfin_spl_bas - X0) / (nf - g_inpset->spline_pol_degree);
    double knots[m];
    for (i = 0; i < m; i++) {
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }

    // Arrays to store the results from the ODE solver
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][kNumParam - 2];

    // Allocate memory for basis functions if required by the ODE solver
    double basis[nf];
    // double *basis = malloc(nf * sizeof(double));
    // if (basis == NULL) {
    //     fprintf(stderr, "Memory allocation failed for basis\n");
    //     // Free previously allocated memory
    //     for (i = 0; i < NNEQ; i++) {
    //         int j;
    //         for (j = 0; j < kNumObs; j++) {
    //             free(results_s[i][j]);
    //         }
    //         free(results_s[i]);
    //     }
    //     free(results_s);
    //     for (i = 0; i < NNEQ; i++) {
    //         free(results_y[i]);
    //     }
    //     free(results_y);
    //     free(parameters_odes);
    //     free(knots);
    //     exit(EXIT_FAILURE);
    // }

    // Solve the SIKR model ODEs using the CVODES solver
    SIR_CVODES_Incidence(
        results_y,                   // Output array for state variables
        results_s,                   // Output array for sensitivities
        parameters_odes,             // Input parameters for the ODEs
        g_inpset->spline_pol_degree, // Degree of the spline polynomials
        m,                           // Number of knots
        nf,                          // Number of spline basis functions
        knots,                       // Knot vector
        basis,                       // Spline basis functions
        sensitivities                // Sensitivity settings
    );

    // Calculate the model-predicted incidence (C_I) from the cumulative counts
    double C_I[kNumObs];

    // For the first observation, subtract the initial infected count to get new cases
    C_I[0] = results_y[NNEQ - 1][0] - (num_pop - t_s0);
    for (i = 1; i < kNumObs; i++) {
        // Calculate new cases as the difference in cumulative counts between consecutive times
        C_I[i] = results_y[NNEQ - 1][i] - results_y[NNEQ - 1][i - 1];
    }

    // Variables for under-reporting adjustment
    double u_rep;
    int T0_under_rep, T1_under_rep;
    double U0_under_rep, U1_under_rep;

    // Check if under-reporting is considered
    if (g_inpset->do_under_report == 0) {
        // If not, set under-reporting factors to default values (no under-reporting)
        T0_under_rep = 0;
        T1_under_rep = tfin_spl_bas;
        U0_under_rep = U1_under_rep = 1.0;
    } else {
        // If under-reporting is considered, use the specified settings
        T0_under_rep = g_inpset->T0_under;
        T1_under_rep = g_inpset->T1_under;
        U0_under_rep = g_inpset->U0_under;
        U1_under_rep = g_inpset->U1_under;
    }

    double logL = 0.0; // Initialize the log-likelihood

    // Precompute constants outside the loop
    double log_t_phi = log(t_phi);
    double lgamma_t_phi = lgamma(t_phi);

    // Loop over each observation time to compute the log-likelihood
    for (i = 0; i < kNumObs; i++) {
        // Calculate the under-reporting factor at time i
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        double obs_adj = g_obs[0][i];

        double C_I_i = u_rep * C_I[i];
        if (C_I_i <= 0) {
            fprintf(stderr, "Non-positive or zero model incidence at time %d, with abs value %e \n", i, fabs(C_I_i));
            logL = -INFINITY;
            break;
        }

        double C_I_i_plus_phi = C_I_i + t_phi;

        // Precompute log values
        double log_C_I_i = log(C_I_i);
        double log_C_I_i_plus_phi = log(C_I_i_plus_phi);

        // Precompute lgamma values
        double lgamma_obs_adj_plus_phi = lgamma(obs_adj + t_phi);
        double lgamma_obs_adj_plus_one = lgamma(obs_adj + 1.0);

        logL += lgamma_obs_adj_plus_phi - lgamma_t_phi - lgamma_obs_adj_plus_one;
        logL += obs_adj * (log_C_I_i - log_C_I_i_plus_phi);
        logL += t_phi * (log_t_phi - log_C_I_i_plus_phi);
    }

    #ifdef DEBUG
    // If debugging is enabled, print the log-likelihood value to the log file
    fprintf(lfp, "Likelihood is %f\n", logL);
    #endif

    // Return the total log-likelihood
    return logL;
}
double LoglikelihoodSIKR_Incidence_2(double *H) {
    // Choose the number of components based on the model
    int NNEQ = 3 + g_inpset->num_comp;

    double t_gamma, t_s0, t_i0, t_phi_inv;

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
            t_gamma = exp((state[curr]->pos[0] + H[0])) + gamma0;
        }
    }

    // Transform the initial susceptible count parameter to stay within [0, num_pop]
    t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-(state[curr]->pos[1] + H[1])));

    // Transform the initial infected count parameter to stay within [0, num_pop]
    t_i0 = 0.0 + (num_pop - 0) / (1 + exp(-(state[curr]->pos[2] + H[2])));

    // Transform the phi inverse parameter using an exponential function and a shift
    t_phi_inv = exp((state[curr]->pos[3] + H[3])) + phiInv0;

    // Calculate the dispersion parameter phi as the inverse of t_phi_inv
    double t_phi = 1 / t_phi_inv;

    int i;

    // Set up sensitivities (sensitivities are required here)
    char *sensitivities[] = {"-nosensi"};

    // Prepare the parameters for the ODE solver
    int nf = g_inpset->num_basis_spline; // Number of spline basis functions
    double parameters_odes[kNumParam - 2];

    // Assign gamma and initial infected count to the parameters array
    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_s0;
    parameters_odes[2] = t_i0;
    // Extract spline coefficients from the state vector and assign to parameters array
    for (i = 3; i < kNumParam - 2; i++) {
        parameters_odes[i] = (state[curr]->pos[i + 2] + H[i + 2]);
    }

    // Calculate the number of knots for the spline basis functions
    int m = (nf - g_inpset->spline_pol_degree - 1) + 2 * (g_inpset->spline_pol_degree + 1);

    // Declare variables for spline basis and knots
    double dx = (tfin_spl_bas - X0) / (nf - g_inpset->spline_pol_degree);
    double knots[m];
    for (i = 0; i < m; i++) {
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }

    // Arrays to store the results from the ODE solver
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][kNumParam - 2];

    // Allocate memory for basis functions if required by the ODE solver
    double basis[nf];
    // double *basis = malloc(nf * sizeof(double));
    // if (basis == NULL) {
    //     fprintf(stderr, "Memory allocation failed for basis\n");
    //     // Free previously allocated memory
    //     for (i = 0; i < NNEQ; i++) {
    //         int j;
    //         for (j = 0; j < kNumObs; j++) {
    //             free(results_s[i][j]);
    //         }
    //         free(results_s[i]);
    //     }
    //     free(results_s);
    //     for (i = 0; i < NNEQ; i++) {
    //         free(results_y[i]);
    //     }
    //     free(results_y);
    //     free(parameters_odes);
    //     free(knots);
    //     exit(EXIT_FAILURE);
    // }

    // Solve the SIKR model ODEs using the CVODES solver
    SIR_CVODES_Incidence(
        results_y,                   // Output array for state variables
        results_s,                   // Output array for sensitivities
        parameters_odes,             // Input parameters for the ODEs
        g_inpset->spline_pol_degree, // Degree of the spline polynomials
        m,                           // Number of knots
        nf,                          // Number of spline basis functions
        knots,                       // Knot vector
        basis,                       // Spline basis functions
        sensitivities                // Sensitivity settings
    );

    // Calculate the model-predicted incidence (C_I) from the cumulative counts
    double C_I[kNumObs];

    // For the first observation, subtract the initial infected count to get new cases
    C_I[0] = results_y[NNEQ - 1][0] - (num_pop - t_s0);
    for (i = 1; i < kNumObs; i++) {
        // Calculate new cases as the difference in cumulative counts between consecutive times
        C_I[i] = results_y[NNEQ - 1][i] - results_y[NNEQ - 1][i - 1];
    }

    // Variables for under-reporting adjustment
    double u_rep;
    int T0_under_rep, T1_under_rep;
    double U0_under_rep, U1_under_rep;

    // Check if under-reporting is considered
    if (g_inpset->do_under_report == 0) {
        // If not, set under-reporting factors to default values (no under-reporting)
        T0_under_rep = 0;
        T1_under_rep = tfin_spl_bas;
        U0_under_rep = U1_under_rep = 1.0;
    } else {
        // If under-reporting is considered, use the specified settings
        T0_under_rep = g_inpset->T0_under;
        T1_under_rep = g_inpset->T1_under;
        U0_under_rep = g_inpset->U0_under;
        U1_under_rep = g_inpset->U1_under;
    }

    double logL = 0.0; // Initialize the log-likelihood

    // Precompute constants outside the loop
    double log_t_phi = log(t_phi);
    double lgamma_t_phi = lgamma(t_phi);

    // Loop over each observation time to compute the log-likelihood
    for (i = 0; i < kNumObs; i++) {
        // Calculate the under-reporting factor at time i
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        double obs_adj = g_obs[0][i];

        double C_I_i = u_rep * C_I[i];
        if (C_I_i <= 0) {
            fprintf(stderr, "Non-positive or zero model incidence at time %d, with abs value %e \n", i, fabs(C_I_i));
            logL = -INFINITY;
            break;
        }

        double C_I_i_plus_phi = C_I_i + t_phi;

        // Precompute log values
        double log_C_I_i = log(C_I_i);
        double log_C_I_i_plus_phi = log(C_I_i_plus_phi);

        // Precompute lgamma values
        double lgamma_obs_adj_plus_phi = lgamma(obs_adj + t_phi);
        double lgamma_obs_adj_plus_one = lgamma(obs_adj + 1.0);

        logL += lgamma_obs_adj_plus_phi - lgamma_t_phi - lgamma_obs_adj_plus_one;
        logL += obs_adj * (log_C_I_i - log_C_I_i_plus_phi);
        logL += t_phi * (log_t_phi - log_C_I_i_plus_phi);
    }

    #ifdef DEBUG
    // If debugging is enabled, print the log-likelihood value to the log file
    fprintf(lfp, "Likelihood is %f\n", logL);
    #endif

    // Return the total log-likelihood
    return logL;
}

void GradientSIKR_Incidence(TypeState *st) {
    // Choose the number of components based on the model
    int NNEQ = 3 + g_inpset->num_comp;

    int i, j;
    int nf = g_inpset->num_basis_spline;  // Number of spline basis functions

    // Extract spline coefficients from the state vector
    double spline_coordinates[nf];
    for (i = 0; i < nf; i++) {
        spline_coordinates[i] = st->pos[i + 5];
    }

    // Transform the gamma parameter
    double t_gamma;
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
    double t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-st->pos[1]));

    // Transform the initial infected count parameter to stay within [0, num_pop]
    double t_i0 = 0.0 + (num_pop - 0) / (1 + exp(-st->pos[2]));

    // Transform the phi inverse parameter using an exponential function and a shift, [phiInv0, Infty] for stability
    double t_phi_inv = exp(st->pos[3]) + phiInv0;

    // Calculate the dispersion parameter phi as the inverse of t_phi_inv
    double t_phi = 1 / t_phi_inv;

    // Transform the tau parameter using an exponential function [0, Infty]
    double t_tau = exp(st->pos[4]);

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
    double exp_neg_pos1 = exp(-st->pos[1]);
    double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
    Jacobian[1] = (num_pop - 1) * exp_neg_pos1 / denom1;
    DJacobian[1] = (1 - exp(st->pos[1])) / (1 + exp(st->pos[1]));

    // Calculate Jacobian and its derivative for t_i0
    exp_neg_pos1 = exp(-st->pos[2]);
    denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
    Jacobian[2] = (num_pop - 1) * exp_neg_pos1 / denom1;
    DJacobian[2] = (1 - exp(st->pos[2])) / (1 + exp(st->pos[2]));

    // Calculate Jacobian and its derivative for t_phi_inv
    Jacobian[3] = t_phi_inv - phiInv0;  // Since t_phi_inv = exp(st->pos[2]) + phiInv0
    DJacobian[3] = 1.0;

    // Calculate Jacobian and its derivative for t_tau
    Jacobian[4] = t_tau;  // Since t_tau = exp(st->pos[3])
    DJacobian[4] = 1.0;

    // For spline coefficients, Jacobian is 1 and derivative is 0
    for (i = 5; i < kNumParam; i++) {
        Jacobian[i] = 1.0;
        DJacobian[i] = 0.0;
    }

    // Set up sensitivities (sensitivities are required here)
    char *sensitivities[] = {"-sensi", "stg", "t"};

    // Prepare the parameters for the ODE solver
    double parameters_odes[kNumParam - 2];
    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_s0;
    parameters_odes[2] = t_i0;
    // printf("%f, %f, ", parameters_odes[0], parameters_odes[1]);
    for (i = 3; i < kNumParam - 2; i++) {
        parameters_odes[i] = st->pos[i + 2];
        // printf("%f, ", parameters_odes[i]);
    }
    // printf("\n");

    // Calculate the number of knots for the spline basis functions
    int m = (nf - g_inpset->spline_pol_degree - 1) + 2 * (g_inpset->spline_pol_degree + 1);

    // Declare variables for spline basis and knots
    double dx = (tfin_spl_bas - X0) / (nf - g_inpset->spline_pol_degree);
    double knots[m];
    for (i = 0; i < m; i++) {
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }

    // Arrays to store the results from the ODE solver
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][kNumParam - 2];

    // Allocate memory for basis functions if required by the ODE solver
    double basis[nf];
    // double *basis = malloc(nf * sizeof(double));
    // if (basis == NULL) {
    //     fprintf(stderr, "Memory allocation failed for basis\n");
    //     // Free previously allocated memory
    //     for (i = 0; i < NNEQ; i++) {
    //         int j;
    //         for (j = 0; j < kNumObs; j++) {
    //             free(results_s[i][j]);
    //         }
    //         free(results_s[i]);
    //     }
    //     free(results_s);
    //     for (i = 0; i < NNEQ; i++) {
    //         free(results_y[i]);
    //     }
    //     free(results_y);
    //     free(parameters_odes);
    //     free(knots);
    //     exit(EXIT_FAILURE);
    // }

    // Solve the SIR model ODEs using the CVODES solver
    SIR_CVODES_Incidence(
        results_y,                   // Output array for state variables
        results_s,                   // Output array for sensitivities
        parameters_odes,             // Input parameters for the ODEs
        g_inpset->spline_pol_degree, // Degree of the spline polynomials
        m,                           // Number of knots
        nf,                          // Number of spline basis functions
        knots,                       // Knot vector
        basis,                        // Spline basis functions
        sensitivities                // Sensitivity settings
    );

    // printf("-------------------------------------------\n");
    // printf("Delta sensitivity:\n");
    // printf("iteration 0, value %f | ", results_s[NNEQ - 1][0][kNumParam - 3]);
    // for (i = 1; i < kNumObs; i++) {
    //     printf("iteration %d, value %f | ", i, results_s[NNEQ - 1][i][kNumParam - 3] - results_s[NNEQ - 1][i - 1][kNumParam - 3]);
    // }
    // printf("\n");
    // Calculate the model-predicted incidence (C_I) from the cumulative counts
    double C_I[kNumObs];
    C_I[0] = results_y[NNEQ - 1][0] - (num_pop - t_s0);
    // printf("-------------------------------------------\n");
    // printf("Commulative incidence:\n");
    // printf("Iteration 0, value %f | ", C_I[0]);
    for (i = 1; i < kNumObs; i++) {
        C_I[i] = results_y[NNEQ - 1][i] - results_y[NNEQ - 1][i - 1];
        // printf("iteration %d, value %f | ", i, C_I[i]);
    }
    // printf("\n");

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
        T0_under_rep = g_inpset->T0_under;
        T1_under_rep = g_inpset->T1_under;
        U0_under_rep = g_inpset->U0_under;
        U1_under_rep = g_inpset->U1_under;
    }

    // Initialize the gradient vector
    for (i = 0; i < kNumParam; i++) {
        st->grad[i] = 0.0;
    }

    // Initialize prior_settings array dynamically
    PriorSetting prior_settings[3] = {
        // Initialize t_gamma prior setting
        initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2),

        // // Initialize t_s0 prior setting
        // initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2),

        // // Initialize t_i0 prior setting
        // initialize_prior(g_inpset->I0_prior, g_inpset->I0_param_1, g_inpset->I0_param_2),

        // Initialize t_phi_inv prior setting
        initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2),

        // Initialize t_tau prior setting
        initialize_prior(g_inpset->tau_prior, g_inpset->tau_param_1, g_inpset->tau_param_2)
    };

    // printf("%f, %f\n", prior_settings[0].param1, prior_settings[0].param2);

    // Compute the gradient for t_gamma
    if (g_inpset->gammaFixed != 1) {
        j = 0;  // Index for t_gamma
        st->grad[j] = 0.0;

        for (i = 0; i < kNumObs; i++) {
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            double obs_adj = g_obs[0][i] * u_rep;
            double C_I_i = u_rep * C_I[i];
            double denom = C_I_i + t_phi;
            // printf("%f, ", C_I_i);

            double term = (obs_adj / C_I_i) - (obs_adj + t_phi * u_rep) / denom;
            // printf("%f", term);
            double delta_C_I = (i == 0) ? results_s[NNEQ - 1][0][j] : results_s[NNEQ - 1][i][j] - results_s[NNEQ - 1][i - 1][j]; // condition ? expression_if_true : expression_if_false;
            st->grad[j] += delta_C_I * term;
        }
        // printf("\n");

        // printf("%f\n", st->grad[j]);

        // Add the derivative of the log-prior density
        double prior_derivative = log_prior_density_derivative(parameters_odes[j], &prior_settings[j]);

        // Adjust the gradient with Jacobian and derivative
        st->grad[j] += prior_derivative;
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }

    // Compute the gradient for t_s0
    j = 1;  // Index for t_s0
    st->grad[j] = 0.0;

    for (i = 0; i < kNumObs; i++) {
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        double obs_adj = g_obs[0][i] * u_rep;
        double C_I_i = u_rep * C_I[i];
        double denom = C_I_i + t_phi;

        double term = (obs_adj / C_I_i) - (obs_adj + t_phi * u_rep) / denom;
        double delta_C_I;
        if (i == 0) {
            delta_C_I = results_s[NNEQ - 1][0][j] + 1.0;  // Adjust for initial infected count
        } else {
            delta_C_I = results_s[NNEQ - 1][i][j] - results_s[NNEQ - 1][i - 1][j];
        }
        st->grad[j] += delta_C_I * term;
    }

    // // Add the derivative of the log-prior density
    // double prior_derivative = log_prior_density_derivative(num_pop - parameters_odes[j]  - t_i0, &prior_settings[j]);
    // Add the derivative of the log-prior density
    double dirich_alpha[3];
    dirich_alpha[0] = g_inpset->dirich_alpha_0;
    dirich_alpha[1] = g_inpset->dirich_alpha_1;
    dirich_alpha[2] = g_inpset->dirich_alpha_2;
    double prior_derivative = dirichlet_logdensity_3_param_derivative(t_s0 / num_pop, t_i0 / num_pop, dirich_alpha, 0) / num_pop;

    // Adjust the gradient with Jacobian and derivative
    st->grad[j] += prior_derivative;
    st->grad[j] *= Jacobian[j];
    st->grad[j] += DJacobian[j];

    // Compute the gradient for t_i0
    j = 2;  // Index for t_i0
    st->grad[j] = 0.0;

    for (i = 0; i < kNumObs; i++) {
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        double obs_adj = g_obs[0][i] * u_rep;
        double C_I_i = u_rep * C_I[i];
        double denom = C_I_i + t_phi;

        double term = (obs_adj / C_I_i) - (obs_adj + t_phi * u_rep) / denom;
        double delta_C_I;
        if (i == 0) {
            delta_C_I = results_s[NNEQ - 1][0][j];  
        } else {
            delta_C_I = results_s[NNEQ - 1][i][j] - results_s[NNEQ - 1][i - 1][j];
        }
        st->grad[j] += delta_C_I * term;
    }

    // Add the derivative of the log-prior density
    // prior_derivative = log_prior_density_derivative(parameters_odes[j], &prior_settings[j]);
    prior_derivative = dirichlet_logdensity_3_param_derivative(t_s0 / num_pop, t_i0 / num_pop, dirich_alpha, 1) / num_pop;

    // Adjust the gradient with Jacobian and derivative
    st->grad[j] += prior_derivative;
    st->grad[j] *= Jacobian[j];
    st->grad[j] += DJacobian[j];

    // Compute the gradient for t_phi_inv
    j = 3;  // Index for t_phi_inv
    st->grad[j] = 0.0;

    for (i = 0; i < kNumObs; i++) {
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        double obs_adj = g_obs[0][i];
        double C_I_i = u_rep * C_I[i];
        double denom = C_I_i + t_phi;

        st->grad[j] += gsl_sf_psi(obs_adj + t_phi) - gsl_sf_psi(t_phi);
        st->grad[j] += (C_I_i - obs_adj) / denom;
        st->grad[j] += log(t_phi) - log(denom);
    }

    // Multiply by derivative of t_phi with respect to t_phi_inv
    st->grad[j] *= -1 / (t_phi_inv * t_phi_inv);

    // Add the derivative of the log-prior density
    prior_derivative = log_prior_density_derivative(t_phi_inv, &prior_settings[j - 2]);

    // Adjust the gradient with Jacobian and derivative
    st->grad[j] += prior_derivative;
    st->grad[j] *= Jacobian[j];
    st->grad[j] += DJacobian[j];

    // Compute the gradient for t_tau
    j = 4;  // Index for t_tau
    // Penalty matrix computations (similar to prior function)
    int penalty_order = g_inpset->spline_pol_degree - 1;  
    int penalty_size = nf - penalty_order;
    double **KK = malloc(penalty_size * sizeof(double *));
    for (i = 0; i < penalty_size; i++) {
        KK[i] = calloc(nf, sizeof(double));
    }

    // Construct the finite difference penalty matrix KK of specified order
    for (i = 0; i < penalty_size; i++) {
        int k;
        double sign = 1.0;
        for (k = 0; k <= penalty_order; k++) {
            // Compute binomial coefficients with alternating signs
            KK[i][i + k] = sign * binomial_coefficient(penalty_order, k);
            sign *= -1.0;
        }
    }

    // Compute KKT (transpose of KK) and result_matrix (KKT * KK)
    double **KKT = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++) {
        KKT[i] = calloc(penalty_size, sizeof(double));
    }
    double **result_matrix = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++) {
        result_matrix[i] = calloc(nf, sizeof(double));
    }

    // Transpose KK to get KKT
    for (i = 0; i < nf; i++) {
        int k;
        for (k = 0; k < penalty_size; k++) {
            KKT[i][k] = KK[k][i];
        }
    }
    

    // Compute result_matrix = KKT * KK
    matrix_matrix_mult(KKT, KK, result_matrix, nf, nf, penalty_size);

    // Compute result_mv = result_matrix * spline_coordinates
    double result_mv[nf];
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);

    // Compute the gradient for t_tau
    st->grad[j] = 0.0;
    double prior_derivative_tau = log_prior_density_derivative(t_tau, &prior_settings[j - 2]);

    // Gradient contribution from the spline penalty
    double penalty_grad = vectors_dot_prod(spline_coordinates, result_mv, nf) / (t_tau * t_tau * t_tau);

    // Total gradient for t_tau
    st->grad[j] += prior_derivative_tau + penalty_grad;

    // Adjust the gradient with Jacobian and derivative
    st->grad[j] *= Jacobian[j];
    st->grad[j] += DJacobian[j];

    // Compute the gradient for spline coefficients
    for (j = 5; j < kNumParam; j++) {
        st->grad[j] = 0.0;

        for (i = 0; i < kNumObs; i++) {
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            double obs_adj = g_obs[0][i] * u_rep;
            double C_I_i = u_rep * C_I[i];
            double denom = C_I_i + t_phi;

            double term = (obs_adj / C_I_i) - (obs_adj + t_phi * u_rep) / denom;
            double delta_C_I;
            if (i == 0) {
                delta_C_I = results_s[NNEQ - 1][0][j - 2];
            } else {
                delta_C_I = results_s[NNEQ - 1][i][j - 2] - results_s[NNEQ - 1][i - 1][j - 2];
            }
            st->grad[j] += delta_C_I * term;
            // if (j == (kNumParam - 1)){
            //     printf("Iteration %d, term %e, delta_C_IValue %e | ", i, term, delta_C_I);
            // }     
        }
        // if(j == (kNumParam - 1)){
        //     printf("Term from C_I: %e\n", st->grad[j]);
        // }

        // Add the gradient contribution from the spline penalty
        // printf("result_mv[%d]: %f\n", j - 5,  result_mv[j - 5]);
        st->grad[j] += - (1 / (t_tau * t_tau)) * result_mv[j - 5];

        // if(j == (kNumParam - 1)){
        //     printf("Term from smoothness: %e\n", st->grad[j]);
        // }

        // Adjust the gradient with Jacobian and derivative
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
        // if(j == (kNumParam - 1)){
        //     printf("Term from jacobian, djacobian: %e, %e\n", Jacobian[j], DJacobian[j]);
        // }     
    }

    // Free allocated memory
    for (i = 0; i < penalty_size; i++) {
        free(KK[i]);
    }
    free(KK);
    for (i = 0; i < nf; i++) {
        free(KKT[i]);
    }
    free(KKT);
    for (i = 0; i < nf; i++) {
        free(result_matrix[i]);
    }
    free(result_matrix);

    // // (Optional) Print the gradient for debugging
    
    // printf("------------------------------\n");
    // printf("Gradient: \n");
    // for (i = 0; i < kNumParam; i++){
    //     printf("%f ", st->grad[i]);
    // }
    // printf("\n");
    // printf("------------------------------\n");
    
}

double LogJacobianSIKR_Incidence() {
    // This function computes the log of the absolute value of the determinant of the Jacobian matrix
    // for the transformation from the parameters in the state vector st->pos[] to the transformed parameters.
    // It accounts for different options for t_gamma (fixed, bounded, unbounded) and other parameters.

    double Jacobian[kNumParam];
    double t_gamma, t_phi_inv;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
        // Since t_gamma is fixed, its derivative with respect to st->pos[0] is zero
        // To avoid affecting the Jacobian determinant, we set Jacobian[0] = 1.0
        Jacobian[0] = 1.0;  // No change in volume for fixed parameter
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma using a logistic function to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[0]));
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
            // Transform gamma using an exponential function and a shift
            t_gamma = exp(state[curr]->pos[0]) + gamma0;
            // Derivative of t_gamma with respect to st->pos[0]
            Jacobian[0] = exp(state[curr]->pos[0]);
        }
    }

    // Transform the initial infected count parameter t_s0
    double t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-state[curr]->pos[1]));
    // Derivative of t_s0 with respect to st->pos[1]
    double exp_neg_pos1 = exp(-state[curr]->pos[1]);
    double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
    double J1 = (num_pop - 0) * exp_neg_pos1 / denom1;
    // Check for NaN or Inf
    if (isnan(J1) || isinf(J1)) {
        J1 = 0.0;
    }
    Jacobian[1] = J1;

    // Transform the initial infected count parameter t_i0
    double t_i0 = 0.0 + (num_pop - 0) / (1 + exp(-state[curr]->pos[2]));
    // Derivative of t_i0 with respect to st->pos[1]
    exp_neg_pos1 = exp(-state[curr]->pos[2]);
    denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
    double J2 = (num_pop - 0) * exp_neg_pos1 / denom1;
    // Check for NaN or Inf
    if (isnan(J2) || isinf(J2)) {
        J2 = 0.0;
    }
    Jacobian[2] = J2;

    // Transform the phi inverse parameter t_phi_inv
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;
    // Derivative of t_phi_inv with respect to st->pos[2]
    Jacobian[3] = exp(state[curr]->pos[3]);

    // Transform the tau parameter t_tau
    double t_tau = exp(state[curr]->pos[4]);
    // Derivative of t_tau with respect to st->pos[3]
    Jacobian[4] = t_tau;

    // For spline coefficients, the transformation is identity, so the Jacobian is 1
    int i;
    for (i = 5; i < kNumParam; i++) {
        Jacobian[i] = 1.0;
    }

    // Compute the log of the absolute value of the determinant of the Jacobian matrix
    // Since the parameters are transformed independently, the total Jacobian determinant is the product of the individual derivatives
    double logJ = 0.0;
    // Multiply Jacobian contributions for transformed parameters
    for (i = 0; i < 5; i++) {
        // Since log(1) = 0, we can skip Jacobian[i] == 1
        if (Jacobian[i] != 1.0) {
            if (Jacobian[i] <= 0.0) {
                // Logarithm of non-positive number is undefined
                // Return -INFINITY to indicate that the Jacobian determinant is zero or invalid
                return -INFINITY;
            }
            logJ += log(Jacobian[i]);
        }
    }

    // Return the total log-Jacobian
    return logJ;
}

double LogJacobianSIKR_Incidence_2(double *H) {
    // This function computes the log of the absolute value of the determinant of the Jacobian matrix
    // for the transformation from the parameters in the state vector st->pos[] to the transformed parameters.
    // It accounts for different options for t_gamma (fixed, bounded, unbounded) and other parameters.

    double Jacobian[kNumParam];
    double t_gamma, t_phi_inv;

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, use the fixed value
        t_gamma = g_inpset->gamma_fixed_value;
        // Since t_gamma is fixed, its derivative with respect to st->pos[0] is zero
        // To avoid affecting the Jacobian determinant, we set Jacobian[0] = 1.0
        Jacobian[0] = 1.0;  // No change in volume for fixed parameter
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma using a logistic function to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-(state[curr]->pos[0] + H[0])));
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
            // Transform gamma using an exponential function and a shift
            t_gamma = exp((state[curr]->pos[0] + H[0])) + gamma0;
            // Derivative of t_gamma with respect to st->pos[0]
            Jacobian[0] = exp((state[curr]->pos[0] + H[0]));
        }
    }

    // // Transform the initial infected count parameter t_s0
    // double t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-(state[curr]->pos[1] + H[1])));
    // Derivative of t_s0 with respect to st->pos[1]
    double exp_neg_pos1 = exp(-(state[curr]->pos[1] + H[1]));
    double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
    double J1 = (num_pop - 0) * exp_neg_pos1 / denom1;
    // Check for NaN or Inf
    if (isnan(J1) || isinf(J1)) {
        J1 = 0.0;
    }
    Jacobian[1] = J1;

    // Transform the initial infected count parameter t_i0
    // double t_i0 = 0.0 + (num_pop - 0) / (1 + exp(-(state[curr]->pos[2] + H[2])));
    // Derivative of t_i0 with respect to st->pos[1]
    exp_neg_pos1 = exp(-(state[curr]->pos[2] + H[2]));
    denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
    double J2 = (num_pop - 0) * exp_neg_pos1 / denom1;
    // Check for NaN or Inf
    if (isnan(J2) || isinf(J2)) {
        J2 = 0.0;
    }
    Jacobian[2] = J2;

    // Transform the phi inverse parameter t_phi_inv
    t_phi_inv = exp((state[curr]->pos[3] + H[3])) + phiInv0;
    // Derivative of t_phi_inv with respect to st->pos[2]
    Jacobian[3] = exp((state[curr]->pos[3] + H[3]));

    // Transform the tau parameter t_tau
    double t_tau = exp((state[curr]->pos[4] + H[4]));
    // Derivative of t_tau with respect to st->pos[3]
    Jacobian[4] = t_tau;

    // For spline coefficients, the transformation is identity, so the Jacobian is 1
    int i;
    for (i = 5; i < kNumParam; i++) {
        Jacobian[i] = 1.0;
    }

    // Compute the log of the absolute value of the determinant of the Jacobian matrix
    // Since the parameters are transformed independently, the total Jacobian determinant is the product of the individual derivatives
    double logJ = 0.0;
    // Multiply Jacobian contributions for transformed parameters
    for (i = 0; i < 5; i++) {
        // Since log(1) = 0, we can skip Jacobian[i] == 1
        if (Jacobian[i] != 1.0) {
            if (Jacobian[i] <= 0.0) {
                // Logarithm of non-positive number is undefined
                // Return -INFINITY to indicate that the Jacobian determinant is zero or invalid
                return -INFINITY;
            }
            logJ += log(Jacobian[i]);
        }
    }

    // Return the total log-Jacobian
    return logJ;
}

void JacobianMatSIKR_Incidence(double *jac_mat) {
    // This function computes the inverse of the Jacobian matrix entries for the transformation
    // from the parameters in the state vector st->pos[] to the transformed parameters.
    // The inverse Jacobian entries are stored in the array jac_mat[].
    // It accounts for different options for t_gamma (fixed, bounded, unbounded) and other parameters.

    double t_gamma, t_phi_inv;

    // Compute the Jacobian entries for each parameter based on their transformations
    // and store their inverses in jac_mat[]

    // Transform the gamma parameter
    if (g_inpset->gammaFixed == 1){
        // If gamma is fixed, the Jacobian entry is 1
        jac_mat[0] = 1.0;  // No change in volume for fixed parameter
    }
    else{
        if (g_inpset->gammaBounded == 1){
            // Transform gamma using a logistic function to ensure it stays within [g_inpset->gammaLower, g_inpset->gammaUpper]
            t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[0]));
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
            // Transform gamma using an exponential function and a shift
            t_gamma = exp(state[curr]->pos[0]) + gamma0;
            // Derivative of t_gamma with respect to st->pos[0]
            double J0 = exp(state[curr]->pos[0]);
            if (isnan(J0) || isinf(J0) || J0 == 0.0) {
                jac_mat[0] = 0.0;
            } else {
                jac_mat[0] = 1.0 / J0;
            }
        }
    }

    // Transform the initial infected count parameter t_s0
    double t_s0 = 0.0 + (num_pop - 0) / (1 + exp(-state[curr]->pos[1]));
    // Derivative of t_s0 with respect to st->pos[1]
    double exp_neg_pos1 = exp(-state[curr]->pos[1]);
    double denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
    double J1 = (num_pop - 0) * exp_neg_pos1 / denom1;
    // Check for NaN or Inf
    if (isnan(J1) || isinf(J1) || J1 == 0.0) {
        jac_mat[1] = 0.0;  // Avoid division by zero
    } else {
        jac_mat[1] = 1.0 / J1;
    }

    // Transform the initial infected count parameter t_i0
    double t_i0 = 0.0 + (num_pop - 0) / (1 + exp(-state[curr]->pos[2]));
    // Derivative of t_i0 with respect to st->pos[1]
    exp_neg_pos1 = exp(-state[curr]->pos[2]);
    denom1 = (1 + exp_neg_pos1) * (1 + exp_neg_pos1);
    double J2 = (num_pop - 0) * exp_neg_pos1 / denom1;
    // Check for NaN or Inf
    if (isnan(J2) || isinf(J2) || J2 == 0.0) {
        jac_mat[2] = 0.0;  // Avoid division by zero
    } else {
        jac_mat[2] = 1.0 / J2;
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

    // Transform the tau parameter t_tau
    double t_tau = exp(state[curr]->pos[4]);
    // Derivative of t_tau with respect to st->pos[3]
    double J4 = t_tau;  // Since derivative of exp(x) is exp(x)
    if (isnan(J4) || isinf(J4) || J4 == 0.0) {
        jac_mat[4] = 0.0;
    } else {
        jac_mat[4] = 1.0 / J4;
    }

    // For spline coefficients and other parameters without transformation, the Jacobian is 1
    int i;
    for (i = 5; i < kNumParam; i++) {
        jac_mat[i] = 1.0;
    }
}

void TransfSIKR_Incidence() {
    // This function applies transformations to the parameters to map them to the real line.
    // It accounts for different options for t_gamma (fixed, bounded, unbounded).

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

    // Transform initial infected count parameter t_s0 (st->pos[1])
    // Original transformation: t_s0 = 0.0 + (num_pop - 0.0) / (1 + exp(-st->pos[1]))
    // Inverse transformation: st->pos[1] = log((t_s0 - 0.0) / (num_pop - t_s0))
    double t_s0 = state[curr]->pos[1];  // Current value of t_s0
    state[curr]->pos[1] = log((t_s0 - 0.0) / (num_pop - t_s0));

    // Transform initial infected count parameter t_i0 (st->pos[2])
    // Original transformation: t_i0 = 0.0 + (num_pop - 0.0) / (1 + exp(-st->pos[2]))
    // Inverse transformation: st->pos[1] = log((t_i0 - 0.0) / (num_pop - t_i0))
    double t_i0 = state[curr]->pos[2];  // Current value of t_i0
    state[curr]->pos[2] = log((t_i0 - 0.0) / (num_pop - t_i0));

    // Transform phi inverse parameter t_phi_inv (st->pos[3])
    // Original transformation: t_phi_inv = exp(st->pos[3]) + phiInv0
    // Inverse transformation: st->pos[3] = log(t_phi_inv - phiInv0)
    double t_phi_inv = state[curr]->pos[3];  // Current value of t_phi_inv
    state[curr]->pos[3] = log(t_phi_inv - phiInv0);

    // Transform tau parameter t_tau (st->pos[4])
    // Original transformation: t_tau = exp(st->pos[4])
    // Inverse transformation: st->pos[4] = log(t_tau)
    double t_tau = state[curr]->pos[4];  // Current value of t_tau
    state[curr]->pos[4] = log(t_tau);
}

void InvTransfSIKR_Incidence() {
    // This function applies inverse transformations to the parameters to map them back to their original space.
    // It accounts for different options for t_gamma (fixed, bounded, unbounded).

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
    state[curr]->pos[2] = 0.0 + (num_pop - 0.0) / (1.0 + exp(-state[curr]->pos[2]));

    // Transform phi inverse parameter t_phi_inv (st->pos[3])
    // Transformation: t_phi_inv = exp(st->pos[3]) + phiInv0
    state[curr]->pos[3] = exp(state[curr]->pos[3]) + phiInv0;

    // Transform tau parameter t_tau (st->pos[4])
    // Transformation: t_tau = exp(st->pos[4])
    state[curr]->pos[4] = exp(state[curr]->pos[4]);
}

void HessianSIKR_Incidence(){
    int iii, kk;
    double h = 1e-8;

    for (kk = 0; kk < kNumParam; kk++){
        for (iii = 0; iii < kNumParam; iii++){
            if (iii != kk){
                temporalStatePlusH[curr]->pos[iii] = state[curr]->pos[iii];
                temporalStateMinusH[curr]->pos[iii] = state[curr]->pos[iii];
            } else {
                temporalStatePlusH[curr]->pos[iii] = (state[curr]->pos[iii] + h);
                temporalStateMinusH[curr]->pos[iii] = (state[curr]->pos[iii] - h);
            }
        }

        GradientSIKR_Incidence(temporalStatePlusH[curr]);
        GradientSIKR_Incidence(temporalStateMinusH[curr]);

        CSM_STORE(state[curr]->hes, kNumParam, kk, kk, temporalStatePlusH[curr]->grad[kk] / (2 * h) - temporalStateMinusH[curr]->grad[kk] / (2 * h));
        for (iii = kk + 1; iii < kNumParam; iii++){
            CSM_STORE(state[curr]->hes, kNumParam, kk, iii, temporalStatePlusH[curr]->grad[iii] / (2 * h) - temporalStateMinusH[curr]->grad[iii] / (2 * h));
        }
    }
}

void GradientDescentSIKR_Incidence() { // double **positions
    int i, j, ii, jj;
    int init_attempts=0, init_max_attempts=100;
    int l_rate_attempts=0, l_rate_attempts_max = 20;
    int save_counter = 0;  // Counter for saved positions
    char filename[100];
    double logPosterior;
    double st_temp[kNumParam];
    PriorSetting prior_settings[kNumParam - 2];
    printf("Starting Gradient Descent\n");
    printf("------------------------------\n");
    printf("LogPrior: %f\n", PriorSIKR_Incidence());
    printf("LogJacobian: %f\n", LogJacobianSIKR_Incidence());
    printf("LogLIkelihood: %f\n", LoglikelihoodSIKR_Incidence());

    logPosterior = PriorSIKR_Incidence() + LoglikelihoodSIKR_Incidence() + LogJacobianSIKR_Incidence();
    // double h = 1e-7;
    // double H_values[1 + kNumParam][kNumParam];
    // double LogPosterior_values[1 + kNumParam];
    // double Num_Derivative[kNumParam];

    // for (i = 0; i < (1 + kNumParam); i++){
    //     if (i == 0){
    //         for (j = 0; j < (kNumParam); j++){
    //             H_values[i][j] = 0;
    //         }
    //     } else{
    //         for (j = 0; j < (i - 1); j++){
    //             H_values[i][j] = 0;
    //         }
    //         H_values[i][i - 1] = h;
    //         for (j = i; j < (kNumParam); j++){
    //             H_values[i][j] = 0;
    //         }
    //     }
    // }
    
    while(!isfinite(logPosterior)){
        if(init_attempts < init_max_attempts){
            init_attempts += 1;
            printf("Initial point is not computable (We will try another initialization)\n");
            
            // Initialize t_gamma prior setting
            prior_settings[0] =initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
            // // Initialize t_s0 prior setting
            // prior_settings[1] =initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2);
            // // Initialize t_i0 prior setting
            // prior_settings[2] =initialize_prior(g_inpset->I0_prior, g_inpset->I0_param_1, g_inpset->I0_param_2);
            // Initialize t_phi_inv prior setting
            prior_settings[1] =initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
            // Initialize t_tau prior setting
            prior_settings[2] =initialize_prior(g_inpset->tau_prior, g_inpset->tau_param_1, g_inpset->tau_param_2);
            


            // For spline coefficients (indices 4 onwards)
            for (i = 0; i < g_inpset->num_basis_spline; i++) {
                prior_settings[3 + i].dist_type = PRIOR_UNIFORM; // Example: Uniform distribution
                prior_settings[3 + i].param1 = 0.5;       // a
                prior_settings[3 + i].param2 = 2.5;    // b
            }

            // Initialize the global random number generator (if not already initialized)
            // if (g_rng == NULL) {
            //     const gsl_rng_type * T;
            //     gsl_rng_env_setup();
            //     T = gsl_rng_default;
            //     g_rng = gsl_rng_alloc(T);
            //     gsl_rng_set(g_rng, seed); // Seed with current time or a fixed seed
            // }

            // Sample initial values for gamma, I_0, phi_inv, and tau from their priors
            // For gamma (state[curr]->pos[0])
            if (g_inpset->gammaFixed == 1) {
                // If gamma is fixed, use the fixed value
                state[curr]->pos[0] = g_inpset->gamma_fixed_value;
            } else {
                state[curr]->pos[0] = sample_from_prior(&prior_settings[0]);  // gamma
            }

            // // For I_0 (state[curr]->pos[2])
            // state[curr]->pos[2] = sample_from_prior(&prior_settings[2]);  // I_0

            // // For S_0 (state[curr]->pos[1])
            // state[curr]->pos[1] = num_pop - (state[curr]->pos[2] + sample_from_prior(&prior_settings[1]));  // S_0
            double dirich_alpha[3], dirich_x[3];
            dirich_alpha[0] = g_inpset->dirich_alpha_0;
            dirich_alpha[1] = g_inpset->dirich_alpha_1;
            dirich_alpha[2] = g_inpset->dirich_alpha_2;

            sample_dirichlet_3_param(dirich_alpha, dirich_x);
            state[curr]->pos[1] = num_pop * dirich_x[0];
            state[curr]->pos[2] = num_pop * dirich_x[1];

            // For phi_inv (state[curr]->pos[2])
            // state[curr]->pos[3] = sample_from_prior(&prior_settings[1]);  // phi_inv
            state[curr]->pos[3] = 0.005;

            // For tau (state[curr]->pos[3])
            state[curr]->pos[4] = sample_from_prior(&prior_settings[2]);  // tau

            double beta_0 = sample_from_prior(&prior_settings[3]);

            // Initialize spline coefficients (state[curr]->pos[5] onwards)
            for (i = 0; i < g_inpset->num_basis_spline; i++) {
                // Assuming the prior for spline coefficients is defined in prior_settings[5 + i]
                state[curr]->pos[5 + i] = -beta_0;
                // state[curr]->pos[5 + i] = sample_from_prior(&prior_settings[5 + i]);
            }
            printf("New Random Initial Starting Point (%d/%d)\n", init_attempts, init_max_attempts);
            for(i=0;i<kNumParam;i++){
                printf("%f ", state[curr]->pos[i]);
            }
            printf("\n");
            printf("S0, I0, R0: %f, %f, %f\n", state[curr]->pos[1], state[curr]->pos[2], num_pop - state[curr]->pos[1] - state[curr]->pos[2]);

            // Transform parameters to the real line
            TransfSIKR_Incidence();

            // printf("------------------------------\n");
            // printf("LogPrior: %f\n", PriorSIKR_Incidence());
            // printf("LogJacobian: %f\n", LogJacobianSIKR_Incidence());
            // printf("LogLIkelihood: %f\n", LoglikelihoodSIKR_Incidence());

            logPosterior = PriorSIKR_Incidence() + LoglikelihoodSIKR_Incidence() + LogJacobianSIKR_Incidence();
            } else{
                printf("Initial point is not computable (All initialization attempts failed). Try changing the priors. Exiting...\n");
                exit(EXIT_FAILURE);
            }  
    }
    
    printf("Initial point is computable.\n");
    
    // printf("%d%%\n", 100 * 0 / g_inpset->max_iter);
    // printf("------------------------------\n");
    // printf("Log-Posterior: %.8f \n", logPosterior);
    // printf("------------------------------\n");

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
        GradientSIKR_Incidence(state[curr]);

        // for(j = 0; j < (1 + kNumParam); j++){
        //     LogPosterior_values[j] = PriorSIKR_Incidence_2(H_values[j]) + LoglikelihoodSIKR_Incidence_2(H_values[j]) + LogJacobianSIKR_Incidence_2(H_values[j]);
        // }


        // for(j = 1; j < (1 + kNumParam); j++){
        //     Num_Derivative[j - 1] = (LogPosterior_values[j] - LogPosterior_values[0]) / h;
        // }
        // printf("Numerical Derivative\n");
        // for (j = 0; j < (kNumParam); j++){
        //     printf("%f ", Num_Derivative[j]);
        // }
        // printf("\n");
        // printf("------------------------------\n");
        // printf("gradient: \n");
        // for (j = 0; j < (kNumParam); j++){
        //     printf("%f ", state[curr]->grad[j]);
        // }
        // printf("\n");
        // printf("------------------------------\n");

        // Update each parameter by moving in the direction of the gradient
        for (j = 0; j < kNumParam; j++) {
            st_temp[j] = state[curr]->pos[j];
            state[curr]->pos[j] += g_inpset->learning_rate * state[curr]->grad[j];
            // state[curr]->pos[j] += g_inpset->learning_rate * Num_Derivative[j];

        }

        // Limit the gradient step to fulfill S + I + R = N. If R = N - S  - I < 0, we truncate S so R = 0. 
        InvTransfSIKR_Incidence();
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

        TransfSIKR_Incidence();

        // printf("------------------------------\n");
        // printf("LogPrior: %f\n", PriorSIKR_Incidence());
        // printf("LogJacobian: %f\n", LogJacobianSIKR_Incidence());
        // printf("LogLIkelihood: %f\n", LoglikelihoodSIKR_Incidence());
        logPosterior = PriorSIKR_Incidence() + LoglikelihoodSIKR_Incidence() + LogJacobianSIKR_Incidence();
        if (!isfinite(logPosterior)){
            l_rate_attempts+=1;
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
                fprintf(stderr, "ReadModelSIKR_Incidence - Failed to open starting point file '%s'\n", filename);
                exit(EXIT_FAILURE);
            }
            // Transform parameters from the real line
            InvTransfSIKR_Incidence();
            for (ii = 0; ii < kNumParam; ii++) {
                fprintf(fp_startingpoint, "%lf\n", state[curr]->pos[ii]);
            }
            fclose(fp_startingpoint);
            save_counter++;
            // Transform parameters to the real line
            TransfSIKR_Incidence();
        }
        if (g_inpset->max_iter <= 10){
            printf("LogPrior: %f\n", PriorSIKR_Incidence());
            printf("LogJacobian: %f\n", LogJacobianSIKR_Incidence());
            printf("LogLIkelihood: %f\n", LoglikelihoodSIKR_Incidence());
            // logPosterior = PriorSIKR_Incidence() + LoglikelihoodSIKR_Incidence() + LogJacobianSIKR_Incidence();
            printf("%d\n", i);
            printf("------------------------------\n");
            printf("Log-Posterior: %.8f \n", logPosterior);
            printf("------------------------------\n");
        } else{
            if (i % (g_inpset->max_iter / 10) == 0 && i > 0){
                // logPosterior = PriorSIKR_Incidence() + LoglikelihoodSIKR_Incidence() + LogJacobianSIKR_Incidence();
                printf("%d%%\n", 100 * i / g_inpset->max_iter);
                printf("------------------------------\n");
                printf("Log-Posterior: %.8f \n", logPosterior);
                InvTransfSIKR_Incidence();
                printf("Current proposal\n");
                for(jj=0;jj<kNumParam;jj++){
                    printf("%f ", state[curr]->pos[jj]);
                }
                printf("\n");
                printf("------------------------------\n");
                TransfSIKR_Incidence();
            }
        }  
    }
}

int ReadModelSIKR_Incidence() {

    char filename[100];
    FILE *fptr_obs;
    FILE *fptr_initialpoint;
    FILE *fptr_splinebasis;
    int i = 0, ii = 0;
    int symbol;
    double num;

    kNumObs = 0;
    
    //open file with dataset
    sprintf(filename, "./benchmarks/epidemics/%s_data.txt", g_inpset->data);

    fptr_obs = fopen(filename, "r");
    if((fptr_obs) == NULL){
        fprintf(stderr, "sir_model_synthetic.c - Failed to open the dataset file '%s'", filename);
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
    kNumParam = 5 + g_inpset->num_basis_spline;
    Prior = PriorSIKR_Incidence;
    Loglikelihood = LoglikelihoodSIKR_Incidence;
    Gradient = GradientSIKR_Incidence;
    Transf = TransfSIKR_Incidence;
    InvTransf = InvTransfSIKR_Incidence;
    LogJacobian = LogJacobianSIKR_Incidence;
    JacobianMat = JacobianMatSIKR_Incidence;
    Hessian = HessianSIKR_Incidence;

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

    //check if spline basis file exists
    int bool_splinebasis_exists = -1;
    sprintf(filename, "./output/%s/splinebasis_%s.txt", id_global, id_global);
    fptr_splinebasis = fopen(filename, "r");
    if ((fptr_splinebasis) == NULL) {
        bool_splinebasis_exists = 0;
    } else{
        bool_splinebasis_exists = 1;
    }

    fprintf(red_lfp, "bool_initialpointfile_exists is %i\n", bool_initialpointfile_exists);
    printf("bool_initialpointfile_exists is %i\n", bool_initialpointfile_exists);

    fprintf(red_lfp, "bool_splinebasis_exists is %i\n", bool_splinebasis_exists);
    printf("bool_splinebasis_exists is %i\n", bool_splinebasis_exists);

    //depending on inputfile generate a random starting point or read in from file
    if((bool_initialpointfile_exists == 0) || (bool_splinebasis_exists == 0)){
        // fprintf(stderr, "Error in SIR model input: initialpoints and/or splinebasis files not found. \n");
        // exit(EXIT_FAILURE);
        // If files do not exist, use default initial values based on the number of components
        // If files do not exist, sample initial values from priors

        // Define prior settings for each parameter
        PriorSetting prior_settings[kNumParam - 2];

        // Initialize prior_settings array dynamically
        
        // Initialize t_gamma prior setting
        prior_settings[0] =initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
        // // Initialize t_s0 prior setting
        // prior_settings[1] =initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2);
        // // Initialize t_i0 prior setting
        // prior_settings[2] =initialize_prior(g_inpset->I0_prior, g_inpset->I0_param_1, g_inpset->I0_param_2);

        // Initialize t_phi_inv prior setting
        prior_settings[1] =initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
        // Initialize t_tau prior setting
        prior_settings[2] =initialize_prior(g_inpset->tau_prior, g_inpset->tau_param_1, g_inpset->tau_param_2);
        


        // For spline coefficients (indices 5 onwards)
        for (i = 0; i < g_inpset->num_basis_spline; i++) {
            prior_settings[3 + i].dist_type = PRIOR_UNIFORM; // Example: Uniform distribution
            prior_settings[3 + i].param1 = 0.5;       // a
            prior_settings[3 + i].param2 = 2.5;    // b
        }

        // Initialize the global random number generator (if not already initialized)
        // if (g_rng == NULL) {
        //     const gsl_rng_type * T;
        //     gsl_rng_env_setup();
        //     T = gsl_rng_default;
        //     g_rng = gsl_rng_alloc(T);
        //     gsl_rng_set(g_rng, seed); // Seed with current time or a fixed seed
        // }

        // Sample initial values for gamma, S_0, I_0, phi_inv, and tau from their priors
        // For gamma (state[curr]->pos[0])
        if (g_inpset->gammaFixed == 1) {
            // If gamma is fixed, use the fixed value
            state[curr]->pos[0] = g_inpset->gamma_fixed_value;
        } else {
            state[curr]->pos[0] = sample_from_prior(&prior_settings[0]);  // gamma
        }

  
        // // For I_0 (state[curr]->pos[2])
        // state[curr]->pos[2] = sample_from_prior(&prior_settings[2]);  // I_0

        // // For S_0 (state[curr]->pos[1])
        // state[curr]->pos[1] = num_pop - (state[curr]->pos[2] + sample_from_prior(&prior_settings[1]));  // S_0

        double dirich_alpha[3], dirich_x[3];
        dirich_alpha[0] = g_inpset->dirich_alpha_0;
        dirich_alpha[1] = g_inpset->dirich_alpha_1;
        dirich_alpha[2] = g_inpset->dirich_alpha_2;

        sample_dirichlet_3_param(dirich_alpha, dirich_x);
        state[curr]->pos[1] = num_pop * dirich_x[0];
        state[curr]->pos[2] = num_pop * dirich_x[1];

        // For phi_inv (state[curr]->pos[3])
        // state[curr]->pos[3] = sample_from_prior(&prior_settings[1]);  // phi_inv
        state[curr]->pos[3] = 0.005;

        // For tau (state[curr]->pos[4])
        state[curr]->pos[4] = sample_from_prior(&prior_settings[2]);  // tau

        double beta_0 = sample_from_prior(&prior_settings[3]);
        printf("Beta_0: %f\n", -beta_0);

        // Initialize spline coefficients (state[curr]->pos[4] onwards)
        for (i = 0; i < g_inpset->num_basis_spline; i++) {
            // Assuming the prior for spline coefficients is defined in prior_settings[5 + i]
            state[curr]->pos[5 + i] = -beta_0;
            // state[curr]->pos[5 + i] = sample_from_prior(&prior_settings[5 + i]);
        }
	    printf("Random Initial Starting Point\n");
        for(i=0;i<kNumParam;i++){
            printf("%f ", state[curr]->pos[i]);
        }
        printf("\n");
        printf("S0, I0, R0: %f, %f, %f\n", state[curr]->pos[1], state[curr]->pos[2], num_pop - state[curr]->pos[1] - state[curr]->pos[2]);
    } else if((bool_initialpointfile_exists == 1) && (bool_splinebasis_exists == 1)){ //read in from file

        // count number of spaces in the first line
        int numSizeInitialPoint = 0;
        symbol = 1; //initialise symbol
        while(!feof(fptr_initialpoint) && symbol != '\n'){
            symbol = fgetc(fptr_initialpoint);
            if(isspace(symbol)){ numSizeInitialPoint++; }
        }

        fprintf(red_lfp, "numSizeInitialPoint is %i\n", numSizeInitialPoint);
        printf("numSizeInitialPoint is %i\n", numSizeInitialPoint);

        int numSizeSplineBasis = 0;
        symbol = 1; //initialise symbol
        while(!feof(fptr_splinebasis) && symbol != '\n'){
            symbol = fgetc(fptr_splinebasis);
            if(isspace(symbol)){ numSizeSplineBasis++; }
        }

        fprintf(red_lfp, "numSizeSplineBasis is %i\n", numSizeSplineBasis);
        printf("numSizeSplineBasis is %i\n", numSizeSplineBasis);

        // if(g_inpset->num_basis_spline != numSizeSplineBasis){
        //     fprintf(stderr, "sir_model_synthetic.c - num_basis_spline (%d) is not equal to numSizeSplineBasis (%i) in sir_model_synthetic.c.\n", g_inpset->num_basis_spline, numSizeSplineBasis);
        //     // exit(EXIT_FAILURE);
        // }

        if(kNumParam != numSizeInitialPoint + numSizeSplineBasis){
            fprintf(stderr, "sir_model_synthetic.c - kNumParam (%d) is not equal to numSizeInitialPoint + numSizeSplineBasis(%i) in sir_model_synthetic.c.\n", kNumParam, numSizeInitialPoint + numSizeSplineBasis);
            fprintf(stderr, "sir_model_synthetic.c - a random initialization point was drawn from the priors.\n");
            // exit(EXIT_FAILURE);

            // Define prior settings for each parameter
            PriorSetting prior_settings[kNumParam - 2];

            // Initialize prior_settings array dynamically
            
            // Initialize t_gamma prior setting
            prior_settings[0] =initialize_prior(g_inpset->gamma_prior, g_inpset->gamma_param_1, g_inpset->gamma_param_2);
            // // Initialize t_s0 prior setting
            // prior_settings[1] =initialize_prior(g_inpset->S0_prior, g_inpset->S0_param_1, g_inpset->S0_param_2);
            // // Initialize t_i0 prior setting
            // prior_settings[2] =initialize_prior(g_inpset->I0_prior, g_inpset->I0_param_1, g_inpset->I0_param_2);

            // Initialize t_phi_inv prior setting
            prior_settings[1] =initialize_prior(g_inpset->phi_inv_prior, g_inpset->phi_inv_param_1, g_inpset->phi_inv_param_2);
            // Initialize t_tau prior setting
            prior_settings[2] =initialize_prior(g_inpset->tau_prior, g_inpset->tau_param_1, g_inpset->tau_param_2);
            


            // For spline coefficients (indices 5 onwards)
            for (i = 0; i < g_inpset->num_basis_spline; i++) {
                prior_settings[3 + i].dist_type = PRIOR_UNIFORM; // Example: Uniform distribution
                prior_settings[3 + i].param1 = 0.5;       // a
                prior_settings[3 + i].param2 = 2.5;    // b
            }

            // Initialize the global random number generator (if not already initialized)
            // if (g_rng == NULL) {
            //     const gsl_rng_type * T;
            //     gsl_rng_env_setup();
            //     T = gsl_rng_default;
            //     g_rng = gsl_rng_alloc(T);
            //     gsl_rng_set(g_rng, seed); // Seed with current time or a fixed seed
            // }

            // Sample initial values for gamma, S_0, I_0, phi_inv, and tau from their priors
            // For gamma (state[curr]->pos[0])
            if (g_inpset->gammaFixed == 1) {
                // If gamma is fixed, use the fixed value
                state[curr]->pos[0] = g_inpset->gamma_fixed_value;
            } else {
                state[curr]->pos[0] = sample_from_prior(&prior_settings[0]);  // gamma
            }

    
            // // For I_0 (state[curr]->pos[2])
            // state[curr]->pos[2] = sample_from_prior(&prior_settings[2]);  // I_0

            // // For S_0 (state[curr]->pos[1])
            // state[curr]->pos[1] = num_pop - (state[curr]->pos[2] + sample_from_prior(&prior_settings[1]));  // S_0

            double dirich_alpha[3], dirich_x[3];
            dirich_alpha[0] = g_inpset->dirich_alpha_0;
            dirich_alpha[1] = g_inpset->dirich_alpha_1;
            dirich_alpha[2] = g_inpset->dirich_alpha_2;

            sample_dirichlet_3_param(dirich_alpha, dirich_x);
            state[curr]->pos[1] = num_pop * dirich_x[0];
            state[curr]->pos[2] = num_pop * dirich_x[1];

            // For phi_inv (state[curr]->pos[3])
            // state[curr]->pos[3] = sample_from_prior(&prior_settings[1]);  // phi_inv
            state[curr]->pos[3] = 0.005;

            // For tau (state[curr]->pos[4])
            state[curr]->pos[4] = sample_from_prior(&prior_settings[2]);  // tau

            double beta_0 = sample_from_prior(&prior_settings[3]);

            // Initialize spline coefficients (state[curr]->pos[4] onwards)
            for (i = 0; i < g_inpset->num_basis_spline; i++) {
                // Assuming the prior for spline coefficients is defined in prior_settings[5 + i]
                state[curr]->pos[5 + i] = -beta_0;
                // state[curr]->pos[5 + i] = sample_from_prior(&prior_settings[5 + i]);
            }
            printf("Random Initial Starting Point\n");
            for(i=0;i<kNumParam;i++){
                printf("%f ", state[curr]->pos[i]);
            }
            printf("\n");
            printf("S0, I0, R0: %f, %f, %f\n", state[curr]->pos[1], state[curr]->pos[2], num_pop - state[curr]->pos[1] - state[curr]->pos[2]);
        } else{
            //go to beginning of file
            rewind(fptr_initialpoint);
            rewind(fptr_splinebasis);

            //read initial points in one by one
            for (ii = 0; ii < numSizeInitialPoint; ii++) {
                if (fscanf(fptr_initialpoint, "%lf", &num) != EOF) {
                    state[curr]->pos[ii] = num;
                }
            }

            fclose(fptr_initialpoint);

            //read spline basis in one by one
            for (ii = 0; ii < numSizeSplineBasis; ii++) {
                if (fscanf(fptr_splinebasis, "%lf", &num) != EOF) {
                    state[curr]->pos[numSizeInitialPoint + ii] = num;
                }
            }

            fprintf(red_lfp, "initialpoints are ");
            for (i = 0; i < numSizeInitialPoint + numSizeSplineBasis; i++) {
                fprintf(red_lfp, "%lf ", state[curr]->pos[i]);
            }
            fprintf(red_lfp, "\n");

            fclose(fptr_splinebasis);

            printf("initialpoints are ");
            for (i = 0; i < numSizeInitialPoint + numSizeSplineBasis; i++) {
                printf("%lf ", state[curr]->pos[i]);
            }
            printf("\n");
        }
    } else{
        fprintf(stderr, "sir_model_synthetic.c - If-else not matched for initial position\n");
        exit(EXIT_FAILURE);
    }
    // end initialize position variables
    printf("SUNDIALS version: %d.%d.%d\n", SUNDIALS_VERSION_MAJOR, SUNDIALS_VERSION_MINOR, SUNDIALS_VERSION_PATCH);

    // Perform gradient descent if enabled
    if (g_inpset->do_grad_desc == 1) {
        // Transform parameters to the real line before gradient descent
        TransfSIKR_Incidence();

        // Perform gradient descent and get positions
        GradientDescentSIKR_Incidence(); // positions

        // Transform parameters back from the real line
        InvTransfSIKR_Incidence();

        // Save the final position from gradient descent to the starting point file
        sprintf(filename, "./output/%s/startingpoint_%s.txt", id_global, id_global);
        FILE *fp_startingpoint = fopen(filename, "w");
        if (fp_startingpoint == NULL) {
            fprintf(stderr, "ReadModelSIKR_Incidence - Failed to open starting point file '%s'\n", filename);
            exit(EXIT_FAILURE);
        }
        for (ii = 0; ii < kNumParam; ii++) {
            fprintf(fp_startingpoint, "%lf\n", state[curr]->pos[ii]);
        }
        fclose(fp_startingpoint);
    }

    /*** initialize momenta vectors from a gaussian distr ****/
    if (state[curr]->mom == NULL) {
        printf("sir_model_synthetic.c - Error in allocating memory (state mom)");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < kNumParam; i++) {
        state[curr]->mom[i] = gsl_ran_gaussian_ziggurat(g_rng, 1.);
    }

    /***** END initialization ******/
    Ham = (TypeHamiltonian *)malloc(sizeof(TypeHamiltonian));

    return 0;
}
