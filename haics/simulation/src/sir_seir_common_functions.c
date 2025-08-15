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
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */
#include <sundials/sundials_config.h>

//#include <legendre/legendre_polynomial.h>
/* Accessor macros */

#include "Globals.h"
#include "utils.h"
#include "sir_seir_common_functions.h"
#include "read_input.h"

#define M_PI 3.14159265358979323846

/*
 *--------------------------------------------------------------------
 * SPLINE Functions
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
 */

int mindex(int i, int j, int nrow) {
    return (j - 1) * nrow + i;
}

double bs(int nknots, int nspline, int updegree, double x, double * knots) {
    double y, y1, y2, temp1, temp2;
    
    // Base case: B-spline of degree 0 (order 1)
    if (updegree == 1) {
        // Check if x is within the support of the basis function
        if ((x >= knots[nspline - 1]) && (x < knots[nspline]))
            y = 1.0; // Basis function is 1 in this interval
        else
            y = 0.0; // Basis function is 0 outside this interval
    }
    else {
        // Initialize coefficients
        temp1 = 0.0;
        temp2 = 0.0;

        // Compute the left fraction coefficient
        if ((knots[nspline + updegree - 2] - knots[nspline - 1]) > 0)
            temp1 = (x - knots[nspline - 1]) / (knots[nspline + updegree - 2] - knots[nspline - 1]);

        // Compute the right fraction coefficient
        if ((knots[nspline + updegree - 1] - knots[nspline]) > 0)
            temp2 = (knots[nspline + updegree - 1] - x) / (knots[nspline + updegree - 1] - knots[nspline]);

        // Recursive calls for lower-degree B-spline basis functions
        y1 = bs(nknots, nspline, updegree - 1, x, knots);
        y2 = bs(nknots, nspline + 1, updegree - 1, x, knots);

        // Combine the results using the recursion formula
        y =  temp1 * y1 + temp2 * y2;
    }
    return y;
}

void splinebasis(int *d, int *n, int *m, double * x, double * knots, double * basis) {
    int mm = *m; // Total number of knots
    int dd = *d; // Degree of the spline
    int nn = *n; // Number of data points
    int k = mm - dd - 1; // Number of basis functions (order)
    int i, j, ir, jr;

    // Loop over each data point
    for (i = 0; i < nn; i++) {
        ir = i + 1; // Adjust index for 1-based systems if needed

        // Special case: if x[i] equals the last knot
        if (x[i] == knots[mm - 1]) {
            // Set the last basis function to 1
            basis[mindex(ir, k, nn) - 1] = 1.0;

            // Set all other basis functions to 0
            for (j = 0; j < (k - 1); j++) {
                jr = j + 1;
                basis[mindex(ir, jr, nn) - 1] = 0.0;
            }
        } else {
            // Compute all basis functions at x[i]
            for (j = 0; j < k; j++) {
                jr = j + 1;
                // Evaluate the j-th B-spline basis function at x[i]
                basis[mindex(ir, jr, nn) - 1] = bs(mm, jr, dd + 1, x[i], knots);
            }
        }
    }
}


double splinefunction(double * beta, int *d, int *n, int *m, double * x, double * knots, double * basis, int *nf){
    double splineValue = 0.0;
    int ii, nnf = *nf;

    // Compute the spline basis functions at the data points x
    splinebasis(d, n, m, x, knots, basis);

    // Compute the spline function value as the dot product of beta and basis
    for (ii = 0; ii < nnf; ii++){
        splineValue += beta[ii] * basis[ii];
    }
    return splineValue;
}


/*
 *--------------------------------------------------------------------
 * General Functions
 *--------------------------------------------------------------------
 */

double underReport(int t, double u0, double u1 , int t0, int t1){
    double res = 0.0;
    if (t <= t0) {
        res = u0;
    } else {
	if (t >= t1){
	    res = u1;
	} else {
	    res = u0 + ((u1 - u0) / (t1 - t0)) * (t - t0);
	}
    }
    return res;
}

double vectors_dot_prod(double *x, double *y, int n){
    double res = 0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

void matrix_vector_mult(double **mat, double *vec, double *result, int rows, int cols){ // in matrix form: result = mat * vec;
    int i;
    for (i = 0; i < rows; i++)
    {
        result[i] = vectors_dot_prod(mat[i], vec, cols);
    }
}

void matrix_matrix_mult(double **matrix_left, double **matrix_right, double **result_matrix, int rows, int cols, int prod_dim){ // in matrix form: result = mat * vec;
    int i,j,k;
    double result;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++){
            result = 0;
            for(k = 0; k <  prod_dim; k++){
                result += matrix_left[i][k] * matrix_right[k][j];
            }
            result_matrix[i][j] = result;
        }

    }
}


/*--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Process and verify arguments to cvsfwddenx.
 */

void ProcessArgs(char *argv[], booleantype *sensi, int *sensi_meth, booleantype *err_con)
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
    WrongArgs("SIR_CVODES_synth");

  if (*sensi) {

    if (strcmp(argv[1],"sim") == 0)
      *sensi_meth = CV_SIMULTANEOUS;
    else if (strcmp(argv[1],"stg") == 0)
      *sensi_meth = CV_STAGGERED;
    else if (strcmp(argv[1],"stg1") == 0)
      *sensi_meth = CV_STAGGERED1;
    else
      WrongArgs("SIR_CVODES_synth");

    if (strcmp(argv[2],"t") == 0)
      *err_con = SUNTRUE;
    else if (strcmp(argv[2],"f") == 0)
      *err_con = SUNFALSE;
    else
      WrongArgs("SIR_CVODES_synth");
  }

}

void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = t or f\n");

    exit(0);
}

int check_retval(void *returnvalue, const char *funcname, int opt)
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

// Function to compute the natural logarithm of the beta function
double lbeta(double a, double b) {
    return lgamma(a) + lgamma(b) - lgamma(a + b);
}

// Function to compute the binomial coefficient C(n, k)
double binomial_coefficient(int n, int k) {
    double result = 1.0;
    int i;
    for (i = 1; i <= k; i++) {
        result *= (double)(n - (k - i)) / (double)i;
    }
    return result;
}

// Function to convert a string to PriorDistributionType
PriorDistributionType get_prior_type_from_string(const char* prior_type_str) {
    if (strcmp(prior_type_str, "PRIOR_NORMAL") == 0) return PRIOR_NORMAL;
    if (strcmp(prior_type_str, "PRIOR_GAMMA") == 0) return PRIOR_GAMMA;
    if (strcmp(prior_type_str, "PRIOR_EXPONENTIAL") == 0) return PRIOR_EXPONENTIAL;
    if (strcmp(prior_type_str, "PRIOR_INVGAMMA") == 0) return PRIOR_INVGAMMA;
    if (strcmp(prior_type_str, "PRIOR_UNIFORM") == 0) return PRIOR_UNIFORM;
    if (strcmp(prior_type_str, "PRIOR_BETA") == 0) return PRIOR_BETA;
    if (strcmp(prior_type_str, "PRIOR_LOGNORMAL") == 0) return PRIOR_LOGNORMAL;
    if (strcmp(prior_type_str, "PRIOR_WEIBULL") == 0) return PRIOR_WEIBULL;
    if (strcmp(prior_type_str, "PRIOR_KUMARASWAMY") == 0) return PRIOR_KUMARASWAMY;
    if (strcmp(prior_type_str, "PRIOR_LORENTZ") == 0) return PRIOR_LORENTZ;
    if (strcmp(prior_type_str, "PRIOR_CHISQUARED") == 0) return PRIOR_CHISQUARED;
    if (strcmp(prior_type_str, "PRIOR_CONSTANT") == 0) return PRIOR_CONSTANT;
    // Handle unexpected cases
    fprintf(stderr, "Error: Unknown prior type string %s\n", prior_type_str);
    exit(EXIT_FAILURE);
}

// Modified initialize_prior function using PriorDistributionType
PriorSetting initialize_prior(const char* prior_type_str, double param1, double param2) {
    PriorDistributionType prior_type = get_prior_type_from_string(prior_type_str);

    switch (prior_type) {
        case PRIOR_NORMAL:
            return PRIOR_SETTING_NORMAL(param1, param2);
        case PRIOR_GAMMA:
            return PRIOR_SETTING_GAMMA(param1, param2);
        case PRIOR_EXPONENTIAL:
            return PRIOR_SETTING_EXPONENTIAL(param1);
        case PRIOR_INVGAMMA:
            return PRIOR_SETTING_INVGAMMA(param1, param2);
        case PRIOR_UNIFORM:
            return PRIOR_SETTING_UNIFORM(param1, param2);
        case PRIOR_BETA:
            return PRIOR_SETTING_BETA(param1, param2);
        case PRIOR_LOGNORMAL:
            return PRIOR_SETTING_LOGNORMAL(param1, param2);
        case PRIOR_WEIBULL:
            return PRIOR_SETTING_WEIBULL(param1, param2);
        case PRIOR_KUMARASWAMY:
            return PRIOR_SETTING_KUMARASWAMY(param1, param2);
        case PRIOR_LORENTZ:
            return PRIOR_SETTING_LORENTZ(param1, param2);
        case PRIOR_CHISQUARED:
            return PRIOR_SETTING_CHISQUARED(param1);
        case PRIOR_CONSTANT:
            return PRIOR_SETTING_CONSTANT(param1);
        default:
            fprintf(stderr, "Error: Unsupported PriorDistributionType %d\n", prior_type);
            exit(EXIT_FAILURE);
    }
}

double dirichlet_logdensity_4_param(double x0, double x1, double x2, const double alpha[4]) {
    // Calculate x3 as 1 - x0 - x1 - x2
    double x3 = 1.0 - x0 - x1 - x2;

    // Sum of the alpha parameters
    double sum_alpha = alpha[0] + alpha[1] + alpha[2] + alpha[3];

    // Sum of log(x) * (alpha - 1) terms
    double sum_log_x;

    if (fabs(alpha[3] - 1) < 1e-8) {
        sum_log_x = (alpha[0] - 1.0) * log(x0) +
                       (alpha[1] - 1.0) * log(x1) +
                       (alpha[2] - 1.0) * log(x2);
    } else{
       sum_log_x = (alpha[0] - 1.0) * log(x0) +
                       (alpha[1] - 1.0) * log(x1) +
                       (alpha[2] - 1.0) * log(x2) +
                       (alpha[3] - 1.0) * log(x3);
    }                   

    // Normalization term using the Gamma function
    double log_gamma_sum_alpha = lgamma(sum_alpha);
    double log_gamma_alpha1 = lgamma(alpha[0]);
    double log_gamma_alpha2 = lgamma(alpha[1]);
    double log_gamma_alpha3 = lgamma(alpha[2]);
    double log_gamma_alpha4 = lgamma(alpha[3]);

    double normalization_term = log_gamma_sum_alpha - (log_gamma_alpha1 + log_gamma_alpha2 + log_gamma_alpha3 + log_gamma_alpha4);

    // Return the log-density
    return sum_log_x + normalization_term;
}

double dirichlet_logdensity_3_param(double x0, double x1, const double alpha[3]) {
    // Calculate x2 as 1 - x0 - x1
    double x2 = 1.0 - x0 - x1;

    // Sum of the alpha parameters
    double sum_alpha = alpha[0] + alpha[1] + alpha[2];

    // Sum of log(x) * (alpha - 1) terms
    double sum_log_x;

    if (fabs(alpha[2] - 1) < 1e-8) {
        sum_log_x = (alpha[0] - 1.0) * log(x0) +
                       (alpha[1] - 1.0) * log(x1);
    } else{
       sum_log_x = (alpha[0] - 1.0) * log(x0) +
                       (alpha[1] - 1.0) * log(x1) +
                       (alpha[2] - 1.0) * log(x2);
    }                 


    // Normalization term using the Gamma function
    double log_gamma_sum_alpha = lgamma(sum_alpha);
    double log_gamma_alpha1 = lgamma(alpha[0]);
    double log_gamma_alpha2 = lgamma(alpha[1]);
    double log_gamma_alpha3 = lgamma(alpha[2]);

    double normalization_term = log_gamma_sum_alpha - (log_gamma_alpha1 + log_gamma_alpha2 + log_gamma_alpha3);

    // Return the log-density
    return sum_log_x + normalization_term;
}

double dirichlet_logdensity_4_param_derivative(double x0, double x1, double x2, const double alpha[4], int i) {
    // Calculate x3 as 1 - x0 - x1 - x2
    double x3 = 1.0 - x0 - x1 - x2;

    // Ensure i is within bounds (0, 1, or 2)
    if (i < 0 || i > 2) {
        return 0.0;  // Return 0 if index is out of bounds
    }

    // Select x[i] and alpha[i] based on input
    double xi = (i == 0) ? x0 : (i == 1) ? x1 : x2;
    double alphai = alpha[i];

    // Calculate the derivative term for x[i]
    double term1 = (alphai - 1.0) / xi;
    double term2;
    if (fabs(alpha[3] - 1) < 1e-8){
        term2 = 0;
    } else {
        term2= (alpha[3] - 1.0) / x3;
    }
    // printf("Derivative logdirichdensity component %d: term1 %e, term2 %e\n", i, term1, term2);

    // Return the derivative with respect to x[i]
    return term1 - term2;
}

double dirichlet_logdensity_3_param_derivative(double x0, double x1, const double alpha[3], int i) {
    // Calculate x2 as 1 - x0 - x1
    double x2 = 1.0 - x0 - x1;

    // Ensure i is within bounds (0 or 1)
    if (i < 0 || i > 1) {
        return 0.0;  // Return 0 if index is out of bounds
    }

    // Select x[i] and alpha[i] based on input
    double xi = (i == 0) ? x0 : x1;
    double alphai = alpha[i];

    // Calculate the derivative term for x[i]
    double term1 = (alphai - 1.0) / xi;
    double term2;
    if (fabs(alpha[2] - 1) < 1e-8){
        term2 = 0;
    } else {
        term2= (alpha[2] - 1.0) / x2;
    }

    // Return the derivative with respect to x[i]
    return term1 - term2;
}

// Function to compute the log-prior density given a value and prior settings
double log_prior_density(double x, PriorSetting *prior) {
    switch (prior->dist_type) {
        case PRIOR_NORMAL: {
            // Normal distribution: param1 = mean (mu), param2 = standard deviation (sigma)
            double mu = prior->param1;
            double sigma = prior->param2;
            double diff = x - mu;
            return - (diff * diff) / (2 * sigma * sigma) - log(sqrt(2 * M_PI) * sigma);
        }
        case PRIOR_GAMMA: {
            // Gamma distribution: param1 = shape (k), param2 = scale (theta)
            double k = prior->param1;
            double theta = prior->param2;
            if (x <= 0) return -INFINITY;  // Support: x > 0
            return (k - 1) * log(x) - x / theta - lgamma(k) - k * log(theta);
        }
        case PRIOR_INVGAMMA: {
            // Inverse gamma distribution: param1 = shape (alpha), param2 = scale (beta)
            double alpha = prior->param1;
            double beta = prior->param2;
            if (x <= 0) return -INFINITY;  // Support: x > 0
            return alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - beta / x;
        }
        case PRIOR_BETA: {
            // Beta distribution: param1 = alpha, param2 = beta
            double a = prior->param1;
            double b = prior->param2;
            if (x <= 0 || x >= 1) return -INFINITY;  // Support: 0 < x < 1
            return (a - 1) * log(x) + (b - 1) * log(1 - x) - lbeta(a, b);
        }
        case PRIOR_LOGNORMAL: {
            // Log-normal distribution: param1 = mean of log(x), param2 = sd of log(x)
            double mu = prior->param1;
            double sigma = prior->param2;
            if (x <= 0) return -INFINITY;  // Support: x > 0
            double ln_x = log(x);
            double diff = ln_x - mu;
            return - (diff * diff) / (2 * sigma * sigma) - ln_x - log(sqrt(2 * M_PI) * sigma);
        }
        case PRIOR_WEIBULL: {
            // Weibull distribution: param1 = shape (k), param2 = scale (lambda)
            double k = prior->param1;
            double lambda = prior->param2;
            if (x <= 0) return -INFINITY;  // Support: x > 0
            return log(k / lambda) + (k - 1) * log(x / lambda) - pow(x / lambda, k);
        }
        case PRIOR_KUMARASWAMY: {
            // Kumaraswamy distribution: param1 = a, param2 = b
            double a = prior->param1;
            double b = prior->param2;
            if (x <= 0 || x >= 1) return -INFINITY;  // Support: 0 < x < 1
            return log(a * b) + (a - 1) * log(x) + (b - 1) * log(1 - pow(x, a));
        }
        case PRIOR_LORENTZ: {
            // Lorentz (Cauchy) distribution: param1 = location (x0), param2 = scale (gamma)
            double x0 = prior->param1;
            double gamma = prior->param2;
            double denom = M_PI * gamma * (1 + ((x - x0) / gamma) * ((x - x0) / gamma));
            return -log(denom);
        }
        case PRIOR_EXPONENTIAL: {
            // Exponential distribution: param1 = rate (lambda)
            double lambda = prior->param1;
            double tolerance = 1e-8;
            // printf("lambda, tolerance, x, (-tolerance <= x), return: %f, %e, %f, %d, %f\n", lambda, tolerance, x, (-tolerance <= x), (log(lambda) - lambda * x));
            if (x < -tolerance) return -INFINITY;  // Support: x >= -tolerance
            if ((-tolerance <= x) && (x < 0)) return log(lambda); // we take it as x = 0;
            return log(lambda) - lambda * x;
        }
        case PRIOR_CHISQUARED: {
            // Chi-squared distribution: param1 = degrees of freedom (k)
            double k = prior->param1;
            if (x <= 0) return -INFINITY;  // Support: x > 0
            return (k / 2 - 1) * log(x) - x / 2 - (k / 2) * log(2) - lgamma(k / 2);
        }
        case PRIOR_UNIFORM: {
            // Uniform distribution: param1 = lower bound (a), param2 = upper bound (b)
            double a = prior->param1;
            double b = prior->param2;
            double tolerance = 1e-8;
            if ((x < a - tolerance) || (x > b + tolerance)) return -INFINITY;  // Support: a - tolerance <= x <= b + tolerance
            if ((a - tolerance <= x) &&  (x <= b + tolerance)) return -log(b - a);
        }
        case PRIOR_CONSTANT: {
            // A fixed value for a parameter
            return 0;
        }
        default:
            // Unknown distribution type
            fprintf(stderr, "Unknown prior distribution type.\n");
            exit(EXIT_FAILURE);
    }
}

double log_prior_density_derivative(double x, PriorSetting *prior)
{
    switch (prior->dist_type) {
        case PRIOR_NORMAL: {
            // Normal distribution: derivative is -(x - mu) / sigma^2
            double mu = prior->param1;
            double sigma = prior->param2;
            return - (x - mu) / (sigma * sigma);
        }
        case PRIOR_GAMMA: {
            // Gamma distribution: derivative is (k - 1)/x - 1/θ
            double k = prior->param1;
            double theta = prior->param2;
            if (x <= 0) return 0.0;  // Derivative is undefined; return zero or handle appropriately
            return (k - 1) / x - 1 / theta;
        }
        case PRIOR_INVGAMMA: {
            // Inverse gamma: derivative is - (α + 1)/x + β / x^2
            double alpha = prior->param1;
            double beta = prior->param2;
            if (x <= 0) return 0.0;
            return - (alpha + 1) / x + beta / (x * x);
        }
        case PRIOR_BETA: {
            // Beta distribution: derivative is (α - 1)/x - (β - 1)/(1 - x)
            double a = prior->param1;
            double b = prior->param2;
            if (x <= 0 || x >= 1) return 0.0;
            return (a - 1) / x - (b - 1) / (1 - x);
        }
        case PRIOR_LOGNORMAL: {
            // Log-normal distribution: derivative is [ - (ln(x) - μ) / (σ^2 * x) ] - 1 / x
            double mu = prior->param1;
            double sigma = prior->param2;
            if (x <= 0) return 0.0;
            double ln_x = log(x);
            return - (ln_x - mu) / (sigma * sigma * x) - 1 / x;
        }
        case PRIOR_WEIBULL: {
            // Weibull distribution: derivative is [ (k - 1)/x - k x^{k - 1} / λ^k ] / x
            double k = prior->param1;
            double lambda = prior->param2;
            if (x <= 0) return 0.0;
            return ((k - 1) / x - k * pow(x, k - 1) / pow(lambda, k)) / x;
        }
        case PRIOR_KUMARASWAMY: {
            // Kumaraswamy's distribution derivative is (a - 1)/x - a b x^{a - 1} / (1 - x^{a})
            double a = prior->param1;
            double b = prior->param2;
            if (x <= 0 || x >= 1) return 0.0;
            return (a - 1) / x - (a * b * pow(x, a - 1)) / (1 - pow(x, a));
        }
        case PRIOR_LORENTZ: {
            // Lorentz (Cauchy) distribution: derivative is -2(x - x0) / [γ ( (x - x0)^2 + γ^2 ) ]
            double x0 = prior->param1;
            double gamma = prior->param2;
            double denom = (x - x0) * (x - x0) + gamma * gamma;
            return -2 * (x - x0) / (gamma * denom);
        }
        case PRIOR_EXPONENTIAL: {
            // Exponential distribution: derivative is -λ
            double lambda = prior->param1;
            double tolerance = 1e-8;
            if (x < -tolerance) return 0.0;
            return -lambda;
        }
        case PRIOR_CHISQUARED: {
            // Chi-squared distribution: derivative is (k/2 - 1)/x - 1/2
            double k = prior->param1;
            if (x <= 0) return 0.0;
            return (k / 2 - 1) / x - 0.5;
        }
        case PRIOR_UNIFORM: {
            // Uniform distribution: derivative is zero within support
            double a = prior->param1;
            double b = prior->param2;
            if (x < a || x > b) return 0.0;
            return 0.0;
        }
        case PRIOR_CONSTANT: {
            // A fixed value for a parameter
            return 0.0;
        }
        default:
            // Unknown distribution type
            fprintf(stderr, "Unknown prior distribution type in derivative.\n");
            exit(EXIT_FAILURE);
    }
}


// Function to sample from prior distribution
double sample_from_prior(PriorSetting *prior) {
    switch (prior->dist_type) {
        case PRIOR_NORMAL: {
            // Normal distribution: param1 = mean (mu), param2 = standard deviation (sigma)
            double mu = prior->param1;
            double sigma = prior->param2;
            return gsl_ran_gaussian(g_rng, sigma) + mu;
        }
        case PRIOR_GAMMA: {
            // Gamma distribution: param1 = shape (k), param2 = scale (theta)
            double k = prior->param1;
            double theta = prior->param2;
            return gsl_ran_gamma(g_rng, k, theta);
        }
        case PRIOR_INVGAMMA: {
            // Inverse gamma distribution: param1 = shape (alpha), param2 = scale (beta)
            // Note: GSL does not have an inverse gamma distribution function directly
            // We can sample from the gamma distribution and take the reciprocal
            double alpha = prior->param1;
            double beta = prior->param2;
            double gamma_sample = gsl_ran_gamma(g_rng, alpha, 1.0 / beta);
            return 1.0 / gamma_sample;
        }
        case PRIOR_BETA: {
            // Beta distribution: param1 = alpha, param2 = beta
            double a = prior->param1;
            double b = prior->param2;
            return gsl_ran_beta(g_rng, a, b);
        }
        case PRIOR_LOGNORMAL: {
            // Log-normal distribution: param1 = mean of log(x), param2 = sd of log(x)
            double z = gsl_ran_gaussian(g_rng, prior->param2) + prior->param1;
            return exp(z);
        }
        case PRIOR_WEIBULL: {
            // Weibull distribution: param1 = scale (lambda), param2 = shape (k)
            double lambda = prior->param1;
            double k = prior->param2;
            return gsl_ran_weibull(g_rng, lambda, k);
        }
        case PRIOR_EXPONENTIAL: {
            // Exponential distribution: param1 = rate (lambda)
            double lambda = prior->param1;
            return gsl_ran_exponential(g_rng, 1.0 / lambda);
        }
        case PRIOR_CHISQUARED: {
            // Chi-squared distribution: param1 = degrees of freedom (k)
            double k = prior->param1;
            return gsl_ran_chisq(g_rng, k);
        }
        case PRIOR_UNIFORM: {
            // Uniform distribution: param1 = lower bound (a), param2 = upper bound (b)
            double a = prior->param1;
            double b = prior->param2;
            return gsl_ran_flat(g_rng, a, b);
        }
        case PRIOR_CONSTANT: {
            // A constant prior, we just return the constant itself
            return prior->param1;
        }
        default:
            // Unknown distribution type
            fprintf(stderr, "Unknown prior distribution type in sample_from_prior.\n");
            exit(EXIT_FAILURE);
    }
}

// Function to generate a random sample from a four-parameter Dirichlet distribution
void sample_dirichlet_4_param(const double alpha[4], double x[4]) {
    double gamma_samples[4];
    double sum = 0.0;
    int i;

    // Generate a Gamma sample for each alpha parameter
    for (i = 0; i < 4; i++) {
        gamma_samples[i] = gsl_ran_gamma(g_rng, alpha[i], 1.0);  // scale parameter theta = 1
        sum += gamma_samples[i];
    }

    // Normalize the Gamma samples to obtain Dirichlet proportions
    for (i = 0; i < 4; i++) {
        x[i] = gamma_samples[i] / sum;
    }
}

void sample_dirichlet_3_param(const double alpha[3], double x[3]) {
    double gamma_samples[3];
    double sum = 0.0;
    int i;

    // Generate a Gamma sample for each alpha parameter
    for (i = 0; i < 3; i++) {
        gamma_samples[i] = gsl_ran_gamma(g_rng, alpha[i], 1.0);  // scale parameter theta = 1
        sum += gamma_samples[i];
    }

    // Normalize the Gamma samples to obtain Dirichlet proportions
    for (i = 0; i < 3; i++) {
        x[i] = gamma_samples[i] / sum;
    }
}
