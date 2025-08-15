#!/bin/bash

###############################################################################
#Script Description
###############################################################################
#
#Description---------------------------
# Shell script that produces GUIs to ease the choice of parameters for run_local.sh
#
#Argument/s----------------------------
# none
#
#Output--------------------------------
# none directly, but calls simulation/src/bin/haics which does produce output files
#
#Author/s------------------------------
# Felix Muller
# Jorge Perez Heredia
#
#------------------------------------------------------------------------------


###############################################################################
#Global stuff
###############################################################################

#set colors for special printing
colRed='\e[1;37;41m'
colNormal='\033[0m'
colYellow='\033[1;37;43m'

function fnExitError {
  echo -e "\n${colRed}ERROR in ${0##*/}${colNormal}: ${1} \nExiting...\n" && exit 1
}

function fnPrintWarning {
  echo -e "\n${colYellow}WARNING in ${0##*/}${colNormal}: ${1} \n"
}

title="Haics--J&F-Edition"

dirHaicsSim=$(pwd)/../../simulation
dirLocalOutput=../../output
dirLocalTrashbin=${dirLocalOutput}/LOCAL_HAICS_TRASH
mkdir -p ${dirLocalTrashbin}

###############################################################################
#Reading argument (tuning or sampling)
###############################################################################

method=${1}

model=${2}

data=${3}

ID=${4}

REUSE=${5}

RERUN=${6}

RunParallel=${7}

ATune=${8}

seed=0

echo "Entering set_up_local_inputfile.sh"

if [ ${RERUN} -eq 1 ]; then

  #move previous files to trashbin
  ID_trashbin=${ID}_${RANDOM}
  echo -e "Moving previous folder to ${dirLocalTrashbin}/${ID_trashbin}"
  cp -r "${dirLocalOutput}/${ID}" "${dirLocalTrashbin}/${ID_trashbin}"
  
  dirReRunSimulation=${dirLocalOutput}/${ID}

  #find the number of chains
  cd ${dirReRunSimulation}
  INT_RUNS=$(ls -dq *${ID}* | wc -l)

  #clean all the folders of the old simulation but the input one
  find ./ -mindepth 1 -maxdepth 1 -type d ! -name "input" -exec rm -r {} +

  echo -e "\nThere are ${INT_RUNS} chains."

  # launch run_local.sh
  (cd ../../scripts/running/ && ./run_local.sh ${ID} ${INT_RUNS} ${REUSE} ${RunParallel})

  echo -e "Reached the end of ${0##*/}"
  exit

fi

#Assign variables that are method- and/or model-dependent

Phi=1
t_Phi=0
num_comp_E=0
AR_target=0.92
delta_AR_target=0.01
iter_tune=200

if [ "$ATune" == "yes" ]; then

    fitting_factor_approach="Gupta"
    integrator="3sAIA-HMC"
    opt_param_setting=1
    L=6
    t_L=2
    stepsize=0.0001
    t_stepsize=0
    stepsize_delta=0.00001
    scaling_value=1

    #yad window: number of iterations for burn-in and production
    yad_iter_ATune=$(yad --title "Number of iterations for burn-in and production" \
    --center --width=500 --height=100 --form --separator=' ' \
    --field="Number of burn-in iterations:" "1000" \
    --field="Number of sampling iterations:" "5000" \
    --field="Number of chains:" "4" \
    --field="thinning" "1" \
    --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

    #Read yad GUI
    yad_iter_sampling=(${yad_iter_ATune}) #this makes it a ZERO-INDEXED array

    iter_burn_in="${yad_iter_sampling[0]}"
    iter_sampling="${yad_iter_sampling[1]}"
    INT_RUNS="${yad_iter_sampling[2]}"
    thinning="${yad_iter_sampling[3]}"

    #mkdir -p ../../output/${ID}/input/
    #(cd ../sampling_optimal_setting/ && ./sampling_optimal_setting.sh ${iter_sampling} ${ID} ${REUSE} ${RunParallel})

    #echo -e "Reached the end of ${0##*/}"
    #exit

elif [ "$ATune" == "no" ]; then

    opt_param_setting=0

    #open yad GUI
    config=$(yad --title "HMC general setting" \
    --center --width=500 --height=350 --form --separator=' ' \
    --field="Number of burn-in iterations:" "1000" \
    --field="Number of sampling iterations:" "5000" \
    --field="integrator:CB" "1sVV!2sVV!2sBCSS-HMC!2sMinError-HMC!2sAIA-HMC!3sBCSS-HMC!3sMinError-HMC!^3sAIA-HMC" \
    --field="trajectory length L" "1" \
    --field="type of L:CB" "^constant!unif(1,L)" \
    --field="stepsize" "0.001" \
    --field="type of stepsize:CB" "^constant!N(stepsize,0.0025*stepsize^2)!unif(stepsize-stepsize_delta,stepsize)!unif(stepsize,stepsize+stepsize_delta)!unif(stepsize-stepsize_delta,stepsize+stepsize_delta)" \
    --field="thinning" "1" \
    --field="stepsize delta" "0.0001"  \
    --field="fitting factor (for s-AIA)" "1" \
    --field="Number of chains:" "4" \
    --button="gtk-cancel:1" \
    --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

    #Read yad GUI
    config=(${config}) #this makes it a ZERO-INDEXED array

    integrator="${config[0]}"
    L="${config[1]}"
    t_L="${config[2]}"
    stepsize="${config[3]}"
    t_stepsize="${config[4]}"
    thinning="${config[5]}"
    stepsize_delta="${config[6]}"
    scaling_value="${config[7]}"

    case $t_L in
        "constant")	t_L=0;;
        "unif(1,L)")	t_L=1;;
        *) fnExitError "No match found in case statement for t_L (which is '${t_L}')" ;;
    esac

    case $t_stepsize in
        "constant")	t_stepsize=0;;
        "N(stepsize,0.0025*stepsize^2)") t_stepsize=1;;
        "unif(stepsize-stepsize_delta,stepsize)")	t_stepsize=20;;
        "unif(stepsize,stepsize+stepsize_delta)") t_stepsize=21;;
        "unif(stepsize-stepsize_delta,stepsize+stepsize_delta)") t_stepsize=3;;
        *) fnExitError "No match found in case statement for t_stepsize (which is '${t_stepsize}')" ;;
    esac

    if [ $method == "GHMC" ]; then
        #open yad GUI
        GHMC_params=$(yad --title "Extra parameters for $method" \
        --center --width=500 --height=150 --form --separator=' ' \
        --field="noise parameter (Phi) (0,1]" "0.5" \
        --field="Type of Phi:CB" "^constant!unif(0,Phi)!OptimalPhi" \
        --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

        #read yad GUI
        GHMC_params=(${GHMC_params}) #this makes it a ZERO-INDEXED array

        Phi="${GHMC_params[0]}"
        t_Phi="${GHMC_params[1]}"

        #check that we got as many elements as we asked
        if [ ! -z ${extra_item} ]; then
            fnExitError "Problem with the number of items introduced. Remember that white spaces are not allowed."
        fi

        case $t_Phi in
            "constant")	t_Phi=0;;
            "unif(0,Phi)")	t_Phi=1;;
            "OptimalPhi") t_Phi=2;;
            *) fnExitError "No match found in case statement for t_Phi (which is '${t_Phi}')" ;;
        esac
    else
        Phi=1
        t_Phi=0
    fi #end of GHMC method
fi #end of ATune if-elif

# Helper function to set the label for the first field based on prior selection
function get_prior_field_label_1() {
    case "$1" in
        "Normal") echo "Mean";;
        "Gamma") echo "k";;
        "Exponential") echo "Lambda (1/Mean)";;
        "Inverse_Gamma") echo "Shape";; # Ensure correct match for "Inverse Gamma"
        "Uniform") echo "a";;
        "Beta") echo "a";;
        "Kumaraswamy") echo "a";;
        "Weibull") echo "k";;
        "Lorentz") echo "x0";;
        "Chi_Squared") echo "k";;
        *) echo "Field 1";;
    esac
}

# Helper function to set the label for the second field based on prior selection
function get_prior_field_label_2() {
    case "$1" in
        "Normal") echo "Standard Deviation";;
        "Gamma") echo "theta";;
        "Exponential") echo " ";; # No second field for Exponential
        "Inverse_Gamma") echo "Scale";; # Ensure correct match for "Inverse Gamma"
        "Uniform") echo "b";;
        "Beta") echo "b";;
        "Kumaraswamy") echo "b";;
        "Weibull") echo "lambda";;
        "Lorentz") echo "gamma";;
        "Chi_Squared") echo " ";;
        *) echo "Field 2";;
    esac
}

working_directory=$(pwd)
# Count the number of lines in the file
line_count=$(wc -l < "${working_directory%/*/*}/simulation/benchmarks/"$data"_data.txt")
line_count=$((line_count - 1))

# # Print the result
# echo "data variable: $data"
# echo "Number of rows: $line_count"

initialCountsFixed=0
S0_fixed_value=2185595
E0_fixed_value=10
I0_fixed_value=0
dirich_alpha_0=999994.424608
dirich_alpha_1=4.575392
dirich_alpha_2=1
dirich_alpha_3=1

if [ "$model" = "SIKR_Incidence" ] || [ "$model" = "SEMIKR_Incidence" ]; then

    ################################################################
    # The windows to load the parameters for the compartmental model
    ################################################################

    # First window to select the prior distributions for I0, Gamma, Phi Inv, and Tau
    function select_priors() {
        prior_selections=$(yad --title "Select Priors for Parameters" \
        --center --width=500 --height=350 --form --separator=' ' \
        --field="Number of spline basis" "15" \
        --field="Degree of spline polynomial" "3" \
        --field="Final time for the spline basis" "100" \
        --field="Time points returned by CVODES" "100" \
        --field="Number of infectious compartments" "1" \
        --field="Is gamma fixed?:CB" "No!Yes" \
        --field="Is gamma bounded?:CB" "No!Yes" \
        --field="Fixed value for Initial number of Susceptibles (S0)" "2185595" \
        --field="Fixed value for Initial number of Recovered (R0):CB" "10" \
        --field="Prior for gamma:CB" "Normal!Gamma!Exponential!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --field="Prior for overdispersion (phi inv):CB" "Exponential!Normal!Gamma!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --field="Prior for spline roughness (tau):CB" "Inverse_Gamma!Normal!Gamma!Exponential!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

        echo "$prior_selections"
    } # A prior on R0 is the same as a prior on S0 since S0 = N0 - (E0 + I0 + R0). The prior on S0 is dependent on the values of E0 and I0, in order to fulfill S0 >= E0 + I0.

    # Second window to request parameters based on the selected priors
    function set_prior_parameters() {
        local S0_fixed_value="$1"
        local prior_I0="$2"
        local prior_gamma="$3"
        local prior_phi="$4"
        local prior_tau="$5"
        local gamma_fixed="$6"
        local gamma_bounded="$7"

        # Dynamically build the form based on the gamma logic
        if [ "$gamma_fixed" == "Yes" ]; then
            # If gamma is fixed, show only one field for its value
            parameter_entries=$(yad --title "Set Prior Parameters" \
            --center --width=500 --height=350 --form --separator=' ' \
            --field="R0--$prior_S0--$(get_prior_field_label_1 "$prior_S0")" "10" \
            --field="R0--$prior_S0--$(get_prior_field_label_2 "$prior_S0")" "1" \
            --field="I0--$prior_I0--$(get_prior_field_label_1 "$prior_I0")" "21.88017" \
            --field="I0--$prior_I0--$(get_prior_field_label_2 "$prior_I0")" "7.29339" \
            --field="Fixed Gamma Value" "0.1" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_1 "$prior_phi")" "10" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_2 "$prior_phi")" "0" \
            --field="tau--$prior_tau--$(get_prior_field_label_1 "$prior_tau")" "1.0" \
            --field="tau--$prior_tau--$(get_prior_field_label_2 "$prior_tau")" "0.005" \
            --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."
        elif [ "$gamma_bounded" == "Yes" ]; then
            # If gamma is bounded, show two fields for bounds and a prior
            parameter_entries=$(yad --title "Set Prior Parameters" \
            --center --width=500 --height=350 --form --separator=' ' \
            --field="R0--$prior_S0--$(get_prior_field_label_1 "$prior_S0")" "10" \
            --field="R0--$prior_S0--$(get_prior_field_label_2 "$prior_S0")" "1" \
            --field="I0--$prior_I0--$(get_prior_field_label_1 "$prior_I0")" "21.88017" \
            --field="I0--$prior_I0--$(get_prior_field_label_2 "$prior_I0")" "7.29339" \
            --field="Gamma Lower Bound" "0.03333333333333" \
            --field="Gamma Upper Bound" "1.0" \
            --field="gamma--$prior_gamma--$(get_prior_field_label_1 "$prior_gamma")" "0.1" \
            --field="gamma--$prior_gamma--$(get_prior_field_label_2 "$prior_gamma")" "0.05" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_1 "$prior_phi")" "10" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_2 "$prior_phi")" "0" \
            --field="tau--$prior_tau--$(get_prior_field_label_1 "$prior_tau")" "1.0" \
            --field="tau--$prior_tau--$(get_prior_field_label_2 "$prior_tau")" "0.005" \
            --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."
        else
            # If gamma is not fixed and not bounded, show only the prior
            parameter_entries=$(yad --title "Set Prior Parameters" \
            --center --width=500 --height=350 --form --separator=' ' \
            --field="R0--$prior_S0--$(get_prior_field_label_1 "$prior_S0")" "10" \
            --field="R0--$prior_S0--$(get_prior_field_label_2 "$prior_S0")" "1" \
            --field="I0--$prior_I0--$(get_prior_field_label_1 "$prior_I0")" "21.88017" \
            --field="I0--$prior_I0--$(get_prior_field_label_2 "$prior_I0")" "7.29339" \
            --field="gamma--$prior_gamma--$(get_prior_field_label_1 "$prior_gamma")" "0.1" \
            --field="gamma--$prior_gamma--$(get_prior_field_label_2 "$prior_gamma")" "0.05" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_1 "$prior_phi")" "10" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_2 "$prior_phi")" "0" \
            --field="tau--$prior_tau--$(get_prior_field_label_1 "$prior_tau")" "1.0" \
            --field="tau--$prior_tau--$(get_prior_field_label_2 "$prior_tau")" "0.005" \
            --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."
        fi

        echo "$parameter_entries"
    }

    # Main script execution
    prior_result=$(select_priors)

    # Read the selected priors into an array
    priors=(${prior_result})

    # Extract the values for gamma logic
    gamma_fixed="${priors[5]}"
    gamma_bounded="${priors[6]}"

    # Display second window based on prior selections and gamma logic
    parameters_result=$(set_prior_parameters "${priors[7]}" "${priors[8]}" "${priors[9]}" "${priors[10]}" "${priors[11]}" "$gamma_fixed" "$gamma_bounded")

    # Read the parameter entries into an array
    parameters=(${parameters_result})

    # Display final output for debugging
    echo "Final selection of priors and parameters:"
    echo "R0 Prior: ${priors[7]}, Param 1: ${parameters[0]}, Param 2: ${parameters[1]}"
    echo "I0 Prior: ${priors[8]}, Param 1: ${parameters[2]}, Param 2: ${parameters[3]}"
    if [ "$gamma_fixed" == "Yes" ]; then
        echo "Gamma is fixed, value: ${parameters[4]}"
        echo "Phi Inv Prior: ${priors[10]}, Param 1: ${parameters[5]}, Param 2: ${parameters[6]}"
        echo "Tau Prior: ${priors[11]}, Param 1: ${parameters[7]}, Param 2: ${parameters[8]}"
    elif [ "$gamma_bounded" == "Yes" ]; then
        echo "Gamma Prior: ${priors[9]}, Lower Bound: ${parameters[4]}, Upper Bound: ${parameters[5]}"
        echo "Gamma Param 1: ${parameters[6]}, Param 2: ${parameters[7]}"
        echo "Phi Inv Prior: ${priors[10]}, Param 1: ${parameters[8]}, Param 2: ${parameters[9]}"
        echo "Tau Prior: ${priors[11]}, Param 1: ${parameters[10]}, Param 2: ${parameters[11]}"
    else
        echo "Gamma Prior: ${priors[9]}, Param 1: ${parameters[4]}, Param 2: ${parameters[5]}"
        echo "Phi Inv Prior: ${priors[10]}, Param 1: ${parameters[6]}, Param 2: ${parameters[7]}"
        echo "Tau Prior: ${priors[11]}, Param 1: ${parameters[8]}, Param 2: ${parameters[9]}"
    fi

    S0_param_1="${parameters[0]}"
    S0_param_2="${parameters[1]}"
    I0_param_1="${parameters[2]}"
    I0_param_2="${parameters[3]}"
    if [ "$gamma_fixed" == "Yes" ]; then
        is_gamma_fixed=1
        gamma_fixed_value="${parameters[4]}"
        gamma_param_1=0
        gamma_param_2=0
        is_gamma_bounded=1
        gammaLower="${parameters[4]}"
        gammaUpper="${parameters[4]}"
        phi_inv_param_1="${parameters[5]}"
        phi_inv_param_2="${parameters[6]}"
        tau_param_1="${parameters[7]}"
        tau_param_2="${parameters[8]}"
    else
        gamma_fixed_value=0
        is_gamma_fixed=0
        if [ "$gamma_bounded" == "Yes" ]; then
            is_gamma_bounded=1
            gammaLower="${parameters[4]}"
            gammaUpper="${parameters[5]}"
            gamma_param_1="${parameters[6]}"
            gamma_param_2="${parameters[7]}"
            phi_inv_param_1="${parameters[8]}"
            phi_inv_param_2="${parameters[9]}"
            tau_param_1="${parameters[10]}"
            tau_param_2="${parameters[11]}"
        else
            gammaLower=0
            gammaUpper=10**50
            is_gamma_bounded=0
            gamma_param_1="${parameters[4]}"
            gamma_param_2="${parameters[5]}"
            phi_inv_param_1="${parameters[6]}"
            phi_inv_param_2="${parameters[7]}"
            tau_param_1="${parameters[8]}"
            tau_param_2="${parameters[9]}"
        fi
    fi

    if [ "${priors[9]}" == "Normal" ]; then
        gamma_prior="PRIOR_NORMAL"
    elif [ "${priors[9]}" == "Gamma" ]; then
        gamma_prior="PRIOR_GAMMA"
    elif [ "${priors[9]}" == "Exponential" ]; then
        gamma_prior="PRIOR_EXPONENTIAL"
    elif [ "${priors[9]}" == "Inverse_Gamma" ]; then
        gamma_prior="PRIOR_INVGAMMA"
    elif [ "${priors[9]}" == "Uniform" ]; then
        gamma_prior="PRIOR_UNIFORM"
    elif [ "${priors[9]}" == "Beta" ]; then
        gamma_prior="PRIOR_BETA"
    elif [ "${priors[9]}" == "Kumaraswamy" ]; then
        gamma_prior="PRIOR_KUMARASWAMY"
    elif [ "${priors[9]}" == "Weibull" ]; then
        gamma_prior="PRIOR_WEIBULL"
    elif [ "${priors[9]}" == "Lorentz" ]; then
        gamma_prior="PRIOR_LORENTZ"
    elif [ "${priors[9]}" == "Chi_Squared" ]; then
        gamma_prior="PRIOR_CHISQUARED"
    fi

    if [ "${priors[10]}" == "Normal" ]; then
        phi_inv_prior="PRIOR_NORMAL"
    elif [ "${priors[10]}" == "Gamma" ]; then
        phi_inv_prior="PRIOR_GAMMA"
    elif [ "${priors[10]}" == "Exponential" ]; then
        phi_inv_prior="PRIOR_EXPONENTIAL"
    elif [ "${priors[10]}" == "Inverse_Gamma" ]; then
        phi_inv_prior="PRIOR_INVGAMMA"
    elif [ "${priors[10]}" == "Uniform" ]; then
        phi_inv_prior="PRIOR_UNIFORM"
    elif [ "${priors[10]}" == "Beta" ]; then
        phi_inv_prior="PRIOR_BETA"
    elif [ "${priors[10]}" == "Kumaraswamy" ]; then
        phi_inv_prior="PRIOR_KUMARASWAMY"
    elif [ "${priors[10]}" == "Weibull" ]; then
        phi_inv_prior="PRIOR_WEIBULL"
    elif [ "${priors[10]}" == "Lorentz" ]; then
        phi_inv_prior="PRIOR_LORENTZ"
    elif [ "${priors[10]}" == "Chi_Squared" ]; then
        phi_inv_prior="PRIOR_CHISQUARED"
    fi

    if [ "${priors[11]}" == "Normal" ]; then
        tau_prior="PRIOR_NORMAL"
    elif [ "${priors[11]}" == "Gamma" ]; then
        tau_prior="PRIOR_GAMMA"
    elif [ "${priors[11]}" == "Exponential" ]; then
        tau_prior="PRIOR_EXPONENTIAL"
    elif [ "${priors[11]}" == "Inverse_Gamma" ]; then
        tau_prior="PRIOR_INVGAMMA"
    elif [ "${priors[11]}" == "Uniform" ]; then
        tau_prior="PRIOR_UNIFORM"
    elif [ "${priors[11]}" == "Beta" ]; then
        tau_prior="PRIOR_BETA"
    elif [ "${priors[11]}" == "Kumaraswamy" ]; then
        tau_prior="PRIOR_KUMARASWAMY"
    elif [ "${priors[11]}" == "Weibull" ]; then
        tau_prior="PRIOR_WEIBULL"
    elif [ "${priors[11]}" == "Lorentz" ]; then
        tau_prior="PRIOR_LORENTZ"
    elif [ "${priors[11]}" == "Chi_Squared" ]; then
        tau_prior="PRIOR_CHISQUARED"
    fi

    num_basis_spline="${priors[0]}"
    spline_pol_degree="${priors[1]}"
    tfin_spl_bas="${priors[2]}"
    num_out_cvodes="${priors[3]}"
    num_comp="${priors[4]}"

    function select_priors_SEMIKR() {
        prior_selections=$(yad --title "Select Priors for Parameters" \
        --center --width=500 --height=350 --form --separator=' ' \
        --field="Number of exposed compartments" "1" \
        --field="Prior for alpha:CB" "Normal!Gamma!Exponential!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --field="Prior for Initial number of exposed (E0):CB" "Normal!Gamma!Exponential!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

        echo "$prior_selections"
    }

    # Second window to request parameters based on the selected priors
    function set_prior_parameters_SEMIKR() {
        local prior_alpha="$1"
        local prior_E0="$2"
            
        parameter_entries=$(yad --title "Set Prior Parameters" \
        --center --width=500 --height=350 --form --separator=' ' \
        --field="alpha--$prior_alpha--$(get_prior_field_label_1 "$prior_alpha")" "0.5" \
        --field="alpha--$prior_alpha--$(get_prior_field_label_2 "$prior_alpha")" "0.05" \
        --field="E0--$prior_E0--$(get_prior_field_label_1 "$prior_E0")" "21.88017" \
        --field="E0--$prior_E0--$(get_prior_field_label_2 "$prior_E0")" "7.29339" \
        --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

        echo "$parameter_entries"
    }

    if [ "$model" == "SEMIKR_Incidence" ]; then
        # Main script execution
        prior_result_SEMIKR=$(select_priors_SEMIKR)

        # Read the selected priors into an array
        priors_SEMIKR=(${prior_result_SEMIKR})

        # Display second window based on prior selections and gamma logic
        parameters_result_SEMIKR=$(set_prior_parameters_SEMIKR "${priors_SEMIKR[1]}" "${priors_SEMIKR[2]}")

        # Read the parameter entries into an array
        parameters_SEMIKR=(${parameters_result_SEMIKR})

        # Display final output for debugging
        echo "Final selection of priors and parameters:"
        echo "Alpha Prior: ${priors_SEMIKR[1]}, Param 1: ${parameters_SEMIKR[0]}, Param 2: ${parameters_SEMIKR[1]}"
        echo "E0 Prior: ${priors_SEMIKR[2]}, Param 1: ${parameters_SEMIKR[2]}, Param 2: ${parameters_SEMIKR[3]}"


        if [ "${priors_SEMIKR[1]}" == "Normal" ]; then
            alpha_prior="PRIOR_NORMAL"
        elif [ "${priors_SEMIKR[1]}" == "Gamma" ]; then
            alpha_prior="PRIOR_GAMMA"
        elif [ "${priors_SEMIKR[1]}" == "Exponential" ]; then
            alpha_prior="PRIOR_EXPONENTIAL"
        elif [ "${priors_SEMIKR[1]}" == "Inverse_Gamma" ]; then
            alpha_prior="PRIOR_INVGAMMA"
        elif [ "${priors_SEMIKR[1]}" == "Uniform" ]; then
            alpha_prior="PRIOR_UNIFORM"
        elif [ "${priors_SEMIKR[1]}" == "Beta" ]; then
            alpha_prior="PRIOR_BETA"
        elif [ "${priors_SEMIKR[1]}" == "Kumaraswamy" ]; then
            alpha_prior="PRIOR_KUMARASWAMY"
        elif [ "${priors_SEMIKR[1]}" == "Weibull" ]; then
            alpha_prior="PRIOR_WEIBULL"
        elif [ "${priors_SEMIKR[1]}" == "Lorentz" ]; then
            alpha_prior="PRIOR_LORENTZ"
        elif [ "${priors_SEMIKR[1]}" == "Chi_Squared" ]; then
            alpha_prior="PRIOR_CHISQUARED"
        fi

        alpha_param_1="${parameters_SEMIKR[0]}"
        alpha_param_2="${parameters_SEMIKR[1]}"
        E0_param_1="${parameters_SEMIKR[2]}"
        E0_param_2="${parameters_SEMIKR[3]}"
        num_comp_E="${priors_SEMIKR[0]}"
    else
        # Not relevant for SIR spline based models
        alpha_prior="PRIOR_NORMAL"
        alpha_param_1=0
        alpha_param_2=0
        num_comp_E=0
    fi #end of SEMIKR_Incidence model

    # Not relevant for spline based models
    beta_prior="PRIOR_NORMAL"
    beta_param_1=0
    beta_param_2=0
fi

if [ "$model" == "SIR_standard_Incidence" ]; then
    # First window to select the prior distributions for I0, Gamma, Phi Inv, and Tau
    function select_priors_stand_SIR() {
        prior_selections=$(yad --title "Select Priors for Parameters" \
        --center --width=500 --height=350 --form --separator=' ' \
        --field="Final time for the spline basis" "100" \
        --field="Time points returned by CVODES" "100" \
        --field="Is gamma fixed?:CB" "No!Yes" \
        --field="Is gamma bounded?:CB" "No!Yes" \
        --field="Prior for Initial number of susceptible (S0):CB" "Normal!Gamma!Exponential!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --field="Prior for Initial number of infected (I0):CB" "Normal!Gamma!Exponential!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --field="Prior for gamma:CB" "Normal!Gamma!Exponential!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --field="Prior for overdispersion (phi inv):CB" "Exponential!Normal!Gamma!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --field="Prior for beta:CB" "Normal!Gamma!Exponential!Inverse_Gamma!Uniform!Beta!LogNormal!Weibull!Kumaraswamy!Lorentz!Chi_Squared" \
        --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

        echo "$prior_selections"
    }

    # Second window to request parameters based on the selected priors
    function set_prior_parameters_stand_SIR() {
        local prior_S0="$1"
        local prior_I0="$2"
        local prior_gamma="$3"
        local prior_phi="$4"
        local prior_beta="$5"
        local gamma_fixed="$6"
        local gamma_bounded="$7"

        # Dynamically build the form based on the gamma logic
        if [ "$gamma_fixed" == "Yes" ]; then
            # If gamma is fixed, show only one field for its value
            parameter_entries=$(yad --title "Set Prior Parameters" \
            --center --width=500 --height=350 --form --separator=' ' \
            --field="S0--$prior_S0--$(get_prior_field_label_1 "$prior_S0")" "10" \
            --field="S0--$prior_S0--$(get_prior_field_label_2 "$prior_S0")" "1" \
            --field="I0--$prior_I0--$(get_prior_field_label_1 "$prior_I0")" "21.88017" \
            --field="I0--$prior_I0--$(get_prior_field_label_2 "$prior_I0")" "7.29339" \
            --field="Fixed Gamma Value" "1.0" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_1 "$prior_phi")" "10" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_2 "$prior_phi")" "0" \
            --field="beta--$prior_beta--$(get_prior_field_label_1 "$prior_beta")" "0.5" \
            --field="beta--$prior_beta--$(get_prior_field_label_2 "$prior_beta")" "0.05" \
            --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."
        elif [ "$gamma_bounded" == "Yes" ]; then
            # If gamma is bounded, show two fields for bounds and a prior
            parameter_entries=$(yad --title "Set Prior Parameters" \
            --center --width=500 --height=350 --form --separator=' ' \
            --field="S0--$prior_S0--$(get_prior_field_label_1 "$prior_S0")" "10" \
            --field="S0--$prior_S0--$(get_prior_field_label_2 "$prior_S0")" "1" \
            --field="I0--$prior_I0--$(get_prior_field_label_1 "$prior_I0")" "21.88017" \
            --field="I0--$prior_I0--$(get_prior_field_label_2 "$prior_I0")" "7.29339" \
            --field="Gamma Lower Bound" "0.03333333333333" \
            --field="Gamma Upper Bound" "1.0" \
            --field="gamma--$prior_gamma--$(get_prior_field_label_1 "$prior_gamma")" "0.1" \
            --field="gamma--$prior_gamma--$(get_prior_field_label_2 "$prior_gamma")" "0.05" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_1 "$prior_phi")" "10" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_2 "$prior_phi")" "0" \
            --field="beta--$prior_beta--$(get_prior_field_label_1 "$prior_beta")" "0.5" \
            --field="beta--$prior_beta--$(get_prior_field_label_2 "$prior_beta")" "0.05" \
            --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."
        else
            # If gamma is not fixed and not bounded, show only the prior
            parameter_entries=$(yad --title "Set Prior Parameters" \
            --center --width=500 --height=350 --form --separator=' ' \
            --field="S0--$prior_S0--$(get_prior_field_label_1 "$prior_S0")" "10" \
            --field="S0--$prior_S0--$(get_prior_field_label_2 "$prior_S0")" "1" \
            --field="I0--$prior_I0--$(get_prior_field_label_1 "$prior_I0")" "21.88017" \
            --field="I0--$prior_I0--$(get_prior_field_label_2 "$prior_I0")" "7.29339" \
            --field="gamma--$prior_gamma--$(get_prior_field_label_1 "$prior_gamma")" "0.1" \
            --field="gamma--$prior_gamma--$(get_prior_field_label_2 "$prior_gamma")" "0.05" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_1 "$prior_phi")" "10" \
            --field="phiInv--$prior_phi--$(get_prior_field_label_2 "$prior_phi")" "0" \
            --field="beta--$prior_beta--$(get_prior_field_label_1 "$prior_beta")" "0.5" \
            --field="beta--$prior_beta--$(get_prior_field_label_2 "$prior_beta")" "0.05" \
            --button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."
        fi

        echo "$parameter_entries"
    }

    # Main script execution
    prior_result_stand_SIR=$(select_priors_stand_SIR)

    # Read the selected priors into an array
    priors_stand_SIR=(${prior_result_stand_SIR})

    # Extract the values for gamma logic
    gamma_fixed="${priors_stand_SIR[2]}"
    gamma_bounded="${priors_stand_SIR[3]}"

    # Display second window based on prior selections and gamma logic
    parameters_result_stand_SIR=$(set_prior_parameters_stand_SIR "${priors_stand_SIR[4]}" "${priors_stand_SIR[5]}" "${priors_stand_SIR[6]}" "${priors_stand_SIR[7]}" "${priors_stand_SIR[8]}" "$gamma_fixed" "$gamma_bounded")

    # Read the parameter entries into an array
    parameters_stand_SIR=(${parameters_result_stand_SIR})

    if [ "$gamma_fixed" == "Yes" ]; then
        is_gamma_fixed=1
        gamma_fixed_value="${parameters_stand_SIR[4]}"
        gamma_param_1=0
        gamma_param_2=0
        is_gamma_bounded=1
        gammaLower="${parameters_stand_SIR[4]}"
        gammaUpper="${parameters_stand_SIR[4]}"
        phi_inv_param_1="${parameters_stand_SIR[5]}"
        phi_inv_param_2="${parameters_stand_SIR[6]}"
        beta_param_1="${parameters_stand_SIR[7]}"
        beta_param_2="${parameters_stand_SIR[8]}"
        tau_param_1=0
        tau_param_2=0
    else
        gamma_fixed_value=0
        is_gamma_fixed=0
        if [ "$gamma_bounded" == "Yes" ]; then
            is_gamma_bounded=1
            gammaLower="${parameters_stand_SIR[4]}"
            gammaUpper="${parameters_stand_SIR[5]}"
            gamma_param_1="${parameters_stand_SIR[6]}"
            gamma_param_2="${parameters_stand_SIR[7]}"
            phi_inv_param_1="${parameters_stand_SIR[8]}"
            phi_inv_param_2="${parameters_stand_SIR[9]}"
            beta_param_1="${parameters_stand_SIR[10]}"
            beta_param_2="${parameters_stand_SIR[11]}"
            tau_param_1=0
            tau_param_2=0
        else
            gammaLower=0
            gammaUpper=10**50
            is_gamma_bounded=0
            gamma_param_1="${parameters_stand_SIR[4]}"
            gamma_param_2="${parameters_stand_SIR[5]}"
            phi_inv_param_1="${parameters_stand_SIR[6]}"
            phi_inv_param_2="${parameters_stand_SIR[7]}"
            beta_param_1="${parameters_stand_SIR[8]}"
            beta_param_2="${parameters_stand_SIR[9]}"
            tau_param_1=0
            tau_param_2=0
        fi
    fi

    if [ "${priors_stand_SIR[6]}" == "Normal" ]; then
        gamma_prior="PRIOR_NORMAL"
    elif [ "${priors_stand_SIR[6]}" == "Gamma" ]; then
        gamma_prior="PRIOR_GAMMA"
    elif [ "${priors_stand_SIR[6]}" == "Exponential" ]; then
        gamma_prior="PRIOR_EXPONENTIAL"
    elif [ "${priors_stand_SIR[6]}" == "Inverse_Gamma" ]; then
        gamma_prior="PRIOR_INVGAMMA"
    elif [ "${priors_stand_SIR[6]}" == "Uniform" ]; then
        gamma_prior="PRIOR_UNIFORM"
    elif [ "${priors_stand_SIR[6]}" == "Beta" ]; then
        gamma_prior="PRIOR_BETA"
    elif [ "${priors_stand_SIR[6]}" == "Kumaraswamy" ]; then
        gamma_prior="PRIOR_KUMARASWAMY"
    elif [ "${priors_stand_SIR[6]}" == "Weibull" ]; then
        gamma_prior="PRIOR_WEIBULL"
    elif [ "${priors_stand_SIR[6]}" == "Lorentz" ]; then
        gamma_prior="PRIOR_LORENTZ"
    elif [ "${priors_stand_SIR[6]}" == "Chi_Squared" ]; then
        gamma_prior="PRIOR_CHISQUARED"
    fi

    if [ "${priors_stand_SIR[7]}" == "Normal" ]; then
        phi_inv_prior="PRIOR_NORMAL"
    elif [ "${priors_stand_SIR[7]}" == "Gamma" ]; then
        phi_inv_prior="PRIOR_GAMMA"
    elif [ "${priors_stand_SIR[7]}" == "Exponential" ]; then
        phi_inv_prior="PRIOR_EXPONENTIAL"
    elif [ "${priors_stand_SIR[7]}" == "Inverse_Gamma" ]; then
        phi_inv_prior="PRIOR_INVGAMMA"
    elif [ "${priors_stand_SIR[7]}" == "Uniform" ]; then
        phi_inv_prior="PRIOR_UNIFORM"
    elif [ "${priors_stand_SIR[7]}" == "Beta" ]; then
        phi_inv_prior="PRIOR_BETA"
    elif [ "${priors_stand_SIR[7]}" == "Kumaraswamy" ]; then
        phi_inv_prior="PRIOR_KUMARASWAMY"
    elif [ "${priors_stand_SIR[7]}" == "Weibull" ]; then
        phi_inv_prior="PRIOR_WEIBULL"
    elif [ "${priors_stand_SIR[7]}" == "Lorentz" ]; then
        phi_inv_prior="PRIOR_LORENTZ"
    elif [ "${priors_stand_SIR[7]}" == "Chi_Squared" ]; then
        phi_inv_prior="PRIOR_CHISQUARED"
    fi

    if [ "${priors_stand_SIR[8]}" == "Normal" ]; then
        beta_prior="PRIOR_NORMAL"
    elif [ "${priors_stand_SIR[8]}" == "Gamma" ]; then
        beta_prior="PRIOR_GAMMA"
    elif [ "${priors_stand_SIR[8]}" == "Exponential" ]; then
        beta_prior="PRIOR_EXPONENTIAL"
    elif [ "${priors_stand_SIR[8]}" == "Inverse_Gamma" ]; then
        beta_prior="PRIOR_INVGAMMA"
    elif [ "${priors_stand_SIR[8]}" == "Uniform" ]; then
        beta_prior="PRIOR_UNIFORM"
    elif [ "${priors_stand_SIR[8]}" == "Beta" ]; then
        beta_prior="PRIOR_BETA"
    elif [ "${priors_stand_SIR[8]}" == "Kumaraswamy" ]; then
        beta_prior="PRIOR_KUMARASWAMY"
    elif [ "${priors_stand_SIR[8]}" == "Weibull" ]; then
        beta_prior="PRIOR_WEIBULL"
    elif [ "${priors_stand_SIR[8]}" == "Lorentz" ]; then
        beta_prior="PRIOR_LORENTZ"
    elif [ "${priors_stand_SIR[8]}" == "Chi_Squared" ]; then
        beta_prior="PRIOR_CHISQUARED"
    fi

    num_out_cvodes="${priors_stand_SIR[1]}"
    # Not relevan in the Standard SIR model
    alpha_param_1=0
    alpha_param_2=0
    alpha_prior="PRIOR_INVGAMMA"
    tau_prior="PRIOR_INVGAMMA"
    num_basis_spline=0
    spline_pol_degree=0
    tfin_spl_bas=0
    num_comp=1
fi #end of SIR_standard_Incidence model

#enter Under Report
SIR_underReport=$(yad --title "Accounting for undereport" \
--center --width=500 --height=250 --form --separator=' ' \
--field="Account for undereport? 1 = Yes, 0 = No" "1" \
--field="Time of 1st measurement of under report" "92" \
--field="1st value of undereport (fraction of detected cases)" "0.15" \
--field="Time of 2nd measurement of under report" "281" \
--field="2nd value of undereport (fraction of detected cases)" "0.54" \
--button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

#read yad GUI
SIR_underReport=(${SIR_underReport}) #this makes it a ZERO-INDEXED array

do_under_report="${SIR_underReport[0]}"
T0_under="${SIR_underReport[1]}"
U0_under="${SIR_underReport[2]}"
T1_under="${SIR_underReport[3]}"
U1_under="${SIR_underReport[4]}"

#enter Gradient DESCENT
SIR_gradDesc=$(yad --title "Gradient descent parameters" \
--center --width=500 --height=150 --form --separator=' ' \
--field="Perform Gradient Descent? 1 = Yes, 0 = No" "0" \
--field="Learning rate" "1e-6" \
--field="Max number of iterations" "3000" \
--button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

#read yad GUI
SIR_gradDesc=(${SIR_gradDesc}) #this makes it a ZERO-INDEXED array

do_grad_desc="${SIR_gradDesc[0]}"
learning_rate="${SIR_gradDesc[1]}"
max_iter="${SIR_gradDesc[2]}"

#copy initialpoint and splinebasis files
SIR_initpoint=$(yad --title "Initialpoint and splinebasis filed" \
--center --width=500 --height=150 --form --separator=' ' \
--field="Initialpoint file direction" "../../simulation/input/example-initialpoint_SIR.txt" \
--field="Splinebasis file direction" "../../simulation/input/example-splinebasis_SIR.txt" \
--button="gtk-cancel:1" --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

#read yad GUI
SIR_initpoint=(${SIR_initpoint}) #this makes it a ZERO-INDEXED array

initialpoint_file="${SIR_initpoint[0]}"
spline_basis_file="${SIR_initpoint[1]}"

initialpoint_file_path_name="../../output/${ID}/input/initialpoint_${ID}.txt"
rm -rf $initialpoint_file_path_name
touch $initialpoint_file_path_name
cp "${initialpoint_file}" "${initialpoint_file_path_name}"

spline_basis_file_path_name="../../output/${ID}/input/splinebasis_${ID}.txt"
rm -rf $spline_basis_file_path_name
touch $spline_basis_file_path_name
cp "${spline_basis_file}" "${spline_basis_file_path_name}"

###############################################################################
#WRITE INPUTFILE FOR HAICS
###############################################################################
file_path_name="../../output/${ID}/input/inputfile_${ID}.txt"
rm -rf $file_path_name
touch $file_path_name
{
  echo "model $model"
  echo "data $data"
  echo "method $method"
  echo "seed $seed"
  echo "iter_sampling $iter_sampling"
  echo "iter_burn_in $iter_burn_in"
  echo "integrator $integrator"
  echo "t_L $t_L"
  echo "L $L"
  echo "t_stepsize $t_stepsize"
  echo "stepsize $stepsize"
  echo "thinning $thinning"
  echo "t_Phi $t_Phi"
  echo "Phi $Phi"
  echo "scaling_value ${scaling_value}"
  echo "stepsize_delta ${stepsize_delta}"
  echo "iter_tune ${iter_tune}"
  echo "AR_target ${AR_target}"
  echo "delta_AR_target ${delta_AR_target}"
  echo "fitting_factor_approach ${fitting_factor_approach}"
  echo "opt_param_setting ${opt_param_setting}"
  echo "num_basis_spline ${num_basis_spline}"
  echo "spline_pol_degree ${spline_pol_degree}"
  echo "tfin_spl_bas ${tfin_spl_bas}"
  echo "num_out_cvodes ${num_out_cvodes}"
  echo "initialCountsFixed ${initialCountsFixed}"
  echo "S0_fixed_value ${S0_fixed_value}"
  echo "E0_fixed_value ${E0_fixed_value}"
  echo "I0_fixed_value ${I0_fixed_value}"
  echo "dirich_alpha_0 ${dirich_alpha_0}"
  echo "dirich_alpha_1 ${dirich_alpha_1}"
  echo "dirich_alpha_2 ${dirich_alpha_2}"
  echo "dirich_alpha_3 ${dirich_alpha_3}"
  echo "alpha_prior ${alpha_prior}"
  echo "alpha_param_1 ${alpha_param_1}"
  echo "alpha_param_2 ${alpha_param_2}"
  echo "gamma_prior ${gamma_prior}"
  echo "gamma_param_1 ${gamma_param_1}"
  echo "gamma_param_2 ${gamma_param_2}"
  echo "phi_inv_prior ${phi_inv_prior}"
  echo "phi_inv_param_1 ${phi_inv_param_1}"
  echo "phi_inv_param_2 ${phi_inv_param_2}"
  echo "tau_prior ${tau_prior}"
  echo "tau_param_1 ${tau_param_1}"
  echo "tau_param_2 ${tau_param_2}"
  echo "num_comp ${num_comp}"
  echo "num_comp_E ${num_comp_E}"
  echo "do_under_report ${do_under_report}"
  echo "T0_under ${T0_under}"
  echo "U0_under ${U0_under}"
  echo "T1_under ${T1_under}"
  echo "U1_under ${U1_under}"
  echo "do_grad_desc ${do_grad_desc}"
  echo "learning_rate ${learning_rate}"
  echo "max_iter ${max_iter}"
  echo "gammaFixed ${is_gamma_fixed}"
  echo "gamma_fixed_value ${gamma_fixed_value}"
  echo "gammaBounded ${is_gamma_bounded}"
  echo "gammaUpper ${gammaUpper}"
  echo "gammaLower ${gammaLower}"
  echo "beta_prior ${beta_prior}"
  echo "beta_param_1 ${beta_param_1}"
  echo "beta_param_2 ${beta_param_2}"
} >> $file_path_name

###############################################################################
#RUN HAICS
###############################################################################

echo "Launching run_local.sh"
./run_local.sh ${ID} ${INT_RUNS} ${REUSE} ${RunParallel}

###############################################################################
#El fin
###############################################################################
echo -e "Reached the end of ${0##*/}"
