#!/bin/bash

###############################################################################
#Script Description
###############################################################################
#
#Description---------------------------
# Shell script that produces GUIs that walks the user through ALL the features that this version of the code supports.
# See https://gitlab.bcamath.org/mslms/haics/README.md for how to use this file.
#
#Arguments-----------------------------
# none, everything specified on the GUI.
#
#Output--------------------------------
# none directly, but calls other scripts that do produce output files
#
#Author/s------------------------------
# Felix Muller
# Jorge Perez Heredia
# Lorenzo Nagar
#------------------------------------------------------------------------------

echo -e "====================================================="
echo -e "Welcome to HaiCS"
echo -e "====================================================="

###############################################################################
#Global stuff
###############################################################################

echo -e "Initialising..."

#set colors for special printing
colRed='\e[1;37;41m'
colNormal='\033[0m'
colYellow='\033[1;37;43m'
colBlue='\e[1;37;44m'

function fnExitError {
  echo -e "\n${colRed}ERROR${colNormal} in ${0##*/}: ${1} \nExiting...\n" && exit 1
}

function fnPrintWarning {
  echo -e "\n${colYellow}WARNING${colNormal} in ${0##*/}: ${1} \n"
}

#main title for all the graphical windows
strYadTitle="HaiCS GUI by Felix & Jorge"

#get candidate for Hipatia's username
strUser=$(whoami)

#Local directories
dirLocalHaics=$(pwd)
dirLocalOutput=${dirLocalHaics}/output
mkdir -p ${dirLocalOutput}
dirLocalScripts=${dirLocalHaics}/scripts
dirLocalSims=${dirLocalHaics}/simulation
dirLocalBenchmarks=${dirLocalSims}/benchmarks/epidemics
mkdir -p ${dirLocalSims}/bin
mkdir -p ${dirLocalSims}/objects
dirLocalTrashbin=${dirLocalOutput}/LOCAL_HAICS_TRASH
mkdir -p ${dirLocalTrashbin}

echo -e "Initialisation done."

###############################################################################
#Main window of HaiCS
###############################################################################

echo -e "Reading what to do..."

strDataNames=$(find ${dirLocalBenchmarks}/*_data.txt -printf %f\\n | tr '\n' '!' | sed 's/.$//g' | sed 's/\_data.txt//g' | sed -e 's/,/!/g')

#open yad GUI for model, sampler and dataset
config=$(yad --title "HaiCS GUI" \
--center --form --separator=' ' --width=400 --height=200 \
--field="re-run old simulation or new?:CB" "old!^new" \
--field="simulation name:" "test" \
--field="sampler:CB" "HMC!^GHMC" \
--field="ATune parameters setting:CB" "^yes!no" \
--field="compartmental model:CB" "SIR_standard_incidence!SIR_standard_prevalence!^SIKR_spline!SEMIKR_spline" \
--field="data:CB" ${strDataNames} \
--field="Re-use the code?:CB" "yes!^no" \
--field="Run in parallel?:CB" "yes!^no" \
--button="gtk-cancel:1" \
--button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."

#--field="data:CB" "synthetic_data!COVID_Andalucía!COVID_Aragón!COVID_Asturias!COVID_Canarias!COVID_Cantabria!COVID_Castilla_La_Mancha!COVID_Castilla_y_León!COVID_Catalunya!COVID_Euskadi!COVID_Extremadura!COVID_Galicia!COVID_Illes_Balears!COVID_La_Rioja!COVID_Madrid!COVID_Murcia!COVID_Navarra!COVID_Comunitat_Valenciana!COVID_Ceuta!COVID_Melilla" \

#Read yad GUI
config=(${config}) #this makes it a ZERO-INDEXED array

RERUN="${config[0]}"
SimName="${config[1]}"
method="${config[2]}"
ATune="${config[3]}"
model="${config[4]}"
data="${config[5]}"
REUSE="${config[6]}"
RunParallel="${config[7]}"

case $RERUN in
  "old") RERUN=1;;
  "new") RERUN=0;;
  *) fnExitError "No match found in case statement for RERUN (which is '${RERUN}')" ;;
esac

case $model in
  "SIR_standard_incidence") model=SIR_standard_Incidence;;
  "SIR_standard_prevalence") model=SIR_standard_Prevalence;;
  "SIKR_spline") model=SIKR_Incidence;;
  "SEMIKR_spline") model=SEMIKR_Incidence;;
  *) fnExitError "No match found in case statement for model (which is '${model}')" ;;
esac

case $REUSE in
  "yes") REUSE=1;;
  "no") REUSE=0;;
  *) fnExitError "No match found in case statement for REUSE (which is '${REUSE}')" ;;
esac

case $RunParallel in
  "yes") RunParallel=1;;
  "no") RunParallel=0;;
  *) fnExitError "No match found in case statement for RunParallel (which is '${RunParallel}')" ;;
esac

if [ ${RERUN} -eq 1 ]; then

  #find the simulation to re-run
  echo -e "----------------------------------------"
  echo -e "\nChoose the simulation to re-run ..."
  
  #cd output/
  #strIDNames=$(find ${dirLocalOutput}/*/input/inputfile_* -printf %f\\n | tr '\n' '!' | sed 's/.$//g' | sed 's/inputfile_//g' | sed 's/\.txt//g' | sed -e 's/,/!/g')
  strIDNames=$(find output/ -maxdepth 1 -type d -printf %f\\n | sort | tr '\n' '!')
  #cd ../

  #open yad GUI
  yad=$(yad --title "Simulation to re-run:" \
  --text "Choose simulation from the output folder, path: \n ${dirLocalOutput} \n" \
  --center --width=400 --height=200 --form --separator=' ' \
  --field="Simulations:CB" ${strIDNames} \
  --button="gtk-cancel:1" \
  --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status on yad interface" }

  #Read yad GUI
  yad=(${yad}) #this makes it a ZERO-INDEXED array
  SimName="${yad[0]}"

  echo -e "\nThe simulation to re-run is ${SimName}"

elif [ ${RERUN} -eq 0 ]; then

  #Cleaning old simulations with the same ID
  #check if ID exists already locally. If yes, ask for confirmation to remove previous files. If not, go on.
  if [ -d ${dirLocalOutput}/${SimName} ]; then
      ID_trashbin=${SimName}_${RANDOM}

      yad=$(yad --title "Check simulation folder" \
      --center --width=400 --form --separator=' ' \
      --text "The given ID is already in use in the local output folder.\n\nClick Yes:\nPrevious files will be moved to \n'${dirLocalTrashbin}/${ID_trashbin}'.\n\nClick No:\nExit." \
      --button="gtk-cancel:1" \
      --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."
  fi

  #move previous files to trashbin
  echo -e "Moving previous folder to ${dirLocalTrashbin}/${ID_trashbin}"
  mv "${dirLocalOutput}/${SimName}" "${dirLocalTrashbin}/${ID_trashbin}"
  #create new folder
  mkdir -p ${dirLocalOutput}/${SimName}
  mkdir -p ${dirLocalOutput}/${SimName}/input

fi #end of if-esle for RERUN

echo "Launching set_up_local_inputfile.sh"
(cd scripts/running &&./set_up_local_inputfile.sh ${method} ${model} ${data} ${SimName} ${REUSE} ${RERUN} ${RunParallel} ${ATune})

###############################################################################
#El fin
###############################################################################
echo -e "Reached the end of ${0##*/}"
