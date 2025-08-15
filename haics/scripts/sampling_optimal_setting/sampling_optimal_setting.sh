#!/bin/bash

echo -e "======================================================================"
echo -e "Script to set the optimal parameters setting for sampling"
echo -e "======================================================================"

#Initialising------------------------------------------------------------------
echo -e "----------------------------------------"
echo -e "\nInitialising..."

#Reading arguments-------------------------------------------------------------
ITER_SAMPLING=${1}

ID_NEW=${2}

REUSE=${3}

RunParallel=${4}

#set colors for special printing
colRed='\e[0;37;41m'
colNormal='\033[0m'
colYellow='\033[0;37;43m'

function fnExitError {
    echo -e "\n${colRed}ERROR${colNormal} in ${0##*/}: ${1} \nExiting...\n" && exit 1
}

function fnPrintWarning {
  echo -e "\n${colYellow}WARNING${colNormal} in ${0##*/}: ${1} \n"
}

#main title for all the graphical windows
strYadTitle="Optimal parameters setting"

#get candidate for Hipatia's username
strUser=$(whoami)
#dirHaicsLocal=/home/${strUser}/haics
dirHaicsLocal=../..

#Local directories
#echo -e "I'm assumming that you have a copy of the haics repo at:\n\t ${pwd}/${dirHaicsLocal}"

dirOutputLocal=${dirHaicsLocal}/output

#Fetch Stepsizes-----------------------------------------------------------------
echo -e "----------------------------------------"
echo -e "\nReading the metrics from the simulation ..."
    
strIDNames=$(find ${dirOutputLocal}/*/input/inputfile_* -printf %f\\n | tr '\n' '!' | sed 's/.$//g' | sed 's/inputfile_//g' | sed 's/\.txt//g' | sed -e 's/,/!/g')

#open yad GUI
yad=$(yad --title "Burn-in simulation to starting from" \
--text "Choose simulation from the output folder, path: \n ${dirHaicsHipatia} \n" \
--center --width=400 --height=200 --form --separator=' ' \
--field="Simulations:CB" ${strIDNames} \
--button="gtk-cancel:1" \
--button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status on yad interface" }

#Read yad GUI
yad=(${yad}) #this makes it a ZERO-INDEXED array
ID_BATCH="${yad[0]}"

dirSimulation=${dirOutputLocal}/${ID_BATCH}
dirInput=${dirSimulation}/input

#find the number of runs
cd ${dirSimulation}
NumOfRuns=$(ls -dq *${ID_BATCH}* | wc -l)

echo -e "\nThere are ${NumOfRuns} chains."

#go to the scripts/sampling_optimal_setting folder
cd ${dirHaicsLocal}/scripts/sampling_optimal_setting

#Cleaning old data folder with the same ID
#check if data folder with ID_BATCH exists already. If yes, ask for confirmation to remove previous files. If not, go on.
if [ -d "data/${ID_BATCH}" ]; then

  yad=$(yad --title "Check data folder" \
  --center --width=400 --form --separator=' ' \
  --text "The given ID is already in use in the local output folder.\n\nClick Yes:\nPrevious files will be deleted'.\n\nClick No:\nExit." \
  --button="gtk-cancel:1" \
  --button="gtk-ok:0" 2> /dev/null) || fnExitError "non-zero exit status from the yad window."
fi

#create new folder
mkdir -p data/${ID_BATCH}

#take the HSL from the burn-in simulation
fileHsl=data/${ID_BATCH}/hsl-${ID_BATCH}.txt
grep -w 'hsl' ${dirSimulation}/${ID_BATCH}_*/art.txt > ${fileHsl}

#take the HSL_filtAR from the burn-in simulation
fileHslFiltAR=data/${ID_BATCH}/hsl_filtAR-${ID_BATCH}.txt
grep -w 'hsl_filtAR' ${dirSimulation}/${ID_BATCH}_*/art.txt > ${fileHslFiltAR}

#take the HSL_Gupta from the burn-in simulation
fileHslGupta=data/${ID_BATCH}/hsl_Gupta-${ID_BATCH}.txt
grep -w 'hsl_Gupta' ${dirSimulation}/${ID_BATCH}_*/art.txt > ${fileHslGupta}

#take the fitting factor from the burn-in simulation
fileFittingFactor=data/${ID_BATCH}/fitting_factor-${ID_BATCH}.txt
grep -w 'fitting_factor' ${dirSimulation}/${ID_BATCH}_*/art.txt > ${fileFittingFactor}

for i in $(seq 1 1 ${NumOfRuns})
do
  #copy the final points to be the initial points of the sampling stage
  cp ${dirSimulation}/${ID_BATCH}_${i}/finalpoint_burnin.txt data/${ID_BATCH}/finalpoint_${i}.txt
done

#inputfile
fileInputfile=${dirInput}/inputfile_${ID_BATCH}.txt
cp ${dirInput}/inputfile_${ID_BATCH}.txt data/${ID_BATCH}/inputfile_${ID_BATCH}.txt

#Make the inputfile---------------------------------------------------------------
echo -e "----------------------------------------"
echo -e "\nMake the inputfile for the sampling stage..."

Rscript sampling_optimal_setting.R ${ID_BATCH} ${NumOfRuns} ${ITER_SAMPLING} ${ID_NEW} || fnExitError "something wrong while executing sampling_optimal_setting.R"

echo -e "Done, please find the results in $(pwd)/data/${ID_BATCH}"

#Lanching simulations(s)----------------------------------------------------------
echo "Launching run_local.sh"
(cd ../running/ && ./run_local.sh ${ID_NEW} ${NumOfRuns} ${REUSE} ${RunParallel})

#The End------------------------------------------------------------------------
echo -e "\nReached the end of ${0##*/}"