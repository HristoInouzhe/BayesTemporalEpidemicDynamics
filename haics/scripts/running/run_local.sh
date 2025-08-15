#!/bin/bash

###############################################################################
#Script Description
###############################################################################
#
#Description---------------------------
# Shell script that prepares the scene and launches simulation/src/bin/haics for a single case of a batch on Hipatia.
#
#Argument/s----------------------------
#	ID: (string) name of the simulations to run.
#	INT_RUNS: (integer) number of independent runs to simulate
#   REUSE: (boolean) where the binary should be compiled again
#
#Output--------------------------------
#	none directly, but calls simulation/src/bin/haics which does produce output files
#
#Author/s------------------------------
# Felix Muller
# Jorge Perez Heredia
# Lorenzo Nagar
#------------------------------------------------------------------------------

echo "======================================="
echo "HaiCS Simulation on a Local Machine"
echo "======================================="

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


###############################################################################
#Get running configuration
###############################################################################

echo -e "Reading arguments from command line..."

ID=${1}
if [ -z "${1}" ]
then
	fnExitError "ID needs to be supplied as the first argument."
fi

INT_RUNS=${2}

REUSE=${3}

RunParallel=${4}

USER=$(whoami)
DIR_HMCS=$(pwd)/../..
DIR_SIM=${DIR_HMCS}/simulation
DIR_OUT=${DIR_HMCS}/output

###############################################################################
#Build binary (if specified)
###############################################################################

if [ ${REUSE} -eq 0 ]
then
	echo -e "\nBuilding the binary..."
	rm -rf ${DIR_SIM}/bin/haics
	(cd ${DIR_SIM}/src/ && make clean && make)
fi

###############################################################################
#Launch simulations
###############################################################################

echo -e "Launching Simulations..."
START_TIME=$(date +%s)

echo -e "RunParallel ${RunParallel}"

if [ ${RunParallel} -eq 0 ]; then
	echo -e "----------------------------------------"
	echo -e "Running ${INT_RUNS} chains sequentially"
	echo -e "----------------------------------------"

	 for run in $(seq 1 1 ${INT_RUNS})
	 	do
		echo -e "----------------------------------------"
		echo -e "Run ${run} out of ${INT_RUNS}"
		echo -e "----------------------------------------"

		ID_run=${ID}_${run}
		DIR_OUT_HAICS=${DIR_SIM}/output/${ID_run}
		echo -e "creating HaiCS output directory at ${DIR_OUT_HAICS}"
		mkdir -p ${DIR_OUT_HAICS}

		#copy the inputfile (if it exists)
		echo -e "copying inputfile to ${DIR_OUT_HAICS}"
		if [ -f "${DIR_OUT}/${ID}/input/inputfile_${ID}.txt" ]; then
			echo -e "copying the inputfile to ${DIR_OUT_HAICS}"
			cp "${DIR_OUT}/${ID}/input/inputfile_${ID}.txt" "${DIR_OUT_HAICS}/inputfile_${ID_run}.txt"
		elif [ -f "${DIR_OUT}/${ID}/input/inputfile_${ID_run}.txt" ]; then
			echo -e "copying the inputfile to ${DIR_OUT_HAICS}"
			cp "${DIR_OUT}/${ID}/input/inputfile_${ID_run}.txt" "${DIR_OUT_HAICS}/inputfile_${ID_run}.txt"
		fi	

		#copy the initialpoint (if it exists)
		echo -e "copying initialpoint file to ${DIR_OUT_HAICS}"
		if [ -f "${DIR_OUT}/${ID}/input/initialpoint_${ID}.txt" ]; then
			echo -e "copying the initialpoint to ${DIR_OUT_HAICS}"
			cp "${DIR_OUT}/${ID}/input/initialpoint_${ID}.txt" "${DIR_OUT_HAICS}/initialpoint_${ID_run}.txt"
		elif [ -f "${DIR_OUT}/${ID}/input/initialpoint_${ID_run}.txt" ]; then
			echo -e "copying the initialpoint to ${DIR_OUT_HAICS}"
			cp "${DIR_OUT}/${ID}/input/initialpoint_${ID_run}.txt" "${DIR_OUT_HAICS}/initialpoint_${ID_run}.txt"
		fi	

		#copy the spline basis file (if it exists)
		echo -e "copying splinebasis file to ${DIR_OUT_HAICS}"
		if [ -f "${DIR_OUT}/${ID}/input/splinebasis_${ID}.txt" ]; then
			echo -e "copying the spline basis file to ${DIR_OUT_HAICS}"
			cp "${DIR_OUT}/${ID}/input/splinebasis_${ID}.txt" "${DIR_OUT_HAICS}/splinebasis_${ID_run}.txt"
		elif [ -f "${DIR_OUT}/${ID}/input/splinebasis_${ID_run}.txt" ]; then
			echo -e "copying the spline basis file to ${DIR_OUT_HAICS}"
			cp "${DIR_OUT}/${ID}/input/splinebasis_${ID_run}.txt" "${DIR_OUT_HAICS}/splinebasis_${ID_run}.txt"
		fi

		echo -e "HAICS Job ID=${ID_run} started..."
		(cd ${DIR_SIM} && bin/haics ${ID_run})
		#the cd is necessary because haics needs to be excuted like "bin/haics <ID> for the paths"

		#remove the empty files after running bin/haics
		#find ${DIR_OUT_HAICS} -type f -empty -delete
		done
elif [ ${RunParallel} -eq 1 ]; then
	export INT_RUNS
	export ID
	export DIR_SIM
	export DIR_OUT

	echo -e "----------------------------------------"
	echo -e "Running ${INT_RUNS} chains in parallel"
	echo -e "----------------------------------------"

	run_haics() {
    	run=$1
    	echo -e "----------------------------------------"
    	echo -e "Run ${run} out of ${INT_RUNS}"
    	echo -e "----------------------------------------"

    	ID_run=${ID}_${run}
    	DIR_OUT_HAICS=${DIR_SIM}/output/${ID_run}
    	echo -e "creating HaiCS output directory at ${DIR_OUT_HAICS}"
    	mkdir -p ${DIR_OUT_HAICS}

    	# Copy the input file (if it exists)
    	echo -e "copying inputfile to ${DIR_OUT_HAICS}"
    	if [ -f "${DIR_OUT}/${ID}/input/inputfile_${ID}.txt" ]; then
        	cp "${DIR_OUT}/${ID}/input/inputfile_${ID}.txt" "${DIR_OUT_HAICS}/inputfile_${ID_run}.txt"
    	elif [ -f "${DIR_OUT}/${ID}/input/inputfile_${ID_run}.txt" ]; then
        	cp "${DIR_OUT}/${ID}/input/inputfile_${ID_run}.txt" "${DIR_OUT_HAICS}/inputfile_${ID_run}.txt"
    	fi

    	# Copy the initial point file (if it exists)
    	echo -e "copying initialpoint file to ${DIR_OUT_HAICS}"
    	if [ -f "${DIR_OUT}/${ID}/input/initialpoint_${ID}.txt" ]; then
        	cp "${DIR_OUT}/${ID}/input/initialpoint_${ID}.txt" "${DIR_OUT_HAICS}/initialpoint_${ID_run}.txt"
    	elif [ -f "${DIR_OUT}/${ID}/input/initialpoint_${ID_run}.txt" ]; then
        	cp "${DIR_OUT}/${ID}/input/initialpoint_${ID_run}.txt" "${DIR_OUT_HAICS}/initialpoint_${ID_run}.txt"
    	fi

    	# Copy the spline basis file (if it exists)
    	echo -e "copying splinebasis file to ${DIR_OUT_HAICS}"
    	if [ -f "${DIR_OUT}/${ID}/input/splinebasis_${ID}.txt" ]; then
        	cp "${DIR_OUT}/${ID}/input/splinebasis_${ID}.txt" "${DIR_OUT_HAICS}/splinebasis_${ID_run}.txt"
    	elif [ -f "${DIR_OUT}/${ID}/input/splinebasis_${ID_run}.txt" ]; then
        	cp "${DIR_OUT}/${ID}/input/splinebasis_${ID_run}.txt" "${DIR_OUT_HAICS}/splinebasis_${ID_run}.txt"
    	fi

    	echo -e "HAICS Job ID=${ID_run} started..."
    	(cd ${DIR_SIM} && bin/haics ${ID_run})

    	# Remove empty files after running bin/haics
    	find ${DIR_OUT_HAICS} -type f -empty -delete
	}

	export -f run_haics

	# Set number of jobs to number of CPU cores (Linux/MacOS detection)
	NUM_JOBS=$(nproc)  # For Linux
	# NUM_JOBS=$(sysctl -n hw.ncpu)  # For macOS

	# Run simulations in parallel
	# seq 1 "${INT_RUNS}" | parallel -j "${NUM_JOBS}" run_haics {}
	seq 1 "${INT_RUNS}" | parallel -j "${NUM_JOBS}" "run_haics {} > /home/hinouzhe/haics_GitHub/logs/job_{}.log 2>&1"
	
fi

###############################################################################
#Sort output folders differently
###############################################################################

mkdir -p ${DIR_OUT}/${ID}
#collect output from haics directory back to DIR_OUT
mv ${DIR_SIM}/output/${ID}_* ${DIR_OUT}/${ID}/.

###############################################################################
#Measure running time
###############################################################################

echo "======================================="
echo -e "The output for HaiCS Job ID=${ID} can be found at: \n haics/output/${ID}"
END_TIME=$(date +%s)
echo -e "TOTAL_RUN_TIME (hours)     = "`echo "$START_TIME $END_TIME" | awk '{printf("%.4f",($2-$1)/60.0/60.0)}'`
echo "======================================="
echo

###############################################################################
#El fin
###############################################################################
echo -e "Reached the end of ${0##*/}"
