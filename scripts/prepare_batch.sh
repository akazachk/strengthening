#!/usr/bin/env bash
# Usage example:
#   prepare_batch.sh /path/to/instance/list.instances /path/to/results/dir [(mode) test / str / gmic / disjset / splits / 6splits ] [(depth) -dX] 

if [ -z "$PROJ_DIR" ]
then
  if [ ! -z "${REPOS_DIR}" ]
  then
    echo "Please define PROJ_DIR (the root project dir, possibly ${REPOS_DIR}/strengthening):"
  else
    echo "Please define PROJ_DIR (the root project dir):"
  fi
  read PROJ_DIR
  if [ -z "$PROJ_DIR" ]
    then echo "Need to define PROJ_DIR. Exiting."
    exit 1
  fi
fi

# Constants
SILENT=1 # set to 0 to print more as the file is being created
MODE="str-a2"
MODE="str"

# Directory options
DEFAULT_DIRS=1     # Set to 0 for custom directories (which you will need to input below)
SAVE_TO_HOME_DIR=1 # SAVE_TO_HOME_DIR = 0: results and instances are relative to ${PROJ_DIR}; = 1: relative to ${HOME}

# Change PROJ_DIR, VPC_DIR to be full path
if [ ! -z "$VPC_DIR" ]; then
  export VPC_DIR=${VPC_DIR}
else
  export VPC_DIR=${PROJ_DIR}/../vpc # For portability to other projects, in which PROJ_DIR and VPC_DIR might be different
fi
if [ "$(uname)" == "Darwin" ]; then
  export PROJ_DIR=`realpath ${PROJ_DIR}`
  export VPC_DIR=`realpath ${VPC_DIR}`
else
  export PROJ_DIR=`realpath -s ${PROJ_DIR}`
  export VPC_DIR=`realpath -s ${VPC_DIR}`
fi

# Set relative path for results/instances
if [ ${DEFAULT_DIRS} == 1 ]; then
  if [ ${SAVE_TO_HOME_DIR} == 1 ]; then
    export LOCAL_DIR=${HOME}
  else
    export LOCAL_DIR=${PROJ_DIR}
  fi
else
  # *** Replace with custom dir below
  export LOCAL_DIR=${HOME}
fi

# Results will be sent to ${RESULTS_DIR}/[date]/[mode]/[instance #]
# Overwritten with second command line argument
if [ ${DEFAULT_DIRS} == 1 ]; then
  export RESULTS_DIR=${LOCAL_DIR}/results
else
  # *** Replace with custom dir below
  export RESULTS_DIR=${PROJ_DIR}/results
fi

# File with IP objective values 
if [ ${DEFAULT_DIRS} == 1 ]; then
  export OPTFILE="${VPC_DIR}/data/ip_obj.csv"
else
  # *** Replace with custom dir below
  export OPTFILE="${VPC_DIR}/data/ip_obj.csv"
fi

# Directory with instances (instance list will provide relative paths from this directory)
if [ ${DEFAULT_DIRS} == 1 ]; then
  export INSTANCE_DIR=${LOCAL_DIR}/instances
else
  # *** Replace with custom dir below
  export INSTANCE_DIR=${LOCAL_DIR}/instances
fi

# Where to find (or save) IP solutions
if [ ${DEFAULT_DIRS} == 1 ]; then
  export SOL_DIR=${LOCAL_DIR}/solutions
else
  # *** Replace with custom dir below
  export SOL_DIR=${LOCAL_DIR}/solutions
fi

# Directory with instance lists for default instances
# Overwritten by first command line argument
if [ ${DEFAULT_DIRS} == 1 ]; then
  export INSTANCE_LIST_DIR=${PROJ_DIR}/data/experiments
else
  # *** Replace with custom dir below
  export INSTANCE_LIST_DIR=${PROJ_DIR}/data/experiments
fi

# Set default instance lists
if [ $MODE == preprocess ]; then
  INSTANCE_LIST="${INSTANCE_LIST_DIR}/original.test"
elif [ $MODE == gmic ]; then
  INSTANCE_LIST="${INSTANCE_LIST_DIR}/gmic.test"
else
  INSTANCE_LIST="${INSTANCE_LIST_DIR}/presolved.test"
fi

# HiPerGator
#export INSTANCE_DIR=/blue/akazachkov/${USER}/instances/vpc
#export RESULTS_DIR=/blue/akazachkov/$USER/results
#export INSTANCE_LIST=${SCRIPT_DIR}/presolved.instances

# Accept user options for instance list, results directory, and mode
if [ ! -z $1 ]; then
  INSTANCE_LIST=$1
fi
if [ ! -z $2 ]; then
  RESULTS_DIR=$2
fi
if [ ! -z $3 ]; then
  MODE=$3
fi
if [ -z $4 ]; then
  depthList=(2 4 8 16 32 64)
  #> $JOB_LIST
else
  depthList=($4)
fi

# JOB_LIST is where commands to be run will be saved (in current directory)
JOB_LIST="job_list_${MODE}.txt"

# EXECUTABLE is what file we ultimately run
EXECUTABLE="${PROJ_DIR}/Release/main"

# Constrain to P-cores on Linux
# (I am not sure about this ... maybe can use 0-15? ... but then two processes might end up on the same core.)
if [ $(uname) != "Darwin" ]; then
  EXECUTABLE="taskset -c 0,2,4,6,8,10,12,14 $EXECUTABLE"
fi

# Set parameters
PARAMS=" --optfile=${OPTFILE}"
PARAMS="$PARAMS -t 3600"
PARAMS="$PARAMS --rounds=1"
PARAMS="$PARAMS --use_all_ones=1"
PARAMS="$PARAMS --use_iter_bilinear=1"
PARAMS="$PARAMS --use_disj_lb=1"
PARAMS="$PARAMS --use_tight_points=0"
PARAMS="$PARAMS --use_tight_rays=0"
PARAMS="$PARAMS --use_unit_vectors=0"
PARAMS="$PARAMS --max_num_hot_start_viol=0"
if [ $MODE = "gmic" ]; then
  depthList=(0)
  PARAMS="$PARAMS --gomory=2"
elif [ $MODE = "str" ]; then
  PARAMS="$PARAMS --strengthen=1"
  PARAMS="$PARAMS --analyze_regularity=1"
  PARAMS="$PARAMS --gomory=-1"
  PARAMS="$PARAMS --bb_runs=0"
  PARAMS="$PARAMS --bb_mode=10"
  PARAMS="$PARAMS --bb_strategy=536"
elif [ $MODE = "str-a2" ]; then
  PARAMS="$PARAMS --strengthen=1"
  PARAMS="$PARAMS --analyze_regularity=2"
  PARAMS="$PARAMS --gomory=-1"
  PARAMS="$PARAMS --bb_runs=0"
  PARAMS="$PARAMS --bb_mode=10"
  PARAMS="$PARAMS --bb_strategy=536"
  PARAMS="$PARAMS --rcvmip_max_iters=1000"
  PARAMS="$PARAMS --rcvmip_total_timelimit=3600"
  PARAMS="$PARAMS --rcvmip_cut_timelimit=600"
  PARAMS="$PARAMS --atilde_compute_rank=0"
elif [ $MODE == "disjset" ] ; then
  depthList=(0)
  PARAMS="$PARAMS --strengthen=1"
  PARAMS="$PARAMS --analyze_regularity=1"
  PARAMS="$PARAMS --gomory=-1"
  PARAMS="$PARAMS --bb_runs=0"
  PARAMS="$PARAMS --bb_mode=10"
  PARAMS="$PARAMS --bb_strategy=536"
  PARAMS="$PARAMS --mode=4"
  PARAMS="$PARAMS --disj_options=\"2:4:8:16:32:64\""
  PARAMS="$PARAMS --cutlimit=-1"
elif [ $MODE == "splits" ] ; then
  depthList=(100000)
  PARAMS="$PARAMS --strengthen=1"
  PARAMS="$PARAMS --analyze_regularity=1"
  PARAMS="$PARAMS --gomory=-1"
  PARAMS="$PARAMS --bb_runs=0"
  PARAMS="$PARAMS --bb_mode=10"
  PARAMS="$PARAMS --bb_strategy=536"
  PARAMS="$PARAMS --mode=1"
  PARAMS="$PARAMS --cutlimit=-1"
elif [ $MODE == "6splits" ] ; then
  depthList=(6)
  PARAMS="$PARAMS --strengthen=1"
  PARAMS="$PARAMS --analyze_regularity=1"
  PARAMS="$PARAMS --gomory=-1"
  PARAMS="$PARAMS --bb_runs=0"
  PARAMS="$PARAMS --bb_mode=10"
  PARAMS="$PARAMS --bb_strategy=536"
  PARAMS="$PARAMS --mode=1"
  PARAMS="$PARAMS --cutlimit=-1"
elif [ $MODE == "test" ]; then
  depthList=(2)
else
  echo "*** ERROR: Option $MODE not recognized"
  exit
fi
PARAMS="$PARAMS --temp=8"
PARAMS="$PARAMS --temp=40" # save IP sol

#PARAMS="$PARAMS --use_all_ones=1"
#PARAMS="$PARAMS --use_iter_bilinear=1"
#PARAMS="$PARAMS --use_disj_lb=1"
#PARAMS="$PARAMS --use_tight_points=0"
#PARAMS="$PARAMS --use_tight_rays=0"
#PARAMS="$PARAMS --use_unit_vectors=0"

# Figure out how many digits to prepend to each folder/case number
# This is done through a couple of calls to bc, which is usually available on most systems
if ! command -v bc &> /dev/null; then
  NUM_DIGITS=3
else
  NUM_JOBS=`< $INSTANCE_LIST wc -l`
  listLength=${#depthList[@]}
  NUM_JOBS=$((NUM_JOBS*listLength))
  LOG_NUM_JOBS=`echo "l(${NUM_JOBS})/l(10)" | bc -l`
  NUM_DIGITS=`echo "${LOG_NUM_JOBS}/1 + 1" | bc`
fi

TASK_ID=0
TOTAL_ERRORS=0
#> $JOB_LIST
for d in ${depthList[*]}; do
  echo "Depth $d"
  while read line; do
    TASK_ID=$((TASK_ID+1))

    # Skip empty lines
    if [ -z "$line" ]; then
      continue
    fi

    # Prepare out directory, based on current date
    CASE_NUM=`printf %0${NUM_DIGITS}d $TASK_ID`
    STUB=`date +%F`
    OUT_DIR=${RESULTS_DIR}/$STUB/${MODE}/${CASE_NUM}

    # Print status (in silent mode, print a ".")
    if [ $SILENT != 1 ]; then
      echo "Preparing command to run instance $line (task $TASK_ID) at `date`"
    else
      echo -n "."
    fi

    # Check if solution exists
    arrIN=(${line//\// })
    arrIN[0]=""
    SOLFILE="$SOL_DIR"
    for entry in ${arrIN[@]}; do
      SOLFILE="${SOLFILE}/${entry}"
    done
    SOLFILE="${SOLFILE}_gurobi.mst.gz"
    if test -f "$SOLFILE"; then
      if [ $SILENT != 1 ]; then
        echo "$SOLFILE exists"
      fi
      SOLPARAM="--solfile=${SOLFILE}"
    else
      if [ $SILENT != 1 ]; then
        echo "*** WARNING: Could not find $SOLFILE"
      fi
      SOLPARAM=""
    fi

    # Check if file exists
    FILE=${INSTANCE_DIR}/$line.mps
    if [ ! -f "$FILE" ] && [ ! -f "$FILE.gz" ] && [ ! -f "$FILE.bz2" ]; then
      echo "*** ERROR: $FILE does not exist; skipping."
      TOTAL_ERRORS=$((TOTAL_ERRORS+1))
    else
      # Finally, write command we will call to a file
      echo -n "mkdir -p ${OUT_DIR}; " >> ${JOB_LIST}
      if [ $(uname) == "Darwin" ]; then
        echo -n "nohup /usr/bin/time " >> ${JOB_LIST}
      else
        echo -n "nohup /usr/bin/time -v " >> ${JOB_LIST}
      fi
      if [ $MODE = "gmic" ]; then
        echo -n "$EXECUTABLE -f ${INSTANCE_DIR}/$line.mps --log=${OUT_DIR}/vpc-$MODE.csv $SOLPARAM $PARAMS --strengthen=0 -d$d >> ${OUT_DIR}/log.out 2>&1; " >> ${JOB_LIST}
        echo "$EXECUTABLE -f ${INSTANCE_DIR}/$line.mps --log=${OUT_DIR}/vpc-$MODE.csv $SOLPARAM $PARAMS --strengthen=1 -d$d >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
      else
        #echo "nohup $EXECUTABLE -f ${INSTANCE_DIR}/$line.mps --log=${OUT_DIR}/vpc-$MODE --optfile=${OPTFILE} $SOLPARAM $PARAMS >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
        echo "$EXECUTABLE -f ${FILE} --log=${OUT_DIR}/vpc-${MODE}.csv $SOLPARAM $PARAMS -d$d >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
      fi
    fi
  done < ${INSTANCE_LIST}
done # loop over depth list

# Shuffle command order to not have dependency in the performance
if [ $(uname) == "Darwin" ]; then
  sort -R ${JOB_LIST} --output=${JOB_LIST}
else
  shuf -o ${JOB_LIST} < ${JOB_LIST}
fi

echo "Done preparing $JOB_LIST. Total errors: $TOTAL_ERRORS."
