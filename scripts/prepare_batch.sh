#!/usr/bin/env bash
# Usage example:
#   prepare_batch.sh /path/to/instance/list.instances /path/to/results/dir [-dX (depth)] [str / gmic (mode)]

if [ ! -z "$VPC_DIR" ]; then
  export VPC_DIR=${VPC_DIR}
fi
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
SILENT=1
MODE="gmic"
MODE="str-a2"
MODE="str"

if [ "$(uname)" == "Darwin" ]; then
  export PROJ_DIR=`realpath ${PROJ_DIR}`
  if [ -z "$VPC_DIR" ]; then
    export VPC_DIR=`realpath ${PROJ_DIR}/../vpc`
  fi
else
  export PROJ_DIR=`realpath -s ${PROJ_DIR}`
  if [ -z "$VPC_DIR" ]; then
    export VPC_DIR=`realpath -s ${PROJ_DIR}/../vpc`
  fi
else
fi

export OPTFILE="${VPC_DIR}/data/ip_obj.csv"
export SCRIPT_DIR=${PROJ_DIR}/scripts

export INSTANCE_DIR=${VPC_DIR}/data/instances
export RESULTS_DIR=${PROJ_DIR}/results
export SOL_DIR=${VPC_DIR}/data/solutions

EXECUTABLE="${PROJ_DIR}/Release/main"

if [ $MODE == preprocess ]; then
  INSTANCE_LIST=${SCRIPT_DIR}/original.instances
else
  INSTANCE_LIST=${SCRIPT_DIR}/presolved.instances
fi

# HiPerGator
#export INSTANCE_DIR=/blue/akazachkov/${USER}/instances/vpc
#export RESULTS_DIR=/blue/akazachkov/$USER/results
#export INSTANCE_LIST=${SCRIPT_DIR}/presolved.instances

if [ "$(uname)" == "Darwin" ]; then
  # MBP19
  export INSTANCE_DIR=${REPOS_DIR}/instances
fi

if [ "$(hostname)" == "ISE-D41L3Q3" ]; then
  # w401
  export LOCAL_DIR=${HOME}
  export INSTANCE_DIR=${LOCAL_DIR}/instances
  export RESULTS_DIR=${LOCAL_DIR}/results
  export SOL_DIR=${INSTANCE_DIR}/solutions
  export INSTANCE_LIST=${VPC_DIR}/scripts/presolved.instances
fi

# Accept user options for instance list, results directory, and mode
if [ ! -z $1 ]; then
  INSTANCE_LIST=$1
fi
if [ ! -z $2 ]; then
  RESULTS_DIR=$2
fi
if [ -z $3 ]; then
  depthList=(2 4 8 16 32 64)
  > $JOB_LIST
else
  depthList=($3)
fi
if [ ! -z $4 ]; then
  MODE=$4
fi
JOB_LIST="job_list_${MODE}.txt"

# Set parameters
PARAMS=" --optfile=${OPTFILE}"
PARAMS="$PARAMS -t 3600"
PARAMS="$PARAMS --rounds=1"
if [ $MODE = "gmic" ]; then
  depthList=(0)
  PARAMS="$PARAMS --gomory=2"
elif [ $MODE = "str" ]; then
  PARAMS="$PARAMS --strengthen=1"
  PARAMS="$PARAMS --gomory=-1"
  PARAMS="$PARAMS --bb_runs=0"
  PARAMS="$PARAMS --bb_mode=10"
  PARAMS="$PARAMS --bb_strategy=536"
  PARAMS="$PARAMS -a1"
elif [ $MODE = "str-a2" ]; then
  PARAMS="$PARAMS --strengthen=1"
  PARAMS="$PARAMS --gomory=-1"
  PARAMS="$PARAMS --bb_runs=0"
  PARAMS="$PARAMS --bb_mode=10"
  PARAMS="$PARAMS --bb_strategy=536"
  PARAMS="$PARAMS -a2"
  PARAMS="$PARAMS --rcvmip_max_iters=1000"
  PARAMS="$PARAMS --rcvmip_total_timelimit=3600"
  PARAMS="$PARAMS --rcvmip_cut_timelimit=600"
  PARAMS="$PARAMS --atilde_compute_rank=0"
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
