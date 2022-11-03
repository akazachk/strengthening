#!/usr/bin/env bash
# Usage example:
#   prepare_batch.sh /path/to/instance/list.instances /path/to/results/dir [-dX (depth)] [str / gmic (mode)]
#   prepare_batch.sh /path/to/instance/list.instances /path/to/results/dir [test / preprocess / bb / bb0]

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
MODE="str"

export PROJ_DIR=`realpath -s ${PROJ_DIR}`
export VPC_DIR=`realpath -s ${PROJ_DIR}/../vpc`
export INSTANCE_DIR=${VPC_DIR}/data/instances
export SOL_DIR=${VPC_DIR}/data/solutions
export SCRIPT_DIR=${PROJ_DIR}/scripts
export INSTANCE_LIST=${SCRIPT_DIR}/small_presolved.instances
export RESULTS_DIR=${PROJ_DIR}/results
export OPTFILE="${VPC_DIR}/data/ip_obj.csv"
JOB_LIST="job_list_strengthen.txt"
EXECUTABLE="${PROJ_DIR}/Release/main"

# HiPerGator
export INSTANCE_DIR=/blue/akazachkov/${USER}/instances/vpc
export RESULTS_DIR=/blue/akazachkov/$USER/results
export INSTANCE_LIST=${SCRIPT_DIR}/presolved.instances

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
else
  echo "*** ERROR: Option $MODE not recognized"
  exit
fi
PARAMS="$PARAMS --temp=8"

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
  while read line; do
    TASK_ID=$((TASK_ID+1))

    # Skip empty lines
    if [ -z "$line" ]; then
      continue
    fi

    # Prepare out directory, based on current date
    CASE_NUM=`printf %0${NUM_DIGITS}d $TASK_ID`
    STUB=`date +%F`
    OUT_DIR="${RESULTS_DIR}/$STUB/str-d$d/${CASE_NUM}"

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
      echo "*** WARNING: Could not find $SOLFILE"
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
      if [ $MODE = "gmic" ]; then
        echo -n "nohup /usr/bin/time -v $EXECUTABLE -f ${INSTANCE_DIR}/$line.mps --log=${OUT_DIR}/vpc-$MODE $SOLPARAM $PARAMS --strengthen=0 -d$d >> ${OUT_DIR}/log.out 2>&1; " >> ${JOB_LIST}
        echo "nohup /usr/bin/time -v $EXECUTABLE -f ${INSTANCE_DIR}/$line.mps --log=${OUT_DIR}/vpc-$MODE $SOLPARAM $PARAMS --strengthen=1 -d$d >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
      else
        #echo "nohup /usr/bin/time -v $EXECUTABLE -f ${INSTANCE_DIR}/$line.mps --log=${OUT_DIR}/vpc-$MODE --optfile=${OPTFILE} $SOLPARAM $PARAMS >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
        echo "nohup /usr/bin/time -v $EXECUTABLE -f ${FILE} --log=${OUT_DIR}/vpc-${MODE}.csv $SOLPARAM $PARAMS -d$d >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
      fi
    fi
  done < ${INSTANCE_LIST}
done # loop over depth list

echo "Done preparing $JOB_LIST. Total errors: $TOTAL_ERRORS."
