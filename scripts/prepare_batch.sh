#!/bin/bash

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

export PROJ_DIR=`realpath -s ${PROJ_DIR}`
export VPC_DIR=`realpath -s ${PROJ_DIR}/../vpc`
export INSTANCE_DIR=${VPC_DIR}/data/instances
export SCRIPT_DIR=${PROJ_DIR}/scripts
export INSTANCE_LIST=${SCRIPT_DIR}/small_presolved.instances
export RESULTS_DIR=${PROJ_DIR}/results
export RESULTS_DIR=/local1/$USER/results

DEPTH="-d2"
PARAMS="$DEPTH"
PARAMS="$PARAMS -t 3600"
PARAMS="$PARAMS --rounds=1"
PARAMS="$PARAMS --strengthen=1"
PARAMS="$PARAMS --gomory=-1"
PARAMS="$PARAMS --bb_runs=0"
PARAMS="$PARAMS --bb_mode=10"
PARAMS="$PARAMS --bb_strategy=536"
#PARAMS="$PARAMS --use_all_ones=1"
#PARAMS="$PARAMS --use_iter_bilinear=1"
#PARAMS="$PARAMS --use_disj_lb=1"
#PARAMS="$PARAMS --use_tight_points=0"
#PARAMS="$PARAMS --use_tight_rays=0"
#PARAMS="$PARAMS --use_unit_vectors=0"
PARAMS="$PARAMS --temp=1"

TASK_ID=0
> job_list_strengthen.txt
while read line; do
  TASK_ID=$((TASK_ID+1))

  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi

  CASE_NUM=`printf %03d $TASK_ID`
  STUB=`date +%F`
  OUT_DIR=${RESULTS_DIR}/$STUB/str$DEPTH/${CASE_NUM}
  echo "Preparing command to run instance $line (task $TASK_ID) at `date`"
  echo "mkdir -p ${OUT_DIR}" >> job_list_strengthen.txt
  echo "nohup /usr/bin/time -v ${PROJ_DIR}/Release/main -f ${INSTANCE_DIR}/$line.mps --log=${OUT_DIR}/vpc-str.csv --optfile=${VPC_DIR}/data/ip_obj.csv $PARAMS 2>&1" >> job_list_strengthen.txt
done < ${INSTANCE_LIST}

