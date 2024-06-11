#/bin/bash

if [ -z ${PROJ_DIR} ]; then
  if [ -z "${REPOS_DIR}" ]; then
    echo "Please define PROJ_DIR (the root strengthening dir, possibly ${REPOS_DIR}/strengthening):"
  else
    echo "Please define PROJ_DIR (the root strengthening dir):"
  fi
  read PROJ_DIR
  echo "Set PROJ_DIR=$PROJ_DIR"
  if [ -z "$PROJ_DIR" ]; then
    echo "Need to define PROJ_DIR. Exiting."
    exit
  fi
fi

VPC_DIR=${PROJ_DIR}/../vpc

# Cbc: --bb_runs=1 --bb_mode=11 --bb_strategy=528
# CPLEX: --bb_runs=1 --bb_mode=11 --bb_strategy=532
# Gurobi: --bb_runs=1 --bb_mode=11 --bb_strategy=536
${PROJ_DIR}/Debug/main -f ${PROJ_DIR}/test/bm23.mps -d 0 --strengthen=1 --gomory=-1 --optfile=${VPC_DIR}/data/ip_obj.csv --disj_options="2:4:8:16:32:64" --analyze_regularity=1 --mode=4 --cutlimit=-6 --temp=40 -v 0 $1 $2 $3 $4 $5 $6 $7 $8
