#/bin/bash

if [ -z ${PROJ_DIR} ]; then
  if [ -z "${REPOS_DIR}" ]; then
    echo "Please define VPC_DIR (the root vpc dir, possibly ${REPOS_DIR}/vpc):"
  else
    echo "Please define VPC_DIR (the root vpc dir):"
  fi
  read VPC_DIR
  echo "Set VPC_DIR=$VPC_DIR"
  if [ -z "$VPC_DIR" ]; then
    echo "Need to define VPC_DIR. Exiting."
    exit
  fi
fi

# Cbc: --bb_runs=1 --bb_mode=11 --bb_strategy=528
# CPLEX: --bb_runs=1 --bb_mode=11 --bb_strategy=532
# Gurobi: --bb_runs=1 --bb_mode=11 --bb_strategy=536
${PROJ_DIR}/Debug/main -f ${PROJ_DIR}/test/bm23.mps --mode=1 --optfile=${PROJ_DIR}/data/ip_obj.csv $1 $2 $3 $4 $5 $6 $7 $8
