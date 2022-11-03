#!/bin/bash
#SBATCH --job-name=vpc-str
#SBATCH --output=slurm_%x_%A-%a.log       # standard output log (%x = job name, %A = job id, %a = array index)
#SBATCH --error=slurm_%x_%A-%a.err        # standard error log (%x = job name, %A = job id, %a = array index)
#SBATCH --ntasks=1                        # run a single task

#SBATCH --mail-user=akazachkov@ufl.edu
#SBATCH --mail-type=BEGIN,FAIL,END        # mail events (NONE, BEGIN, END, FAIL, ALL)

#SBATCH --account=akazachkov
#SBATCH --qos=akazachkov

#SBATCH --time=01:20:00                   # time limit hrs:min:sec
#SBATCH --mem=600M                        # job memory
##SBATCH --array=16,29,44,109,214,217,218,421,431
##SBATCH --array=1-15,17-28,30-43,45-108,110-213,215-216,219-420,422-430,432-435

#SBATCH --array=1-2982
##SBATCH --array=1                         # test

#########################
## To run this script, call (for example)
##     sbatch array_batch.sh [job list file]
## See arguments below
echo "=== START SLURM SCRIPT MESSAGES ==="
pwd; hostname; date


#########################
## Arguments
# Argument 1: CMD_FILE Set command file if given
if [ ! -z $1 ]
then
  CMD_FILE=$1
else
  CMD_FILE=job_list_strengthen.txt
fi


#########################
## Constants
export PROJ_DIR="${REPOS_DIR}/strengthening"
echo "Set PROJ_DIR=$PROJ_DIR"
export CMD_FILE_DIR=${PROJ_DIR}/scripts
echo "Set CMD_FILE_DIR=$CMD_FILE_DIR"


#########################
## Prepare run command
echo ""
echo "Starting ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CMD_FILE_DIR}/${CMD_FILE})

echo ""
echo -e "Calling:\n$CMD"
echo "=== START RUNNING JOB ==="
echo ""


#########################
## RUN COMMAND HERE
eval $CMD


#########################
## Wrap up
echo "=== END JOB ==="
date
#echo "=== END SLURM SCRIPT MESSAGES ==="
