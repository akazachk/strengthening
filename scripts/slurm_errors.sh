#/usr/bin/env bash

if [ -z $1 ]; then
  echo "Job ID needs to be specified via the first argument of the script."
  exit 1
fi
JOB_ID=$1

# pull all out of memory, timeout, failed jobs, 
# then take only the line with ".batch",
# then remove $JOB_ID_ from the start of the line
# then remove the trailing .batch
# then replace newlines with commas
# then remove extra space
sacct -u akazachkov -j $JOB_ID -s to,f,oom --format user,jobid%-30 | grep batch | sed "s/${JOB_ID}_//" | sed 's/\.batch//' | paste -s -d, | sed 's/ //g'
