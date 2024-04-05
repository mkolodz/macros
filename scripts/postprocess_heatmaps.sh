#!/bin/bash
# . /scratch/gccb/software/framework/20230109-4to1-ubuntu20/profile.sh
SLURM_ARRAY_TASK_ID=3
# i=$((i+1))
REAL_SLURM_ARRAY_TASK_ID=$((SLURM_ARRAY_TASK_ID+20))
OUTPUTFILE=_sifi_b$REAL_SLURM_ARRAY_TASK_ID.root
echo $OUTPUTFILE

