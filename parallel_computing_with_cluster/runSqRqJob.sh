#!/bin/bash
#SBATCH --partition=main
#SBATCH --array=1-500
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=slurm.out    # STDOUT output file

L=500	# Parse variable 1
sigma=10	# Parse variable 2
W=90	# Parse variable 3

srun julia Sq_Rq.jl `printf "%d 0.%02d 0.%02d" $L $sigma $W` > `printf "/scratch/mw936/data/L%d/sigma%03d/W%03d/%s_%s.csv" $L $sigma $W $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID`
