#!/usr/bin/env bash

#SBATCH -J dask-worker
#SBATCH -p normal
#SBATCH -A rev
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH -t 01:00:00

JOB_ID=${SLURM_JOB_ID%;*}

/home/twillia2/.conda-envs/weto/bin/python -m distributed.cli.dask_worker tcp://10.60.1.120:41166 --nthreads 2 --nprocs 4 --memory-limit 4.00GB --name name --nanny --death-timeout 60 map_vals.py
