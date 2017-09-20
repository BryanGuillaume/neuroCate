#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=20:00:00
#PBS -o /scratch/users/nus/biebrlg/logs/
#PBS -e /scratch/users/nus/biebrlg/logs/

module add R/3.2.4
R CMD BATCH --no-save --no-restore "--args $PBS_ARRAY_INDEX 2 0" /home/users/nus/biebrlg/scripts/simuCompaCateStrategy.R /scratch/users/nus/biebrlg/logs/simuCompaCateStrategy_2Z_000alpha.$PBS_JOBID.Rout