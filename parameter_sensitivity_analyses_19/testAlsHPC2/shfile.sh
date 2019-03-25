#!/bin/csh
#PBS -l wd
#PBS -q normal
#PBS -l walltime=48:00:00,jobfs=4000Mb
#PBS -l mem=20144MB
#PBS -l software=R
setenv TMPDIR $PBS_JOBFS
module load phyml
module load R/3.2.2
R --vanilla < run.
