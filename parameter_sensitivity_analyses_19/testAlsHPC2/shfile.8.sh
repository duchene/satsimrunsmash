#!/bin/bash
#PBS -P RDS-FSC-Phylogenomics-RW
#PBS -l select=1:ncpus=1:mem=16GB
#PBS -l walltime=240:00:00
#PBS -q defaultQ
cd $PBS_O_WORKDIR
module load phyml
module load R/3.2.2
R --vanilla < run.sattr.8.Rscript
