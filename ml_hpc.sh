#!/bin/bash
#SBATCH -J CP_ml_hpc
#SBATCH -c 100
#SBATCH --mem=992GB
#SBATCH -p cpu-long
#SBATCH -t 96:00:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/ml_hpc_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/ml_hpc_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/ml_hpc.R