#!/bin/bash
#SBATCH -J lasso_selection
#SBATCH -c 64
#SBATCH --mem=64GB
#SBATCH -p cpu
#SBATCH -t 00:01:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/lasso_prediction_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/lasso_prediction_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/lasso_selection.R