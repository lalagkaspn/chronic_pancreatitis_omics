#!/bin/bash
#SBATCH -J glmnet_selection
#SBATCH -c 64
#SBATCH --mem=64GB
#SBATCH -p cpu
#SBATCH -t 00:20:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/glmnet_selection_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/glmnet_selection_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/glmnet_selection.R