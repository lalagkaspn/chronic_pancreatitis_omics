#!/bin/bash
#SBATCH -J linear_svm
#SBATCH -c 64
#SBATCH --mem=176GB
#SBATCH -p cpu
#SBATCH -t 01:00:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/linear_svm_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/linear_svm_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/prediction_models/linear_svm.R