#!/bin/bash
#SBATCH -J L2_svm
#SBATCH -c 64
#SBATCH --mem=124GB
#SBATCH -p cpu
#SBATCH -t 02:30:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/L2_svm_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/L2_svm_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/prediction_models/L2_svm.R