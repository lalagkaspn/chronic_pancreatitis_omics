#!/bin/bash
#SBATCH -J C5.0_ROC
#SBATCH -c 64
#SBATCH --mem=32GB
#SBATCH -p cpu
#SBATCH -t 00:05:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/C5.0_ROC_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/C5.0_ROC_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/prediction_models/C5.0_ROC.R