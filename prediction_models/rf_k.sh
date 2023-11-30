#!/bin/bash
#SBATCH -J rf_k
#SBATCH -c 64
#SBATCH --mem=4GB
#SBATCH -p cpu
#SBATCH -t 00:02:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/rf_k_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/rf_k_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/prediction_models/rf_k.R