#!/bin/bash
#SBATCH -J rbf_svm_model
#SBATCH -c 64
#SBATCH --mem=48GB
#SBATCH -p cpu-long
#SBATCH -t 00:05:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/rbf_svm_model_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/rbf_svm_model_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/prediction_models/rbf_svm_model.R