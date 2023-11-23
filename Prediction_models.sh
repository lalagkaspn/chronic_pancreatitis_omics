#!/bin/bash
#SBATCH -J CP_prediction_models
#SBATCH -c 100
#SBATCH --mem=1024GB
#SBATCH -p cpu-long
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/prediction_model_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/prediction_model_%J.err

# Assign the provided arguments to variables
#arg1="$1" # Subject, GTEX.13O3O
#arg2="$2" # Region, Putamen_basal_ganglia

#Rscript /project/pi_rachel_melamed_uml_edu/Jianfeng/Allen/src/R/CP_gtex_subject_region_multiple.R "$arg1" "$arg2"

Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/Prediction_models.R