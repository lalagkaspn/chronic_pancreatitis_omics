#!/bin/bash
#SBATCH -J dt_boosting
#SBATCH -c 64
#SBATCH --mem=503GB
#SBATCH -p cpu-long
#SBATCH -t 72:00:00
#SBATCH -o /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/dt_boosting_%J.out
#SBATCH -e /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/logs/dt_boosting_%J.err

# load R version specific
module load r/4.2.0.lua

# run R script
Rscript /home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/decision_tree_boosting.R