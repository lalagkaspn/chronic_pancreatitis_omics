#!/usr/bin/env Rscript

# This script is used for the fitting prediction models using 10x10 repeated
# cross-validation and evaluation on the validation set. Then the final selected
# prediction model is evaluated on the microarray data.

# Preparation #####
.libPaths(c("/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2",
            "/usr/local/lib/R/site-library",
            "/usr/local/lib/R/library",
            "/modules/spack/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.4.0/r-4.2.0-3kitpfbxevyhxd2adiznenkjqqdbekzs/rlib/R/library"))
library(stringr, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
# library(ggplot2)
library(openxlsx, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
library(rpart, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2/")
library(C50)
library(MLmetrics, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
library(caret, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
library(ipred, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
library(doParallel, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
library(kernlab, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")

# HPC test
print("Packages successfully loaded!")

# Import data
selected_markers_pre = read.xlsx("/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/genes_for_model_training.xlsx")
train_set = read.xlsx("/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/ML_output/train_set.xlsx")
validation_set = read.xlsx("/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/ML_output/validation_set.xlsx")
lasso_markers = read.xlsx("/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/ML_output/Lasso_coefficients.xlsx")
exprs_microarray = readRDS("/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/DGEA/z_exprs_nonas_microarrays.RDS")
pheno_microarray = read.xlsx("/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/DGEA/Pheno_microarrays.xlsx") %>% 
  dplyr::rename(Sample.ID = GEO_accession)

# Filter train and validation sets for the lasso markers
train_set = train_set[, c("Tissue_type", lasso_markers$Predictor)]
validation_set = validation_set[, c("Tissue_type", lasso_markers$Predictor)]

# Formation of the test set
# Create an "X_Entrez" column for easy joins between datasets
selected_markers_pre$x_entrez = paste0("X_", selected_markers_pre$EntrezGene.ID)
x_entrez_microarray = rownames(exprs_microarray)
str_sub(x_entrez_microarray, 0, 0) = "X_"
rownames(exprs_microarray) = x_entrez_microarray

# Convert to data frame and add gene symbol info
exprs_microarray_df = as.data.frame(exprs_microarray)
exprs_microarray_df$x_entrez = x_entrez_microarray
exprs_microarray_df = exprs_microarray_df %>% dplyr::inner_join(selected_markers_pre,
                                                                by = "x_entrez") %>%
  dplyr::select(-x_entrez, -EntrezGene.ID, -HGNC_Official)
rownames(exprs_microarray_df) = exprs_microarray_df$Gene.Symbol
exprs_microarray_df = exprs_microarray_df[lasso_markers$Predictor, ]

# Transpose to create a samples by genes data frame and add labels
test_set = as.data.frame(t(exprs_microarray_df))
test_set$Sample.ID = rownames(test_set)
test_set = test_set %>% dplyr::inner_join(pheno_microarray, by = "Sample.ID") %>%
  dplyr::select(-Sample.ID, -Study)

# Parallelization #####

# This script requires extensive computational resources and it was submitted
# to a HPC, where 50 cores were employed along with x GB of RAM
# Output is saved as an .xlsx file
# Intermediate objects can be requested from the authors

# Enable parallel programming
system <- Sys.info()['sysname']
print(paste0("Number of cores used: ", detectCores(), collapse = ""))
gc()
cores <- detectCores()

if (system == 'Windows') {
  cl <- makeCluster(getOption('cl.cores', cores), type='PSOCK')
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
} else {
  registerDoParallel(cores)
}

# Fitting Support Vector Machines #####
# Preparing grids and controls
cost_values = seq(from = 0.5, to = 20, by = 0.5)
linear_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                               selectionFunction = "best", allowParallel = TRUE)
linear_svm_tune = expand.grid(C = cost_values)
L2_linear_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                                  selectionFunction = "best", allowParallel = TRUE)
L2_linear_svm_tune = expand.grid(cost = cost_values,
                                 Loss = "L2")
rbf_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                            selectionFunction = "best", allowParallel = TRUE)
rbf_svm_tune = expand.grid(C = cost_values,
                           sigma = c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9))

# Radial Basis Function
RNGversion("4.2.2")
set.seed(123)
print("Fitting radial basis function SVM model!")
rbf_svm_model = train(Tissue_type ~., data = train_set, 
                      method = "svmRadial", 
                      preProcess = NULL,
                      trControl = rbf_svm_ctrl,
                      tuneGrid = rbf_svm_tune)
saveRDS(rbf_svm_model, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/prediction_models/rbf_svm_model.rds")
print("Finished fitting radial basis function SVM model!")
