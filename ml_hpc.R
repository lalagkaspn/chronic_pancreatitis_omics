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

# HPC test
print("NOTE: packages successfully loaded!")

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
library(doParallel, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")

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

# Fitting Decision Tree models #####
# Tuning grids and controls
grid1 = expand.grid(model = c("tree", "rule"),
                    trials = c(10,15,25,50,100),
                    winnow = c(TRUE, FALSE))
ctrl1 = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                     selectionFunction = "best", allowParallel = TRUE)
ctrl2 = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                     selectionFunction = "best", savePredictions = TRUE,
                     classProbs = TRUE, summaryFunction = multiClassSummary,
                     allowParallel = TRUE)
grid2 = expand.grid(mtry=c(2, 4, 8, 10, 20, 30, 40, 50, 60, 67))
boost_grid = expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                         eta = c(0.025, 0.05, 0.1, 0.3),
                         max_depth = c(2, 3, 4, 5, 6),
                         min_child_weight = c(1, 2, 3),
                         colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
                         subsample = c(0.5, 0.75, 1.0),
                         gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0))
boost_control = trainControl(method = "repeatedcv", # cross-validation
                             number = 10, repeats = 10,
                             verboseIter = FALSE, # no training log
                             allowParallel = TRUE # FALSE for reproducible results 
)

# C5.0 - k
RNGversion("4.2.2")
set.seed(123)
print("NOTE: Training the C5.0 - k model! (1/9)")
C50_model = train(Tissue_type ~., data = train_set,
                  method = "C5.0",
                  trControl = ctrl1,
                  metric = "Kappa",
                  tuneGrid = grid1)
saveRDS(C50_model, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/C5.0_k.rds")
print("NOTE: Finished training the C5.0 - k model! (1/9)")

# Boosting
RNGversion("4.2.2")
set.seed(123)
library(xgboost, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
print("NOTE: Training the boosting model! (2/9)")
boost = train(Tissue_type ~ .,
              data = train_set,
              trControl = boost_control,
              tuneGrid = boost_grid,
              method = "xgbTree",
              verbose = FALSE)
saveRDS(boost, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/xgboost.rds")
print("NOTE: Finished training the boosting model! (2/9)")

# Random Forest
RNGversion("4.2.2")
set.seed(123)
library(randomForest, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
print("NOTE: Training the random forest model! (3/9)")
rf = train(Tissue_type ~., data = train_set,
           method = "rf",
           metric = "Kappa", trControl = ctrl1, tuneGrid = grid2)
saveRDS(rf, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/rf_k.rds")
print("NOTE: Finished training the random forest model! (3/9)")

# RF with grid2 and ctrls, ROC
RNGversion("4.2.2")
set.seed(123)
print("NOTE: Training the RF with grid2 and ctrls, ROC model! (4/9)")
comp_rf = train(Tissue_type ~., data = train_set,
                method = "rf",
                metric = "ROC", trControl = ctrl2, tuneGrid = grid2)
saveRDS(comp_rf, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/rf_roc.rds")
print("NOTE: Finished training the RF with grid2 and ctrls, ROC model! (4/9)")

# C5.0 with grid_c50
RNGversion("4.2.2")
set.seed(123)
print("NOTE: Training the C5.0 with grid_c50 model! (5/9)")
m_c50 = train(Tissue_type ~., data = train_set,
              method="C5.0", metric = "ROC",
              trControl = ctrl2, tuneGrid = grid1)
saveRDS(m_c50, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/C5.0_ROC.rds")
print("NOTE: Finished training the C5.0 with grid_c50 model! (5/9)")

# Fitting Support Vector Machines #####
# Preparing grids and controls
print("NOTE: Fitting SVM!")
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

# Linear kernel - simple
RNGversion("4.2.2")
set.seed(123)
print("NOTE: Fitting linear kernel - simple! (6/9)")
library(kernlab, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
linear_svm_model = train(Tissue_type ~., data = train_set, 
                         method = "svmLinear", 
                         preProcess = NULL,
                         trControl = linear_svm_ctrl,
                         tuneGrid = linear_svm_tune)
saveRDS(linear_svm_model, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/linear_svm.rds")
print("NOTE: Finished fitting linear kernel - simple! (6/9)")

# Linear kernel - L2 regularised
RNGversion("4.2.2")
set.seed(123)
print("NOTE: Fitting linear kernel - L2 regularised! (7/9)")
library(LiblineaR, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
L2_linear_svm_model = train(Tissue_type ~., data = train_set, 
                            method = "svmLinear3", 
                            preProcess = NULL,
                            trControl = L2_linear_svm_ctrl,
                            tuneGrid = L2_linear_svm_tune)
saveRDS(L2_linear_svm_model, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/L2_svm.rds")
print("NOTE: Finished fitting linear kernel - L2 regularised! (7/9)")

# Radial Basis Function
RNGversion("4.2.2")
set.seed(123)
print("NOTE: Fitting radial basis function! (8/9)")
rbf_svm_model = train(Tissue_type ~., data = train_set, 
                      method = "svmRadial", 
                      preProcess = NULL,
                      trControl = rbf_svm_ctrl,
                      tuneGrid = rbf_svm_tune)
saveRDS(rbf_svm_model, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/rbf_svm_model.rds")
print("NOTE: Finished fitting radial basis function! (8/9)")

# Fitting k-NN models #####
# Preparing grids and controls
ctrl_knn = trainControl("repeatedcv", number = 10, repeats = 10,
                        allowParallel = TRUE, 
                        selectionFunction = "best")
grid_knn = expand.grid(k = c(1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 40, 50, 100))

# k-NN
RNGversion("4.2.2")
set.seed(123)
print("NOTE: Fitting k-NN model! (9/9)")
kNN_model = train(Tissue_type ~., 
                  data = train_set, 
                  method = "knn", 
                  preProcess = NULL,
                  trControl = ctrl_knn_k,
                  tuneGrid = grid_knn_k,
                  metric = "Accuracy")
saveRDS(kNN_model, "/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/knn_model.rds")
print("NOTE: Finished fitting k-NN model! (9/9)")
