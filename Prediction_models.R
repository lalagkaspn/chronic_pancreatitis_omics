# This script is used for the fitting prediction models using 10x10 repeated
# cross-validation and evaluation on the validation set. Then the final selected
# prediction model is evaluated on the microarray data.

# Preparation #####
setwd("/home/panagiotisnikolaos_lalagkas_student_uml_edu/chronic_pancreatitis_omics/")
## set path to installed packages so the HPC can find them
.libPaths("/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
#library(limma)
#library(org.Hs.eg.db)
library(dplyr)
library(stringr)
#library(tidyr)
#library(EnhancedVolcano)
#library(readr)
library(ggplot2)
#library(ggpubr)
library(openxlsx)
#library(enrichplot)
library(C50)
library(caret)
#library(irr)
#library(ipred)
#library(adabag)
#library(vcd)
#library(randomForest)
#library(pROC)
#library(stringr)
#library(kernlab)
#library(glmnet)
#library(MASS)
#library(LiblineaR)
#library(genefu)

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
gc()

# Formation of the test set
# Create an "X_Entrez" column for easy joins between datasets
selected_markers_pre$x_entrez = paste0("X_", selected_markers_pre$EntrezGene.ID)
x_entrez_microarray = rownames(exprs_microarray)
str_sub(x_entrez_microarray, 0, 0) = "X_"
rownames(exprs_microarray) = x_entrez_microarray

# Convert to data frame and add gene symbol info
exprs_microarray_df = as.data.frame(exprs_microarray)
exprs_microarray_df$x_entrez = x_entrez_microarray
exprs_microarray_df = exprs_microarray_df %>% inner_join(selected_markers_pre,
                                                 by = "x_entrez") %>%
  dplyr::select(-x_entrez, -EntrezGene.ID, -HGNC_Official)
rownames(exprs_microarray_df) = exprs_microarray_df$Gene.Symbol
exprs_microarray_df = exprs_microarray_df[lasso_markers$Predictor, ]

# Transpose to create a samples by genes data frame and add labels
test_set = as.data.frame(t(exprs_microarray_df))
test_set$Sample.ID = rownames(test_set)
test_set = test_set %>% inner_join(pheno_microarray, by = "Sample.ID") %>%
  dplyr::select(-Sample.ID, -Study)

# Parallelization #####

# This script requires extensive computational resources and it was submitted
# to a HPC, where 50 cores were employed along with x GB of RAM
# Output is saved as an .xlsx file
# Intermediate objects can be requested from the authors

# Enable parallel programming
library(doParallel)

system <- Sys.info()['sysname']
cores <- 100 # detectCores()

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
 
DT_Accuracy = data.frame(matrix(ncol=3,0))
colnames(DT_Accuracy) = c("Model", "Accuracy_train", "Accuracy_validation")

# C5.0 - k
RNGversion("4.2.2")
set.seed(123)
C50_model = train(Tissue_type ~., data = train_set,
                  method = "C5.0",
                  trControl = ctrl1,
                  metric = "Kappa",
                  tuneGrid = grid1,
                  na.action = "na.omit")

model_preds_validation = predict(C50_model, validation_set)
model_table_validation = table(model_preds_validation, validation_set$Tissue_type)
model_accuracy_validation = (model_table_validation[1] + model_table_validation[4])/
  sum(model_table_validation)

DT_Accuracy = rbind(DT_Accuracy, c("C5.0 - k", max(C50_model$results$Accuracy),
                                   model_accuracy_validation))

rm(model_preds_validation, model_table_validation, model_accuracy_validation)

# Bagging
RNGversion("4.2.2")
set.seed(123)
model_bag = bagging(Tissue_type ~., data = train_set,
                    nbagg = 100)

bag_pred_train = predict(model_bag, train_set)
bag_train_t = table(bag_pred_train$class, train_set$Tissue_type)
bagging_train_accuracy = 1 - (bag_train_t[1] + bag_train_t[4])/
  sum(bag_train_t)

bag_pred_validation = predict(model_bag, validation_set)
bag_validation_t = table(bag_pred_validation$class, validation_set$Tissue_type)
bagging_accuracy_validation = 1 - (bag_validation_t[1] + bag_validation_t[4])/
  sum(bag_validation_t)

DT_Accuracy = rbind(DT_Accuracy, c("100X Bagging", bagging_train_accuracy, 
                                   bagging_accuracy_validation))

rm(bag_pred_train, bag_train_t, bagging_train_accuracy,
   bag_pred_validation, bag_validation_t, bagging_accuracy_validation)

# Boosting
RNGversion("4.2.2")
set.seed(123)
boost = boosting(Tissue_type ~., data = train_set)

boost_train_pred = predict(boost, train_set)
cm_boost_train = boost_train_pred$confusion
boost_accuracy_train = 1 - (cm_boost_train[1] + cm_boost_train[4])/
  sum(cm_boost_train)

boost_validation_pred = predict(boost, validation_set)
cm_boost_validation = boost_validation_pred$confusion
boost_accuracy_validation = 1 - (cm_boost_validation[1] + cm_boost_validation[4])/
  sum(cm_boost_validation)

DT_Accuracy = rbind(DT_Accuracy, c("Boosting", boost_accuracy_train,
                                   boost_accuracy_validation))

rm(boost_train_pred, cm_boost_train, boost_accuracy_train,
   boost_validation_pred, cm_boost_validation, boost_accuracy_validation)

# Boosting with cross-validation on the whole dataset
# train
RNGversion("4.2.2")
set.seed(123)
adaboost_cv = boosting.cv(Tissue_type ~., data = cbind(rbind(train_set,
                                                          validation_set)),
                          mfinal = 100)

cm_adaboost = adaboost_cv$confusion
adaboost_accuracy = 1 - (cm_adaboost[1] + cm_adaboost[4])/sum(cm_adaboost)

DT_Accuracy = rbind(DT_Accuracy, c("Adaboost", adaboost_accuracy, adaboost_accuracy))

rm(cm_adaboost, adaboost_accuracy)

# Random Forest
RNGversion("4.2.2")
set.seed(123)
rf = train(Tissue_type ~., data = train_set,
           method = "rf",
           metric = "Kappa", trControl = ctrl, tuneGrid = grid2,
           na.action = "na.omit")

rf_preds_validation = predict(rf, validation_set)
rf_table_validation = table(rf_preds_validation, validation_set$Tissue_type)
rf_validation_accuracy = (rf_table_validation[1] + rf_table_validation[4])/
  sum(rf_table_validation)

DT_Accuracy = rbind(DT_Accuracy, c("Random Forest - k", max(rf$results$Accuracy),
                                   rf_validation_accuracy))

rm(rf_preds_validation, rf_table_validation, rf_validation_accuracy)

# RF with grid2 and ctrls, ROC
RNGversion("4.2.2")
set.seed(123)
comp_rf = train(Tissue_type ~., data = train_set,
                method = "rf",
                metric = "ROC", trControl = ctrls, tuneGrid = grid2,
                na.action = "na.omit")

comp_rf_preds_validation = predict(comp_rf, validation_set)
comp_rf_table_validation = table(comp_rf_preds_validation, validation_set$Tissue_type)
comp_rf_accuracy_validation = (comp_rf_table_validation[1] + comp_rf_table_validation[4])/
  sum(comp_rf_table_validation)

DT_Accuracy = rbind(DT_Accuracy, c("Random Forest - ROC", max(comp_rf$results$ROC),
                                   comp_rf_accuracy_validation))

rm(comp_rf_preds_validation, comp_rf_table_validation, comp_rf_accuracy_validation)

# C5.0 with grid_c50
RNGversion("4.2.2")
set.seed(123)
m_c50 = train(Tissue_type ~., data = train_set,
              method="C5.0", metric = "ROC",
              trControl = ctrls, tuneGrid = grid_c50, 
              na.action = "na.omit")

m_c50_preds_validation = predict(m_c50, validation_set)
m_c50_table_validation = table(m_c50_preds_validation, validation_set$Tissue_type)
m_c50_accuracy_validation = (m_c50_table_validation[1] + m_c50_table_validation[4])/
  sum(m_c50_table_validation)

DT_Accuracy = rbind(DT_Accuracy, c("C5.0 - ROC", max(m_c50$results$ROC),
                                   m_c50_accuracy_validation))

rm(m_c50_preds_validation, m_c50_table_validation, m_c50_accuracy_validation)

# Finalising
DT_models = list(C50_model, model_bag, boost, adaboost_cv, 
                 rf, comp_rf, m_c50)
names(DT_models) = c("C5.0 - k", "100X Bagging", "Boosting", "Boosting.cv",
                     "RForest - k", "RForest - ROC", 
                     "C5.0 - ROC")
DT[[i]] = DT_models
DT_Accuracy = DT_Accuracy[2:nrow(DT_Accuracy),]
addWorksheet(dtwb, names(Training)[i])
writeData(dtwb, names(Training)[i], DT_Accuracy)
rm(C50_model, model_bag, boost, adaboost_cv, rf, comp_rf, m_c50, 
   DT_Accuracy, DT_models) 
gc()
cat(paste0("Done with ", names(Training)[i], "\n"))

  