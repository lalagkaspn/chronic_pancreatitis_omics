# This script loads all models that were trained and evaluates their performance
# on the validation set in order to pick the optimal model

# Libraries #####
library(dplyr)
library(openxlsx)
library(pROC)

# Importing #####
model_filenames = list(`C5.0 - k` = "C5.0_k", `C5.0 - logLoss` = "C5.0_ROC",
                     `k-NN` = "knn_model", `Linear SVM` = "linear_svm",
                     `L2-reg SVM` = "L2_svm", `RBF SVM` = "rbf_svm_model",
                     `Random Forest - k` = "rf_k", `Random Forest - logLoss` = "rf_roc")

# import loop
models = list()
for (i in 1:length(model_filenames)) {
  models[[i]] = readRDS(paste0("prediction_models/", model_filenames[[i]], ".rds"))
}
names(models) = names(model_filenames)
rm(i, model_filenames); gc()

# Import training, validation and test set
train_set = read.xlsx("ML_output/train_set.xlsx")
validation_set = read.xlsx("ML_output/validation_set.xlsx")
glmnet_markers = read.xlsx("ML_output/glmnet_coefficients100.xlsx")

# Filter train and validation sets for the lasso markers
train_set = train_set[, c("Tissue_type", glmnet_markers$Predictor)]
validation_set = validation_set[, c("Tissue_type", glmnet_markers$Predictor)]
gc()

# Formation of the test set #####
selected_markers_pre = read.xlsx("genes_for_model_training.xlsx")
exprs_microarray = readRDS("DGEA/z_exprs_nonas_microarrays.RDS")
pheno_microarray = read.xlsx("DGEA/Pheno_microarrays.xlsx") %>% 
  dplyr::rename(Sample.ID = GEO_accession)

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
exprs_microarray_df = exprs_microarray_df[glmnet_markers$Predictor, ]

# Transpose to create a samples by genes data frame and add labels
test_set = as.data.frame(t(exprs_microarray_df))
test_set$Sample.ID = rownames(test_set)
test_set = test_set %>% dplyr::inner_join(pheno_microarray, by = "Sample.ID") %>%
  dplyr::select(-Sample.ID, -Study)

test_set[!colnames(test_set) %in% "Tissue_type"] =
  lapply(test_set[!colnames(test_set) %in% "Tissue_type"],
         function(x) as.numeric(x))

rm(selected_markers_pre, x_entrez_microarray, glmnet_markers, exprs_microarray_df)
gc()

# Model evaluations #####
# Set up a table to record results
Accuracy = as.data.frame(matrix(ncol=6, 0))
colnames(Accuracy) = c("model", "max_accuracy_during_training", "performance_train", "performance_val", 
                       "performance_train_%", "performance_val_%")

# C5.0 - k
train_pred = predict(models$`C5.0 - k`, train_set)
agree_train = ifelse(train_pred == train_set$Tissue_type, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(models$`C5.0 - k`, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Tissue_type, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

Accuracy = rbind(Accuracy, c("C5.0 - k", max(models$`C5.0 - k`$results$Accuracy),
                             accuracy_train, accuracy_validation,
                             paste0(round(100*accuracy_train, 2), "%"),
                             paste0(round(100*accuracy_validation, 2), "%")))
Accuracy = Accuracy[-1, ]

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# C5.0 - logLoss
train_pred = predict(models$`C5.0 - logLoss`, train_set)
agree_train = ifelse(train_pred == train_set$Tissue_type, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(models$`C5.0 - logLoss`, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Tissue_type, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

Accuracy = rbind(Accuracy, c("C5.0 - logLoss", max(models$`C5.0 - logLoss`$results$Accuracy),
                             accuracy_train, accuracy_validation,
                             paste0(round(100*accuracy_train, 2), "%"),
                             paste0(round(100*accuracy_validation, 2), "%")))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# Random Forest - k
train_pred = predict(models$`Random Forest - k`, train_set)
agree_train = ifelse(train_pred == train_set$Tissue_type, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(models$`Random Forest - k`, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Tissue_type, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

Accuracy = rbind(Accuracy, c("Random Forest - k", max(models$`Random Forest - k`$results$Accuracy),
                             accuracy_train, accuracy_validation,
                             paste0(round(100*accuracy_train, 2), "%"),
                             paste0(round(100*accuracy_validation, 2), "%")))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# Random Forest - logLoss
train_pred = predict(models$`Random Forest - logLoss`, train_set)
agree_train = ifelse(train_pred == train_set$Tissue_type, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(models$`Random Forest - logLoss`, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Tissue_type, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

Accuracy = rbind(Accuracy, c("Random Forest - logLoss", max(models$`Random Forest - logLoss`$results$Accuracy),
                             accuracy_train, accuracy_validation,
                             paste0(round(100*accuracy_train, 2), "%"),
                             paste0(round(100*accuracy_validation, 2), "%")))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# Linear SVM
train_pred = predict(models$`Linear SVM`, train_set)
agree_train = ifelse(train_pred == train_set$Tissue_type, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(models$`Linear SVM`, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Tissue_type, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

Accuracy = rbind(Accuracy, c("Linear SVM", max(models$`Linear SVM`$results$Accuracy),
                             accuracy_train, accuracy_validation,
                                     paste0(round(100*accuracy_train, 2), "%"),
                                     paste0(round(100*accuracy_validation, 2), "%")))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# L2-Linear SVM
train_pred = predict(models$`L2-reg SVM`, train_set)
agree_train = ifelse(train_pred == train_set$Tissue_type, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(models$`L2-reg SVM`, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Tissue_type, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

Accuracy = rbind(Accuracy, c("L2-reg SVM", max(models$`L2-reg SVM`$results$Accuracy),
                             accuracy_train, accuracy_validation,
                             paste0(round(100*accuracy_train, 2), "%"),
                             paste0(round(100*accuracy_validation, 2), "%")))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# RBF SVM
train_pred = predict(models$`RBF SVM`, train_set)
agree_train = ifelse(train_pred == train_set$Tissue_type, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(models$`RBF SVM`, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Tissue_type, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

Accuracy = rbind(Accuracy, c("RBF SVM", max(models$`RBF SVM`$results$Accuracy),
                             accuracy_train, accuracy_validation,
                             paste0(round(100*accuracy_train, 2), "%"),
                             paste0(round(100*accuracy_validation, 2), "%")))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# k-NN
train_pred = predict(models$`k-NN`, train_set)
agree_train = ifelse(train_pred == train_set$Tissue_type, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(models$`k-NN`, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Tissue_type, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

Accuracy = rbind(Accuracy, c("k-NN", max(models$`k-NN`$results$Accuracy),
                             accuracy_train, accuracy_validation,
                             paste0(round(100*accuracy_train, 2), "%"),
                             paste0(round(100*accuracy_validation, 2), "%")))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# AUC metric and ROC curves #####

# SVM models cannot produce predicted probabilities so they are excluded from 
# this step of the evaluation

# Get predicted probabilities
prob_preds_val_c50_k = predict(models$`C5.0 - k`, validation_set, type = "prob")
prob_preds_val_c50_logLoss = predict(models$`C5.0 - logLoss`, validation_set, type = "prob")
prob_preds_val_rf_k = predict(models$`Random Forest - k`, validation_set, type = "prob")
prob_preds_val_rf_logLoss = predict(models$`Random Forest - logLoss`, validation_set, type = "prob")
prob_preds_val_knn = predict(models$`k-NN`, validation_set, type = "prob")

# Joining
join_val_c50_k = cbind(prob_preds_val_c50_k, validation_set)
join_val_c50_k = join_val_c50_k[order(join_val_c50_k$chronic_pancreatitis),]
join_val_c50_logLoss = cbind(prob_preds_val_c50_logLoss, validation_set)
join_val_c50_logLoss = join_val_c50_logLoss[order(join_val_c50_logLoss$chronic_pancreatitis),]
join_val_rf_k = cbind(prob_preds_val_rf_k, validation_set)
join_val_rf_k = join_val_rf_k[order(join_val_rf_k$chronic_pancreatitis),]
join_val_rf_logLoss = cbind(prob_preds_val_rf_logLoss, validation_set)
join_val_rf_logLoss = join_val_rf_logLoss[order(join_val_rf_logLoss$chronic_pancreatitis),]
join_val_knn = cbind(prob_preds_val_knn, validation_set)
join_val_knn = join_val_knn[order(join_val_knn$chronic_pancreatitis),]

# ROC
model_roc_val_c50_k = multiclass.roc(predictor = join_val_c50_k$chronic_pancreatitis, 
                          response = as.character(join_val_c50_k$Tissue_type))
model_roc_val_c50_logLoss = multiclass.roc(predictor = join_val_c50_logLoss$chronic_pancreatitis, 
                                response = as.character(join_val_c50_logLoss$Tissue_type))
model_roc_val_rf_k = multiclass.roc(predictor = join_val_rf_k$chronic_pancreatitis, 
                         response = as.character(join_val_rf_k$Tissue_type))
model_roc_val_rf_logLoss = multiclass.roc(predictor = join_val_rf_logLoss$chronic_pancreatitis, 
                               response = as.character(join_val_rf_logLoss$Tissue_type))
model_roc_val_knn = multiclass.roc(predictor = join_val_knn$chronic_pancreatitis, 
                        response = as.character(join_val_knn$Tissue_type))

# AUC
AUC_val_c50_k = round(auc(model_roc_val_c50_k), 3) # 0.79
AUC_val_c50_logLoss = round(auc(model_roc_val_c50_logLoss), 3) # 0.79
AUC_val_rf_k = round(auc(model_roc_val_rf_k), 3) # 0.797
AUC_val_rf_logLoss = round(auc(model_roc_val_rf_logLoss), 3) # 0.797
AUC_val_knn = round(auc(model_roc_val_knn), 3) # 0.603
cat("multiclass AUC values:\n", paste("C5.0-k:", AUC_val_c50_k, "\n"),
    paste("C5.0-logLoss:", AUC_val_c50_logLoss, "\n"), 
    paste("Random Forest-k:", AUC_val_rf_k, "\n"),
    paste("Random Forest-logLoss:", AUC_val_rf_logLoss, "\n"),
    paste("k-NN:", AUC_val_knn, "\n"))

AUCs = c(0.79, 0.79, 0.797, 0.797, NA, NA, NA, 0.603)
Accuracy = cbind(Accuracy, AUCs)

# Write out an Excel sheet with the recorded performance metrics
write.xlsx(Accuracy, "ML_output/Model_performance.xlsx")

# Evaluate best model on test set #####

# Random Forests had the best classification accuracy and ROC in the validation
# set. We therefore pick the logLoss-optimized Random Forest as our best model
# to evaluate on the test set

final_model = models$`Random Forest - logLoss`
test_pred = predict(final_model, test_set)
agree_test = ifelse(test_pred == test_set$Tissue_type, 1, 0)
accuracy_test= sum(agree_test) / nrow(test_set)
prob_preds_test_final_model = predict(final_model, test_set, type = "prob")
join_test_final_model = cbind(prob_preds_test_final_model, test_set)
join_test_final_model = join_test_final_model[order(join_test_final_model$chronic_pancreatitis),]
model_roc_test_final_model = multiclass.roc(predictor = join_test_final_model$chronic_pancreatitis, 
                          response = as.character(join_test_final_model$Tissue_type))
AUC_test_final_model = round(auc(model_roc_test_final_model), 3) # 0.79

cat(paste0("The final model (Random Forest, logLoss-optimized) achieved ", 
             round(accuracy_test*100, 2), "% accuracy and AUC: ", 
             AUC_test_final_model, ", in the test set."))

x = data.frame(actual = test_set$Tissue_type, predicted = predict(final_model, test_set))
table(x)

rm(test_pred, agree_test)
