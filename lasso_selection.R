# Preparation #####
library(dplyr)
library(stringr)
library(openxlsx)
library(caret)
library(glmnet)

# This script requires extensive computational resources and it was submitted
# to a HPC, where 50 cores were employed along with x GB of RAM
# Output is saved as an .xlsx file
# Intermediate objects can be requested from the authors

# Enable parallel programming
library(doParallel)

system <- Sys.info()['sysname']
cores <- 50 # detectCores()

if (system == 'Windows') {
  cl <- makeCluster(getOption('cl.cores', cores), type='PSOCK')
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
} else {
  registerDoParallel(cores)
}

# Import data
selected_markers_pre = read.xlsx("genes_for_model_training.xlsx")
pheno_rnaseq = read.xlsx("DGEA/Pheno_RNAseq.xlsx") %>% dplyr::rename(Sample.ID = GEO_accession)
exprs_rnaseq = readRDS("DGEA/z_exprs_nonas_RNAseq.RDS")

# Training & Validation set splitting #####

# The RNA-seq data will be split into a training and a validation set

# Create an "X_Entrez" column for easy joins between datasets
selected_markers_pre$x_entrez = paste0("X_", selected_markers_pre$EntrezGene.ID)
x_entrez_rnaseq = rownames(exprs_rnaseq)
str_sub(x_entrez_rnaseq, 0, 0) = "X_"
rownames(exprs_rnaseq) = x_entrez_rnaseq

# Convert to data frame and add gene symbol info
exprs_rnaseq_df = as.data.frame(exprs_rnaseq)
exprs_rnaseq_df$x_entrez = x_entrez_rnaseq
exprs_rnaseq_df = exprs_rnaseq_df %>% inner_join(selected_markers_pre,
                                                 by = "x_entrez") %>%
  dplyr::select(-x_entrez, -EntrezGene.ID, -HGNC_Official)
rownames(exprs_rnaseq_df) = exprs_rnaseq_df$Gene.Symbol

# Transpose to create a samples by genes data frame and add labels
ml_input = as.data.frame(t(exprs_rnaseq_df))
ml_input$Sample.ID = rownames(ml_input)
ml_input = ml_input %>% inner_join(pheno_rnaseq, by = "Sample.ID") %>%
  dplyr::select(-Sample.ID)

# Stratified split in a 70%-30% manner
RNGversion("4.2.2")
set.seed(123)
tobesplit = ml_input %>%
  mutate(nr = row_number()) %>%
  dplyr::select(nr, everything()) %>%
  as.data.frame()
RNGversion("4.2.2")
set.seed(123)
train_set = tobesplit %>%
  group_by(Tissue_type, Study) %>%
  sample_frac(0.7) %>%
  as.data.frame()
validation_set = anti_join(tobesplit, train_set)

# Ungroup
train_set = train_set %>% dplyr::ungroup()
validation_set = validation_set %>% dplyr::ungroup()

# Convert outcome to factor type and remove study variable
train_set$Tissue_type = factor(train_set$Tissue_type, levels = c("chronic_pancreatitis", "normal", "tumor"),
                               labels = c("chronic_pancreatitis", "normal", "tumor"))
validation_set$Tissue_type = factor(validation_set$Tissue_type, levels = c("chronic_pancreatitis", "normal", "tumor"),
                                    labels = c("chronic_pancreatitis", "normal", "tumor"))
train_set = train_set %>% dplyr::select(-Study, -nr)
validation_set = validation_set %>% dplyr::select(-Study, -nr)

# Write out
write.xlsx(train_set, "ML_output/train_set.xlsx")
write.xlsx(validation_set, "ML_output/validation_set.xlsx")

# Logistic Regression: backward and regularised models #####

# Define control for repeated cross-validation
train_control <- trainControl(method = "repeatedcv", 
                              number = 10, 
                              repeats = 10,
                              # search = "grid",
                              allowParallel = TRUE
)

# Define the grid for tuning alpha and lambda
alpha_values <- seq(0, 1, length = 10)
lambda_values <- exp(seq(log(0.01), log(1), length = 10))

tuning_grid <- expand.grid(alpha = alpha_values, lambda = lambda_values)
tuning_grid_lasso <- expand.grid(alpha = 1, lambda = lambda_values)

# Train the model (computationally intensive - you can load the RDS, see below)
# Replace dataset and y with your dataset and response variable names

RNGversion("4.2.2")
set.seed(123)
model <- train(Tissue_type ~ ., 
               data = train_set, 
               method = "glmnet", 
               trControl = train_control,
               tuneGrid = tuning_grid,
               family = "multinomial",
               type.multinomial = "grouped")

saveRDS(model, "model.RDS")

# Load the RDS
# model = readRDS("model.rds")

# Print the basic model info
print(model)

# Print the best tune values
print(model$bestTune)
#     alpha   lambda
# 26 0.3333333 0.129155

# Print the basic model info for lasso only
lasso_results = model$results[model$results$alpha == 1,] %>%
  dplyr::arrange(desc(Accuracy))
print(lasso_results)

# Best lasso model
#     alpha   lambda
#      1   0.07742637

# The classification accuracy of the best lasso model and that of the
# best elastic net model overall, are almost equal:

# Model  alpha  lambda      accuracy    kappa
# Lasso: 1      0.07742637  0.5958757  0.2937924
# ENet:  0.33   0.129155    0.5985210  0.3160064577

# We can therefore use lasso for feature selection
# Fit lasso using glmnet, no cross validation and set lambda equal to the optimal lambda from caret
optimal_lasso_glmnet = glmnet(x = train_set %>% dplyr::select(-Tissue_type),
                              y = train_set$Tissue_type,
                              family = "multinomial", 
                              type.multinomial = "grouped",
                              alpha = 1,
                              lambda = model_lasso$bestTune$lambda)

# Gather coefficients in a data frame and keep the non-zero ones
coef_lasso_optimal = as.data.frame(as.matrix(optimal_lasso_glmnet$beta$chronic_pancreatitis))
coef_lasso_optimal$Predictor = rownames(coef_lasso_optimal)
coef_lasso_optimal$PDAC_coef = as.numeric(as.matrix(optimal_lasso_glmnet$beta$tumor)[,1])
coef_lasso_optimal$normal_coef = as.numeric(as.matrix(optimal_lasso_glmnet$beta$normal)[,1])
coef_lasso_optimal = coef_lasso_optimal %>% dplyr::select(Predictor, everything())
colnames(coef_lasso_optimal)[2] = "CP_coef"
coef_lasso_optimal = coef_lasso_optimal %>% 
  dplyr::arrange(desc(abs(CP_coef))) %>% # rank by CP coefficient
  dplyr::filter(abs(CP_coef) > 0)

# Coefficients for all three levels sum to zero
# coef_lasso_optimal$sum_coef = apply(coef_lasso_optimal[, c(2:4)], 1, function(x) sum(x))

write.xlsx(coef_lasso_optimal, "ML_output/Lasso_coefficients.xlsx")
