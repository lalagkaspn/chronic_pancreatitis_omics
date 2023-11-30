#!/usr/bin/env Rscript

# Preparation #####
.libPaths(c("/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2",
            "/usr/local/lib/R/site-library",
            "/usr/local/lib/R/library",
            "/modules/spack/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.4.0/r-4.2.0-3kitpfbxevyhxd2adiznenkjqqdbekzs/rlib/R/library"))
library(dplyr)
library(stringr, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
library(openxlsx, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
library(caret, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")
library(glmnet)
library(doParallel, lib.loc = "/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2")

# HPC test
print("Packages successfully loaded!")

# This script requires extensive computational resources and it was submitted
# to a HPC, where 50 cores were employed along with x GB of RAM
# Output is saved as an .xlsx file
# Intermediate objects can be requested from the authors

# Enable parallel programming
system <- Sys.info()['sysname']
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(cores)

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

# Convert gene columns to the numeric type
train_set[!colnames(train_set) %in% "Tissue_type"] =
  lapply(train_set[!colnames(train_set) %in% "Tissue_type"],
         function(x) as.numeric(x))
validation_set[!colnames(validation_set) %in% "Tissue_type"] =
  lapply(validation_set[!colnames(validation_set) %in% "Tissue_type"],
         function(x) as.numeric(x))

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

# # Load the RDS
# # model = readRDS("model.rds")

# Print the basic model info
print(model)

# Print the best tune values (a ridge regression model with 67.3% accuracy)
print(model$bestTune)
#     alpha   lambda
# 10      0        1

# Print the basic model info for lasso only
lasso_results = model$results[model$results$alpha == 1,] %>%
  dplyr::arrange(desc(Accuracy))
print(lasso_results)

# Best lasso model (54% accuracy)
#     alpha   lambda
#      1         0.1

# Fit glmnet, no cross validation and set lambda equal to the optimal lambda from caret
optimal_glmnet = glmnet(x = train_set %>% dplyr::select(-Tissue_type),
                              y = train_set$Tissue_type,
                              family = "multinomial", 
                              type.multinomial = "grouped",
                              alpha = model$bestTune$alpha,
                              lambda = model$bestTune$lambda)

# Gather coefficients in a data frame and keep the non-zero ones
coef_optimal = as.data.frame(as.matrix(optimal_glmnet$beta$chronic_pancreatitis))
coef_optimal$Predictor = rownames(coef_optimal)
coef_optimal$PDAC_coef = as.numeric(as.matrix(optimal_glmnet$beta$tumor)[,1])
coef_optimal$normal_coef = as.numeric(as.matrix(optimal_glmnet$beta$normal)[,1])
coef_optimal = coef_optimal %>% dplyr::select(Predictor, everything())
colnames(coef_optimal)[2] = "CP_coef"
coef_optimal = coef_optimal %>% 
  dplyr::arrange(desc(abs(CP_coef))) %>% # rank by CP coefficient
  dplyr::filter(abs(CP_coef) > 0)

# Coefficients for all three levels sum to zero
coef_optimal$sum_coef = apply(coef_optimal[, c(2:4)], 1, function(x) sum(x))
write.xlsx(coef_optimal, "ML_output/glmnet_coefficients.xlsx")

# Pick the top 100 features (in terms of CP coefficient magnitude)
coef_optimal_100 = coef_optimal[1:100, ]
write.xlsx(coef_optimal_100, "ML_output/glmnet_coefficients100.xlsx")
