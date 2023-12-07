## This script is used to download and pre-process series (microarrays) from GEO for 
## patients with chronic pancreatitis (normal and chronic pancreatitis tissues)

library(dplyr)
library(matrixStats)
library(GEOquery) # updated
library(org.Hs.eg.db)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(tidyr)
library(limma)
library(openxlsx)
library(EnhancedVolcano)
library(impute)
library(statmod)

## Set type for image compression based on operating system
## For macOS, X11 installation is required (link: https://www.xquartz.org/)
# if Windows OS
if (Sys.info()['sysname'] == "Windows") {
  type_compression = "windows"
} 
# if Linux
if (Sys.info()['sysname'] == "Linux") {
  type_compression = "windows"
}
# if macOS
if (Sys.info()['sysname'] == "Darwin") {
  type_compression = "cairo"
} 

##### Downloading data #####
# Vector with datasets to download from GEO
datasets = c("GSE143754", "GSE61166", "GSE71989", "GSE101462", "GSE77858")

# Download data
GEOsets = list()
for (i in 1:length(datasets)){
  GEOsets[[i]] = getGEO(datasets[i])
}; rm(i)
GEOsets = unlist(GEOsets)
names(GEOsets) = datasets; rm(datasets)

## --------- ##
## pData
## --------- ##

## Isolate pData
pdata = list()
for(i in 1:length(GEOsets)) {
  pdata[[i]] = pData(GEOsets[[i]])
}
names(pdata) = names(GEOsets); rm(i)

### Keep only necessary pdata
## - Sample names
## - Sample type (chronic pancreatitis, normal/non_tumor)
## - Available demographics (gender, age,)
filt_pdata = list()

## -- GSE143754 -- ##
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143754
# Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7648960/
# Comments:
#   - 26 Pancreatic headmass tissue samples
#     - 6 Chronic Pancreatitis samples
#     - 11 PDAC samples
#     - 9 Adjacent Normal samples

# Select informative data only
filt_pdata[["GSE143754"]] = pdata$GSE143754 %>%
  dplyr::select(GEO_accession = geo_accession,
                Patient_ID = title,
                Platform = platform_id,
                Gender = "Sex:ch1",
                Age = "age:ch1",
                Tissue_type = "disease state:ch1")

# Change Tissue_type to chronic pancreatitis and normal
filt_pdata$GSE143754$Tissue_type <- case_when(
  filt_pdata$GSE143754$Tissue_type == "Tumor" ~ "tumor",
  filt_pdata$GSE143754$Tissue_type == "Adjacent Normal" ~ "normal",
  filt_pdata$GSE143754$Tissue_type == "Chronic Pancreatitis" ~ "chronic_pancreatitis",
  TRUE ~ filt_pdata$GSE143754$Tissue_type
)

## NOTE: We keep all samples (including PDAC) because we want to adjust for them in the DGEA between normal and CPs

# Clear patient ID
filt_pdata[["GSE143754"]]$Patient_ID = gsub("Benign Tissue, Biological Replicate ", "", filt_pdata[["GSE143754"]]$Patient_ID)
filt_pdata[["GSE143754"]]$Patient_ID = gsub("Malignant Tissue, Biological Replicate ", "", filt_pdata[["GSE143754"]]$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE143754$Tissue_type = factor(x = filt_pdata$GSE143754$Tissue_type,
                                          levels = c("chronic_pancreatitis", "normal", "tumor"),
                                          labels = c("chronic_pancreatitis", "normal", "tumor"))
filt_pdata$GSE143754$Gender = factor(x = filt_pdata$GSE143754$Gender,
                                     levels = c("Female","Male"),
                                     labels = c("female","male"))
# Transform age to numeric
filt_pdata$GSE143754$Age = as.numeric(filt_pdata$GSE143754$Age)

## -- GSE61166 -- ##
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61166
# Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4742134/ (5 in total)
# Comments:
#   - 12 samples (12 patients)
#     - 4 Pancreatic tumors with type3C Diabetes Mellitus samples (A1-A4)
#     - 4 Pancreatic tumors without type3C Diabetes Mellitus samples (B1-B4)
#     - 4 Pancreatitis samples
# NOTE: I keep only the pancreatitis samples
#### samples are pancreatitis, not Chronic Pancreatitis

# Select informative data only
filt_pdata[["GSE61166"]] = pdata$GSE61166 %>%
  dplyr::select(GEO_accession = geo_accession,
                Patient_ID = title,
                Platform = platform_id,
                Gender = "gender:ch1",
                Age = "age:ch1",
                Tissue_type = "disease status:ch1")

# Change Tissue_type to chronic pancreatitis and normal
filt_pdata$GSE61166$Tissue_type <- case_when(
  filt_pdata$GSE61166$Tissue_type == "pancreatitis" ~ "chronic_pancreatitis",
  filt_pdata$GSE61166$Tissue_type == "pancreatic tumor" ~ "tumor",
  TRUE ~ filt_pdata$GSE61166$Tissue_type
)

# Clear patient ID
filt_pdata$GSE61166$Patient_ID = gsub("Pancreatitis_tissues_of_Patient", "", filt_pdata$GSE61166$Patient_ID)
filt_pdata$GSE61166$Patient_ID = gsub("Pancreatic_tumors_of_Patient", "", filt_pdata$GSE61166$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE61166$Tissue_type = factor(x = filt_pdata$GSE61166$Tissue_type,
                                         levels = c("chronic_pancreatitis", "normal", "tumor"),
                                         labels = c("chronic_pancreatitis", "normal", "tumor"))
filt_pdata$GSE61166$Gender = factor(x = filt_pdata$GSE61166$Gender,
                                    levels = c("female","male"),
                                    labels = c("female","male"))

# Transform age to numeric
filt_pdata$GSE61166$Age = as.numeric(filt_pdata$GSE61166$Age)

## -- GSE71989 -- ##
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71989
# Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5288176/
# Comments:
#   - 22 samples (22 patients)
#     - 13 PDAC samples
#     - 1 chronic pancreatitis sample
#     - 8 normal pancreatic samples

# Select informative data only
# No AGE and GENDER information provided
filt_pdata[["GSE71989"]] = pdata$GSE71989 %>%
  dplyr::select(GEO_accession = geo_accession,
                Patient_ID = title,
                Platform = platform_id,
                Tissue_type = "tissue subtype:ch1")

# Change Tissue_type to chronic pancreatitis and normal
filt_pdata$GSE71989$Tissue_type <- case_when(
  filt_pdata$GSE71989$Tissue_type == "CP" ~ "chronic_pancreatitis",
  filt_pdata$GSE71989$Tissue_type == "normal pancreatic tissue" ~ "normal",
  filt_pdata$GSE71989$Tissue_type == "PDAC" ~ "tumor",
  TRUE ~ filt_pdata$GSE71989$Tissue_type
)

# Clear patient ID
filt_pdata$GSE71989$Patient_ID = gsub(", human normal pancreatic tissue", "", filt_pdata$GSE71989$Patient_ID)
filt_pdata$GSE71989$Patient_ID = gsub(", human Chronic Pancreatitis tissue", "", filt_pdata$GSE71989$Patient_ID)
filt_pdata$GSE71989$Patient_ID = gsub(", human PDAC tissue", "", filt_pdata$GSE71989$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE71989$Tissue_type = factor(x = filt_pdata$GSE71989$Tissue_type,
                                         levels = c("chronic_pancreatitis", "normal", "tumor"),
                                         labels = c("chronic_pancreatitis", "normal", "tumor"))
# NOTE: No age or gender information was provided

## -- GSE101462 -- ##
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101462
# Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5895731/
# Comments:
#   - 24 samples
#     - 3 fresh-frozen (FF) normal samples
#     - 5 formalin-fixed paraffin embedded (FFPE) normal sampls
#     - 3 fresh-frozen (FF) PDAC samples
#     - 3 formalin-fixed paraffin embedded (FFPE) PDAC samples
#     - 10 formalin-fixed paraffin embedded (FFPE) pancreatitis samples

#### NOTE: samples are pancreatitis, not Chronic Pancreatitis
#### Also, they used different tissue_storage (FF or FFPE) --> how does this affect our analysis?

# Select informative data only
filt_pdata[["GSE101462"]] = pdata$GSE101462 %>%
  dplyr::select(GEO_accession = geo_accession,
                Patient_ID = title,
                Platform = platform_id,
                Tissue_type = "tissue type:ch1",
                Tissue_storage = "tissue storage:ch1")

# Change Tissue_type to chronic pancreatitis and normal
filt_pdata$GSE101462$Tissue_type = case_when(
  filt_pdata$GSE101462$Tissue_type == "pancreatitis" ~ "chronic_pancreatitis",
  filt_pdata$GSE101462$Tissue_type == "normal" ~ "normal",
  filt_pdata$GSE101462$Tissue_type == "PDAC" ~ "tumor",
  TRUE ~ filt_pdata$GSE101462$Tissue_type
)

# Change Tissue_storage
filt_pdata$GSE101462$Tissue_storage <- case_when(
  filt_pdata$GSE101462$Tissue_storage == "fresh frozen (FF)" ~ "fresh_frozen",
  filt_pdata$GSE101462$Tissue_storage == "formalin-fixed paraffin embedded (FFPE)" ~ "formalin_fixed_paraffin_embedded",
  TRUE ~ filt_pdata$GSE101462$Tissue_storage
)

# Transform to factors with consistent universal levels
filt_pdata$GSE101462$Tissue_type = factor(x = filt_pdata$GSE101462$Tissue_type,
                                          levels = c("chronic_pancreatitis", "normal", "tumor"),
                                          labels = c("chronic_pancreatitis", "normal", "tumor"))
# NOTE: no age or gender information was provided

## -- GSE77858 -- ##
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77858
# Paper: -
# Comments:
#   - 85 samples
#     - 77 tumors
#     - 3 normal samples
#     - 5 pancraetitis samples

# Select informative data only
filt_pdata[["GSE77858"]] = pdata$GSE77858 %>%
  dplyr::select(GEO_accession = geo_accession,
                Patient_ID = title,
                Platform = platform_id,
                Tissue_type = "morphology:ch2")

# Change Tissue_type to chronic pancreatitis and normal
filt_pdata$GSE77858$Tissue_type <- case_when(
  filt_pdata$GSE77858$Tissue_type == "Normal" ~ "normal",
  filt_pdata$GSE77858$Tissue_type == "Panreatitis" ~ "chronic_pancreatitis",
  filt_pdata$GSE77858$Tissue_type == "Pancreatitis" ~ "chronic_pancreatitis",
  filt_pdata$GSE77858$Tissue_type == "Tumor" ~ "tumor",
  TRUE ~ filt_pdata$GSE77858$Tissue_type
)

# Clear patient ID
filt_pdata[["GSE77858"]]$Patient_ID = gsub("PancTuRef2 vs. ", "", filt_pdata[["GSE77858"]]$Patient_ID)
filt_pdata[["GSE77858"]]$Patient_ID = gsub("PancTuRef vs. ", "", filt_pdata[["GSE77858"]]$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE77858$Tissue_type = factor(x = filt_pdata$GSE77858$Tissue_type,
                                         levels = c("chronic_pancreatitis", "normal", "tumor"),
                                         labels = c("chronic_pancreatitis", "normal", "tumor"))

## -- full_pdata -- ##
# Keep only information for Study, GEO_accession and Tissue_type
# Useful for QC analysis
pdata143754 = filt_pdata$GSE143754 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE143754")
pdata61166 = filt_pdata$GSE61166 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE61166")
pdata71989 = filt_pdata$GSE71989 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE71989")
pdata101462 = filt_pdata$GSE101462 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE101462")
pdata77858 = filt_pdata$GSE77858 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE77858")

full_pdata = rbind(pdata143754, pdata61166, pdata71989, pdata101462, pdata77858)
rownames(full_pdata) = full_pdata$GEO_accession
rm(pdata143754, pdata61166, pdata71989, pdata101462, pdata77858)
table(full_pdata$Tissue_type)

## Save pdata
openxlsx::write.xlsx(full_pdata, "DGEA/Pheno_microarrays.xlsx", overwrite = TRUE)

## -------------------------------
## Splitting into Train, Val, Test
## -------------------------------

# Stratified splitting into training, validation and test #####
RNGversion("4.2.2")
set.seed(123)
tobesplit = full_pdata %>%
  mutate(nr = row_number()) %>%
  dplyr::select(nr, everything()) %>%
  as.data.frame()
RNGversion("4.2.2")
set.seed(123)
train_set = tobesplit %>%
  group_by(Tissue_type, Study) %>%
  sample_frac(0.7) %>%
  as.data.frame()
RNGversion("4.2.2")
set.seed(123)
validation_set = anti_join(tobesplit, train_set) %>%
  group_by(Tissue_type, Study) %>%
  sample_frac(0.667) %>%
  as.data.frame()
test_set = anti_join(tobesplit, as.data.frame(rbind(train_set, validation_set)))

train_set = train_set %>% dplyr::ungroup()
validation_set = validation_set %>% dplyr::ungroup()
test_set = test_set %>% dplyr::ungroup()

train_samples = train_set$GEO_accession
validation_samples = validation_set$GEO_accession
test_samples = test_set$GEO_accession
trainval_samples = c(train_samples, validation_samples)

# Create lists of GEO sets for train-val and test sets
GEOsets_trainval = list(); GEOsets_test = list()
for (i in 1:length(GEOsets)){
  GEOsets_trainval[[i]] = GEOsets[[i]][, intersect(trainval_samples, colnames(GEOsets[[i]]))]
  GEOsets_test[[i]] = GEOsets[[i]][, intersect(test_samples, colnames(GEOsets[[i]]))]
}
names(GEOsets_test) = names(GEOsets_trainval) = names(GEOsets)
rm(i); gc()

## ---------------------------
## Expression data
## ---------------------------

esets = list()
for (i in 1:length(GEOsets)) {
  esets[[i]] = exprs(GEOsets[[i]])
  esets[[i]] = as.data.frame(esets[[i]])
}; rm(i)
names(esets) = names(GEOsets)

esets_trainval = list()
for (i in 1:length(GEOsets_trainval)) {
  esets_trainval[[i]] = exprs(GEOsets_trainval[[i]])
  esets_trainval[[i]] = as.data.frame(esets_trainval[[i]])
}; rm(i)
names(esets_trainval) = names(GEOsets_trainval)

esets_test = list()
for (i in 1:length(GEOsets_test)) {
  esets_test[[i]] = exprs(GEOsets_test[[i]])
  esets_test[[i]] = as.data.frame(esets_test[[i]])
}; rm(i)
names(esets_test) = names(GEOsets_test)

## Calculate NAs values
na_esets = c()
for(i in 1:length(esets)) {
  na_esets[i] = sum(is.na(esets[[i]]))
}; rm(i)
names(na_esets) = names(GEOsets)
na_esets 

na_esets_trainval = c()
for(i in 1:length(esets_trainval)) {
  na_esets_trainval[i] = sum(is.na(esets_trainval[[i]]))
}; rm(i)
names(na_esets_trainval) = names(GEOsets_trainval)
na_esets_trainval 

na_esets_test = c()
for(i in 1:length(esets_test)) {
  na_esets_test[i] = sum(is.na(esets_test[[i]]))
}; rm(i)
names(na_esets_test) = names(GEOsets_test)
na_esets_test 

# Missing values in 2 studies 
# GSE101462: 0 (full sets)
# GSE77858: 0 (full sets)
# GSE101462: 0 (train-val sets)
# GSE77858: 0 (train-val sets)
# GSE101462: 1 (test sets)
# GSE77858: 325 (test sets)

## We will impute missing values using k-NN imputation. 
## However, GSE77858 has a lot of missing values. 
## We remove rows with more than 25% missing values because they will negatively affect imputation
GEOsets[["GSE101462"]] = GEOsets[["GSE101462"]][rowSums(is.na(GEOsets[["GSE101462"]]@assayData[["exprs"]]))/
                                                  length(colnames(GEOsets[["GSE101462"]]@assayData[["exprs"]])) < 0.25, ]
GEOsets[["GSE77858"]] = GEOsets[["GSE77858"]][rowSums(is.na(GEOsets[["GSE77858"]]@assayData[["exprs"]]))/
                                                length(colnames(GEOsets[["GSE77858"]]@assayData[["exprs"]])) < 0.25, ]

# For consistency, we remove these rows from the train-val subset and the test subset
GEOsets_trainval[["GSE101462"]] = GEOsets_trainval[["GSE101462"]][rownames(GEOsets[["GSE101462"]]), ]
GEOsets_test[["GSE101462"]] = GEOsets_test[["GSE101462"]][rownames(GEOsets[["GSE101462"]]), ]
GEOsets_trainval[["GSE77858"]] = GEOsets_trainval[["GSE77858"]][rownames(GEOsets[["GSE77858"]]), ]
GEOsets_test[["GSE77858"]] = GEOsets_test[["GSE77858"]][rownames(GEOsets[["GSE77858"]]), ]

## -- k-NN imputation -- ##
# Firstly, we perform k-NN imputation on the train-val subsets
# GSE101462
RNGversion("4.2.2")
eset101462 = GEOsets_trainval[["GSE101462"]]@assayData[["exprs"]]
eset101462 = impute.knn(eset101462, k = 10, maxp = nrow(eset101462),
                        rng.seed = 123)
eset101462 = eset101462[["data"]]

# GSE102238
RNGversion("4.2.2")
eset77858 = GEOsets_trainval[["GSE77858"]]@assayData[["exprs"]]
eset77858 = impute.knn(eset77858, k = 10, maxp = nrow(eset77858),
                       rng.seed = 123)
eset77858 = eset77858[["data"]]

# Update esets object
esets_trainval[["GSE77858"]] = as.data.frame(eset77858); rm(eset77858)
esets_trainval[["GSE101462"]] = as.data.frame(eset101462); rm(eset101462)
esets_test[["GSE77858"]] = esets_test[["GSE77858"]][rownames(esets_trainval[["GSE77858"]]), ]
esets_test[["GSE101462"]] = esets_test[["GSE101462"]][rownames(esets_trainval[["GSE101462"]]), ]

# In order to impute the test samples for each study, we do the following in order:
# 1. find the 10 nearest neighbors for every gene in our imputed training and validation
#    samples of a study.
# 2. for each sample, impute the missing value in a gene using 10-NN averaging of the 
#    values in that gene's neighbours, as these were determined in step 1.
# 3. repeat for all studies

library(data.table)

neighbors_list <- lapply(esets_trainval, function(eset) {
  dists <- as.matrix(dist(eset, method = "euclidean"))
  neighbors <- apply(dists, 1, function(row) {
    head(order(row), 11)[-1] # Excluding the gene itself, take first 10
  })
  colnames(neighbors) <- rownames(dists)
  as.data.frame(neighbors)
})
names(neighbors_list) <- names(esets_trainval)
gc()

# Imputing the sets that correspond to neighbors_list:
esets_test_new <- list()

for (i in 1:length(esets_test)) {
  eset <- esets_test[[i]]
  imp_map <- neighbors_list[[names(esets_test)[i]]]
  imputed_data <- apply(eset, 2, function(sample) {
    sapply(seq_along(sample), function(k) {
      if (is.na(sample[k])) {
        neighbors <- imp_map[, rownames(eset)[k]]
        neighbor_values <- sample[neighbors]
        if (all(is.na(neighbor_values))) {
          mean(sample, na.rm = TRUE)
        } else {
          mean(neighbor_values, na.rm = TRUE)
        }
      } else {
        sample[k]
      }
    })
  })
  esets_test_new[[i]] <- imputed_data
  colnames(esets_test_new[[i]]) <- colnames(eset)
  rownames(esets_test_new[[i]]) <- rownames(eset)
}

rm(i); gc()
names(esets_test_new) <- names(esets_test)

# Convert to data frames
esets_test_new <- lapply(esets_test_new, as.data.frame)

##### Annotation esets with Entrez ID's #####
# GSE143754
# Platform [HTA-2_0] Affymetrix Human Transcriptome Array 2.0 [transcript (gene) version] -->
# gene data for each probe were downloaded from: https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_legacy/unconverted_chips/HTA_2_0.chip
hta2_0_probe_to_gene = data.table::fread("platforms_probes_to_genes/HTA_2_0.chip.txt")
fdata143754 = fData(GEOsets_trainval$GSE143754) %>%
  dplyr::select(probeset_id) %>%
  filter(probeset_id %in% hta2_0_probe_to_gene$`Probe Set ID`) %>%
  left_join(hta2_0_probe_to_gene, by = c("probeset_id" = "Probe Set ID")) %>%
  na.omit() %>%
  dplyr::select(ID = probeset_id, ENTREZ_GENE_ID = "Entrez Gene") %>%
  distinct()

## For probes that match to the same EntrezID, calculate variance across all 
## samples and keep the probe with max variance

# train-val samples
esets_trainval[["GSE143754"]]$variance = rowVars(as.matrix(esets_trainval[["GSE143754"]]))
esets_trainval[["GSE143754"]] = esets_trainval[["GSE143754"]] %>%
  mutate(ID = rownames(esets_trainval[["GSE143754"]])) %>%
  inner_join(fdata143754) %>%
  # dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)

# test samples
esets_test_new[["GSE143754"]] = esets_test_new[["GSE143754"]] %>%
  mutate(ID = rownames(esets_test_new[["GSE143754"]])) %>%
  inner_join(esets_trainval[["GSE143754"]] %>% dplyr::select(ID, ENTREZ_GENE_ID),
             by = "ID") %>%
  dplyr::select(-ID)
esets_trainval[["GSE143754"]] = esets_trainval[["GSE143754"]] %>%
  dplyr::select(-ID)
rm(fdata143754, hta2_0_probe_to_gene)

# GSE61166 
# need to match gene sequence to gene ID
# No information for corresponding genes provided in GEOset
# Three separate .fasta files were prepared based on the probe sequences available at the GPL file.
# BLASTN alignment (against Nucleotide sequences, "nt") was performed to annotate probes with RefSeq gene IDs: https://blast.ncbi.nlm.nih.gov/Blast.cgi

# Prepare files for BLASTN
gse61166_sequences = fData(GEOsets_trainval$GSE61166) %>% dplyr::select(SEQUENCE) %>% distinct()
gse61166_sequences = gse61166_sequences %>% mutate(seq_id = paste0(">sequence_", seq(1:nrow(gse61166_sequences))), .before = SEQUENCE)
gse61166_sequences = do.call(rbind, lapply(seq(nrow(gse61166_sequences)), function(i) t(gse61166_sequences[i, ])))
gse61166_sequences = data.frame(gse61166_sequences)
gse61166_sequences = data.frame(sequences = gse61166_sequences$ASHGA5P000001)
gse61166_sequences_1 = data.frame(sequences = gse61166_sequences[1:30000, ]) # the number of rows included in each file is based on the capacity of BLASTN for a run
gse61166_sequences_2 = data.frame(sequences = gse61166_sequences[30001:60000, ])
gse61166_sequences_3 = data.frame(sequences = gse61166_sequences[60001:91178, ])
# save files
data.table::fwrite(gse61166_sequences_1, "data/GSE61166_sequences_1.fasta", sep = "\t", row.names = FALSE, col.names = FALSE) ; rm(gse61166_sequences_1)
data.table::fwrite(gse61166_sequences_2, "data/GSE61166_sequences_2.fasta", sep = "\t", row.names = FALSE, col.names = FALSE) ; rm(gse61166_sequences_2)
data.table::fwrite(gse61166_sequences_3, "data/GSE61166_sequences_3.fasta", sep = "\t", row.names = FALSE, col.names = FALSE) ; rm(gse61166_sequences_3)

## RefSeq to EntrezID reference data frame
ref = org.Hs.egREFSEQ2EG
mapped_genes_official = mappedkeys(ref)
ref_df = as.data.frame(ref[mapped_genes_official]) ; rm(ref, mapped_genes_official)
ref_df = ref_df %>% dplyr::rename(ENTREZ_GENE_ID = gene_id, RefSeq = accession)
length(unique(ref_df$RefSeq)) == length(ref_df$RefSeq) # TRUE --> No duplicates in RefSeq

# Read csv results from BLASTN
blastn_1 = read.csv(file = "data/Sequences_1.csv", header = FALSE)
blastn_2 = read.csv(file = "data/Sequences_2.csv", header = FALSE)
blastn_3 = read.csv(file = "data/Sequences_3.csv", header = FALSE)
blastn = rbind(blastn_1, blastn_2, blastn_3)
blastn$V2 = gsub("\\..", "", blastn$V2) # RefSeq column: keep main nomenclature - ignore variants
blastn = blastn %>%
  dplyr::select(V1, V2, V3, V11) %>%
  dplyr::filter(V3 == 100.00) %>% # keep rows with 100% confidence aligned score
  distinct()
rm(blastn_1, blastn_2, blastn_3)
colnames(blastn) = c("probe", "RefSeq", "Alignment_perc", "E_value")
blastn$Alignment_perc = as.numeric(blastn$Alignment_perc)
blastn$E_value = as.numeric(blastn$E_value)
blastn = blastn[order(blastn$probe, 1/blastn$E_value),]

# Annotate sequences with probe IDs
gse61166_sequences = fData(GEOsets_trainval$GSE61166) %>% dplyr::select(SEQUENCE) %>% distinct()
gse61166_sequences = gse61166_sequences %>% mutate(seq_id = paste0("sequence_", seq(1:nrow(gse61166_sequences))), .before = SEQUENCE)
gse61166_sequences$probe_ID = rownames(gse61166_sequences)
gse61166_sequences = gse61166_sequences %>% dplyr::select(seq_id, probe_ID) %>% distinct()
blastn = left_join(blastn, gse61166_sequences, by = c("probe" = "seq_id"))
rm(gse61166_sequences)

# Map through org.Hs.eg.db
mapped_blastn = inner_join(blastn, ref_df, by = "RefSeq") %>%
  dplyr::select(probe = probe_ID, RefSeq, ENTREZ_GENE_ID) %>%
  distinct() %>%
  dplyr::select(probe, ENTREZ_GENE_ID) # we do not use distinct() here for filtering purposes
rm(ref_df)

# Probes which match to multiple ID's: 
# Check if >50% of Entrez ID's mapping to a probe are actually a unique Entrez ID (aka a probe matches to one Entrez ID)
# If yes, map the probe to that Entrez ID. If not, discard the probe
# Process described here and suggested by Ensembl: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431719/

mapped_blastn = mapped_blastn[order(mapped_blastn$probe),]
mapped_blastn$unique_Entrez_perc = NA

# Add a column to filter as described previously
mapped_blastn = mapped_blastn %>% 
  group_by(probe) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(probe, ENTREZ_GENE_ID) %>%
  mutate(matching_entrez = n()) %>%
  mutate(unique_Entrez_perc = (matching_entrez / total) * 100) %>%
  ungroup() %>%
  dplyr::select(-total, -matching_entrez)

# We now keep everything with unique Entrez mapping percentage over 50%.
# That is guaranteed to keep one Entrez ID for each probe and discard probes 
# for which only lower mapping percentages exist.

mapped_blastn_filt = mapped_blastn %>%
  dplyr::filter(unique_Entrez_perc > 50) %>%
  dplyr::select(probe, ENTREZ_GENE_ID) %>%
  distinct()

length(which(duplicated(mapped_blastn_filt$probe)))
# 0
length(which(duplicated(mapped_blastn_filt$ENTREZ_GENE_ID)))
# 11,869

# The numbers above mean each probe is mapped to a unique Entrez ID, but multiple probes
# may map to the same ID. For each of these probes we calculate the variance and keep the one 
# with the max variance, as we did previously.

# GSE61166
# train-val samples
esets_trainval[["GSE61166"]]$variance = rowVars(as.matrix(esets_trainval[["GSE61166"]]))
esets_trainval[["GSE61166"]] = esets_trainval[["GSE61166"]] %>%
  mutate(probe = rownames(esets_trainval[["GSE61166"]])) %>%
  inner_join(mapped_blastn_filt) %>%
  # dplyr::select(-probe) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)
esets_trainval[["GSE61166"]]$ENTREZ_GENE_ID = as.integer(esets_trainval[["GSE61166"]]$ENTREZ_GENE_ID)

# test samples
esets_test_new[["GSE61166"]] = esets_test_new[["GSE61166"]] %>%
  mutate(probe = rownames(esets_test_new[["GSE61166"]]))%>%
  inner_join(esets_trainval[["GSE61166"]] %>% dplyr::select(probe, ENTREZ_GENE_ID),
             by = "probe") %>%
  dplyr::select(-probe)
esets_trainval[["GSE61166"]] = esets_trainval[["GSE61166"]] %>%
  dplyr::select(-probe)
rm(blastn, mapped_blastn, mapped_blastn_filt)
esets_test_new[["GSE61166"]]$ENTREZ_GENE_ID = as.integer(esets_test_new[["GSE61166"]]$ENTREZ_GENE_ID)

# GSE71989
fdata71989 = fData(GEOsets_trainval$GSE71989) %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)

# train-val samples
esets_trainval[["GSE71989"]]$variance = rowVars(as.matrix(esets_trainval[["GSE71989"]]))
esets_trainval[["GSE71989"]] = esets_trainval[["GSE71989"]] %>%
  mutate(ID = rownames(esets_trainval[["GSE71989"]])) %>%
  inner_join(fdata71989) %>%
  # dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)
esets_trainval[["GSE71989"]]$ENTREZ_GENE_ID = as.integer(esets_trainval[["GSE71989"]]$ENTREZ_GENE_ID)

# test samples
esets_test_new[["GSE71989"]] = esets_test_new[["GSE71989"]] %>%
  mutate(ID = rownames(esets_test_new[["GSE71989"]]))%>%
  inner_join(esets_trainval[["GSE71989"]] %>% dplyr::select(ID, ENTREZ_GENE_ID),
             by = "ID") %>%
  dplyr::select(-ID)
esets_trainval[["GSE71989"]] = esets_trainval[["GSE71989"]] %>%
  dplyr::select(-ID)
rm(fdata71989)
esets_test_new[["GSE71989"]]$ENTREZ_GENE_ID = as.integer(esets_test_new[["GSE71989"]]$ENTREZ_GENE_ID)

# GSE101462
fdata101462 = fData(GEOsets$GSE101462) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)

# train-val samples
esets_trainval[["GSE101462"]]$variance = rowVars(as.matrix(esets_trainval[["GSE101462"]]))
esets_trainval[["GSE101462"]] = esets_trainval[["GSE101462"]] %>%
  mutate(ID = rownames(esets_trainval[["GSE101462"]])) %>%
  inner_join(fdata101462) %>%
  # dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)

# test samples
esets_test_new[["GSE101462"]] = esets_test_new[["GSE101462"]] %>%
  mutate(ID = rownames(esets_test_new[["GSE101462"]]))%>%
  inner_join(esets_trainval[["GSE101462"]] %>% dplyr::select(ID, ENTREZ_GENE_ID),
             by = "ID") %>%
  dplyr::select(-ID)
esets_trainval[["GSE101462"]] = esets_trainval[["GSE101462"]] %>%
  dplyr::select(-ID)
rm(fdata101462)
esets_test_new[["GSE101462"]]$ENTREZ_GENE_ID = as.integer(esets_test_new[["GSE101462"]]$ENTREZ_GENE_ID)

# GSE77858
fdata77858 = fData(GEOsets$GSE77858) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = GENE) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)

# train-val samples
esets_trainval[["GSE77858"]]$variance = rowVars(as.matrix(esets_trainval[["GSE77858"]]))
esets_trainval[["GSE77858"]] = esets_trainval[["GSE77858"]] %>%
  mutate(ID = rownames(esets_trainval[["GSE77858"]])) %>%
  inner_join(fdata77858) %>%
  # dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)

# test samples
esets_test_new[["GSE77858"]] = esets_test_new[["GSE77858"]] %>%
  mutate(ID = rownames(esets_test_new[["GSE77858"]]))%>%
  inner_join(esets_trainval[["GSE77858"]] %>% dplyr::select(ID, ENTREZ_GENE_ID),
             by = "ID") %>%
  dplyr::select(-ID)
esets_trainval[["GSE77858"]] = esets_trainval[["GSE77858"]] %>%
  dplyr::select(-ID)
rm(fdata77858)
esets_test_new[["GSE77858"]]$ENTREZ_GENE_ID = as.integer(esets_test_new[["GSE77858"]]$ENTREZ_GENE_ID)

##### Calculate NAs values in the train-val and test sets #####
na_esets_trainval = c()
for(i in 1:length(esets_trainval)) {
  na_esets_trainval[i] = sum(is.na(esets_trainval[[i]]))
}; rm(i)
names(na_esets_trainval) = names(GEOsets_trainval)
na_esets_trainval # No NA values exist in the expression files (train-val)
rm(na_esets_trainval, na_esets)

na_esets_test_new = c()
for(i in 1:length(esets_test_new)) {
  na_esets_test_new[i] = sum(is.na(esets_test_new[[i]]))
}; rm(i)
names(na_esets_test_new) = names(GEOsets_test)
na_esets_test_new # No NA values exist in the expression files (test)
rm(na_esets_test_new)

##### z-score-transformation #####
# KBZ transformation method (https:://www.biostars.org/p/283083/)

calculate_z <- function(eset) {
  df <- as.data.frame(eset[, -1]) # first column is ENTREZ_GENE_ID
  t_df <- t(df)
  means <- colMeans(t_df, na.rm = TRUE)
  sds <- apply(t_df, 2, sd, na.rm = TRUE)
  z_t <- sweep(sweep(t_df, 2, means, FUN = "-"), 2, sds, FUN = "/")
  z <- t(z_t)
  rownames(z) <- eset$ENTREZ_GENE_ID
  colnames(z) <- colnames(df)
  z_df <- as.data.frame(z)
  z_df$ENTREZ_GENE_ID <- eset$ENTREZ_GENE_ID
  list(means = means, sd = sds, z = z_df)
}

results <- lapply(esets_trainval, calculate_z)
means <- lapply(results, `[[`, "means") # Extract element from sublist
sd <- lapply(results, `[[`, "sd")
z_trainval <- lapply(results, `[[`, "z")
names(means) <- names(sd) <- names(z_trainval) <- names(esets_trainval)

# To apply normalisation consistently on the test set we will use the saved 
# means and standard deviations in the means and sd lists that we created:

normalize_test_set <- function(test_set, i) {
  # data frame with only one numeric column
  if (ncol(test_set) == 2) {
    # Extract the numeric column (assuming it's the first column)
    numeric_col <- test_set[, 1, drop = FALSE]
    
    # Normalize the numeric column
    z_normalized <- (numeric_col - means[[i]]) / sd[[i]]
    
    # Create the result dataframe
    z_df <- as.data.frame(z_normalized)
    rownames(z_df) <- test_set$ENTREZ_GENE_ID
    
    # Adding back the ENTREZ_GENE_ID
    z_df$ENTREZ_GENE_ID <- test_set$ENTREZ_GENE_ID
    
    z_df
  } else {
    df <- as.data.frame(test_set[, -ncol(test_set)]) # last column is ENTREZ_GENE_ID
    t_df <- t(df)
    z_t <- sweep(sweep(t_df, 2, means[[i]], FUN = "-"), 2, sd[[i]], FUN = "/")
    z <- t(z_t)
    rownames(z) <- test_set$ENTREZ_GENE_ID
    colnames(z) <- colnames(df)
    z_df <- as.data.frame(z)
    z_df$ENTREZ_GENE_ID <- test_set$ENTREZ_GENE_ID
    z_df
  }
}

z_test <- Map(normalize_test_set, esets_test_new, seq_along(esets_test_new))
names(z_test) <- names(z_trainval)

##### Quality Control #####

# Joining in one expression matrix: original version
for (i in 1:length(esets_trainval)){
  esets_trainval[[i]]$ENTREZ_GENE_ID = as.character(esets_trainval[[i]]$ENTREZ_GENE_ID)
  esets_trainval[[i]] = esets_trainval[[i]][,c("ENTREZ_GENE_ID", 
                                               intersect(colnames(esets_trainval[[i]]),
                                                         filt_pdata[[i]]$GEO_accession))]
}
original_exprs = esets_trainval[[1]] %>% inner_join(esets_trainval[[2]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets_trainval[[3]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets_trainval[[4]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets_trainval[[5]], by = "ENTREZ_GENE_ID") %>%
  dplyr::select(ENTREZ_GENE_ID, everything())

rows = original_exprs$ENTREZ_GENE_ID
original_exprs = as.matrix(original_exprs %>% dplyr::select(-ENTREZ_GENE_ID))
rownames(original_exprs) = rows; rm(rows) # 13,749 x 149
# Making sure we do not have NAs in any row
original_exprs = original_exprs[rowSums(is.na(original_exprs)) != ncol(original_exprs), ]
original_exprs_nonas = na.omit(original_exprs) # 13,749 x 149
# Keeping all samples in our full_pdata object
original_exprs_nonas = original_exprs_nonas[, intersect(full_pdata$GEO_accession,
                                                      colnames(original_exprs_nonas))] # 13,749 x 149

# Joining in one expression matrix
# Test (normalised)
for (i in 1:length(z_test)){
  z_test[[i]]$ENTREZ_GENE_ID = as.character(z_test[[i]]$ENTREZ_GENE_ID)
  z_test[[i]] = z_test[[i]][,c("ENTREZ_GENE_ID", 
                               intersect(colnames(z_test[[i]]),
                                         filt_pdata[[i]]$GEO_accession))]
}

z_test_exprs = z_test[[1]] %>% inner_join(z_test[[2]], by = "ENTREZ_GENE_ID") %>%
  inner_join(z_test[[3]], by = "ENTREZ_GENE_ID") %>%
  inner_join(z_test[[4]], by = "ENTREZ_GENE_ID") %>%
  inner_join(z_test[[5]], by = "ENTREZ_GENE_ID") %>%
  dplyr::select(ENTREZ_GENE_ID, everything())

rownames(z_test_exprs) = z_test_exprs$ENTREZ_GENE_ID
z_test_exprs = as.matrix(z_test_exprs %>% dplyr::select(-ENTREZ_GENE_ID)) # 13749 x 16
# Making sure we do not have NAs in any row
z_test_exprs = z_test_exprs[rowSums(is.na(z_test_exprs)) != ncol(z_test_exprs), ]
z_test_exprs_nonas = na.omit(z_test_exprs) # 13749 x 16
z_test_exprs_nonas = z_test_exprs_nonas[, intersect(full_pdata$GEO_accession,
                                          colnames(z_test_exprs_nonas))] # 13749 x 16

# Joining in one expression matrix: z-score normalized version (train-val)
for (i in 1:length(z_trainval)){
  z_trainval[[i]]$ENTREZ_GENE_ID = as.character(z_trainval[[i]]$ENTREZ_GENE_ID)
  z_trainval[[i]] = z_trainval[[i]][,c("ENTREZ_GENE_ID", 
                                       intersect(colnames(z_trainval[[i]]),
                                                 filt_pdata[[i]]$GEO_accession))]
}
z_exprs = z_trainval[[1]] %>% inner_join(z_trainval[[2]], by = "ENTREZ_GENE_ID") %>%
  inner_join(z_trainval[[3]], by = "ENTREZ_GENE_ID") %>%
  inner_join(z_trainval[[4]], by = "ENTREZ_GENE_ID") %>%
  inner_join(z_trainval[[5]], by = "ENTREZ_GENE_ID") %>%
  dplyr::select(ENTREZ_GENE_ID, everything())

rownames(z_exprs) = z_exprs$ENTREZ_GENE_ID
z_exprs = as.matrix(z_exprs %>% dplyr::select(-ENTREZ_GENE_ID)) # 13749   149
# Making sure we do not have NAs in any row
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
z_exprs_nonas = na.omit(z_exprs) # 13749   149
z_exprs_nonas = z_exprs_nonas[, intersect(full_pdata$GEO_accession,
                                          colnames(z_exprs_nonas))] # 13749   149
# save
saveRDS(z_exprs_nonas, "TVT_split_scripts/z_exprs_nonas_microarrays_trainval.rds")

# Multidimensional scaling plot: original matrix #####
original_mds = plotMDS(original_exprs_nonas)
original_pca = data.frame(cbind(original_mds$x, original_mds$y, 
                                as.character(full_pdata$Study), 
                                full_pdata$GEO_accession, 
                                as.character(full_pdata$Tissue_type)))
colnames(original_pca) = c("X1", "X2", "Study", "GSM", "Type")
original_pca$Study = factor(original_pca$Study)
original_pca$Type = factor(original_pca$Type)
original_pca$X1 = as.numeric(original_pca$X1)
original_pca$X2 = as.numeric(original_pca$X2)

original_MDS = ggplot(original_pca, aes(X1, X2, color = Study, shape = Type)) +
  geom_point(size = 0.2) +
  scale_color_brewer(palette = "Dark2") +
  theme(plot.title = element_text(face = "bold", size = 5, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 margin = margin(t = 1, unit = "cm"),
                                 size = 2.5),
        axis.title = element_text(angle = 0, hjust = 0.5, face = "bold", 
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 3.5),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "right",
        legend.key.size = unit(1, units = "mm"),
        legend.key.height = unit(1.5, "mm"),
        legend.text = element_text(size = 2.5),
        legend.title = element_text(face = "bold", size = 3),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(1, units = "mm"))+
  labs(title = "Multidimensional Scaling Plot",
       # x = paste0("\nPC1 (", round(100*original_mds$var.explained[1],2), "% of variance)"),
       # y = paste0("PC2 (", round(100*original_mds$var.explained[2],2), "% of variance)\n")
       x = "MDS1", y = "MDS2")
original_MDS
ggsave(filename = "Original_MDS.tiff",
       path = "QC/GSE_microarrays", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# Multidimensional scaling plot: z-score normalised matrix
z_mds = plotMDS(z_exprs_nonas)
z_pca = data.frame(cbind(z_mds$x, z_mds$y, 
                         as.character(full_pdata$Study), full_pdata$GEO_accession, 
                         as.character(full_pdata$Tissue_type)))
colnames(z_pca) = c("X1", "X2", "Study", "GSM", "Type")
z_pca$Study = factor(z_pca$Study)
z_pca$Type = factor(z_pca$Type)
z_pca$X1 = as.numeric(z_pca$X1)
z_pca$X2 = as.numeric(z_pca$X2)

KBZ_MDS_plot = ggplot(z_pca, aes(X1, X2, color = Study, shape = Type)) +
  geom_point(size = 0.2) +
  scale_color_brewer(palette = "Dark2") +
  theme(plot.title = element_text(face = "bold", size = 5, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 margin = margin(t = 1, unit = "cm"),
                                 size = 2.5),
        axis.title = element_text(angle = 0, hjust = 0.5, face = "bold", 
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 3.5),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "right",
        legend.key.size = unit(1, units = "mm"),
        legend.key.height = unit(1.5, "mm"),
        legend.text = element_text(size = 2.5),
        legend.title = element_text(face = "bold", size = 3),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(1, units = "mm"))+
  labs(title = "Multidimensional Scaling Plot: normalized data",
       # x = paste0("\nPC1 (", round(100*original_mds$var.explained[1],2), "% of variance)"),
       # y = paste0("PC2 (", round(100*original_mds$var.explained[2],2), "% of variance)\n")
       x = "MDS1", y = "MDS2")
KBZ_MDS_plot
ggsave(filename = "z_MDS.tiff",
       path = "QC/GSE_microarrays", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# Defining the multiplot function

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# End of multiplot function

if (Sys.info()['sysname'] == "Windows") {
  type_compression = "windows"
} 
# if Linux
if (Sys.info()['sysname'] == "Linux") {
  type_compression = "windows"
}
# if macOS
if (Sys.info()['sysname'] == "Darwin") {
  type_compression = "cairo"
} 

tiff("QC/GSE_microarrays/MDS_multiplot.tiff", 
     width = 1920, height = 2160, res = 700, compression = "lzw")
multiplot(original_MDS, KBZ_MDS_plot, cols = 1)
m = ggplot(multiplot(original_MDS, KBZ_MDS_plot, cols = 1))
dev.off(); rm(m)

# Global expression boxplot: original matrix
trainval_pheno_ordered = full_pdata[full_pdata$GEO_accession %in% trainval_samples,]
trainval_pheno_ordered = trainval_pheno_ordered[match(trainval_samples,
                                                      trainval_pheno_ordered$GEO_accession), ]
original_eset = as.data.frame(original_exprs_nonas[, trainval_samples])
data_boxplot = reshape2::melt(original_eset[, trainval_samples])
data_boxplot = left_join(data_boxplot, full_pdata[, c(1,3)], by = c("variable" = "GEO_accession")) %>%
  arrange(Study, variable)
original_boxplot = ggplot(data_boxplot, aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.01, outlier.shape = 20, linewidth = 0.02,
               fill = c(rep("cyan", 18),
                        rep("chartreuse", 23),
                        rep("orange", 11),
                        rep("red", 20),
                        rep("grey", 77))) +
  scale_y_continuous("Expression",
                     limits = c(-3,
                                round(max(reshape2::melt(original_eset)$value, na.rm = TRUE)+1)), 
                     breaks = seq(-3,
                                  round(max(reshape2::melt(original_eset)$value, na.rm = TRUE)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, #margin = margin(t = 1, unit = "cm"),
                                   size = 3),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5), 
        #margin = margin(t = .05, unit = "cm") ),
        axis.title = element_text(angle = 0, hjust = 0.5, #margin = margin(t = 1, unit = "cm"),
                                  size = 10, face = "bold"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_line(linewidth = 0.1))+
  labs(title = "Boxplot of expression: original values",
       x = "\nSamples",
       y = "Expression\n")
rm(data_boxplot)

original_boxplot
ggsave(filename = "Original_boxplot.tiff",
       path = "QC/GSE_microarrays/", 
       width = 7680, height = 3240, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# Global expression boxplot: z-score normalised matrix
z_eset = as.data.frame(z_exprs_nonas[, trainval_samples])
data_boxplot = reshape2::melt(z_eset[, trainval_samples])
data_boxplot = left_join(data_boxplot, full_pdata[, c(1,3)], by = c("variable" = "GEO_accession")) %>%
  arrange(Study, variable)
KBZ_boxplot = ggplot(data_boxplot, aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.01, outlier.shape = 20, linewidth = 0.02,
               fill = c(rep("cyan", 18),
                        rep("chartreuse", 23),
                        rep("orange", 11),
                        rep("red", 20),
                        rep("grey", 77))) +
  scale_y_continuous("Expression",
                     limits = c(round(min(reshape2::melt(z_eset)$value, na.rm = TRUE)-1),
                                round(max(reshape2::melt(z_eset)$value, na.rm = TRUE)+1)), 
                     breaks = seq(round(min(reshape2::melt(z_eset)$value, na.rm = TRUE)-1),
                                  round(max(reshape2::melt(z_eset)$value, na.rm = TRUE)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, #margin = margin(t = 1, unit = "cm"),
                                   size = 3),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5), 
        #margin = margin(t = .05, unit = "cm") ),
        axis.title = element_text(angle = 0, hjust = 0.5, #margin = margin(t = 1, unit = "cm"),
                                  size = 10, face = "bold"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_line(linewidth = 0.1))+
  labs(title = "Boxplot of expression: z-score-normalised data",
       x = "\nSamples",
       y = "Expression\n")
rm(data_boxplot)

KBZ_boxplot
ggsave(filename = "z_boxplot.tiff",
       path = "QC/GSE_microarrays/", 
       width = 7680, height = 3240, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

tiff("QC/GSE_microarrays/Boxplot_multiplot.tiff", 
     width = 7680, height = 6480, res = 700, compression = "lzw")
multiplot(original_boxplot, KBZ_boxplot, cols = 1)
m = ggplot(multiplot(original_boxplot, KBZ_boxplot, cols = 1))
dev.off(); rm(m)

# Heatmaps
save_pheatmap_png <- function(x, filename, width=2600*2, height=1800*2, res = 700) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_tiff <- function(x, filename, width=2600*2, height=1800*2, res = 700) {
  tiff(filename, width = width, height = height, res = res, compression = "lzw")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Original
annotation_for_heatmap = trainval_pheno_ordered[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = trainval_pheno_ordered$GEO_accession

original_dists = as.matrix(dist(t(original_exprs_nonas), method = "manhattan"))

rownames(original_dists) = trainval_pheno_ordered$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(original_dists) <- NULL
diag(original_dists) <- NA

ann_colors <- list(
  Tissue_type = c(chronic_pancreatitis = "deeppink4", normal = "dodgerblue4", tumor = "gray"),
  Study = c(GSE143754 = "darkseagreen", GSE61166 = "darkorange",
            GSE71989 = "darkcyan", GSE101462 = "darkred", GSE77858 = "darkmagenta")
)

original_heatmap = pheatmap(t(original_dists), col = hmcol,
                            annotation_col = annotation_for_heatmap,
                            annotation_colors = ann_colors,
                            legend = TRUE,
                            show_rownames = F,
                            show_colnames = F,
                            treeheight_col = 0,
                            fontsize = 5,
                            legend_breaks = c(min(original_dists, na.rm = TRUE), 
                                              max(original_dists, na.rm = TRUE)), 
                            legend_labels = (c("small distance", "large distance")),
                            main = "Original heatmap")
save_pheatmap_tiff(original_heatmap, "QC/GSE_microarrays/original_heatmap.tiff")
dev.off()

# Z-score version
annotation_for_heatmap = trainval_pheno_ordered[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = trainval_pheno_ordered$GEO_accession

z_dists = as.matrix(dist(t(z_exprs_nonas), method = "manhattan"))

rownames(z_dists) = trainval_pheno_ordered$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

ann_colors <- list(
  Tissue_type = c(chronic_pancreatitis = "deeppink4", normal = "dodgerblue4", tumor = "gray"),
  Study = c(GSE143754 = "darkseagreen", GSE61166 = "darkorange",
            GSE71989 = "darkcyan", GSE101462 = "darkred", GSE77858 = "darkmagenta")
)

z_heatmap = pheatmap(t(z_dists), col = hmcol,
                     annotation_col = annotation_for_heatmap,
                     annotation_colors = ann_colors,
                     legend = TRUE,
                     show_rownames = F,
                     show_colnames = F,
                     treeheight_col = 0,
                     fontsize = 5,
                     legend_breaks = c(min(z_dists, na.rm = TRUE), 
                                       max(z_dists, na.rm = TRUE)), 
                     legend_labels = (c("small distance", "large distance")),
                     main = "Normalized data heatmap")
save_pheatmap_tiff(z_heatmap, "QC/GSE_microarrays/KBZ_heatmap.tiff")
dev.off()

#
# ##### Differential Gene Expression (DGEA) #####
# 

# Annotation with official gene symbols
# official_df, Aliases, ID_Map
official = org.Hs.egSYMBOL
mapped_genes_official = mappedkeys(official)
official_df = as.data.frame(official[mapped_genes_official])
official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
official_df$HGNC_Official = "Yes"
official_df = official_df[-which(duplicated(official_df$Gene.Symbol)==T),]
official_df = distinct(official_df)

alias = org.Hs.egALIAS2EG
mapped_genes_alias = mappedkeys(alias)
alias_df = as.data.frame(alias[mapped_genes_alias])
alias_df = alias_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = alias_symbol)
alias_df = alias_df[-which(alias_df$Gene.Symbol %in% official_df$Gene.Symbol),]
alias_df$HGNC_Official = "No"

ID_Map = rbind(official_df, alias_df) %>% distinct()
ID_Map$EntrezGene.ID = as.numeric(ID_Map$EntrezGene.ID)
ID_Map = ID_Map[order(ID_Map$EntrezGene.ID),] %>%
  dplyr::rename(probe=Gene.Symbol) %>%
  dplyr::select(probe, EntrezGene.ID, HGNC_Official)

# Aliases
aliases_for_join = alias_df %>% dplyr::rename(Alias = Gene.Symbol)
Aliases = official_df %>% inner_join(aliases_for_join,
                                     by = "EntrezGene.ID") %>%
  dplyr::select(Alias, Gene.Symbol, EntrezGene.ID) %>%
  dplyr::rename(probe = Alias, HGNC_Symbol = Gene.Symbol,
                Entrez = EntrezGene.ID) %>%
  distinct()

ID_Map$EntrezGene.ID = as.character(ID_Map$EntrezGene.ID)
ID_Map = ID_Map %>% dplyr::rename(Gene.Symbol = probe)
rm(alias, alias_df, aliases_for_join, official,
   mapped_genes_alias, mapped_genes_official)

# Create design and contrast matrix
design = model.matrix(~0 + trainval_pheno_ordered$Tissue_type + 
                        trainval_pheno_ordered$Study)
colnames(design) = c("chronic_pancreatitis", "normal", "tumor", 
                     "GSE143754", "GSE61166", "GSE71989",
                     "GSE77858") 
rownames(design) = colnames(z_exprs_nonas)
cont.matrix = makeContrasts(CPvsNormal = chronic_pancreatitis - normal,
                            CPvsTumor = chronic_pancreatitis - tumor,
                            TumorvsNormal = tumor - normal,
                            levels = design)

# Limma
z_fit = lmFit(z_exprs_nonas, design)
z_fit2 = contrasts.fit(z_fit, cont.matrix)
z_fit2 = eBayes(z_fit2, robust = TRUE)
z_results = decideTests(z_fit2)
results = as.data.frame(cbind(summary(z_results), rownames(summary(z_results))))
colnames(results)[ncol(results)] = "direction"
results = results %>% dplyr::select(direction, everything())

DGEA_topTables = createWorkbook()
addWorksheet(DGEA_topTables, "summary")
writeData(DGEA_topTables, "summary", results)

DE_maps = list()

for (i in 1:ncol(z_fit2)){
  z_DE = as.data.frame(topTable(z_fit2, adjust.method="BH", 
                                number = Inf, coef = colnames(z_fit2)[i]))
  z_DE$EntrezGene.ID = rownames(z_DE)
  
  # Annotation with official gene symbols
  z_DE_mapped = z_DE %>% left_join(ID_Map, by = "EntrezGene.ID", multiple = "all")
  z_DE_mapped$Filter = NA
  unmapped = which(is.na(z_DE_mapped$HGNC_Official))
  z_DE_mapped$HGNC_Official[unmapped] = "unmapped"
  for(j in 1:nrow(z_DE_mapped)){
    if(z_DE_mapped$HGNC_Official[j] == "Yes"){
      z_DE_mapped$Filter[j] = "Keep"
    } else if(length(unique(z_DE_mapped$HGNC_Official[z_DE_mapped$EntrezGene.ID ==
                                                      z_DE_mapped$EntrezGene.ID[j]])) > 1 &&
              z_DE_mapped$HGNC_Official[j] == "No"){
      z_DE_mapped$Filter[j] = "Discard"
    } else if(unique(z_DE_mapped$HGNC_Official[z_DE_mapped$EntrezGene.ID ==
                                               z_DE_mapped$EntrezGene.ID[j]]) == "No"){
      z_DE_mapped$Filter[j] = "Keep"
      z_DE_mapped$Gene.Symbol[j] = z_DE_mapped$EntrezGene.ID[j]
    } else if(z_DE_mapped$HGNC_Official[j] == "unmapped"){
      z_DE_mapped$Gene.Symbol[j] = z_DE_mapped$EntrezGene.ID[j]
      z_DE_mapped$Filter[j] = "Keep"
    }
  }
  
  z_DE_mapped = z_DE_mapped %>% 
    dplyr::filter(Filter == "Keep") %>%
    dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
    dplyr::select(-Filter) %>%
    distinct()
  z_DE_mapped = z_DE_mapped[order(z_DE_mapped$adj.P.Val),]
  rownames(z_DE_mapped) = z_DE_mapped$EntrezGene.ID
  addWorksheet(DGEA_topTables, colnames(z_fit2)[i])
  writeData(DGEA_topTables, colnames(z_fit2)[i], z_DE_mapped)
  DE_maps[[i]] = z_DE_mapped
}

saveWorkbook(DGEA_topTables, "DGEA/GSE_microarrays/DGEA_results.xlsx", overwrite = TRUE)
names(DE_maps) = colnames(z_fit2)

##### Volcano plots #####
# create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
keyvals.colours = list()
for (i in 1:length(DE_maps)) {
  tab = DE_maps[[i]]
  keyvals.colour <- ifelse(
    tab$logFC < -1 & tab$adj.P.Val < 0.05, 'royalblue',
    ifelse(tab$logFC > 1 & tab$adj.P.Val < 0.05, 'red4',
           ifelse(abs(tab$logFC) < 1 & tab$adj.P.Val < 0.05, 'pink', 
                  'grey')))
  # keyvals.colour[is.na(keyvals.colour)] <- 'black'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
  names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
  names(keyvals.colour)[keyvals.colour == 'pink'] <- '|DE| < 1'
  names(keyvals.colour)[keyvals.colour == 'grey'] <- 'p.adj > 0.05'
  keyvals.colours[[i]] = keyvals.colour
}
names(keyvals.colours) = names(DE_maps)

# CP vs. Normal samples
volcano_CPvsNormal = EnhancedVolcano(DE_maps[["CPvsNormal"]],
                                     lab = DE_maps[["CPvsNormal"]][, "Gene.Symbol"],
                                     caption = NULL,
                                     x = 'logFC',
                                     y = 'adj.P.Val',
                                     title = "Chronic Pancreatitis vs. Normal",
                                     pCutoff = 0.05,
                                     cutoffLineType = "dashed",
                                     cutoffLineWidth = 0.3,
                                     cutoffLineCol = "black",
                                     FCcutoff = 1,
                                     colCustom = keyvals.colours[["CPvsNormal"]],
                                     colAlpha = 0.7,
                                     xlim = c(-5, 5),
                                     ylab = bquote(bold(-log[10]("BH adj. p-value"))),
                                     xlab = "\nDifferential expression",
                                     pointSize = 1.5,
                                     axisLabSize = 7,
                                     subtitle = NULL,
                                     labSize = 2,
                                     selectLab = DE_maps[["CPvsNormal"]][1:20, "Gene.Symbol"], # top 20 genes
                                     legendLabSize = 6,
                                     legendIconSize = 4,
                                     labFace = "bold",
                                     boxedLabels = TRUE,
                                     drawConnectors = TRUE,
                                     typeConnectors = "closed",
                                     arrowheads = FALSE,
                                     widthConnectors = 0.3,
                                     max.overlaps = Inf,
                                     legendLabels = c("NS", "|DE| > 1 s.d.", 
                                                      "p.adj < 0.05", 
                                                      "p.adj < 0.05 & |DE| > 1 s.d."))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.4),
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold", size = 10),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.ticks.length = unit(1, units = "mm"),
        legend.position = "bottom",
        #legend.text = element_text(size = 8),
        #legend.title = element_blank(),
        #legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        #legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(0.3, units = "mm")#,
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"))
  )
volcano_CPvsNormal
ggsave(filename = "CPvsNormal_Volcano.tiff",
       path = "DGEA/GSE_microarrays/", 
       width = 100, height = 142, device = 'tiff', units = "mm",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# CP vs. Tumor samples
volcano_CPvsTumor = EnhancedVolcano(DE_maps[["CPvsTumor"]],
                                    lab = DE_maps[["CPvsTumor"]][, "Gene.Symbol"],
                                    caption = NULL,
                                    x = 'logFC',
                                    y = 'adj.P.Val',
                                    title = "Chronic Pancreatitis vs. PDAC",
                                    pCutoff = 0.05,
                                    cutoffLineType = "dashed",
                                    cutoffLineWidth = 0.3,
                                    cutoffLineCol = "black",
                                    FCcutoff = 1,
                                    colCustom = keyvals.colours[["CPvsTumor"]],
                                    colAlpha = 0.7,
                                    xlim = c(-5, 5),
                                    ylab = bquote(bold(-log[10]("BH adj. p-value"))),
                                    xlab = "\nDifferential expression",
                                    pointSize = 1.5,
                                    axisLabSize = 7,
                                    subtitle = NULL,
                                    labSize = 2,
                                    selectLab = DE_maps[["CPvsTumor"]][1:20, "Gene.Symbol"], # top 20 genes
                                    legendLabSize = 6,
                                    legendIconSize = 4,
                                    labFace = "bold",
                                    boxedLabels = TRUE,
                                    drawConnectors = TRUE,
                                    typeConnectors = "closed",
                                    arrowheads = FALSE,
                                    widthConnectors = 0.3,
                                    max.overlaps = Inf,
                                    legendLabels = c("NS", "|DE| > 1 s.d.", 
                                                     "p.adj < 0.05", 
                                                     "p.adj < 0.05 & |DE| > 1 s.d."))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.4),
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold", size = 10),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.ticks.length = unit(1, units = "mm"),
        legend.position = "bottom",
        #legend.text = element_text(size = 8),
        #legend.title = element_blank(),
        #legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        #legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(0.3, units = "mm")#,
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"))
  )
volcano_CPvsTumor
ggsave(filename = "CPvsTumor_Volcano.tiff",
       path = "DGEA/GSE_microarrays/", 
       width = 100, height = 142, device = 'tiff', units = "mm",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# Tumor vs. Normal samples
volcano_TumorvsNormal = EnhancedVolcano(DE_maps[["TumorvsNormal"]],
                                        lab = DE_maps[["TumorvsNormal"]][, "Gene.Symbol"],
                                        caption = NULL,
                                        x = 'logFC',
                                        y = 'adj.P.Val',
                                        title = "PDAC vs. Normal",
                                        pCutoff = 0.05,
                                        cutoffLineType = "dashed",
                                        cutoffLineWidth = 0.3,
                                        cutoffLineCol = "black",
                                        FCcutoff = 1,
                                        colCustom = keyvals.colours[["TumorvsNormal"]],
                                        colAlpha = 0.7,
                                        xlim = c(-5, 5),
                                        ylab = bquote(bold(-log[10]("BH adj. p-value"))),
                                        xlab = "\nDifferential expression",
                                        pointSize = 1.5,
                                        axisLabSize = 7,
                                        subtitle = NULL,
                                        labSize = 2,
                                        selectLab = DE_maps[["TumorvsNormal"]][1:20, "Gene.Symbol"], # top 20 genes
                                        legendLabSize = 6,
                                        legendIconSize = 4,
                                        labFace = "bold",
                                        boxedLabels = TRUE,
                                        drawConnectors = TRUE,
                                        typeConnectors = "closed",
                                        arrowheads = FALSE,
                                        widthConnectors = 0.3,
                                        max.overlaps = Inf,
                                        legendLabels = c("NS", "|DE| > 1 s.d.", 
                                                         "p.adj < 0.05", 
                                                         "p.adj < 0.05 & |DE| > 1 s.d."))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.4),
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold", size = 10),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.ticks.length = unit(1, units = "mm"),
        legend.position = "bottom",
        #legend.text = element_text(size = 8),
        #legend.title = element_blank(),
        #legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        #legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(0.3, units = "mm")#,
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"))
  )
volcano_TumorvsNormal
ggsave(filename = "TumorvsNormal_Volcano.tiff",
       path = "DGEA/GSE_microarrays/", 
       width = 100, height = 142, device = 'tiff', units = "mm",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# Write out trainval and test subsets #####
dir.create("TVT_split_scripts/microarray_objects")
saveRDS(z_eset, "TVT_split_scripts/microarray_objects/microarray_z_trainval.rds")
z_test_eset = as.data.frame(z_test_exprs_nonas[, test_samples])
test_pheno_ordered = full_pdata[full_pdata$GEO_accession %in% test_samples,]
test_pheno_ordered = test_pheno_ordered[match(test_samples,
                                                      test_pheno_ordered$GEO_accession), ]
trainval_pheno_ordered$Sample_type = "Microarray"
test_pheno_ordered$Sample_type = "Microarray"
saveRDS(z_test_eset, "TVT_split_scripts/microarray_objects/microarray_z_test.rds")
write.xlsx(trainval_pheno_ordered, "TVT_split_scripts/microarray_objects/microarray_z_trainval_pheno.xlsx")
write.xlsx(test_pheno_ordered, "TVT_split_scripts/microarray_objects/microarray_z_test_pheno.xlsx")

# ---------------------------------------- #
# ---------------------------------------- #
# ---------------------------------------- #