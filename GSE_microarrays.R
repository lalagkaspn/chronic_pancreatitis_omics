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
for (i in 1:length(filt_pdata$GSE143754$Tissue_type)) {
  if (filt_pdata$GSE143754$Tissue_type[i] == "Tumor") {
    filt_pdata$GSE143754$Tissue_type[i] = "tumor"
  } 
  if (filt_pdata$GSE143754$Tissue_type[i] == "Adjacent Normal") {
    filt_pdata$GSE143754$Tissue_type[i] = "non_tumor"
  } 
  if (filt_pdata$GSE143754$Tissue_type[i] == "Chronic Pancreatitis") {
    filt_pdata$GSE143754$Tissue_type[i] = "chronic_pancreatitis"
  } 
}; rm(i)

## NOTE: We keep all samples (including PDAC) because we want to adjust for them in the DGEA between normal and CPs

# Clear patient ID
filt_pdata[["GSE143754"]]$Patient_ID = gsub("Benign Tissue, Biological Replicate ", "", filt_pdata[["GSE143754"]]$Patient_ID)
filt_pdata[["GSE143754"]]$Patient_ID = gsub("Malignant Tissue, Biological Replicate ", "", filt_pdata[["GSE143754"]]$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE143754$Tissue_type = factor(x = filt_pdata$GSE143754$Tissue_type,
                                          levels = c("chronic_pancreatitis", "non_tumor", "tumor"),
                                          labels = c("chronic_pancreatitis", "non_tumor", "tumor"))
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
for (i in 1:length(filt_pdata$GSE61166$Tissue_type)) {
  if (filt_pdata$GSE61166$Tissue_type[i] == "pancreatitis") {
    filt_pdata$GSE61166$Tissue_type[i] = "chronic_pancreatitis"
  } 
  if (filt_pdata$GSE61166$Tissue_type[i] == "pancreatic tumor") {
    filt_pdata$GSE61166$Tissue_type[i] = "tumor"
  } 
}; rm(i)

# Clear patient ID
filt_pdata$GSE61166$Patient_ID = gsub("Pancreatitis_tissues_of_Patient", "", filt_pdata$GSE61166$Patient_ID)
filt_pdata$GSE61166$Patient_ID = gsub("Pancreatic_tumors_of_Patient", "", filt_pdata$GSE61166$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE61166$Tissue_type = factor(x = filt_pdata$GSE61166$Tissue_type,
                                         levels = c("chronic_pancreatitis", "non_tumor", "tumor"),
                                         labels = c("chronic_pancreatitis", "non_tumor", "tumor"))
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
for (i in 1:length(filt_pdata$GSE71989$Tissue_type)) {
  if (filt_pdata$GSE71989$Tissue_type[i] == "CP") {
    filt_pdata$GSE71989$Tissue_type[i] = "chronic_pancreatitis"
  } 
  if (filt_pdata$GSE71989$Tissue_type[i] == "normal pancreatic tissue") {
    filt_pdata$GSE71989$Tissue_type[i] = "non_tumor"
  } 
  if (filt_pdata$GSE71989$Tissue_type[i] == "PDAC") {
    filt_pdata$GSE71989$Tissue_type[i] = "tumor"
  } 
}; rm(i)

# Clear patient ID
filt_pdata$GSE71989$Patient_ID = gsub(", human normal pancreatic tissue", "", filt_pdata$GSE71989$Patient_ID)
filt_pdata$GSE71989$Patient_ID = gsub(", human Chronic Pancreatitis tissue", "", filt_pdata$GSE71989$Patient_ID)
filt_pdata$GSE71989$Patient_ID = gsub(", human PDAC tissue", "", filt_pdata$GSE71989$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE71989$Tissue_type = factor(x = filt_pdata$GSE71989$Tissue_type,
                                         levels = c("chronic_pancreatitis", "non_tumor", "tumor"),
                                         labels = c("chronic_pancreatitis", "non_tumor", "tumor"))
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
for (i in 1:length(filt_pdata$GSE101462$Tissue_type)) {
  if (filt_pdata$GSE101462$Tissue_type[i] == "pancreatitis") {
    filt_pdata$GSE101462$Tissue_type[i] = "chronic_pancreatitis"
  } 
  if (filt_pdata$GSE101462$Tissue_type[i] == "normal") {
    filt_pdata$GSE101462$Tissue_type[i] = "non_tumor"
  } 
  if (filt_pdata$GSE101462$Tissue_type[i] == "PDAC") {
    filt_pdata$GSE101462$Tissue_type[i] = "tumor"
  } 
}; rm(i)

# Change Tissue_storage
for (i in 1:length(filt_pdata$GSE101462$Tissue_storage)) {
  if (filt_pdata$GSE101462$Tissue_storage[i] == "fresh frozen (FF)") {
    filt_pdata$GSE101462$Tissue_storage[i] = "fresh_frozen"
  } 
  if (filt_pdata$GSE101462$Tissue_storage[i] == "formalin-fixed paraffin embedded (FFPE)") {
    filt_pdata$GSE101462$Tissue_storage[i] = "formalin_fixed_paraffin_embedded"
  } 
}; rm(i)

# Transform to factors with consistent universal levels
filt_pdata$GSE101462$Tissue_type = factor(x = filt_pdata$GSE101462$Tissue_type,
                                          levels = c("chronic_pancreatitis", "non_tumor", "tumor"),
                                          labels = c("chronic_pancreatitis", "non_tumor", "tumor"))
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
for (i in 1:length(filt_pdata$GSE77858$Tissue_type)) {
  if (filt_pdata$GSE77858$Tissue_type[i] == "Normal") {
    filt_pdata$GSE77858$Tissue_type[i] = "non_tumor"
  } 
  if (filt_pdata$GSE77858$Tissue_type[i] == "Panreatitis") {
    filt_pdata$GSE77858$Tissue_type[i] = "chronic_pancreatitis"
  } 
  if (filt_pdata$GSE77858$Tissue_type[i] == "Pancreatitis") {
    filt_pdata$GSE77858$Tissue_type[i] = "chronic_pancreatitis"
  } 
  if (filt_pdata$GSE77858$Tissue_type[i] == "Tumor") {
    filt_pdata$GSE77858$Tissue_type[i] = "tumor"
  } 
}; rm(i)

# Clear patient ID
filt_pdata[["GSE77858"]]$Patient_ID = gsub("PancTuRef2 vs. ", "", filt_pdata[["GSE77858"]]$Patient_ID)
filt_pdata[["GSE77858"]]$Patient_ID = gsub("PancTuRef vs. ", "", filt_pdata[["GSE77858"]]$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE77858$Tissue_type = factor(x = filt_pdata$GSE77858$Tissue_type,
                                         levels = c("chronic_pancreatitis", "non_tumor", "tumor"),
                                         labels = c("chronic_pancreatitis", "non_tumor", "tumor"))

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

## ---------------------------
## Expression data
## ---------------------------

esets = list()
for (i in 1:length(GEOsets)) {
  esets[[i]] = exprs(GEOsets[[i]])
  esets[[i]] = as.data.frame(esets[[i]])
}; rm(i)
names(esets) = names(GEOsets)

## Calculate NAs values
na_esets = c()
for(i in 1:length(esets)) {
  na_esets[i] = sum(is.na(esets[[i]]))
}; rm(i)
names(na_esets) = names(GEOsets)
na_esets 
# Missing values in 2 studies
# GSE101462: 22
# GSE77858: 5907

## We will impute missing values using KNN imputation. 
## However, GSE77858 has a lot of missing values. 
## We remove rows with more than 25% missing values because they will negatively affect imputation
GEOsets[["GSE101462"]] = GEOsets[["GSE101462"]][rowSums(is.na(GEOsets[["GSE101462"]]@assayData[["exprs"]]))/
                                                  length(colnames(GEOsets[["GSE101462"]]@assayData[["exprs"]])) < 0.25, ]
GEOsets[["GSE77858"]] = GEOsets[["GSE77858"]][rowSums(is.na(GEOsets[["GSE77858"]]@assayData[["exprs"]]))/
                                                length(colnames(GEOsets[["GSE77858"]]@assayData[["exprs"]])) < 0.25, ]

## -- KNN imputation -- ##
# GSE101462
RNGversion("4.0.2")
eset101462 = GEOsets[["GSE101462"]]@assayData[["exprs"]]
eset101462 = impute.knn(eset101462, k = 10, maxp = nrow(eset101462),
                        rng.seed = 123)
eset101462 = eset101462[["data"]]

# GSE102238
RNGversion("4.0.2")
eset77858 = GEOsets[["GSE77858"]]@assayData[["exprs"]]
eset77858 = impute.knn(eset77858, k = 10, maxp = nrow(eset77858),
                       rng.seed = 123)
eset77858 = eset77858[["data"]]

# Update esets object
esets[["GSE77858"]] = as.data.frame(eset77858); rm(eset77858)
esets[["GSE101462"]] = as.data.frame(eset101462); rm(eset101462)

##### Annotation esets with Entrez ID's #####
# GSE143754
# Platform [HTA-2_0] Affymetrix Human Transcriptome Array 2.0 [transcript (gene) version] -->
# gene data for each probe were downloaded from: https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_legacy/unconverted_chips/HTA_2_0.chip
hta2_0_probe_to_gene = data.table::fread("platforms_probes_to_genes/HTA_2_0.chip.txt")
fdata143754 = fData(GEOsets$GSE143754) %>%
  dplyr::select(probeset_id) %>%
  filter(probeset_id %in% hta2_0_probe_to_gene$`Probe Set ID`) %>%
  left_join(hta2_0_probe_to_gene, by = c("probeset_id" = "Probe Set ID")) %>%
  na.omit() %>%
  dplyr::select(ID = probeset_id, ENTREZ_GENE_ID = "Entrez Gene") %>%
  distinct()

## For probes that match to the same EntrezID, calculate variance across all samples and keep the probe with max variance
esets[["GSE143754"]]$variance = rowVars(as.matrix(esets[["GSE143754"]]))
esets[["GSE143754"]] = esets[["GSE143754"]] %>%
  mutate(ID = rownames(esets[["GSE143754"]])) %>%
  inner_join(fdata143754) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)
rm(fdata143754, hta2_0_probe_to_gene)

# GSE61166 
# need to match gene sequence to gene ID
# No information for corresponding genes provided in GEOset
# Three separate .fasta files were prepared based on the probe sequences available at the GPL file.
# BLASTN alignment (against Nucleotide sequences, "nt") was performed to annotate probes with RefSeq gene IDs: https://blast.ncbi.nlm.nih.gov/Blast.cgi

# Prepare files for BLASTN
gse61166_sequences = fData(GEOsets$GSE61166) %>% dplyr::select(SEQUENCE) %>% distinct()
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
gse61166_sequences = fData(GEOsets$GSE61166) %>% dplyr::select(SEQUENCE) %>% distinct()
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
for (i in 1:nrow(mapped_blastn)){
  if (length(unique(mapped_blastn$ENTREZ_GENE_ID[mapped_blastn$probe == 
                                                 mapped_blastn$probe[i]])) == 1) {
    mapped_blastn$unique_Entrez_perc[i] = 100
  } else {
    mapped_blastn$unique_Entrez_perc[i] = 100*length(which(mapped_blastn$ENTREZ_GENE_ID[mapped_blastn$probe == 
                                                                                          mapped_blastn$probe[i]] == mapped_blastn$ENTREZ_GENE_ID[i]))/
      nrow(mapped_blastn[mapped_blastn$probe == mapped_blastn$probe[i],])
  }
  cat(i, "\n")
} ; rm(i)

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
# 11,950
# The above numbers mean each probe is mapped to a unique Entrez ID, but multiple probes
# may map to the same ID. For each of these probes we calculate the variance and keep the one 
# with the max variance, as we did previously.

# GSE61166
esets[["GSE61166"]]$variance = rowVars(as.matrix(esets[["GSE61166"]]))
esets[["GSE61166"]] = esets[["GSE61166"]] %>%
  mutate(probe = rownames(esets[["GSE61166"]])) %>%
  inner_join(mapped_blastn_filt) %>%
  dplyr::select(-probe) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)
rm(blastn, mapped_blastn, mapped_blastn_filt)
esets[["GSE61166"]]$ENTREZ_GENE_ID = as.integer(esets[["GSE61166"]]$ENTREZ_GENE_ID)

# GSE71989
fdata71989 = fData(GEOsets$GSE71989) %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE71989"]]$variance = rowVars(as.matrix(esets[["GSE71989"]]))
esets[["GSE71989"]] = esets[["GSE71989"]] %>%
  mutate(ID = rownames(esets[["GSE71989"]])) %>%
  inner_join(fdata71989) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)
rm(fdata71989)
esets[["GSE71989"]]$ENTREZ_GENE_ID = as.integer(esets[["GSE71989"]]$ENTREZ_GENE_ID)

# GSE101462
fdata101462 = fData(GEOsets$GSE101462) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE101462"]]$variance = rowVars(as.matrix(esets[["GSE101462"]]))
esets[["GSE101462"]] = esets[["GSE101462"]] %>%
  mutate(ID = rownames(esets[["GSE101462"]])) %>%
  inner_join(fdata101462) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)
rm(fdata101462)

# GSE77858
fdata77858 = fData(GEOsets$GSE77858) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = GENE) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE77858"]]$variance = rowVars(as.matrix(esets[["GSE77858"]]))
esets[["GSE77858"]] = esets[["GSE77858"]] %>%
  mutate(ID = rownames(esets[["GSE77858"]])) %>%
  inner_join(fdata77858) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)
rm(fdata77858)

##### Calculate NAs values #####
na_esets = c()
for(i in 1:length(esets)) {
  na_esets[i] = sum(is.na(esets[[i]]))
}; rm(i)
names(na_esets) = names(GEOsets)
na_esets # No NA values exist in the expression files
rm(na_esets)

## Save intermediate files ##
saveRDS(esets, "intermediate_files/GSE_microarrays/esets_raw.RDS")
saveRDS(filt_pdata, "intermediate_files/GSE_microarrays/filt_pdata.RDS")
saveRDS(pdata, "intermediate_files/GSE_microarrays/pdata_raw.RDS")
openxlsx::write.xlsx(full_pdata, "DGEA/Pheno.xlsx", overwrite = TRUE)

##### z-score-transformation #####
# KBZ transformation method ( https:://www.biostars.org/p/283083/ )
z = list()
for(i in 1:length(esets)){
  df = as.data.frame(esets[[i]]) %>%
    dplyr::select(-ENTREZ_GENE_ID)
  t = as.data.frame(t(df))
  z_t = sapply(t, function(t) (t-mean(t, na.rm = T))/sd(t, na.rm = T))
  z[[i]] = as.matrix(t(z_t))
  rownames(z[[i]]) = esets[[i]]$ENTREZ_GENE_ID
  colnames(z[[i]]) = colnames(df)
  z[[i]] = as.data.frame(z[[i]])
  z[[i]]$EntrezGene.ID = esets[[i]]$ENTREZ_GENE_ID
  rm(t, z_t, df)
}; rm(i)

##### Quality Control #####

# Joining in one expression matrix: original version
for (i in 1:length(esets)){
  esets[[i]]$ENTREZ_GENE_ID = as.character(esets[[i]]$ENTREZ_GENE_ID)
  esets[[i]] = esets[[i]][,c("ENTREZ_GENE_ID", 
                             intersect(colnames(esets[[i]]),
                                       filt_pdata[[i]]$GEO_accession))]
}
original_exprs = esets[[1]] %>% inner_join(esets[[2]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets[[3]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets[[4]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets[[5]], by = "ENTREZ_GENE_ID") %>%
  dplyr::select(ENTREZ_GENE_ID, everything())

rows = original_exprs$ENTREZ_GENE_ID
original_exprs = as.matrix(original_exprs %>% dplyr::select(-ENTREZ_GENE_ID))
rownames(original_exprs) = rows; rm(rows) # 13,744 x 165
# Making sure we do not have NAs in any row
original_exprs = original_exprs[rowSums(is.na(original_exprs)) != ncol(original_exprs), ]
original_exprs_nonas = na.omit(original_exprs) # 13,744 x 165
# Keeping all samples in our full_pdata_filt object
original_exprs_nonas = original_exprs_nonas[, full_pdata$GEO_accession] # 13,744 x 165

# Joining in one expression matrix: z-score normalized version
for (i in 1:length(z)){
  z[[i]]$EntrezGene.ID = as.character(z[[i]]$EntrezGene.ID)
  z[[i]] = z[[i]][,c("EntrezGene.ID", 
                     intersect(colnames(z[[i]]),
                               filt_pdata[[i]]$GEO_accession))]
}
z_exprs = z[[1]] %>% inner_join(z[[2]], by = "EntrezGene.ID") %>%
  inner_join(z[[3]], by = "EntrezGene.ID") %>%
  inner_join(z[[4]], by = "EntrezGene.ID") %>%
  inner_join(z[[5]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

rownames(z_exprs) = z_exprs$EntrezGene.ID
z_exprs = as.matrix(z_exprs %>% dplyr::select(-EntrezGene.ID)) # 13,744 x 165
# Making sure we do not have NAs in any row
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
z_exprs_nonas = na.omit(z_exprs) # 13,744 x 165
z_exprs_nonas = z_exprs_nonas[, full_pdata$GEO_accession] # 13,744 x 165

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
  geom_point(size = 3, alpha = 1) +
  scale_color_brewer(palette = "Dark2") +
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                 size = 15),
        axis.title = element_text(angle = 0, hjust = 0.5, margin = margin(t = 3, unit = "cm"),
                                  size = 20),
        axis.line = element_line(),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.key.size = unit(1, "cm"))+
  labs(title = "Multidimensional Scaling Plot",
       x = paste0("\nPC1 (", round(100*original_mds$var.explained[1],2), "% of variance)"),
       y = paste0("PC2 (", round(100*original_mds$var.explained[2],2), "% of variance)\n"))
tiff("QC/Original_MDS.tif", width = 1920, height = 1080, res = 100)
original_MDS
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
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2") +
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                 size = 15),
        axis.title = element_text(angle = 0, hjust = 0.5, margin = margin(t = 3, unit = "cm"),
                                  size = 20),
        axis.line = element_line(),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.key.size = unit(1, "cm"))+
  labs(title = "Multidimensional Scaling Plot: z-score-normalised data",
       x = paste0("\nPC1 (", round(100*z_mds$var.explained[1],2), "% of variance)"),
       y = paste0("PC2 (", round(100*z_mds$var.explained[2],2), "% of variance)\n"))
tiff("QC/KBZ_MDS.tif", width = 1920, height = 1080, res = 100)
KBZ_MDS_plot
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

tiff("QC/MDS_multiplot.tif", 
     width = 2160, height = 3840, res = 150)
multiplot(original_MDS, KBZ_MDS_plot, cols = 1)
m = ggplot(multiplot(original_MDS, KBZ_MDS_plot, cols = 1))
dev.off(); rm(m)

# Global expression boxplot: original matrix
original_eset = as.data.frame(original_exprs_nonas)
original_boxplot = ggplot(melt(original_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20, outlier.alpha = 0.1,
               fill = c(rep("cyan", 26),
                        rep("chartreuse", 12),
                        rep("orange", 22),
                        rep("red", 20),
                        rep("grey", 85))) +
  scale_y_continuous("Expression", limits = c(0,round(max(melt(original_eset)$value)+1)), 
                     breaks = seq(0,round(max(melt(original_eset)$value)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                   size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, 
                                   margin = margin(t = .05, unit = "cm") ),
        axis.title = element_text(angle = 0, hjust = 0.5, margin = margin(t = 1, unit = "cm"),
                                  size = 25, face = "bold"),
        axis.line = element_line())+
  labs(title = "Boxplot of expression",
       x = "\nSamples",
       y = "Expression\n")
tiff("QC/Original_boxplot.tif", width = 1920, height = 1080, res = 100)
original_boxplot
dev.off()

# Global expression boxplot: z-score normalised matrix
z_eset = as.data.frame(z_exprs_nonas)
KBZ_boxplot = ggplot(melt(z_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20,  outlier.alpha = 0.1,
               fill = c(rep("cyan", 26),
                        rep("chartreuse", 12),
                        rep("orange", 22),
                        rep("red", 20),
                        rep("grey", 85))) +
  scale_y_continuous("Expression", limits = c(0,round(max(melt(z_eset)$value)+1)), 
                     breaks = seq(0,round(max(melt(z_eset)$value)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                   size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, 
                                   margin = margin(t = .05, unit = "cm") ),
        axis.title = element_text(angle = 0, hjust = 0.5, margin = margin(t = 1, unit = "cm"),
                                  size = 25, face = "bold"),
        axis.line = element_line())+
  labs(title = "Boxplot of expression: z-score-normalised data",
       x = "\nSamples",
       y = "Expression\n")
tiff("QC/KBZ_boxplot.tif", width = 1920, height = 1080, res = 100)
KBZ_boxplot
dev.off()

tiff("QC/Boxplot_multiplot.tif", 
     width = 3840, height = 3840, res = 150)
multiplot(original_boxplot, KBZ_boxplot, cols = 1)
m = ggplot(multiplot(original_boxplot, KBZ_boxplot, cols = 1))
dev.off(); rm(m)

# Heatmaps
save_pheatmap_png <- function(x, filename, width=2600, height=1800, res = 130) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_tiff <- function(x, filename, width=2600, height=1800, res = 130) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Original
annotation_for_heatmap = full_pdata[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata$GEO_accession

original_dists = as.matrix(dist(t(original_exprs_nonas), method = "manhattan"))

rownames(original_dists) = full_pdata$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(original_dists) <- NULL
diag(original_dists) <- NA

ann_colors <- list(
  Tissue_type = c(chronic_pancreatitis = "deeppink4", non_tumor = "dodgerblue4", tumor = "gray"),
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
                            legend_breaks = c(min(original_dists, na.rm = TRUE), 
                                              max(original_dists, na.rm = TRUE)), 
                            legend_labels = (c("small distance", "large distance")),
                            main = "Original heatmap")
save_pheatmap_tiff(original_heatmap, "QC/original_heatmap.tif")

# Z-score version
annotation_for_heatmap = full_pdata[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata$GEO_accession

z_dists = as.matrix(dist(t(z_exprs_nonas), method = "manhattan"))

rownames(z_dists) = full_pdata$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

ann_colors <- list(
  Tissue_type = c(chronic_pancreatitis = "deeppink4", non_tumor = "dodgerblue4", tumor = "gray"),
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
                     legend_breaks = c(min(z_dists, na.rm = TRUE), 
                                       max(z_dists, na.rm = TRUE)), 
                     legend_labels = (c("small distance", "large distance")),
                     main = "Z-score normalisation heatmap")
save_pheatmap_tiff(z_heatmap, "QC/KBZ_heatmap.tif")
# 
# ##### Differential Gene Expression (DGEA) #####
# 
# # Annotation with official gene symbols
# # official_df, Aliases, ID_Map
# official = org.Hs.egSYMBOL
# mapped_genes_official = mappedkeys(official)
# official_df = as.data.frame(official[mapped_genes_official])
# official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
# official_df$HGNC_Official = "Yes"
# official_df = official_df[-which(duplicated(official_df$Gene.Symbol)==T),]
# official_df = distinct(official_df)
# 
# alias = org.Hs.egALIAS2EG
# mapped_genes_alias = mappedkeys(alias)
# alias_df = as.data.frame(alias[mapped_genes_alias])
# alias_df = alias_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = alias_symbol)
# alias_df = alias_df[-which(alias_df$Gene.Symbol %in% official_df$Gene.Symbol),]
# alias_df$HGNC_Official = "No"
# 
# ID_Map = rbind(official_df, alias_df) %>% distinct()
# ID_Map$EntrezGene.ID = as.numeric(ID_Map$EntrezGene.ID)
# ID_Map = ID_Map[order(ID_Map$EntrezGene.ID),] %>%
#   dplyr::rename(probe=Gene.Symbol) %>%
#   dplyr::select(probe, EntrezGene.ID, HGNC_Official)
# 
# # Aliases
# aliases_for_join = alias_df %>% dplyr::rename(Alias = Gene.Symbol)
# Aliases = official_df %>% inner_join(aliases_for_join,
#                                      by = "EntrezGene.ID") %>%
#   dplyr::select(Alias, Gene.Symbol, EntrezGene.ID) %>%
#   dplyr::rename(probe = Alias, HGNC_Symbol = Gene.Symbol,
#                 Entrez = EntrezGene.ID) %>%
#   distinct()
# 
# ID_Map$EntrezGene.ID = as.character(ID_Map$EntrezGene.ID)
# ID_Map = ID_Map %>% dplyr::rename(Gene.Symbol = probe)
# rm(alias, alias_df, aliases_for_join, official,
#    mapped_genes_alias, mapped_genes_official)
# 
# ##### Union #####
# # In this section we perform DGEA on the union of the gene expression matrices,
# # not the intersection, in order to keep all genes from all
# # platforms. NA's will be introduced in this manner, but limma ignores them
# # during model fitting.
# 
# # Generating our union z-score matrix
# union_z_exprs = z[[1]] %>% full_join(z[[2]], by = "EntrezGene.ID") %>%
#   full_join(z[[3]], by = "EntrezGene.ID") %>%
#   full_join(z[[4]], by = "EntrezGene.ID") %>%
#   full_join(z[[5]], by = "EntrezGene.ID") %>%
#   full_join(z[[6]], by = "EntrezGene.ID") %>%
#   dplyr::select(EntrezGene.ID, everything())
# 
# rownames(union_z_exprs) = union_z_exprs$EntrezGene.ID
# union_z_exprs = as.matrix(union_z_exprs %>% dplyr::select(-EntrezGene.ID)) # 37969 x 73
# union_z_exprs_final = union_z_exprs[, full_pdata$GEO_accession] # 37969 x 73
# 
# # Design and contrast matrix
# design_full = model.matrix(~0 + full_pdata$Tissue_type + full_pdata$Study)
# colnames(design_full) = c("chronic_pancreatitis", "non_tumor",
#                           "GSE101462", "GSE143754", "GSE61166", "GSE71989", "GSE77858", "GSE91035")
# rownames(design_full) = colnames(union_z_exprs_final)
# cont.matrix_full = makeContrasts(onevsnormal = (Stage_1a + Stage_1b)/2 - normal, 
#                                  twovsnormal = (Stage_2a + Stage_2b)/2 - normal,
#                                  threevsnormal = Stage_3 - normal, 
#                                  fourvsnormal = Stage_4 - normal, 
#                                  earlyvslate = (Stage_1a + Stage_1b + Stage_2a + Stage_2b)/4 - (Stage_3 + Stage_4)/2,
#                                  onevsfour = (Stage_1a + Stage_1b)/2 - Stage_4, 
#                                  onevstwo = (Stage_1a + Stage_1b)/2 - (Stage_2a + Stage_2b)/2, levels = design_full)
# 
# 
# 
# # --------------------------------------------- #
# # --------------------------------------------- #
# # --------------------------------------------- #
# # --------------------------------------------- #
# # --------------------------------------------- #
# # --------------------------------------------- #
# 
# 
# # Write out the union_z_exprs_final as an RDS file
# saveRDS(union_z_exprs_final, "tumor_expression.rds")
# 
# # Model fitting with limma
# union_z_fit = lmFit(union_z_exprs_final, design_full)
# union_z_fit2 = contrasts.fit(union_z_fit, cont.matrix_full)
# union_z_fit2 = eBayes(union_z_fit2, robust = TRUE)
# union_z_results = decideTests(union_z_fit2)
# results = as.data.frame(cbind(summary(union_z_results), rownames(summary(union_z_results))))
# colnames(results)[8] = "direction"
# results = results %>% dplyr::select(direction, everything())
# 
# DGEA_topTables = createWorkbook()
# addWorksheet(DGEA_topTables, "summary")
# writeData(DGEA_topTables, "summary", results)
# 
# DE_maps = list()
# 
# for (i in 1:ncol(union_z_fit2)){
#   union_z_DE = as.data.frame(topTable(union_z_fit2, adjust.method="BH", 
#                                       number = Inf, coef = colnames(union_z_fit2)[i]))
#   union_z_DE$EntrezGene.ID = rownames(union_z_DE)
#   
#   # Annotation with official gene symbols
#   union_z_DE_mapped = union_z_DE %>% left_join(ID_Map, by = "EntrezGene.ID")
#   union_z_DE_mapped$Filter = NA
#   unmapped = which(is.na(union_z_DE_mapped$HGNC_Official))
#   union_z_DE_mapped$HGNC_Official[unmapped] = "unmapped"
#   for(j in 1:nrow(union_z_DE_mapped)){
#     if(union_z_DE_mapped$HGNC_Official[j] == "Yes"){
#       union_z_DE_mapped$Filter[j] = "Keep"
#     } else if(length(unique(union_z_DE_mapped$HGNC_Official[union_z_DE_mapped$EntrezGene.ID ==
#                                                             union_z_DE_mapped$EntrezGene.ID[j]])) > 1 &&
#               union_z_DE_mapped$HGNC_Official[j] == "No"){
#       union_z_DE_mapped$Filter[j] = "Discard"
#     } else if(unique(union_z_DE_mapped$HGNC_Official[union_z_DE_mapped$EntrezGene.ID ==
#                                                      union_z_DE_mapped$EntrezGene.ID[j]]) == "No"){
#       union_z_DE_mapped$Filter[j] = "Keep"
#       union_z_DE_mapped$Gene.Symbol[j] = union_z_DE_mapped$EntrezGene.ID[j]
#     } else if(union_z_DE_mapped$HGNC_Official[j] == "unmapped"){
#       union_z_DE_mapped$Gene.Symbol[j] = union_z_DE_mapped$EntrezGene.ID[j]
#       union_z_DE_mapped$Filter[j] = "Keep"
#     }
#   }
#   
#   union_z_DE_mapped = union_z_DE_mapped %>% 
#     dplyr::filter(Filter == "Keep") %>%
#     dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
#     dplyr::select(-Filter) %>%
#     distinct()
#   union_z_DE_mapped = union_z_DE_mapped[order(union_z_DE_mapped$adj.P.Val),]
#   rownames(union_z_DE_mapped) = union_z_DE_mapped$EntrezGene.ID
#   addWorksheet(DGEA_topTables, colnames(union_z_fit2)[i])
#   writeData(DGEA_topTables, colnames(union_z_fit2)[i], union_z_DE_mapped)
#   DE_maps[[i]] = union_z_DE_mapped
# }
# 
# saveWorkbook(DGEA_topTables, "DGEA/Union/DGEA_results.xlsx", overwrite = TRUE)
# names(DE_maps) = colnames(union_z_fit2)
# 
# ##### Condcordance and Discordance between stages #####
# # Stage 1 vs. Normal / Stage 2 vs. Normal
# # Save differences between Stage1/Normal and Stage2/Normal DEGs
# # in a variable called "dichotomizers". However, the real dichotomizers are both
# # genes not in common in the two lists as well as genes which are found to be
# # stat. sig. diff. expressed towards opposite directions
# 
# common_genes_1_2 = intersect(DE_maps[["twovsnormal"]]$EntrezGene.ID[DE_maps[["twovsnormal"]]$adj.P.Val<0.05],
#                              DE_maps[["onevsnormal"]]$EntrezGene.ID[DE_maps[["onevsnormal"]]$adj.P.Val<0.05])
# 
# gene_union = union(DE_maps[["twovsnormal"]]$EntrezGene.ID[DE_maps[["twovsnormal"]]$adj.P.Val<0.05],
#                    DE_maps[["onevsnormal"]]$EntrezGene.ID[DE_maps[["onevsnormal"]]$adj.P.Val<0.05])
# 
# diffDEGs_1_2 = as.data.frame(gene_union[!gene_union %in% common_genes_1_2])
# colnames(diffDEGs_1_2) = "EntrezGene.ID"
# diffDEGs_1_2 = diffDEGs_1_2 %>% inner_join(DE_maps[["onevsnormal"]], by = "EntrezGene.ID") %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# 
# # We need to establish which of the overlapping genes are differentially expressed
# # towards the same direction (up-/down-regulated):
# 
# stage_1_subset = DE_maps[["onevsnormal"]][DE_maps[["onevsnormal"]]$EntrezGene.ID
#                                           %in% common_genes_1_2, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_1 = logFC, adj_p_val_stage_1 = adj.P.Val)
# stage_2_subset = DE_maps[["twovsnormal"]][DE_maps[["twovsnormal"]]$EntrezGene.ID
#                                           %in% common_genes_1_2, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_2 = logFC, adj_p_val_stage_2 = adj.P.Val)
# common_set_1_2 = inner_join(stage_1_subset, stage_2_subset, by = c("EntrezGene.ID", "Gene.Symbol"))
# common_set_1_2$concordance = ifelse(common_set_1_2$logFC_stage_1*common_set_1_2$logFC_stage_2 > 0, 1, 0)
# concordant_set_1_2 = common_set_1_2 %>%
#   dplyr::filter(concordance == 1)
# discordant_set_1_2 = common_set_1_2 %>%
#   dplyr::filter(concordance == 0)
# concordant_set_1_2 = concordant_set_1_2[order(concordant_set_1_2$adj_p_val_stage_1, 
#                                               concordant_set_1_2$adj_p_val_stage_2), ]
# discordant_set_1_2 = discordant_set_1_2[order(discordant_set_1_2$adj_p_val_stage_1, 
#                                               discordant_set_1_2$adj_p_val_stage_2), ]
# cwb = createWorkbook()
# addWorksheet(cwb, "Concordance")
# writeData(cwb, "Concordance", concordant_set_1_2)
# addWorksheet(cwb, "Discordance")
# writeData(cwb, "Discordance", discordant_set_1_2)
# saveWorkbook(cwb, file = "DGEA/Stage_1_Stage_2_union_comparison.xlsx",
#              overwrite = TRUE); rm(cwb)
# 
# # 9727 genes are stat. sig. diff. expressed towards the same direction between
# # the two comparisons. We are interested in the remaining genes which were 
# # stat. sig. diff. expressed towards different directions. 0 genes in the discordant set
# 
# discordants_1_2 = discordant_set_1_2 %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# dichotomizers_1_2 = rbind(discordants_1_2, diffDEGs_1_2)
# dichotomizers_1_2 = dichotomizers_1_2 %>%
#   left_join(DE_maps[["onevsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_1 = logFC, adj.P.Val_Stage_1 = adj.P.Val) %>%
#   left_join(DE_maps[["twovsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC_Stage_1, adj.P.Val_Stage_1, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_2 = logFC, adj.P.Val_Stage_2 = adj.P.Val)
# 
# # 6821 genes
# dichotomizers = createWorkbook()
# addWorksheet(dichotomizers, "Stage_1_vs_2")
# writeData(dichotomizers, "Stage_1_vs_2", dichotomizers_1_2)
# 
# # Additional concordance and dichotomizers #####
# 
# # Overall concordance between Stage 2 vs. Normal/3 vs. Normal, Stage 2 vs. Normal/4 vs. Normal,
# # Stage 3 vs. Normal/4 vs. Normal
# 
# # Stage 2 vs. Normal/3 vs. Normal: 11437/12934
# # length(intersect(DE_maps[["threevsnormal"]]$EntrezGene.ID[DE_maps[["threevsnormal"]]$adj.P.Val<0.05],
# # DE_maps[["twovsnormal"]]$EntrezGene.ID[DE_maps[["twovsnormal"]]$adj.P.Val<0.05]))
# 
# # Stage 2 vs. Normal/4 vs. Normal: 10249/11512
# # length(intersect(DE_maps[["fourvsnormal"]]$EntrezGene.ID[DE_maps[["fourvsnormal"]]$adj.P.Val<0.05],
# # DE_maps[["twovsnormal"]]$EntrezGene.ID[DE_maps[["twovsnormal"]]$adj.P.Val<0.05]))
# 
# # Stage 3 vs. Normal/4 vs. Normal: 9537/11512
# # length(intersect(DE_maps[["fourvsnormal"]]$EntrezGene.ID[DE_maps[["fourvsnormal"]]$adj.P.Val<0.05],
# # DE_maps[["threevsnormal"]]$EntrezGene.ID[DE_maps[["threevsnormal"]]$adj.P.Val<0.05]))
# 
# # Stage 1 vs. Normal / Stage 3 vs. Normal
# # Save differences between Stage1/Normal and Stage3/Normal DEGs
# # in a variable called "dichotomizers". However, the real dichotomizers are both
# # genes not in common in the two lists as well as genes which are found to be
# # stat. sig. diff. expressed towards opposite directions
# 
# common_genes_1_3 = intersect(DE_maps[["threevsnormal"]]$EntrezGene.ID[DE_maps[["threevsnormal"]]$adj.P.Val<0.05],
#                              DE_maps[["onevsnormal"]]$EntrezGene.ID[DE_maps[["onevsnormal"]]$adj.P.Val<0.05])
# 
# gene_union = union(DE_maps[["threevsnormal"]]$EntrezGene.ID[DE_maps[["threevsnormal"]]$adj.P.Val<0.05],
#                    DE_maps[["onevsnormal"]]$EntrezGene.ID[DE_maps[["onevsnormal"]]$adj.P.Val<0.05])
# 
# diffDEGs_1_3 = as.data.frame(gene_union[!gene_union %in% common_genes_1_3])
# colnames(diffDEGs_1_3) = "EntrezGene.ID"
# diffDEGs_1_3 = diffDEGs_1_3 %>% inner_join(DE_maps[["onevsnormal"]], by = "EntrezGene.ID") %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# 
# # We need to establish which of the overlapping genes are differentially expressed
# # towards the same direction (up-/down-regulated):
# 
# stage_1_subset = DE_maps[["onevsnormal"]][DE_maps[["onevsnormal"]]$EntrezGene.ID
#                                           %in% common_genes_1_3, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_1 = logFC, adj_p_val_stage_1 = adj.P.Val)
# stage_3_subset = DE_maps[["threevsnormal"]][DE_maps[["threevsnormal"]]$EntrezGene.ID
#                                             %in% common_genes_1_3, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_3 = logFC, adj_p_val_stage_3 = adj.P.Val)
# common_set_1_3 = inner_join(stage_1_subset, stage_3_subset, by = c("EntrezGene.ID", "Gene.Symbol"))
# common_set_1_3$concordance = ifelse(common_set_1_3$logFC_stage_1*common_set_1_3$logFC_stage_3 > 0, 1, 0)
# concordant_set_1_3 = common_set_1_3 %>%
#   dplyr::filter(concordance == 1)
# discordant_set_1_3 = common_set_1_3 %>%
#   dplyr::filter(concordance == 0)
# concordant_set_1_3 = concordant_set_1_3[order(concordant_set_1_3$adj_p_val_stage_1, 
#                                               concordant_set_1_3$adj_p_val_stage_3), ]
# discordant_set_1_3 = discordant_set_1_3[order(discordant_set_1_3$adj_p_val_stage_1, 
#                                               discordant_set_1_3$adj_p_val_stage_3), ]
# cwb = createWorkbook()
# addWorksheet(cwb, "Concordance")
# writeData(cwb, "Concordance", concordant_set_1_3)
# addWorksheet(cwb, "Discordance")
# writeData(cwb, "Discordance", discordant_set_1_3)
# saveWorkbook(cwb, file = "DGEA/Stage_1_Stage_3_union_comparison.xlsx",
#              overwrite = TRUE); rm(cwb)
# 
# # 7607 genes are stat. sig. diff. expressed towards the same direction between
# # the two comparisons. We are interested in the remaining genes which were 
# # stat. sig. diff. expressed towards different directions. 0 discordant genes
# 
# discordants_1_3 = discordant_set_1_3 %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# dichotomizers_1_3 = rbind(discordants_1_3, diffDEGs_1_3)
# dichotomizers_1_3 = dichotomizers_1_3 %>%
#   left_join(DE_maps[["onevsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_1 = logFC, adj.P.Val_Stage_1 = adj.P.Val) %>%
#   left_join(DE_maps[["threevsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC_Stage_1, adj.P.Val_Stage_1, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_3 = logFC, adj.P.Val_Stage_3 = adj.P.Val)
# 
# # 5185 genes
# addWorksheet(dichotomizers, "Stage_1_vs_3")
# writeData(dichotomizers, "Stage_1_vs_3", dichotomizers_1_3)
# 
# # Stage 1 vs. Normal / Stage 4 vs. Normal
# # Save differences between Stage1/Normal and Stage4/Normal DEGs
# # in a variable called "dichotomizers". However, the real dichotomizers are both
# # genes not in common in the two lists as well as genes which are found to be
# # stat. sig. diff. expressed towards opposite directions
# 
# common_genes_1_4 = intersect(DE_maps[["fourvsnormal"]]$EntrezGene.ID[DE_maps[["fourvsnormal"]]$adj.P.Val<0.05],
#                              DE_maps[["onevsnormal"]]$EntrezGene.ID[DE_maps[["onevsnormal"]]$adj.P.Val<0.05])
# 
# gene_union = union(DE_maps[["fourvsnormal"]]$EntrezGene.ID[DE_maps[["fourvsnormal"]]$adj.P.Val<0.05],
#                    DE_maps[["onevsnormal"]]$EntrezGene.ID[DE_maps[["onevsnormal"]]$adj.P.Val<0.05])
# 
# diffDEGs_1_4 = as.data.frame(gene_union[!gene_union %in% common_genes_1_4])
# colnames(diffDEGs_1_4) = "EntrezGene.ID"
# diffDEGs_1_4 = diffDEGs_1_4 %>% inner_join(DE_maps[["onevsnormal"]], by = "EntrezGene.ID") %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# 
# # We need to establish which of the overlapping genes are differentially expressed
# # towards the same direction (up-/down-regulated):
# 
# stage_1_subset = DE_maps[["onevsnormal"]][DE_maps[["onevsnormal"]]$EntrezGene.ID
#                                           %in% common_genes_1_4, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_1 = logFC, adj_p_val_stage_1 = adj.P.Val)
# stage_4_subset = DE_maps[["fourvsnormal"]][DE_maps[["fourvsnormal"]]$EntrezGene.ID
#                                            %in% common_genes_1_4, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_4 = logFC, adj_p_val_stage_4 = adj.P.Val)
# common_set_1_4 = inner_join(stage_1_subset, stage_4_subset, by = c("EntrezGene.ID", "Gene.Symbol"))
# common_set_1_4$concordance = ifelse(common_set_1_4$logFC_stage_1*common_set_1_4$logFC_stage_4 > 0, 1, 0)
# concordant_set_1_4 = common_set_1_4 %>%
#   dplyr::filter(concordance == 1)
# discordant_set_1_4 = common_set_1_4 %>%
#   dplyr::filter(concordance == 0)
# concordant_set_1_4 = concordant_set_1_4[order(concordant_set_1_4$adj_p_val_stage_1, 
#                                               concordant_set_1_4$adj_p_val_stage_4), ]
# discordant_set_1_4 = discordant_set_1_4[order(discordant_set_1_4$adj_p_val_stage_1, 
#                                               discordant_set_1_4$adj_p_val_stage_4), ]
# cwb = createWorkbook()
# addWorksheet(cwb, "Concordance")
# writeData(cwb, "Concordance", concordant_set_1_4)
# addWorksheet(cwb, "Discordance")
# writeData(cwb, "Discordance", discordant_set_1_4)
# saveWorkbook(cwb, file = "DGEA/Stage_1_Stage_4_union_comparison.xlsx",
#              overwrite = TRUE); rm(cwb)
# 
# # 6591 genes are stat. sig. diff. expressed towards the same direction between
# # the two comparisons. We are interested in the remaining genes which were 
# # stat. sig. diff. expressed towards different directions. 0 genes are discordant
# 
# discordants_1_4 = discordant_set_1_4 %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# dichotomizers_1_4 = rbind(discordants_1_4, diffDEGs_1_4)
# dichotomizers_1_4 = dichotomizers_1_4 %>%
#   left_join(DE_maps[["onevsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_1 = logFC, adj.P.Val_Stage_1 = adj.P.Val) %>%
#   left_join(DE_maps[["fourvsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC_Stage_1, adj.P.Val_Stage_1, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_4 = logFC, adj.P.Val_Stage_4 = adj.P.Val)
# 
# # 5105 genes
# addWorksheet(dichotomizers, "Stage_1_vs_4")
# writeData(dichotomizers, "Stage_1_vs_4", dichotomizers_1_4)
# 
# # Stage 2 vs. Normal / Stage 4 vs. Normal
# # Save differences between Stage2/Normal and Stage4/Normal DEGs
# # in a variable called "dichotomizers". However, the real dichotomizers are both
# # genes not in common in the two lists as well as genes which are found to be
# # stat. sig. diff. expressed towards opposite directions
# 
# common_genes_2_4 = intersect(DE_maps[["twovsnormal"]]$EntrezGene.ID[DE_maps[["twovsnormal"]]$adj.P.Val<0.05],
#                              DE_maps[["fourvsnormal"]]$EntrezGene.ID[DE_maps[["fourvsnormal"]]$adj.P.Val<0.05])
# 
# gene_union = union(DE_maps[["twovsnormal"]]$EntrezGene.ID[DE_maps[["twovsnormal"]]$adj.P.Val<0.05],
#                    DE_maps[["fourvsnormal"]]$EntrezGene.ID[DE_maps[["fourvsnormal"]]$adj.P.Val<0.05])
# 
# diffDEGs_2_4 = as.data.frame(gene_union[!gene_union %in% common_genes_2_4])
# colnames(diffDEGs_2_4) = "EntrezGene.ID"
# diffDEGs_2_4 = diffDEGs_2_4 %>% inner_join(DE_maps[["fourvsnormal"]], by = "EntrezGene.ID") %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# 
# # We need to establish which of the overlapping genes are differentially expressed
# # towards the same direction (up-/down-regulated):
# 
# stage_2_subset = DE_maps[["twovsnormal"]][DE_maps[["twovsnormal"]]$EntrezGene.ID
#                                           %in% common_genes_2_4, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_2 = logFC, adj_p_val_stage_2 = adj.P.Val)
# stage_4_subset = DE_maps[["fourvsnormal"]][DE_maps[["fourvsnormal"]]$EntrezGene.ID
#                                            %in% common_genes_2_4, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_4 = logFC, adj_p_val_stage_4 = adj.P.Val)
# common_set_2_4 = inner_join(stage_2_subset, stage_4_subset, by = c("EntrezGene.ID", "Gene.Symbol"))
# common_set_2_4$concordance = ifelse(common_set_2_4$logFC_stage_2*common_set_2_4$logFC_stage_4 > 0, 1, 0)
# concordant_set_2_4 = common_set_2_4 %>%
#   dplyr::filter(concordance == 1)
# discordant_set_2_4 = common_set_2_4 %>%
#   dplyr::filter(concordance == 0)
# concordant_set_2_4 = concordant_set_2_4[order(concordant_set_2_4$adj_p_val_stage_2, 
#                                               concordant_set_2_4$adj_p_val_stage_4), ]
# discordant_set_2_4 = discordant_set_2_4[order(discordant_set_2_4$adj_p_val_stage_2, 
#                                               discordant_set_2_4$adj_p_val_stage_4), ]
# cwb = createWorkbook()
# addWorksheet(cwb, "Concordance")
# writeData(cwb, "Concordance", concordant_set_2_4)
# addWorksheet(cwb, "Discordance")
# writeData(cwb, "Discordance", discordant_set_2_4)
# saveWorkbook(cwb, file = "DGEA/Stage_2_Stage_4_union_comparison.xlsx",
#              overwrite = TRUE); rm(cwb)
# 
# # 8111 genes are stat. sig. diff. expressed towards the same direction between
# # the two comparisons. We are interested in the remaining genes which were 
# # stat. sig. diff. expressed towards different directions. 0 genes are discordant
# 
# discordants_2_4 = discordant_set_2_4 %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# dichotomizers_2_4 = rbind(discordants_2_4, diffDEGs_2_4)
# dichotomizers_2_4 = dichotomizers_2_4 %>%
#   left_join(DE_maps[["twovsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_2 = logFC, adj.P.Val_Stage_2 = adj.P.Val) %>%
#   left_join(DE_maps[["fourvsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC_Stage_2, adj.P.Val_Stage_2, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_4 = logFC, adj.P.Val_Stage_4 = adj.P.Val)
# 
# # 8446 genes
# addWorksheet(dichotomizers, "Stage_2_vs_4")
# writeData(dichotomizers, "Stage_2_vs_4", dichotomizers_2_4)
# 
# # Stage 2 vs. Normal / Stage 3 vs. Normal
# # Save differences between Stage2/Normal and Stage3/Normal DEGs
# # in a variable called "dichotomizers". However, the real dichotomizers are both
# # genes not in common in the two lists as well as genes which are found to be
# # stat. sig. diff. expressed towards opposite directions
# 
# common_genes_2_3 = intersect(DE_maps[["twovsnormal"]]$EntrezGene.ID[DE_maps[["twovsnormal"]]$adj.P.Val<0.05],
#                              DE_maps[["threevsnormal"]]$EntrezGene.ID[DE_maps[["threevsnormal"]]$adj.P.Val<0.05])
# 
# gene_union = union(DE_maps[["twovsnormal"]]$EntrezGene.ID[DE_maps[["twovsnormal"]]$adj.P.Val<0.05],
#                    DE_maps[["threevsnormal"]]$EntrezGene.ID[DE_maps[["threevsnormal"]]$adj.P.Val<0.05])
# 
# diffDEGs_2_3 = as.data.frame(gene_union[!gene_union %in% common_genes_2_3])
# colnames(diffDEGs_2_3) = "EntrezGene.ID"
# diffDEGs_2_3 = diffDEGs_2_3 %>% inner_join(DE_maps[["threevsnormal"]], by = "EntrezGene.ID") %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# 
# # We need to establish which of the overlapping genes are differentially expressed
# # towards the same direction (up-/down-regulated):
# 
# stage_2_subset = DE_maps[["twovsnormal"]][DE_maps[["twovsnormal"]]$EntrezGene.ID
#                                           %in% common_genes_2_3, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_2 = logFC, adj_p_val_stage_2 = adj.P.Val)
# stage_3_subset = DE_maps[["threevsnormal"]][DE_maps[["threevsnormal"]]$EntrezGene.ID
#                                             %in% common_genes_2_3, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_3 = logFC, adj_p_val_stage_3 = adj.P.Val)
# common_set_2_3 = inner_join(stage_2_subset, stage_3_subset, by = c("EntrezGene.ID", "Gene.Symbol"))
# common_set_2_3$concordance = ifelse(common_set_2_3$logFC_stage_2*common_set_2_3$logFC_stage_3 > 0, 1, 0)
# concordant_set_2_3 = common_set_2_3 %>%
#   dplyr::filter(concordance == 1)
# discordant_set_2_3 = common_set_2_3 %>%
#   dplyr::filter(concordance == 0)
# concordant_set_2_3 = concordant_set_2_3[order(concordant_set_2_3$adj_p_val_stage_2, 
#                                               concordant_set_2_3$adj_p_val_stage_3), ]
# discordant_set_2_3 = discordant_set_2_3[order(discordant_set_2_3$adj_p_val_stage_2, 
#                                               discordant_set_2_3$adj_p_val_stage_3), ]
# cwb = createWorkbook()
# addWorksheet(cwb, "Concordance")
# writeData(cwb, "Concordance", concordant_set_2_3)
# addWorksheet(cwb, "Discordance")
# writeData(cwb, "Discordance", discordant_set_2_3)
# saveWorkbook(cwb, file = "DGEA/Stage_2_Stage_3_union_comparison.xlsx",
#              overwrite = TRUE); rm(cwb)
# 
# # 10112 genes are stat. sig. diff. expressed towards the same direction between
# # the two comparisons. We are interested in the remaining genes which were 
# # stat. sig. diff. expressed towards different directions. 0 genes are discordant
# # DNAH17
# 
# discordants_2_3 = discordant_set_2_3 %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# dichotomizers_2_3 = rbind(discordants_2_3, diffDEGs_2_3)
# dichotomizers_2_3 = dichotomizers_2_3 %>%
#   left_join(DE_maps[["twovsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_2 = logFC, adj.P.Val_Stage_2 = adj.P.Val) %>%
#   left_join(DE_maps[["threevsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC_Stage_2, adj.P.Val_Stage_2, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_3 = logFC, adj.P.Val_Stage_3 = adj.P.Val)
# 
# # 6556 genes
# addWorksheet(dichotomizers, "Stage_2_vs_3")
# writeData(dichotomizers, "Stage_2_vs_3", dichotomizers_2_3)
# 
# # Stage 3 vs. Normal / Stage 4 vs. Normal
# # Save differences between Stage3/Normal and Stage4/Normal DEGs
# # in a variable called "dichotomizers". However, the real dichotomizers are both
# # genes not in common in the two lists as well as genes which are found to be
# # stat. sig. diff. expressed towards opposite directions
# 
# common_genes_3_4 = intersect(DE_maps[["threevsnormal"]]$EntrezGene.ID[DE_maps[["threevsnormal"]]$adj.P.Val<0.05],
#                              DE_maps[["fourvsnormal"]]$EntrezGene.ID[DE_maps[["fourvsnormal"]]$adj.P.Val<0.05])
# 
# gene_union = union(DE_maps[["threevsnormal"]]$EntrezGene.ID[DE_maps[["threevsnormal"]]$adj.P.Val<0.05],
#                    DE_maps[["fourvsnormal"]]$EntrezGene.ID[DE_maps[["fourvsnormal"]]$adj.P.Val<0.05])
# 
# diffDEGs_3_4 = as.data.frame(gene_union[!gene_union %in% common_genes_3_4])
# colnames(diffDEGs_3_4) = "EntrezGene.ID"
# diffDEGs_3_4 = diffDEGs_3_4 %>% inner_join(DE_maps[["fourvsnormal"]], by = "EntrezGene.ID") %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# 
# # We need to establish which of the overlapping genes are differentially expressed
# # towards the same direction (up-/down-regulated):
# 
# stage_3_subset = DE_maps[["threevsnormal"]][DE_maps[["threevsnormal"]]$EntrezGene.ID
#                                             %in% common_genes_3_4, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_3 = logFC, adj_p_val_stage_3 = adj.P.Val)
# stage_4_subset = DE_maps[["fourvsnormal"]][DE_maps[["fourvsnormal"]]$EntrezGene.ID
#                                            %in% common_genes_3_4, ] %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_stage_4 = logFC, adj_p_val_stage_4 = adj.P.Val)
# common_set_3_4 = inner_join(stage_3_subset, stage_4_subset, by = c("EntrezGene.ID", "Gene.Symbol"))
# common_set_3_4$concordance = ifelse(common_set_3_4$logFC_stage_3*common_set_3_4$logFC_stage_4 > 0, 1, 0)
# concordant_set_3_4 = common_set_3_4 %>%
#   dplyr::filter(concordance == 1)
# discordant_set_3_4 = common_set_3_4 %>%
#   dplyr::filter(concordance == 0)
# concordant_set_3_4 = concordant_set_3_4[order(concordant_set_3_4$adj_p_val_stage_3, 
#                                               concordant_set_3_4$adj_p_val_stage_4), ]
# discordant_set_3_4 = discordant_set_3_4[order(discordant_set_3_4$adj_p_val_stage_3, 
#                                               discordant_set_3_4$adj_p_val_stage_4), ]
# cwb = createWorkbook()
# addWorksheet(cwb, "Concordance")
# writeData(cwb, "Concordance", concordant_set_3_4)
# addWorksheet(cwb, "Discordance")
# writeData(cwb, "Discordance", discordant_set_3_4)
# saveWorkbook(cwb, file = "DGEA/Stage_3_Stage_4_union_comparison.xlsx",
#              overwrite = TRUE); rm(cwb)
# 
# # 6872 genes are stat. sig. diff. expressed towards the same direction between
# # the two comparisons. We are interested in the remaining genes which were 
# # stat. sig. diff. expressed towards different directions. 0 genes are left!
# 
# discordants_3_4 = discordant_set_3_4 %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol)
# dichotomizers_3_4 = rbind(discordants_3_4, diffDEGs_3_4)
# dichotomizers_3_4 = dichotomizers_3_4 %>%
#   left_join(DE_maps[["threevsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_3 = logFC, adj.P.Val_Stage_3 = adj.P.Val) %>%
#   left_join(DE_maps[["fourvsnormal"]], by = c("EntrezGene.ID", "Gene.Symbol")) %>%
#   dplyr::select(EntrezGene.ID, Gene.Symbol, logFC_Stage_3, adj.P.Val_Stage_3, logFC, adj.P.Val) %>%
#   dplyr::rename(logFC_Stage_4 = logFC, adj.P.Val_Stage_4 = adj.P.Val)
# 
# # 5048 genes
# addWorksheet(dichotomizers, "Stage_3_vs_4")
# writeData(dichotomizers, "Stage_3_vs_4", dichotomizers_3_4)
# saveWorkbook(dichotomizers, "DGEA/Union/Dichotomizers.xlsx", overwrite = TRUE)
# 
# # Union volcanoes #####
# union_stages_volcano = EnhancedVolcano(DE_maps[["earlyvslate"]],
#                                        lab = DE_maps[["earlyvslate"]][, "Gene.Symbol"],
#                                        x = 'logFC',
#                                        y = 'adj.P.Val',
#                                        title = "Stage 1/2 vs. Stage 3/4",
#                                        pCutoff = 0.05,
#                                        FCcutoff = 1,
#                                        col=c('grey', 'pink', 'purple4', 'red4'),
#                                        colAlpha = 0.7,
#                                        xlim = c(-5, 5),
#                                        ylab = bquote(~-Log[10] ~ (italic(adj.p.value))),
#                                        xlab = "\nDifferential expression (units: sd)",
#                                        legendLabels = c("NS", "|DE| > 1 s.d.", "FDR < 0.05", "FDR < 0.05 & |DE| > 1 s.d."))
# tiff("DGEA/Union/Union_Stages_Volcano.tif", width = 1500, height = 1920, res = 100)
# union_stages_volcano
# dev.off()
# 
# union_one_four_stages_volcano = EnhancedVolcano(DE_maps[["onevsfour"]],
#                                                 lab = DE_maps[["onevsfour"]][, "Gene.Symbol"],
#                                                 x = 'logFC',
#                                                 y = 'adj.P.Val',
#                                                 title = "Stage 1 vs. Stage 4",
#                                                 pCutoff = 0.05,
#                                                 FCcutoff = 1,
#                                                 col=c('grey', 'pink', 'purple4', 'red4'),
#                                                 colAlpha = 0.7,
#                                                 xlim = c(-5, 5),
#                                                 ylab = bquote(~-Log[10] ~ (italic(adj.p.value))),
#                                                 xlab = "\nDifferential expression (units: sd)",
#                                                 legendLabels = c("NS", "|DE| > 1 s.d.", "FDR < 0.05", "FDR < 0.05 & |DE| > 1 s.d."))
# tiff("DGEA/Union/Union_One_Four_Stages_Volcano.tif", width = 1500, height = 1920, res = 100)
# union_one_four_stages_volcano
# dev.off()
# 
# union_one_normal_volcano = EnhancedVolcano(DE_maps[["onevsnormal"]],
#                                            lab = DE_maps[["onevsnormal"]][, "Gene.Symbol"],
#                                            x = 'logFC',
#                                            y = 'adj.P.Val',
#                                            title = "Stage 1 vs. Normal",
#                                            pCutoff = 0.05,
#                                            FCcutoff = 1,
#                                            col=c('grey', 'pink', 'purple4', 'red4'),
#                                            colAlpha = 0.7,
#                                            xlim = c(-5, 5),
#                                            ylim = c(0, 10),
#                                            ylab = bquote(~-Log[10] ~ (italic(adj.p.value))),
#                                            xlab = "\nDifferential expression (units: sd)",
#                                            legendLabels = c("NS", "|DE| > 1 s.d.", "FDR < 0.05", "FDR < 0.05 & |DE| > 1 s.d."))
# tiff("DGEA/Union/Union_One_Normal_Volcano.tif", width = 1500, height = 1920, res = 100)
# union_one_normal_volcano
# dev.off()
# 
# union_two_normal_volcano = EnhancedVolcano(DE_maps[["twovsnormal"]],
#                                            lab = DE_maps[["twovsnormal"]][, "Gene.Symbol"],
#                                            x = 'logFC',
#                                            y = 'adj.P.Val',
#                                            title = "Stage 2 vs. Normal",
#                                            pCutoff = 0.05,
#                                            FCcutoff = 1,
#                                            col=c('grey', 'pink', 'purple4', 'red4'),
#                                            colAlpha = 0.7,
#                                            xlim = c(-5, 5),
#                                            ylim = c(0, 35),
#                                            ylab = bquote(~-Log[10] ~ (italic(adj.p.value))),
#                                            xlab = "\nDifferential expression (units: sd)",
#                                            legendLabels = c("NS", "|DE| > 1 s.d.", "FDR < 0.05", "FDR < 0.05 & |DE| > 1 s.d."))
# tiff("DGEA/Union/Union_Two_Normal_Volcano.tif", width = 1500, height = 1920, res = 100)
# union_two_normal_volcano
# dev.off()
# 
# union_three_normal_volcano = EnhancedVolcano(DE_maps[["threevsnormal"]],
#                                              lab = DE_maps[["threevsnormal"]][, "Gene.Symbol"],
#                                              x = 'logFC',
#                                              y = 'adj.P.Val',
#                                              title = "Stage 3 vs. Normal",
#                                              pCutoff = 0.05,
#                                              FCcutoff = 1,
#                                              col=c('grey', 'pink', 'purple4', 'red4'),
#                                              colAlpha = 0.7,
#                                              xlim = c(-5, 5),
#                                              ylim = c(0, 12.5),
#                                              ylab = bquote(~-Log[10] ~ (italic(adj.p.value))),
#                                              xlab = "\nDifferential expression (units: sd)",
#                                              legendLabels = c("NS", "|DE| > 1 s.d.", "FDR < 0.05", "FDR < 0.05 & |DE| > 1 s.d."))
# tiff("DGEA/Union/Union_Three_Normal_Volcano.tif", width = 1500, height = 1920, res = 100)
# union_three_normal_volcano
# dev.off()
# 
# union_four_normal_volcano = EnhancedVolcano(DE_maps[["fourvsnormal"]],
#                                             lab = DE_maps[["fourvsnormal"]][, "Gene.Symbol"],
#                                             x = 'logFC',
#                                             y = 'adj.P.Val',
#                                             title = "Stage 4 vs. Normal",
#                                             pCutoff = 0.05,
#                                             FCcutoff = 1,
#                                             col=c('grey', 'pink', 'purple4', 'red4'),
#                                             colAlpha = 0.7,
#                                             xlim = c(-5, 5),
#                                             ylim = c(0, 15),
#                                             ylab = bquote(~-Log[10] ~ (italic(adj.p.value))),
#                                             xlab = "\nDifferential expression (units: sd)",
#                                             legendLabels = c("NS", "|DE| > 1 s.d.", "FDR < 0.05", "FDR < 0.05 & |DE| > 1 s.d."))
# tiff("DGEA/Union/Union_Four_Normal_Volcano.tif", width = 1500, height = 1920, res = 100)
# union_four_normal_volcano
# dev.off()
# 
# rm(gene_union, i, j)
# 
# tiff("DGEA/Union//Volcano_multiplot.tif", 
#      width = 3000, height = 3840, res = 150)
# multiplot(union_one_normal_volcano, union_two_normal_volcano,
#           union_three_normal_volcano, union_four_normal_volcano, cols = 2)
# m = ggplot(multiplot(union_one_normal_volcano, union_two_normal_volcano,
#                      union_three_normal_volcano, union_four_normal_volcano, cols = 2))
# dev.off(); rm(m)
# 
# # Useful to save the volcano plots only as an .RData file that will be later used
# # to create additional plots after pathway analysis:
# 
# # rm(list=setdiff(ls(), c("union_four_normal_volcano", "union_three_normal_volcano",
# # "union_two_normal_volcano", "union_one_normal_volcano")))
# 
# # Augment that later with the blood samples volcano
