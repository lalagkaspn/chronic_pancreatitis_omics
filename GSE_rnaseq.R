## This script is used to download and pre-process series (RNA-seq) from GEO for 
## patients with Chronic Pancreatitis (normal and CP tissues)

library(dplyr)
library(GEOquery) # updated

##### Downloading data #####
datasets = c("GSE194331", "GSE133684")

GEOsets = list()
for (i in 1:length(datasets)){
  GEOsets[[i]] = getGEO(datasets[i])
}; rm(i)
GEOsets = unlist(GEOsets)
names(GEOsets) = datasets; rm(datasets)

## NCBI recently released a workflow that generates RNA-seq count data for all the publicly available RNA-seq data arcivhed by SRA
## https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html --> PAY ATTENTION TO THE LIMITATIONS MENTIONED
## I downloaded the raw counts for the two datasets from here:
## - GSE194331: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194331 (NCBI-generated data --> Series RNA-seq raw counts matrix)
## - GSE133684: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133684 (NCBI-generated data --> Series RNA-seq raw counts matrix)

## Load raw counts
GSE194331_raw_counts = data.table::fread("/Users/panagiotisnikolaoslalagkas/Desktop/chronic_pancreatitis/data/GSE194331_raw_counts_GRCh38.p13_NCBI.tsv")
GSE133684_raw_counts = data.table::fread("/Users/panagiotisnikolaoslalagkas/Desktop/chronic_pancreatitis/data/GSE133684_raw_counts_GRCh38.p13_NCBI.tsv")

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
filt_pdata = list()

# GSE194331 
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194331
# Paper: https://pubmed.ncbi.nlm.nih.gov/35426393/
# Comments:
#   - 119 participants
#     - 57 mild CP
#     - 20 moderate-severe CP
#     - 10 severe CP
#     - 32 healthy

# Select informative data only
filt_pdata[["GSE194331"]] = pdata$GSE194331 %>%
  dplyr::select(GEO_accession = geo_accession,
                Patient_ID = title,
                Platform = platform_id,
                Tissue_type = "pathology:ch1")

# Change Tissue_type to chronic pancreatitis and normal
for (i in 1:length(filt_pdata$GSE194331$Tissue_type)) {
  if (filt_pdata$GSE194331$Tissue_type[i] != "Healthy control") {
    filt_pdata$GSE194331$Tissue_type[i] = "chronic_pancreatitis"
  } 
  if (filt_pdata$GSE194331$Tissue_type[i] == "Healthy control") {
    filt_pdata$GSE194331$Tissue_type[i] = "non_tumor"
  } 
}; rm(i)

# Transform to factors with consistent universal levels
filt_pdata$GSE194331$Tissue_type = factor(x = filt_pdata$GSE194331$Tissue_type,
                                          levels = c("chronic_pancreatitis", "non_tumor"),
                                          labels = c("chronic_pancreatitis", "non_tumor"))

# Check overlap of GSE samples with the NCBI raw count dataframe
length(intersect(filt_pdata$GSE194331$GEO_accession, colnames(GSE194331_raw_counts)[-1])) # 119 --> 100% overlap

## Filter both pdata and raw counts for common samples
common_samples_gse194331 = intersect(colnames(GSE194331_raw_counts), filt_pdata$GSE194331$GEO_accession)
columns_to_keep = c(1, which(colnames(GSE194331_raw_counts) %in% common_samples_gse194331))
GSE194331_raw_counts = GSE194331_raw_counts[, ..columns_to_keep]
filt_pdata$GSE194331 = filt_pdata$GSE194331 %>% filter(GEO_accession %in% common_samples_gse194331)
rm(common_samples_gse194331, columns_to_keep)

# GSE133684 
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133684
# Paper: https://pubmed.ncbi.nlm.nih.gov/31562239/
# Comments:
#   - 501 samples
#     - 284 PDAC
#     - 100 chronic pancreatitis
#     - 117 healthy

# Select informative data only
filt_pdata[["GSE133684"]] = pdata$GSE133684 %>%
  dplyr::select(GEO_accession = geo_accession,
                Patient_ID = title,
                Platform = platform_id,
                Tissue_type = "disease state:ch1")

# Samples downloaded through GSE contain only PDAC and healthy (284 + 117 = 401)
# Therefore, the samples found in the columns of the raw_counts data frame but not in the GSE data are the chronic pancreatitis
colnames_gse133684 = colnames(GSE133684_raw_counts)[-1]
GSE133684_cp_samples = setdiff(colnames_gse133684, filt_pdata[["GSE133684"]]$GEO_accession) # empty --> that means that NCBI raw count data frame contains only raw counts for PDAC and healthy samples

## Keep only normals
filt_pdata[["GSE133684"]] = filt_pdata[["GSE133684"]] %>%
  filter(Tissue_type == "healthy")

# Change Tissue_type to chronic pancreatitis and normal
for (i in 1:length(filt_pdata$GSE133684$Tissue_type)) {
  if (filt_pdata$GSE133684$Tissue_type[i] != "healthy") {
    filt_pdata$GSE133684$Tissue_type[i] = "chronic_pancreatitis"
  } 
  if (filt_pdata$GSE133684$Tissue_type[i] == "healthy") {
    filt_pdata$GSE133684$Tissue_type[i] = "non_tumor"
  } 
}; rm(i)

## Filter both pdata and raw counts for common samples
common_samples_GSE133684 = intersect(colnames(GSE133684_raw_counts), filt_pdata$GSE133684$GEO_accession)
columns_to_keep = c(1, which(colnames(GSE133684_raw_counts) %in% common_samples_GSE133684))
GSE133684_raw_counts = GSE133684_raw_counts[, ..columns_to_keep]
filt_pdata$GSE133684 = filt_pdata$GSE133684 %>% filter(GEO_accession %in% common_samples_GSE133684)
rm(common_samples_GSE133684, columns_to_keep)


