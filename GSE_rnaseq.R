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

## NCBI recently released a workflow that generates RNA-seq count data for all the publicly available RNA-seq data archived by SRA
## https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html --> PAY ATTENTION TO THE LIMITATIONS MENTIONED
## We downloaded the raw counts for the two datasets from here:
## - GSE194331: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194331 (NCBI-generated data --> Series RNA-seq raw counts matrix)
## - GSE133684: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133684 (NCBI-generated data --> Series RNA-seq raw counts matrix)
## However, chronic pancreatitis raw counts for the GSE133684 are not available --> we download the TPM data from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133684

## Load raw counts
GSE194331_raw_counts = data.table::fread("data/GSE194331_raw_counts_GRCh38.p13_NCBI.tsv")
GSE133684_tpm = data.table::fread("data/GSE133684_TPM_all.txt.gz")

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
                                          levels = c("chronic_pancreatitis", "non_tumor", "tumor"),
                                          labels = c("chronic_pancreatitis", "non_tumor", "tumor"))

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

# Change Tissue_type to chronic pancreatitis and normal
for (i in 1:length(filt_pdata$GSE133684$Tissue_type)) {
  if (filt_pdata$GSE133684$Tissue_type[i] == "PDAC") {
    filt_pdata$GSE133684$Tissue_type[i] = "tumor"
  } 
  if (filt_pdata$GSE133684$Tissue_type[i] == "healthy") {
    filt_pdata$GSE133684$Tissue_type[i] = "non_tumor"
  } 
}; rm(i)

# Add to filt_pdata the CP patient IDs
cp_patients_ids = setdiff(colnames(GSE133684_tpm), filt_pdata$GSE133684$Patient_ID)[-1]
cp_data = data.frame(GEO_accession = "", Patient_ID = cp_patients_ids, Platform = "GPL20795", Tissue_type = "chronic_pancreatitis")
filt_pdata$GSE133684 = rbind(filt_pdata$GSE133684, cp_data) ; rm(cp_patients_ids, cp_data)

# Transform to factors with consistent universal levels
filt_pdata$GSE133684$Tissue_type = factor(x = filt_pdata$GSE133684$Tissue_type,
                                          levels = c("chronic_pancreatitis", "non_tumor", "tumor"),
                                          labels = c("chronic_pancreatitis", "non_tumor", "tumor"))

## ---------------------------
## Expression data
## ---------------------------

## GSE194331: convert raw counts to TPM

# Create a TPM matrix
library(AnnotationHub)

# Create an AnnotationHub object
ah = AnnotationHub()

# Query for the latest Homo sapiens EnsDb object
latest_homo_sapiens = query(ah, c("EnsDb", "Homo sapiens")) 

# Print the latest version
tail(latest_homo_sapiens)

# ID of interest is "AH113665"
edb = ah[["AH113665"]]

# Get gene data
genes_data = genes(edb)
genes_map = data.frame(gene_id = genes_data$gene_id,
                       entrezID = genes_data$entrezid)

# Get transcript data
transcripts_data = transcripts(edb)
transcript_lengths = data.frame(transcript_id = transcripts_data$tx_id,
                                gene_id = transcripts_data$gene_id,
                                tx_length = width(transcripts_data)) %>%
  inner_join(genes_map, by = "gene_id") %>%
  dplyr::filter(!Gene.Symbol == "") %>%
  group_by(Gene.Symbol) %>%
  summarise(longest_transcript_length = max(tx_length, na.rm = TRUE))

# Calculate gene lengths in kilobases
transcript_lengths$transcript_length_kb = transcript_lengths$longest_transcript_length / 1000

# Map gene lengths to rownames of the counts matrix
rownames(transcript_lengths) = transcript_lengths$Gene.Symbol
mapped = transcript_lengths[transcript_lengths$Gene.Symbol %in% 
                              x_filt$genes$Gene.Symbol,]
mapped_lengths = mapped$transcript_length_kb

# Calculate counts per kilobase (CPK)
cpk = x_filt$counts[mapped$Gene.Symbol, ] / mapped_lengths
# Calculate the sum of CPK values for each sample
cpk_sum = colSums(cpk)

# Calculate TPM values
tpm = cpk / cpk_sum * 1e6
log2tpm = log2(tpm + 1)
rm(cpk, cpk_sum, keep.exprs, mapped_lengths, ah, edb)







