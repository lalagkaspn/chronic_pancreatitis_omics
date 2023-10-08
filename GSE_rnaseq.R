## This script is used to download and pre-process series (RNA-seq) from GEO for 
## patients with Chronic Pancreatitis (normal and CP tissues)

library(dplyr)
library(GEOquery) # updated
library(org.Hs.eg.db)

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
filt_pdata$GSE194331 = filt_pdata$GSE194331 %>% dplyr::filter(GEO_accession %in% common_samples_gse194331)
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
entrezid_temp = genes_data$entrezid

genes_map = data.frame(gene_id = NA, entrezID = NA)
for (i in 1:length(entrezid_temp)) {
  genes_map[i, "gene_id"] = names(entrezid_temp)[i]
  genes_map[i, "entrezID"] = paste0(entrezid_temp[[i]], collapse = ",")
  cat(i, "\n")
} ; rm(i)

# Get transcript data
transcripts_data = transcripts(edb)
transcript_lengths = data.frame(transcript_id = transcripts_data$tx_id,
                                gene_id = transcripts_data$gene_id,
                                tx_length = width(transcripts_data)) %>%
  inner_join(genes_map, by = "gene_id") %>%
  dplyr::filter(!entrezID == "") %>%
  group_by(entrezID) %>%
  summarise(longest_transcript_length = max(tx_length, na.rm = TRUE))

# Calculate gene lengths in kilobases
transcript_lengths$transcript_length_kb = transcript_lengths$longest_transcript_length / 1000

# Each row may contain multiple entrezIDs --> split them into separate rows --> find the longest transcript length for duplicated entrezIDs
temp = stringr::str_split_fixed(transcript_lengths$entrezID, ",", n = Inf)
transcript_lengths = cbind(transcript_lengths[, 2], temp)
transcript_lengths = reshape2::melt(transcript_lengths, "longest_transcript_length", colnames(transcript_lengths)[2:ncol(transcript_lengths)])
transcript_lengths = transcript_lengths %>% dplyr::select(entrezID = value, longest_transcript_length) %>% distinct()
transcript_lengths = transcript_lengths[-which(transcript_lengths$entrezID == ""), ] ; rownames(transcript_lengths) = NULL
transcript_lengths = transcript_lengths %>% 
  group_by(entrezID) %>%
  summarise(longest_transcript_length = max(longest_transcript_length, na.rm = TRUE))

# Map gene lengths to rownames of the counts matrix
transcript_lengths = as.data.frame(transcript_lengths)
rownames(transcript_lengths) = transcript_lengths$entrezID
transcript_lengths_in_raw_counts = transcript_lengths %>% dplyr::filter(entrezID %in% GSE194331_raw_counts$GeneID)
transcript_lengths_in_raw_counts$entrezID = as.integer(transcript_lengths_in_raw_counts$entrezID)
transcript_lengths_in_raw_counts = transcript_lengths_in_raw_counts %>% arrange(entrezID)
GSE194331_raw_counts = GSE194331_raw_counts %>% dplyr::filter(GeneID %in% transcript_lengths_in_raw_counts$entrezID)
GSE194331_raw_counts = GSE194331_raw_counts %>% dplyr::arrange(GeneID)

# Calculate counts per kilobase (CPK)
GSE194331_cpk = data.frame(GSE194331_raw_counts)
for (i in 1:nrow(GSE194331_cpk)) {
  gene_temp = GSE194331_cpk[i, "GeneID"]
  transcript_length = transcript_lengths_in_raw_counts %>% dplyr::filter(entrezID == gene_temp)
  GSE194331_cpk[i, 2:ncol(GSE194331_cpk)] = GSE194331_cpk[i, 2:ncol(GSE194331_cpk)] / transcript_length$longest_transcript_length
  cat(i, "\n")
} ; rm(i, gene_temp, transcript_length)
# Calculate the sum of CPK values for each sample
GSE194331_cpk_sum = colSums(GSE194331_cpk)

# Calculate TPM values
GSE194331_tpm = GSE194331_cpk / GSE194331_cpk_sum * 1e6
GSE194331_tpm$GeneID = GSE194331_raw_counts$GeneID

rm(GSE194331_cpk, GSE194331_cpk_sum,  ah, edb, temp, transcript_lengths, transcript_lengths_in_raw_counts, transcript_length_kb, transcripts_data, genes_data, genes_map, entrezid_temp)

## GSE133684_tpm: convert ensembl gene IDs to Entrez Gene IDs

## RefSeq to EntrezID reference data frame
ref = org.Hs.egENSEMBL2EG
mapped_genes_official = mappedkeys(ref)
ref_df = as.data.frame(ref[mapped_genes_official]) ; rm(ref, mapped_genes_official)
ref_df = ref_df %>% dplyr::rename(ENTREZ_GENE_ID = gene_id)
length(unique(ref_df$ensembl_id)) == length(ref_df$ensembl_id) # FALSE --> duplicates in RefSeq --> remove them
# remove ensembl_ids that match to more than one entrezID
duplicated_ensembl = unique(ref_df[which(duplicated(ref_df$ensembl_id)), "ensembl_id"])
ref_df = ref_df %>% filter(!ensembl_id %in% duplicated_ensembl) ; rm(duplicated_ensembl)
GSE133684_tpm = left_join(GSE133684_tpm, ref_df, by = c("V1" = "ensembl_id"))
GSE133684_tpm = GSE133684_tpm %>% dplyr::relocate(ENTREZ_GENE_ID, .before = V1)

