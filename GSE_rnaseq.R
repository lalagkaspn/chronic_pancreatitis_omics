## This script is used to download and pre-process series (RNA-seq) from GEO for 
## patients with Chronic Pancreatitis (normal and CP tissues)

library(dplyr)
library(GEOquery) # updated
library(org.Hs.eg.db)
library(matrixStats)
library(limma)

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

pdata194331 = filt_pdata$GSE194331 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE194331")
pdata133684 = filt_pdata$GSE133684 %>%
  dplyr::select(GEO_accession = Patient_ID, Tissue_type) %>%
  dplyr::mutate(Study = "GSE133684")

full_pdata = rbind(pdata194331, pdata133684)
openxlsx::write.xlsx(full_pdata, "DGEA/Pheno_RNAseq.xlsx", overwrite = TRUE)

## ---------------------------
## Expression data
## ---------------------------

## GSE194331: convert raw counts to TPM

# Summary of the counts table
summary(GSE194331_raw_counts[, 2:5]) # --> needs normalization --> will calculate TPM

## To calculate TPM:
## 1. Calculate RPK (divide counts with transcript lengths in kB) --> reads per kilobase
## 2. Find scaling factor for each sample by summing all the RPK of each sample and divide by 1,000,000
## 3. Multiply each value within each sample with that scaling factor --> transcript per million (TPM)

## -- Add transcript length -- ##
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

sum(transcript_lengths_in_raw_counts$entrezID == GSE194331_raw_counts$GeneID) == nrow(GSE194331_raw_counts) # TRUE --> entrezIDs are in the same order in both files

# Calculate RPK
# Divide the counts by transcript length in kB
GSE194331_rpk = apply(subset(GSE194331_raw_counts, select = c(-GeneID)), 2,
                      function(x) x / (transcript_lengths_in_raw_counts$longest_transcript_length / 1000))
GSE194331_rpk = as.data.frame(GSE194331_rpk)

# Calculate TPM
GSE194331_tpm = apply(GSE194331_rpk, 2,
                      function (x) x / sum(as.numeric(x)) * 10^6)
GSE194331_tpm = as.data.frame(GSE194331_tpm)
GSE194331_tpm$GeneID = GSE194331_raw_counts$GeneID ; GSE194331_tpm = GSE194331_tpm %>% dplyr::relocate(GeneID, .before = GSM5833563) %>% dplyr::rename(ENTREZ_GENE_ID = GeneID)

rm(GSE194331_rpk, ah, edb, temp, transcript_lengths, transcript_lengths_in_raw_counts, transcripts_data, genes_data, genes_map, entrezid_temp, latest_homo_sapiens)

## GSE133684_tpm: convert ensembl gene IDs to Entrez Gene IDs

## RefSeq to EntrezID reference data frame
ref = org.Hs.egENSEMBL2EG
mapped_genes_official = mappedkeys(ref)
ref_df = as.data.frame(ref[mapped_genes_official]) ; rm(ref, mapped_genes_official)
ref_df = ref_df %>% dplyr::rename(ENTREZ_GENE_ID = gene_id)
length(unique(ref_df$ensembl_id)) == length(ref_df$ensembl_id) # FALSE --> duplicates in RefSeq --> remove them
# remove ensembl_ids that match to more than one entrezID
duplicated_ensembl = unique(ref_df[which(duplicated(ref_df$ensembl_id)), "ensembl_id"])
ref_df = ref_df %>% dplyr::filter(!ensembl_id %in% duplicated_ensembl) ; rm(duplicated_ensembl)
GSE133684_tpm = left_join(GSE133684_tpm, ref_df, by = c("V1" = "ensembl_id"))
GSE133684_tpm = GSE133684_tpm %>% dplyr::relocate(ENTREZ_GENE_ID, .before = V1)
GSE133684_tpm$ENTREZ_GENE_ID = as.integer(GSE133684_tpm$ENTREZ_GENE_ID)
# remove NAs, ensemble gene ID column and arrange by entrezID
GSE133684_tpm_nonas = GSE133684_tpm %>% 
  na.omit() %>%
  as.data.frame() %>%
  dplyr::select(-V1) %>% 
  arrange(ENTREZ_GENE_ID)

# there are duplicated EnterzIDs --> keep the rows with max variance
GSE133684_tpm_nonas$variance = rowVars(as.matrix(GSE133684_tpm_nonas))
GSE133684_tpm_nonas = GSE133684_tpm_nonas %>%
  group_by(ENTREZ_GENE_ID) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  dplyr::select(-variance)

# Create list with expression data of both studies
esets = list()
esets[["GSE194331"]] = GSE194331_tpm
esets[["GSE133684"]] = GSE133684_tpm_nonas

##### Quality Control #####

# Joining in one expression matrix: original version
esets[[1]]$ENTREZ_GENE_ID = as.character(esets[[1]]$ENTREZ_GENE_ID)
esets[[1]] = esets[[1]][,c("ENTREZ_GENE_ID", 
                           intersect(colnames(esets[[1]]),
                                     filt_pdata[[1]]$GEO_accession))]
esets[[2]]$ENTREZ_GENE_ID = as.character(esets[[2]]$ENTREZ_GENE_ID)
esets[[2]] = esets[[2]][,c("ENTREZ_GENE_ID", 
                           intersect(colnames(esets[[2]]),
                                     filt_pdata[[2]]$Patient_ID))]

original_exprs = esets[[1]] %>% 
  inner_join(esets[[2]], by = "ENTREZ_GENE_ID") %>%
  dplyr::select(ENTREZ_GENE_ID, everything())

rows = original_exprs$ENTREZ_GENE_ID
original_exprs = as.matrix(original_exprs %>% dplyr::select(-ENTREZ_GENE_ID))
rownames(original_exprs) = rows; rm(rows) # 23722 x 620
# Making sure we do not have NAs in any row
original_exprs = original_exprs[rowSums(is.na(original_exprs)) != ncol(original_exprs), ]
original_exprs_nonas = na.omit(original_exprs) # 23722 x 620
# Keeping all samples in our full_pdata_filt object
original_exprs_nonas = original_exprs_nonas[, full_pdata$GEO_accession] # 23722 x 620

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
tiff("QC/GSE_RNAseq/Original_MDS.tif", width = 1920, height = 1080, res = 100)
original_MDS
dev.off()

# Global expression boxplot: original matrix
original_eset = as.data.frame(original_exprs_nonas)
original_boxplot = ggplot(melt(original_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20, outlier.alpha = 0.1,
               fill = c(rep("cyan", 119),
                        rep("chartreuse", 501))) +
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
  labs(title = "Boxplot of TPMs",
       x = "\nSamples",
       y = "TPM\n")
tiff("QC/GSE_RNAseq/Original_boxplot.tif", width = 1920, height = 1080, res = 100)
original_boxplot
dev.off()

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
  Study = c(GSE194331 = "darkseagreen", GSE133684 = "darkorange")
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
save_pheatmap_tiff(original_heatmap, "QC/GSE_RNAseq/original_heatmap.tif")
