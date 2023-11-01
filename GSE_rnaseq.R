## This script is used to download and pre-process series (RNA-seq) from GEO for 
## patients with Chronic Pancreatitis (normal and CP tissues)

library(dplyr)
library(GEOquery) # updated
library(org.Hs.eg.db)
library(matrixStats)
library(limma)
library(ggplot2)
library(reshape2)
library(pheatmap)

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

# Change Tissue_type to chronic pancreatitis, tumor and normal
table(filt_pdata$GSE194331$Tissue_type)
filt_pdata$GSE194331$Tissue_type = gsub("Healthy control", "normal", filt_pdata$GSE194331$Tissue_type)
filt_pdata$GSE194331$Tissue_type = gsub("Mild AP", "chronic_pancreatitis", filt_pdata$GSE194331$Tissue_type)
filt_pdata$GSE194331$Tissue_type = gsub("Moderately-severe AP", "chronic_pancreatitis",  filt_pdata$GSE194331$Tissue_type)
filt_pdata$GSE194331$Tissue_type = gsub("Severe AP", "chronic_pancreatitis", filt_pdata$GSE194331$Tissue_type)
table(filt_pdata$GSE194331$Tissue_type)

# Transform to factors with consistent universal levels
filt_pdata$GSE194331$Tissue_type = factor(x = filt_pdata$GSE194331$Tissue_type,
                                          levels = c("chronic_pancreatitis", "normal", "tumor"),
                                          labels = c("chronic_pancreatitis", "normal", "tumor"))

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
table(filt_pdata$GSE133684$Tissue_type)
filt_pdata$GSE133684$Tissue_type = gsub("PDAC", "tumor", filt_pdata$GSE133684$Tissue_type)
filt_pdata$GSE133684$Tissue_type = gsub("healthy", "normal", filt_pdata$GSE133684$Tissue_type)
table(filt_pdata$GSE133684$Tissue_type)

# Add to filt_pdata the CP patient IDs
cp_patients_ids = setdiff(colnames(GSE133684_tpm), filt_pdata$GSE133684$Patient_ID)[-1]
cp_data = data.frame(GEO_accession = "not_available", Patient_ID = cp_patients_ids, Platform = "GPL20795", Tissue_type = "chronic_pancreatitis")
filt_pdata$GSE133684 = rbind(filt_pdata$GSE133684, cp_data) ; rm(cp_patients_ids, cp_data)

# Transform to factors with consistent universal levels
filt_pdata$GSE133684$Tissue_type = factor(x = filt_pdata$GSE133684$Tissue_type,
                                          levels = c("chronic_pancreatitis", "normal", "tumor"),
                                          labels = c("chronic_pancreatitis", "normal", "tumor"))

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
transcript_lengths = reshape2::melt(transcript_lengths, 
                                    "longest_transcript_length",
                                    colnames(transcript_lengths)[2:ncol(transcript_lengths)])
transcript_lengths = transcript_lengths %>% dplyr::select(entrezID = value,
                                                          longest_transcript_length) %>% distinct()
transcript_lengths = transcript_lengths[-which(transcript_lengths$entrezID == ""), ]
rownames(transcript_lengths) = NULL
transcript_lengths = transcript_lengths %>% 
  group_by(entrezID) %>%
  summarise(longest_transcript_length = max(longest_transcript_length, na.rm = TRUE))

# Map gene lengths to rownames of the counts matrix
transcript_lengths = as.data.frame(transcript_lengths)
rownames(transcript_lengths) = transcript_lengths$entrezID
transcript_lengths_in_raw_counts = transcript_lengths %>% 
  dplyr::filter(entrezID %in% GSE194331_raw_counts$GeneID)
transcript_lengths_in_raw_counts$entrezID = as.integer(transcript_lengths_in_raw_counts$entrezID)
transcript_lengths_in_raw_counts = transcript_lengths_in_raw_counts %>% arrange(entrezID)
GSE194331_raw_counts = GSE194331_raw_counts %>% 
  dplyr::filter(GeneID %in% transcript_lengths_in_raw_counts$entrezID)
GSE194331_raw_counts = GSE194331_raw_counts %>% dplyr::arrange(GeneID)

sum(transcript_lengths_in_raw_counts$entrezID == 
      GSE194331_raw_counts$GeneID) == 
  nrow(GSE194331_raw_counts) # TRUE --> entrezIDs are in the same order in both files

# Calculate RPK
# Divide the counts by transcript length in kB
GSE194331_rpk = apply(subset(GSE194331_raw_counts, select = c(-GeneID)), 2,
                      function(x) {
                          x / (transcript_lengths_in_raw_counts$longest_transcript_length / 1000)
                        }
                      )
GSE194331_rpk = as.data.frame(GSE194331_rpk)

# Calculate TPM
GSE194331_tpm = apply(GSE194331_rpk, 2,
                      function (x) x / (sum(as.numeric(x)) / 10^6))
GSE194331_tpm = as.data.frame(GSE194331_tpm)
GSE194331_tpm$GeneID = GSE194331_raw_counts$GeneID
GSE194331_tpm = GSE194331_tpm %>% 
  dplyr::select(ENTREZ_GENE_ID = GeneID, everything())

rm(GSE194331_rpk, ah, edb, temp, transcript_lengths, 
   transcript_lengths_in_raw_counts, transcripts_data, genes_data, genes_map, 
   entrezid_temp, latest_homo_sapiens)

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

# log2(TPM + 1) transformation --> +1 is used to avoid -inf values for genes with TPM=0, as log2(0) = -Inf
for (i in 1:length(esets)) {
  cols_to_transform = setdiff(colnames(esets[[i]]), "ENTREZ_GENE_ID")
  esets[[i]][cols_to_transform] = lapply(esets[[i]][cols_to_transform],
                                         function(x) log2(x + 1))
}

original_exprs = esets[[1]] %>% 
  inner_join(esets[[2]], by = "ENTREZ_GENE_ID") %>%
  dplyr::select(ENTREZ_GENE_ID, everything())

rows = original_exprs$ENTREZ_GENE_ID
original_exprs = as.matrix(original_exprs %>% dplyr::select(-ENTREZ_GENE_ID))
rownames(original_exprs) = rows; rm(rows) # 23722 x 620
# Making sure we do not have NAs in any row
original_exprs = original_exprs[rowSums(is.na(original_exprs)) != ncol(original_exprs), ]
original_exprs_nonas = na.omit(original_exprs) # 23722 x 620
# Keeping all samples in our full_pdata object
original_exprs_nonas = original_exprs_nonas[, full_pdata$GEO_accession] # 23722 x 620

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

# Joining in one expression matrix: z-score normalized version
## Separate runs because the 2nd study does not have GSM sample IDs for the chronic_pancreatitis samples --> i use the Patient_ID column for the intersect instead of the GSM
z[[1]]$EntrezGene.ID = as.character(z[[1]]$EntrezGene.ID)
z[[1]] = z[[1]][,c("EntrezGene.ID", 
                   intersect(colnames(z[[1]]),
                             filt_pdata[[1]]$GEO_accession))]
z[[2]]$EntrezGene.ID = as.character(z[[2]]$EntrezGene.ID)
z[[2]] = z[[2]][,c("EntrezGene.ID", 
                   intersect(colnames(z[[2]]),
                             filt_pdata[[2]]$Patient_ID))]
z_exprs = z[[1]] %>% inner_join(z[[2]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

rownames(z_exprs) = z_exprs$EntrezGene.ID
z_exprs = as.matrix(z_exprs %>% dplyr::select(-EntrezGene.ID)) # 23,722 x 620
# Making sure we do not have NAs in any row
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
z_exprs_nonas = na.omit(z_exprs) # 20,550 x 620
z_exprs_nonas = z_exprs_nonas[, full_pdata$GEO_accession] # 20,550 x 620

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
       path = "QC/GSE_RNAseq", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# MDS: normalized matrix
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
       path = "QC/GSE_RNAseq", 
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

tiff("QC/GSE_RNAseq/MDS_multiplot.tiff", 
     width = 1920, height = 2160, res = 700, compression = "lzw", type = type_compression)
multiplot(original_MDS, KBZ_MDS_plot, cols = 1)
m = ggplot(multiplot(original_MDS, KBZ_MDS_plot, cols = 1))
dev.off(); rm(m)

# Global expression boxplot: original matrix
original_eset = as.data.frame(original_exprs_nonas)
original_boxplot = ggplot(melt(original_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.01, outlier.shape = 20, linewidth = 0.02,
               fill = c(rep("cyan", 119),
                        rep("chartreuse", 501))) +
  scale_y_continuous("Expression",
                     limits = c(-1,
                                round(max(reshape2::melt(original_eset)$value, na.rm = TRUE)+1)), 
                     breaks = seq(-1,
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
  labs(title = bquote("Boxplot of expression:" ~ bold(log[2](TPM + "1") ~ .("values"))),
       x = "\nSamples",
       y = bquote("Expression:" ~ bold(log[2](TPM + "1"))))

original_boxplot
ggsave(filename = "Original_boxplot.tiff",
       path = "QC/GSE_RNAseq/", 
       width = 7680, height = 3240, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# Global expression boxplot: z-score normalised matrix
z_eset = as.data.frame(z_exprs_nonas)
KBZ_boxplot = ggplot(melt(z_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.01, outlier.shape = 20, linewidth = 0.02,
               fill = c(rep("cyan", 119),
                        rep("chartreuse", 501))) +
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
  labs(title = bquote("Boxplot of expression: normalized" ~ bold(log[2](TPM + "1") ~ .("values"))),
       x = "\nSamples",
       y = bquote("Expression: normalized" ~ bold(log[2](TPM + "1"))))

KBZ_boxplot
ggsave(filename = "z_boxplot.tiff",
       path = "QC/GSE_RNAseq/", 
       width = 7680, height = 3240, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

tiff("QC/GSE_RNAseq/Boxplot_multiplot.tiff", 
     width = 7680, height = 6480, res = 700, compression = "lzw", type = type_compression)
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
  tiff(filename, width = width, height = height, res = res, compression = "lzw", type = type_compression)
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
  Tissue_type = c(chronic_pancreatitis = "deeppink4", normal = "dodgerblue4", tumor = "gray"),
  Study = c(GSE194331 = "darkseagreen", GSE133684 = "darkorange")
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
save_pheatmap_tiff(original_heatmap, "QC/GSE_RNAseq/original_heatmap.tiff")

# Z-score version
annotation_for_heatmap = full_pdata[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata$GEO_accession

z_dists = as.matrix(dist(t(z_exprs_nonas), method = "manhattan"))

rownames(z_dists) = full_pdata$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

ann_colors <- list(
  Tissue_type = c(chronic_pancreatitis = "deeppink4", normal = "dodgerblue4", tumor = "gray"),
  Study = c(GSE194331 = "darkseagreen", GSE133684 = "darkorange")
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
save_pheatmap_tiff(z_heatmap, "QC/GSE_RNAseq/KBZ_heatmap.tiff")

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
design = model.matrix(~0 + full_pdata$Tissue_type + full_pdata$Study)
colnames(design) = c("chronic_pancreatitis", "non_tumor", "tumor", 
                     "GSE194331") 
rownames(design) = colnames(z_exprs_nonas)
cont.matrix = makeContrasts(CPvsNormal = chronic_pancreatitis - non_tumor,
                            CPvsTumor = chronic_pancreatitis - tumor,
                            TumorvsNormal = tumor - non_tumor, 
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
  z_DE_mapped = z_DE %>% left_join(ID_Map, by = "EntrezGene.ID")
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

saveWorkbook(DGEA_topTables, "DGEA/GSE_RNAseq/DGEA_results.xlsx", overwrite = TRUE)
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
       path = "DGEA/GSE_RNAseq/", 
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
       path = "DGEA/GSE_RNAseq/", 
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
       path = "DGEA/GSE_RNAseq/", 
       width = 100, height = 142, device = 'tiff', units = "mm",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# ---------------------------------------- #
# ---------------------------------------- #
# ---------------------------------------- #