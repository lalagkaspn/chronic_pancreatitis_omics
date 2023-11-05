### Compare the two list of differentially expressed genes obtained by analyzing microarray and RNA-seq data

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

## load libraries
library(dplyr)
library(openxlsx)
library(ggplot2)

## load data
cp_vs_normal_microarrays = read.xlsx("DGEA/GSE_microarrays/DGEA_results.xlsx", sheet = 2)
cp_vs_normal_rnaseq = read.xlsx("DGEA/GSE_RNAseq/DGEA_results.xlsx", sheet = 2)

## description of the data
nr_genes_microarrays = nrow(cp_vs_normal_microarrays)
nr_genes_rnaseq = nrow(cp_vs_normal_rnaseq)
nr_genes_microarrays
nr_genes_rnaseq

## find the overlap of genes in both files
nr_genes_overlap = length(intersect(cp_vs_normal_microarrays$EntrezGene.ID,
                                    cp_vs_normal_rnaseq$EntrezGene.ID))
nr_genes_overlap

## barplot of genes overlap
# create data frame with numbers to be visualized
barplot_geneoverlap = data.frame(method = c("Microarrays", "RNA-seq", "Overlap"),
                                 values = c(nr_genes_microarrays, nr_genes_rnaseq, nr_genes_overlap))
# convert "method" to be a factor --> this will help us visualize the data in the row they appear in the barplot_geneoverlap
barplot_geneoverlap$method = factor(barplot_geneoverlap$method, levels = barplot_geneoverlap$method, labels = barplot_geneoverlap$method)
# create barplot
ggplot(barplot_geneoverlap, aes(x = method, y = values)) +
  geom_col() +
  labs(title = "Genes tested using microarray or RNA-seq data") +
  xlab("") + # empty title for x axis 
  ylab("Number of genes") + # name for y axis
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5), # bold and center title
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, # change font size and position of y axis text
                                   margin = margin(t = 1, unit = "cm"),
                                   size = 14),
        axis.text.x = element_text(face = "bold", color = "black", # change font size and position of x axis text
                                   angle = 0, hjust = 0.5, vjust = 0.5, 
                                   size = 14),
        axis.title = element_text(angle = 0, hjust = 0.5, face = "bold", # change font size and location of axes titles
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 18))
# save plot
ggsave(filename = "Number_of_genes_overlap.tiff",
       path = "DGEA/results_comparison/", 
       width = 2500, height = 1080, device = 'tiff', units = "px",
       dpi = 300, compression = "lzw", type = type_compression)

## create a barplot that shows the percentages of overlapping genes in each method rather than the actual number of genes
# data
barplot_geneoverlap_percentages = data.frame(method = c("Microarrays", "RNA-seq"),
                                             values = c((nr_genes_overlap / nr_genes_microarrays) * 100, # calculate percentages
                                                        (nr_genes_overlap / nr_genes_rnaseq) * 100))
# create barplot
ggplot(barplot_geneoverlap_percentages, aes(x = method, y = values)) +
  geom_col() +
  labs(title = "Percentage of genes overlapping \nin microarray and RNA-seq data") +
  xlab("") + # empty title for x axis 
  ylab("Percentage of genes") + # name for y axis
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) + # set value breaks and limits for y axis
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5), # bold and center title
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, # change font size and position of y axis text
                                   margin = margin(t = 1, unit = "cm"),
                                   size = 14),
        axis.text.x = element_text(face = "bold", color = "black", # change font size and position of x axis text
                                   angle = 0, hjust = 0.5, vjust = 0.5, 
                                   size = 14),
        axis.title = element_text(angle = 0, hjust = 0.5, face = "bold", # change font size and location of axes titles
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 18), aspect.ratio = 1.3)
# save plot
ggsave(filename = "Number_of_genes_overlap_percentages.tiff",
       path = "DGEA/results_comparison/", 
       width = 2500, height = 1600, device = 'tiff', units = "px",
       dpi = 300, compression = "lzw", type = type_compression)

## filter for overlapping genes
overlapping_genes = intersect(cp_vs_normal_microarrays$EntrezGene.ID, cp_vs_normal_rnaseq$EntrezGene.ID)
cp_vs_normal_microarrays_common_genes = cp_vs_normal_microarrays %>% filter(EntrezGene.ID %in% overlapping_genes)
cp_vs_normal_rnaseq_common_genes = cp_vs_normal_rnaseq %>% filter(EntrezGene.ID %in% overlapping_genes)

## filter for significantly dysregulated genes
cp_vs_normal_microarrays_common_genes_sig = cp_vs_normal_microarrays_common_genes %>% filter(adj.P.Val < 0.05)
cp_vs_normal_rnaseq_common_genes_sig = cp_vs_normal_rnaseq_common_genes %>% filter(adj.P.Val < 0.05)

## number of significantly dysregulated genes
nr_dysregulated_genes_microarrays = nrow(cp_vs_normal_microarrays_common_genes_sig)
nr_dysregulated_genes_microarrays
nr_dysregulated_genes_rnaseq = nrow(cp_vs_normal_rnaseq_common_genes_sig)
nr_dysregulated_genes_rnaseq

## overlap of significantly dysregulated genes
overlap_sig_dys_genes = intersect(cp_vs_normal_microarrays_common_genes_sig$EntrezGene.ID, cp_vs_normal_rnaseq_common_genes_sig$EntrezGene.ID)
length(overlap_sig_dys_genes)

## check the significance of the overlap
# the background genes are the genes tested in both microarry and RNA-seq data (n=13,029)
# create 2x2 table
contingency_table = matrix(ncol = 2, nrow = 2)
colnames(contingency_table) = c("in RNA-seq", "not in RNA-seq")
rownames(contingency_table) = c("in microarray", "not in microarray")
# populate the 2x2 table
# genes found to be significantly dysregulated in both data
contingency_table[1, 1] = length(overlap_sig_dys_genes) 
# genes found to be significantly dysregulated in RNA-seq but not in microarray
contingency_table[2, 1] = length(setdiff(cp_vs_normal_rnaseq_common_genes_sig$EntrezGene.ID, cp_vs_normal_microarrays_common_genes_sig$EntrezGene.ID)) 
# genes found to be significantly dysregulated in microarray but not in RNA-seq
contingency_table[1, 2] = length(setdiff(cp_vs_normal_microarrays_common_genes_sig$EntrezGene.ID, cp_vs_normal_rnaseq_common_genes_sig$EntrezGene.ID)) 
# genes not found to be significantly dysregulated in any data
contingency_table[2, 2] = length(overlapping_genes) - contingency_table[1, 1] - contingency_table[2, 1] - contingency_table[1, 2]
# assess the significance of the overlap using the Fisher's exact test (hypergeometric test)
fisher.test(contingency_table, alternative = "greater")

## barplot to summarize the overlap of dysregulated genes in both data
# data
barplot_dys_geneoverlap = data.frame(method = c("Microarrays", "RNA-seq", "Overlap"),
                                     values = c(nr_dysregulated_genes_microarrays, nr_dysregulated_genes_rnaseq, length(overlap_sig_dys_genes)))
barplot_dys_geneoverlap$method = factor(barplot_dys_geneoverlap$method, levels = barplot_dys_geneoverlap$method, 
                                        labels = barplot_dys_geneoverlap$method)
# create barplot
ggplot(barplot_dys_geneoverlap, aes(x = method, y = values)) +
  geom_col() +
  labs(title = "Significantly differentially expressed genes \nin microarray and RNA-seq data") +
  xlab("") + # empty title for x axis 
  ylab("Number of genes") + # name for y axis
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5), # bold and center title
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, # change font size and position of y axis text
                                   margin = margin(t = 1, unit = "cm"),
                                   size = 14),
        axis.text.x = element_text(face = "bold", color = "black", # change font size and position of x axis text
                                   angle = 0, hjust = 0.5, vjust = 0.5, 
                                   size = 14),
        axis.title = element_text(angle = 0, hjust = 0.5, face = "bold", # change font size and location of axes titles
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 18), aspect.ratio = 1.3)
# save plot
ggsave(filename = "Number_of_sig_dysregulated_genes_overlap.tiff",
       path = "DGEA/results_comparison/", 
       width = 3200, height = 2200, device = 'tiff', units = "px",
       dpi = 300, compression = "lzw", type = type_compression)

## filter for overlapping significantly dysregualted genes
cp_vs_normal_microarrays_common_genes_sig_ovelap = cp_vs_normal_microarrays_common_genes_sig %>% filter(EntrezGene.ID %in% overlap_sig_dys_genes)
cp_vs_normal_rnaseq_common_genes_sig_overlap = cp_vs_normal_rnaseq_common_genes_sig %>% filter(EntrezGene.ID %in% overlap_sig_dys_genes)

## arrange the two data frames based on EntrezGene.ID
cp_vs_normal_microarrays_common_genes_sig_ovelap = cp_vs_normal_microarrays_common_genes_sig_ovelap %>% arrange(EntrezGene.ID)
cp_vs_normal_rnaseq_common_genes_sig_overlap = cp_vs_normal_rnaseq_common_genes_sig_overlap %>% arrange(EntrezGene.ID)

## create a data frame that contains the logFC from each data type
logfc = data.frame(EntrezGene.ID = cp_vs_normal_microarrays_common_genes_sig_ovelap$EntrezGene.ID, # first column contains the gene IDs
                   microarrays_logFC = cp_vs_normal_microarrays_common_genes_sig_ovelap$logFC, # second column contains the logFC from microarrays
                   rnaseq_logFC = cp_vs_normal_rnaseq_common_genes_sig_overlap$logFC) # third column contains the logFC from RNA-seq

## create scatterplot
ggplot(logfc, aes(x = microarrays_logFC, y = rnaseq_logFC)) +
  geom_point(size = 0.5) +
  labs(title = "logFC of overlapping significantly dysregulated genes \n in microarry and RNA-seq data") +
  xlab("Microarray (logFC)") + # empty title for x axis 
  ylab("RNA-seq (logFC)") + # name for y axis
  geom_vline(xintercept = 0) + # add vertical line
  geom_hline(yintercept = 0) + # add horizontal line
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), # bold and center title
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, # change font size and position of y axis text
                                   margin = margin(t = 1, unit = "cm"),
                                   size = 14),
        axis.text.x = element_text(face = "bold", color = "black", # change font size and position of x axis text
                                   angle = 0, hjust = 0.5, vjust = 0.5, 
                                   size = 10),
        axis.title = element_text(angle = 0, hjust = 0.5, face = "bold", # change font size and location of axes titles
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 12), aspect.ratio = 1.3)
# save plot
ggsave(filename = "logFC_of_sig_dysregulated_genes_overlap.tiff",
       path = "DGEA/results_comparison/", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 300, compression = "lzw", type = type_compression)

## Pearson correlation analysis for logFC in microarray and RNA-seq data
pearson_correlation = cor.test(logfc$microarrays_logFC, logfc$rnaseq_logFC, method = "pearson")
pearson_correlation
  
  

