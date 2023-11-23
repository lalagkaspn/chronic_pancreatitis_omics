## Preparing list of genes for training model

library(openxlsx)
library(dplyr)

## load RNA-seq differentially expressed genes
rnaseq_dgea_cp_vs_normal = read.xlsx("DGEA/GSE_RNAseq/DGEA_results.xlsx", sheet = 2)
rnaseq_dgea_cp_vs_tumor = read.xlsx("DGEA/GSE_RNAseq/DGEA_results.xlsx", sheet = 3)
rnaseq_dgea_tumor_vs_normal = read.xlsx("DGEA/GSE_RNAseq/DGEA_results.xlsx", sheet = 4)

## filter for significantly DEG
rnaseq_dgea_cp_vs_normal_sig = rnaseq_dgea_cp_vs_normal %>% filter(adj.P.Val < 0.05)
rnaseq_dgea_cp_vs_tumor = rnaseq_dgea_cp_vs_tumor %>% filter(adj.P.Val < 0.05)
rnaseq_dgea_tumor_vs_normal = rnaseq_dgea_tumor_vs_normal %>% filter(adj.P.Val < 0.05)

## keep genes for the training model
## - unique DEG genes CP vs normal
## - unique DEG genes tumor vs normal
## - unique DEG genes CP vs tumor
## - common DEG genes between CP vs normal and tumor vs normal with opposite direction in CP and tumor

## unique DEG genes CP vs normal
unique_deg_cp_vs_normal = rnaseq_dgea_cp_vs_normal_sig %>% 
  filter(!EntrezGene.ID %in% c(rnaseq_dgea_cp_vs_tumor$EntrezGene.ID, rnaseq_dgea_tumor_vs_normal$EntrezGene.ID))

## unique DEG genes tumor vs normal
unique_deg_tumor_vs_normal = rnaseq_dgea_tumor_vs_normal %>%
  filter(!EntrezGene.ID %in% c(rnaseq_dgea_cp_vs_normal_sig$EntrezGene.ID, rnaseq_dgea_cp_vs_tumor$EntrezGene.ID))

## - unique DEG genes CP vs tumor
unique_deg_cp_vs_tumor = rnaseq_dgea_cp_vs_tumor %>%
  filter(!EntrezGene.ID %in% c(rnaseq_dgea_cp_vs_normal_sig$EntrezGene.ID, rnaseq_dgea_tumor_vs_normal$EntrezGene.ID))

## - common DEG genes between CP vs normal and tumor vs normal with opposite direction in CP and tumor
common_deg_cpVSnormal_tumorVSnormal = inner_join(rnaseq_dgea_cp_vs_normal_sig, rnaseq_dgea_tumor_vs_normal, by = c("EntrezGene.ID", "Gene.Symbol", "HGNC_Official")) %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol,
                logFC_CPvsNormal = logFC.x, AveExpr_CPvsNormal = AveExpr.x, t_CPvsNormal = t.x, P.Value_CPvsNormal = P.Value.x, adj.P.Val_CPvsNormal = adj.P.Val.x, B_CPvsNormal = B.x,
                logFC_TumorvsNormal = logFC.y, AveExpr_TumorvsNormal = AveExpr.y, t_TumorvsNormal = t.y, P.Value_TumorvsNormal = P.Value.y, adj.P.Val_TumorvsNormal = adj.P.Val.y, B_TumorvsNormal = B.y,
                HGNC_Official)
# These are 1,012 genes:
# - 321 are common between CPvsNormal and TumorVSnormal only
# - 691 are commono between CPvsNormal, TumorVSnormal and CPvsTumor
# We keep only the 321 and we further filter them for those that have opposite direction in CPvsNormal vs TumorVSnormal
common_deg_cpVSnormal_tumorVSnormal = common_deg_cpVSnormal_tumorVSnormal %>% 
  filter(!EntrezGene.ID %in% intersect(rnaseq_dgea_cp_vs_normal_sig$EntrezGene.ID, rnaseq_dgea_cp_vs_tumor$EntrezGene.ID)) %>%
  mutate(logFC_combined = logFC_CPvsNormal * logFC_TumorvsNormal) %>%
  filter(logFC_combined < 0) # DEG with opposite logFC direction
nrow(common_deg_cpVSnormal_tumorVSnormal) # all the genes had the same logFC diretction!

## combine genes
genes_for_model_training = data.frame(EntrezGene.ID = c(unique_deg_cp_vs_normal$EntrezGene.ID, 
                                                        unique_deg_tumor_vs_normal$EntrezGene.ID,
                                                        unique_deg_cp_vs_tumor$EntrezGene.ID,
                                                        common_deg_cpVSnormal_tumorVSnormal$EntrezGene.ID),
                                      Gene.Symbol = c(unique_deg_cp_vs_normal$Gene.Symbol, 
                                                      unique_deg_tumor_vs_normal$Gene.Symbol,
                                                      unique_deg_cp_vs_tumor$Gene.Symbol,
                                                      common_deg_cpVSnormal_tumorVSnormal$Gene.Symbol),
                                      HGNC_Official = c(unique_deg_cp_vs_normal$HGNC_Official, 
                                                        unique_deg_tumor_vs_normal$HGNC_Official,
                                                        unique_deg_cp_vs_tumor$HGNC_Official,
                                                        common_deg_cpVSnormal_tumorVSnormal$HGNC_Official))
genes_for_model_training = unique(genes_for_model_training) # 3,453 genes

## filter those genes for the ones also tested in microarrays
# load microarrays
microarray_dgea_cp_vs_normal = read.xlsx("DGEA/GSE_microarrays/DGEA_results.xlsx", sheet = 2) # I only need one as the same genes were tested in each comparison
# filtering the RNAseq genes for model training for those also tested in microarrays
genes_for_model_training = genes_for_model_training %>%
  filter(EntrezGene.ID %in% microarray_dgea_cp_vs_normal$EntrezGene.ID) # 2,340 / 3,453 genes remain

## save
write.xlsx(genes_for_model_training, "genes_for_model_training.xlsx", overwrite = TRUE)
  
