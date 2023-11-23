
## -- Installation of CRAN and BioConductor R packages -- ##

# dplyr
install.packages("dplyr")

# matrixStats
install.packages("matrixStats")

# reshape2
install.packages("reshape2")

# openxlsx
install.packages("openxlsx")

# GEOquery
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

# impute
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")

# org.Hs.eg.db
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

# limma
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
  
# ggplot2
install.packages("ggplot2")

# pheatmap
install.packages("pheatmap")

# EnhancedVolcano
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")

