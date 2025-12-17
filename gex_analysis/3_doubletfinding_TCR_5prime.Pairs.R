# Goal: Run doublet finding

working_dir <- "./sc_analysis/"
setwd(working_dir)
source("Resources/functions_sc_analysis.R")

library(Seurat)
library(knitr)
library(devtools)
library(dplyr)
library(parallel)
library(data.table)

load("./objects/seurat_obj.csfpb.Pairs.UMIregress.TCR.RData")

# Add meta data
tcr_metadata <- "./metadata/metadata.csv"
seurat.obj.csfpbmc.pairs <- addRNAseqMetaData(seurat.obj.csfpbmc.pairs,metadata_sheet = tcr_metadata)

# Doublet Finding
seurat.obj.csfpbmc.pairs <- find_doublets(seurat.obj.csfpbmc.pairs)

# SAVE
fh <- "./objects/seurat_obj.csfpb.Pairs.UMIregress.TCR.postDF.RData"
save(seurat.obj.csfpbmc.pairs, file = fh)

cat(paste0("\n\n>>> DONE <<<\n\n"))

