library(Seurat)
library(knitr)
library(SingleR)
library(dplyr)
library(parallel)
library(data.table)

setwd("./sc_analysis/")
source("Resources/functions_sc_analysis.R")

# Load most updated 5 prime Seurat Object
cat(paste0(">> Loading Seurat Object\n\n"))
load("./objects/seurat_obj.csfpb.Pairs.UMIregress.TCR.RData")

# Choose Reference
cat(paste0(">> Loading Cell References\n\n"))
# Monaco Reference
monaco            <- MonacoImmuneData()
blueprint_encode  <- BlueprintEncodeData()

run_bencode   <- T
run_monaco    <- T

# Run SingleR --------------------
singler_results_dir <- "./objects/pairs_singler_results/"
test <- as.SingleCellExperiment(seurat.obj.csfpbmc.pairs)

### BlueprintENCODE reference ---------------
if (run_bencode) {
  cat(paste0(">> Running SingleR: BlueprintENCODE\n\n"))
  immun_ref <- blueprint_encode
  singler.obj.csfpbmc.bencodeRef.sub <- SingleR(test,ref = immun_ref,labels = immun_ref$label.fine,clusters = seurat.obj.csfpbmc.pairs$seurat_clusters)
  save(singler.obj.csfpbmc.bencodeRef.sub, file = paste0(singler_results_dir,"singler_obj.csfpb.BlueprintEncodeRef.SUB.5prime.RData"))
  singler.obj.csfpbmc.bencodeRef.main <- SingleR(test,ref = immun_ref,labels = immun_ref$label.main,clusters = seurat.obj.csfpbmc.pairs$seurat_clusters)
  save(singler.obj.csfpbmc.bencodeRef.main, file = paste0(singler_results_dir,"singler_obj.csfpb.BlueprintEncodeRef.MAIN.5prime.RData"))
}


### Monaco Immune reference ------
if (run_monaco) {
  cat(paste0(">> Running SingleR: Monaco\n\n"))
  immun_ref <- monaco
  singler.obj.csfpbmc.monacoRef.sub <- SingleR(test,ref = immun_ref,labels = immun_ref$label.fine,clusters = seurat.obj.csfpbmc.pairs$seurat_clusters)
  save(singler.obj.csfpbmc.monacoRef.sub, file = paste0(singler_results_dir,"singler_obj.csfpb.MonacoImmuneRef.5prime.RData"))
  singler.obj.csfpbmc.monacoRef.main <- SingleR(test,ref = immun_ref,labels = immun_ref$label.main,clusters = seurat.obj.csfpbmc.pairs$seurat_clusters)
  save(singler.obj.csfpbmc.monacoRef.main, file = paste0(singler_results_dir,"singler_obj.csfpb.MonacoImmuneRef.MAIN.5prime.RData"))
}
