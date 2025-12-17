# CREATE SEURAT OBJECT FOR TCR

# PURPOSE: Create a single cell RNA-Seq object with paired Bulk CSF & PBMC

# DEFAULT conditions for creating Seurat Object Currently:
# min.genes = 200,
# min.cells = 2,
# regress.out = 'nUMI',
# npca = 10,
# resolution=0.8,

setwd("./sc_analysis/")
source("Resources/functions_sc_analysis.R")

library(Seurat)
library(knitr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(SingleR)
library(kableExtra)
library(devtools)
library(dplyr)
library(parallel)
library(data.table)

### INPUT and OUTPUT Directories
dataset_loc_pairs <- "./gex_data/"

# RData Objects to save
fh_raw_seurat_obj_pairs <- "./objects/seurat_obj.csfpb.Pairs.UMIregress.TCR.RData"
MAKE_PAIRS <- TRUE

# Increase size of ENV
options(future.globals.maxSize= 891289600)


### 1. CREATE PAIRED OBJECT ----------------------------------------------------------------------------
# a. Create vector of patient IDs and vector of sample names
samples.vec.pairs <- dir(dataset_loc_pairs)
patients.vec      <- unique(unlist(lapply(samples.vec.pairs, function(x) strsplit(x,'_')[[1]][1])))

# b. create Seurat Object
if (MAKE_PAIRS) {
  cat(paste0(">>> CREATING PAIRS OBJECT\n\n"))
  dataset_loc <- dataset_loc_pairs
  samples.vec <- samples.vec.pairs
  patient.combined.matrix <- generatePatientCombinedMatrix(dataset_loc, samples.vec,THREADS = 10)
  
  # Only keep cells with total number of genes >400 and <2500
  # Testing mingenes = 700
  num_genes_per_cell <- (colSums(patient.combined.matrix !=0))
  target_bool <- (num_genes_per_cell > 700 & num_genes_per_cell < 2500)
  cells_to_keep <- target_bool[target_bool == T]
  patient.combined.matrix <- patient.combined.matrix[,colnames(patient.combined.matrix) %in% names(cells_to_keep)]
  
  seurat.obj.csfpbmc.pairs <- SingleR.CreateSeurat.fixed(
    sc.data      = patient.combined.matrix,
    project.name = "PAIRS_CSF_PBMC_5prime_results",
    npca         = 20, min.genes = 700
  )
  
  # Add sample.name to meta data with full sample names for each cell
  cell.id.vec <- names(seurat.obj.csfpbmc.pairs$orig.ident)
  sample.ident.vec <- sapply(strsplit(cell.id.vec,'-'), '[',2)
  seurat.obj.csfpbmc.pairs$sample.name <- sample.ident.vec
  names(seurat.obj.csfpbmc.pairs$sample.name) <- cell.id.vec
  
  save(seurat.obj.csfpbmc.pairs, file = fh_raw_seurat_obj_pairs)
  rm(seurat.obj.csfpbmc.pairs,patient.combined.matrix)
}
