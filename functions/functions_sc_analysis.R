## Support Functions for Single Cell analysis
library(DoubletFinder)
library(data.table)
library(ggrepel)
library(ggplot2)

# OBJECT GENERATION ---------------------------------------

# Takes samples from one patient and creates a Seurat Object
generatePatientCombinedMatrix <- function(dataset_loc, samples.vec, THREADS = 4) {
  # Access More clusters
  cl <- makeCluster(getOption("cl.cores", THREADS))
  clusterExport(cl,c('Read10X', "dataset_loc", "samples.vec"))
  start.time <- Sys.time()
  
  # Load sample Matricies
  d10x.data <- sapply(samples.vec, function(i){
    d10x <- Read10X(file.path(dataset_loc,i,"/"))
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    d10x
  })
  
  # Stop running mult-threads
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  stopCluster(cl)
  rm(cl)
  
  # Combine all patient matricies
  experiment.data <- do.call("cbind", d10x.data)
  
  return(experiment.data)
}


# SingleR combine, which combines a list of singleR objects into one
SingleR.Combine.fixed <- function (singler.list, order = NULL, clusters = NULL, expr = NULL, 
                                   xy = NULL) 
{
  singler = c()
  singler$singler = singler.list[[1]]$singler
  for (j in 1:length(singler.list[[1]]$singler)) {
    singler$singler[[j]]$SingleR.cluster = c()
    singler$singler[[j]]$SingleR.cluster.main = c()
    singler$singler[[j]]$SingleR.single$clusters = c()
  }
  singler$meta.data = singler.list[[1]]$meta.data
  singler$meta.data$clusters = c()
  singler$meta.data$xy = c()
  singler$meta.data$data.sets = rep(singler$meta.data$project.name, 
                                    length(singler$meta.data$orig.ident))
  for (i in 2:length(singler.list)) {
    for (j in 1:length(singler$singler)) {
      if (singler.list[[i]]$singler[[j]]$about$RefData != 
          singler.list[[1]]$singler[[j]]$about$RefData) {
        stop("The objects are not ordered by the same reference data.")
      }
      singler$singler[[j]]$about$Organism = c(singler$singler[[j]]$about$Organism, 
                                              singler.list[[i]]$singler[[j]]$about$Organism)
      singler$singler[[j]]$about$Citation = c(singler$singler[[j]]$about$Citation, 
                                              singler.list[[i]]$singler[[j]]$about$Citation)
      singler$singler[[j]]$about$Technology = c(singler$singler[[j]]$about$Technology, 
                                                singler.list[[i]]$singler[[j]]$about$Technology)
      singler$singler[[j]]$SingleR.single$labels = rbind(singler$singler[[j]]$SingleR.single$labels, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$labels)
      if (!is.null(singler$singler[[j]]$SingleR.single$labels1)) {
        singler$singler[[j]]$SingleR.single$labels1 = rbind(singler$singler[[j]]$SingleR.single$labels1, 
                                                            singler.list[[i]]$singler[[j]]$SingleR.single$labels1)
      }
      singler$singler[[j]]$SingleR.single$scores = rbind(singler$singler[[j]]$SingleR.single$scores, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$scores)
      singler$singler[[j]]$SingleR.single.main$labels = rbind(singler$singler[[j]]$SingleR.single.main$labels, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$labels)
      if (!is.null(singler$singler[[j]]$SingleR.single.main$labels1)) {
        singler$singler[[j]]$SingleR.single.main$labels1 = rbind(singler$singler[[j]]$SingleR.single.main$labels1, 
                                                                 singler.list[[i]]$singler[[j]]$SingleR.single.main$labels1)
      }
      singler$singler[[j]]$SingleR.single.main$scores = rbind(singler$singler[[j]]$SingleR.single.main$scores, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$scores)
      singler$singler[[j]]$SingleR.single$cell.names = c(singler$singler[[j]]$SingleR.single$cell.names, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$cell.names)
      singler$singler[[j]]$SingleR.single.main$cell.names = c(singler$singler[[j]]$SingleR.single.main$cell.names, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$cell.names)
      if (!is.null(singler$singler[[j]]$SingleR.single.main$pval)) {
        singler$singler[[j]]$SingleR.single.main$pval = c(singler$singler[[j]]$SingleR.single.main$pval, 
                                                          singler.list[[i]]$singler[[j]]$SingleR.single.main$pval)
      }
      if (!is.null(singler$singler[[j]]$SingleR.single$pval)) {
        singler$singler[[j]]$SingleR.single$pval = c(singler$singler[[j]]$SingleR.single$pval, 
                                                     singler.list[[i]]$singler[[j]]$SingleR.single$pval)
      }
    }
    singler$meta.data$project.name = paste(singler$meta.data$project.name, 
                                           singler.list[[i]]$meta.data$project.name, sep = "+")
    singler$meta.data$orig.ident = c(singler$meta.data$orig.ident, 
                                     singler.list[[i]]$meta.data$orig.ident)
    singler$meta.data$data.sets = c(singler$meta.data$data.sets, 
                                    rep(singler.list[[i]]$meta.data$project.name, length(singler.list[[i]]$meta.data$orig.ident)))
  }
  
  ## DEBUG filter order to match singler cell barcodes
  cells_singler_obj <- singler$singler[[j]]$SingleR.single$labels
  order <- order[order %in% row.names(cells_singler_obj)]
  ## 
  for (j in 1:length(singler$singler)) {
    if (!is.null(order)) {
      singler$singler[[j]]$SingleR.single$labels = singler$singler[[j]]$SingleR.single$labels[order,] 
      if (!is.null(singler$singler[[j]]$SingleR.single$labels1)) {
        singler$singler[[j]]$SingleR.single$labels1 = singler$singler[[j]]$SingleR.single$labels1[order, 
        ]
      }
      singler$singler[[j]]$SingleR.single$scores = singler$singler[[j]]$SingleR.single$scores[order, 
      ]
      singler$singler[[j]]$SingleR.single$cell.names = singler$singler[[j]]$SingleR.single$cell.names[order]
      singler$singler[[j]]$SingleR.single.main$labels = singler$singler[[j]]$SingleR.single.main$labels[order, 
      ]
      if (!is.null(singler$singler[[j]]$SingleR.single.main$labels1)) {
        singler$singler[[j]]$SingleR.single.main$labels1 = singler$singler[[j]]$SingleR.single.main$labels1[order, 
        ]
      }
      singler$singler[[j]]$SingleR.single.main$scores = singler$singler[[j]]$SingleR.single.main$scores[order, 
      ]
      singler$singler[[j]]$SingleR.single.main$cell.names = singler$singler[[j]]$SingleR.single.main$cell.names[order]
      if (!is.null(singler$singler[[j]]$SingleR.single$pval)) {
        singler$singler[[j]]$SingleR.single$pval = singler$singler[[j]]$SingleR.single$pval[order]
        singler$singler[[j]]$SingleR.single.main$pval = singler$singler[[j]]$SingleR.single.main$pval[order]
      }
    }
  }
  if (!is.null(clusters) && !is.null(expr)) {
    for (j in 1:length(singler$singler)) {
      if (is.character(singler$singler[[j]]$about$RefData)) {
        ref = get(singler$singler[[j]]$about$RefData)
      }
      else {
        ref = singler$singler[[j]]$about$RefData
      }
      singler$singler[[j]]$SingleR.clusters = SingleR("cluster", 
                                                      expr, ref$data, types = ref$types, clusters = factor(clusters), 
                                                      sd.thres = ref$sd.thres, genes = "de", fine.tune = T)
      singler$singler[[j]]$SingleR.clusters.main = SingleR("cluster", 
                                                           expr, ref$data, types = ref$main_types, clusters = factor(clusters), 
                                                           sd.thres = ref$sd.thres, genes = "de", fine.tune = T)
    }
    singler$meta.data$clusters = clusters
    if (!is.null(xy)) {
      singler$meta.data$xy = xy
    }
  }
  singler
}


# SingleR create Seurat: a automated workflow for generating a Seurat object with all the PC and regression done automatically
SingleR.CreateSeurat.fixed <-  function (project.name, sc.data, min.genes = 200, min.cells = 2, 
                                         regress.out = "nCount_RNA", npca = 10, resolution = 0.8, temp.dir = NULL) 
{
  mtgenes = "^mt-"
  if (packageVersion("Seurat") >= 3) {
    sc = CreateSeuratObject(sc.data, min.cells = min.cells, 
                            min.features = min.genes, project = project.name)
    percent.mito <- PercentageFeatureSet(object = sc, pattern = "^(?i)mt-")
    sc <- AddMetaData(object = sc, metadata = percent.mito, 
                      col.name = "percent.mito")
  }
  else {
    sc = CreateSeuratObject(sc.data, min.cells = min.cells, 
                            min.genes = min.genes, project = project.name)
    mito.genes <- grep(pattern = mtgenes, x = rownames(x = sc@data), 
                       value = TRUE, ignore.case = TRUE)
    percent.mito <- colSums((sc.data[mito.genes, ]))/colSums(sc.data)
    sc <- AddMetaData(object = sc, metadata = percent.mito, 
                      col.name = "percent.mito")
    sc <- NormalizeData(object = sc, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  }
  if (packageVersion("Seurat") >= 3) {
    sc <- SCTransform(object = sc, verbose = FALSE, do.correct.umi = T, vars.to.regress = "nCount_RNA")
    sc <- RunPCA(object = sc, verbose = FALSE)
    sc <- FindNeighbors(object = sc, dims = 1:30)
    sc <- FindClusters(object = sc)
    if (ncol(sc@assays$RNA@data) < 100) {
      sc <- RunTSNE(sc, perplexity = 10, dims = 1:npca)
    }
    else {
      sc <- RunTSNE(sc, dims = 1:30)
    }
    sc <- RunUMAP(sc, dims = 1:30, verbose = FALSE)
  }
  else {
    sc <- FindVariableGenes(object = sc, mean.function = ExpMean, 
                            dispersion.function = LogVMR, x.low.cutoff = 0.0125, 
                            x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F, 
                            do.plot = F)
    if (!is.null(regress.out)) {
      sc <- ScaleData(object = sc, vars.to.regress = regress.out)
    }
    else {
      sc <- ScaleData(object = sc)
    }
    sc <- RunPCA(object = sc, pc.genes = sc@var.genes, do.print = FALSE)
    sc <- FindClusters(object = sc, reduction.type = "pca", 
                       dims.use = 1:npca, resolution = resolution, print.output = 0, 
                       save.SNN = F, temp.file.location = temp.dir)
    if (ncol(sc@data) < 100) {
      sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T, 
                    perplexity = 10)
    }
    else {
      sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T, 
                    check_duplicates = FALSE)
    }
  }
  sc
}




# Adding Meta Data -----------------------------------------

# Add metadata from table
addRNAseqMetaData <- function(seurat_obj, gex_column = 2,metadata_sheet = "./metadata/metadata.csv"){
  # Read in Meta Data table
  metadata.df <- read.csv(metadata_sheet)
  # Filter out columns you do not want
  cols_to_exclude <- c("CSF_WBC", "CSF_WBC_Lymph.", "Notes", "MRI")
  metadata.df <- metadata.df[,!names(metadata.df) %in% cols_to_exclude]
  cols_metadata <- names(metadata.df)
  
  # Omit entries with no RNA-seq data
  metadata.df <- metadata.df[!metadata.df[,gex_column] %in% c(""),]
  
  sample_name_vec        <- seurat_obj$sample.name
  patient_vec            <- as.character(unique(metadata.df$Patient_ID))
  
  for(i in 3:length(cols_metadata)){
    categ <- cols_metadata[i]
    categ_vec <- sample_name_vec
    
    for(patient in patient_vec){
      # Skip if patient is not in seurat object
      if (length(sample_name_vec[grep(patient,sample_name_vec)]) == 0){next}
      
      patient_subset <- metadata.df[metadata.df$Patient_ID == patient,]
      for (r in 1:nrow(patient_subset)){
        sample <- as.character(patient_subset[r,gex_column])
        categ_value <- as.character(patient_subset[r,i])
        if (sample %in% categ_vec){
          categ_vec[grep(sample, categ_vec)] <- categ_value
        }else{
          next
        }
      }
    }
    # Add categ_vec to Meta Data of Seurat object
    seurat_obj <- AddMetaData(object = seurat_obj,metadata = categ_vec,col.name = categ)
  }
  
  return(seurat_obj)
}


# RNA-Seq ANALYSIS -----------------------------------------
addCellAnnotationsToObject <- function(singler_seurat_obj){
  subtype_cellannot_vec <- singler_seurat_obj$singler[[2]]$SingleR.single$labels
  main_cellannot_vec <- singler_seurat_obj$singler[[2]]$SingleR.single.main$labels
  
  singler_seurat_obj$seurat$subtype.cell.annot <- subtype_cellannot_vec
  singler_seurat_obj$seurat$main.cell.annot    <- main_cellannot_vec
  
  return(singler_seurat_obj)
}

# Only works with results from SingleR
addSingleRtoSeurat <- function(seurat,singler.main,singler.fine,main.name="main.cell.annot",fine.name="subtype.cell.annot"){
  # Add BlueprintENCODE annotations
  subtype_annot        <- singler.fine$labels
  names(subtype_annot) <- row.names(singler.fine)
  subtype_annot <- subtype_annot[colnames(seurat)]
  seurat[[fine.name]] <- subtype_annot
  
  main_annot        <- singler.main$labels
  names(main_annot) <- row.names(singler.main)
  main_annot <- main_annot[colnames(seurat)]
  seurat[[main.name]] <- main_annot
  
  return(seurat)
}




# add 10X clonotype results to seurat object
#  - Read in formatted table of all single cell clonotype results
#  - <create_igseq_results_table_5prime.pl> was used to create this table
add5primeClonotypeResults <- function(seurat_obj, igseq_results = "./Resources/igseq_results_tables/MERGED_filtered_contigs_ALLpatients.forSeurat.csv"){
  igseq.results.df  <- read.csv(igseq_results)
  col_igseq_results <- names(igseq.results.df)
  cols_to_keep      <- c("clonotype_ID","consensus_ID","flag","CDR3","cluster_id_Heavy","expanded_clonotype_status","shared_CDR3_within_1patient","shared_CDR3_across_patients","flag_contam_cell_entry","SUSET_EXTENDED")
  target_cols       <- col_igseq_results[col_igseq_results %in% cols_to_keep]
  
  sample_name_vec        <- seurat_obj$sample.name
  patient_vec            <- tstrsplit(names(sample_name_vec),"_")[[1]]
  patient_vec            <- unique(tstrsplit(patient_vec,"-")[[2]])

  for(i in 1:length(col_igseq_results)){
    categ <- col_igseq_results[i]
    if (! categ %in% target_cols){next}
    categ_vec <- sample_name_vec
    
    present_cells <- c()
    for(patient in patient_vec){
      # Skip if patient is not in seurat object
      if (length(sample_name_vec[grep(patient,sample_name_vec)]) == 0){next}
      
      patient_subset <- igseq.results.df[igseq.results.df$Patient_ID == patient,]
      for (r in 1:nrow(patient_subset)){
        sample_temp <- as.character(patient_subset[r,3]) # Sample name for GEX is in column 3
        categ_value_temp <- as.character(patient_subset[r,i])
        if (is.na(categ_value_temp)){categ_value_temp <- ""}
        chain            <- as.character(patient_subset[r,10])
        cell_barcode     <- as.character(patient_subset[r,9])
        
        # for downstream clonotype results (cluster_id_Heavy) -- skip all Light chain entries
        if (chain != "IGH"){next}
        
        # Need a way to account for the fact that the same cell
        # could have 2 rows of results for Heavy and Light chain
        
        if (!is.na(categ_vec[cell_barcode]) && length(categ_vec[cell_barcode]) != 0){
          categ_vec[cell_barcode] <- categ_value_temp
          present_cells <- c(present_cells,cell_barcode)
        }
      }
    }
    cells_with_results <- unique(present_cells)
    
    # For cells not in the 10X VDJ results, leave empty or put NA
    categ_vec[!names(categ_vec) %in% cells_with_results] <- ""
    # Add categ_vec to Meta Data of Seurat object
    seurat_obj <- AddMetaData(object = seurat_obj,metadata = categ_vec,col.name = categ)
  }
  
  return(seurat_obj)
}



# Adding Immcantation result into the Seurat Object
addImmcantationResults <- function(seurat_obj, TYPE = "TCR",gex_sample_col=5,seurat_cell_barcode_col=6,igseq_results = "./Resources/igseq_results_tables/tables_v3.1/FORMATTED_immcantation_results_table.ALL.v3.1.csv"){
  igseq.results.df  <- read.csv(igseq_results)
  
  # Subet Heavy chains (currently only Heavy chain meta data supported)
  if (TYPE == "BCR"){
    igseq.results.df <- igseq.results.df[igseq.results.df$chain %in% c("IGH"),]
  }
  if (TYPE == "TCR"){
    igseq.results.df <- igseq.results.df[igseq.results.df$chain %in% c("TRB"),]
  }
  
  col_igseq_results <- names(igseq.results.df)
  cols_to_keep      <- c("cgene_imm","clonotypeID_imm","clonotype_cell_count_imm","clonotype_CSF_cell_count_imm","clonotype_PB_cell_count_imm","expanded_clonotype_status_imm","MERGED_clonotypeID_imm")
  target_cols       <- col_igseq_results[col_igseq_results %in% cols_to_keep]
  
  sample_name_vec        <- seurat_obj$sample.name
  patient_vec            <- tstrsplit(names(sample_name_vec),"_")[[1]]
  patient_vec            <- unique(tstrsplit(patient_vec,"-")[[2]])
  
  for(i in 1:length(col_igseq_results)){
    categ <- col_igseq_results[i]
    if (! categ %in% target_cols){next}
    categ_vec <- sample_name_vec
    
    present_cells <- c()
    for(patient in patient_vec){

      # Skip if patient is not in seurat object
      if (length(sample_name_vec[grep(patient,sample_name_vec)]) == 0){next}
      # Skip if patient not in Ig-Seq results (happens if there are patients in the Seurat object, but not the BCR results)
      if (length(igseq.results.df$Patient_ID[igseq.results.df$Patient_ID == patient]) == 0){next}
      
      patient_subset <- igseq.results.df[igseq.results.df$Patient_ID == patient,]
      for (r in 1:nrow(patient_subset)){
        sample_temp <- as.character(patient_subset[r,gex_sample_col])      # Sample name for GEX is in column 4
        categ_value_temp <- as.character(patient_subset[r,i])
        if (is.na(categ_value_temp)){categ_value_temp <- ""}
        cell_barcode     <- as.character(patient_subset[r,seurat_cell_barcode_col]) # Seurat Cell Barcode
        
        # These results only have Heavy Chain
        if (!is.na(categ_vec[cell_barcode]) && length(categ_vec[cell_barcode]) != 0){
          categ_vec[cell_barcode] <- categ_value_temp
          present_cells <- c(present_cells,cell_barcode)
        }
      }
    }
    cells_with_results <- unique(present_cells)
    
    # For cells not in the 10X VDJ results, leave empty or put NA
    categ_vec[!names(categ_vec) %in% cells_with_results] <- ""
    # Add categ_vec to Meta Data of Seurat object
    seurat_obj <- AddMetaData(object = seurat_obj,metadata = categ_vec,col.name = categ)
  }
  
  return(seurat_obj)
}


# Filter out cells with extremely poor CD79B and XBP1 expression
filter_out_cells_using_raw_gene_expression <- function(seurat_obj, threshCD79B = 0){
  cell_ids   <- colnames(seurat_obj)
  raw_counts <- GetAssay(object = seurat_obj,slot="counts")
  CD79B_counts <- raw_counts[row.names(raw_counts) %in% c("CD79B"),]

  valid_cells <- c()
  for(cell in colnames(raw_counts)){
    CD79B  <- CD79B_counts[,cell]

    if (CD79B == "."){CD79B <- 0}
    else             {CD79B <- as.numeric(CD79B)}
    
    if (CD79B > 0){
      if (CD79B > threshCD79B){valid_cells <- c(valid_cells,cell)}
    }
  }
  
  # Subset for valid cells
  seurat_obj <- seurat_obj[,valid_cells]

  return(seurat_obj)
}

recluster_seurat_obj_v3 <- function(seurat_obj,dims = 1:30,npca = 10) {
  seurat_obj <- RunPCA(object = seurat_obj, verbose = FALSE)
  seurat_obj <- FindNeighbors(object = seurat_obj, dims = dims)
  seurat_obj <- FindClusters(object = seurat_obj)
  
  if (ncol(seurat_obj@assays$RNA@data) < 100) {
    seurat_obj <- RunTSNE(seurat_obj, perplexity = 10, dims = 1:npca)
  }
  else {
    seurat_obj <- RunTSNE(seurat_obj, dims = dims)
  }
  seurat_obj <- RunUMAP(seurat_obj, dims = dims, verbose = FALSE)
  
  return(seurat_obj)
}



# Re-cluster B cell Subset
recluster_Bcells <- function(bcells,use.pcs=1:3,dims.use=1:3,resolution=seq(0.5),perplexity=10){
  ## 1. Run PCA
  bcells <- RunPCA(
    object = bcells,
    pc.genes = bcells@var.genes,
    do.print = TRUE,
    pcs.print = 1:5,
    genes.print = 5,
    pcs.compute = 40,
    maxit = 500)
  
  PrintPCAParams(bcells)
  
  # Plot PCs 1 to 40 to determine what range of PCs to use
  PCElbowPlot(
    bcells,
    num.pc = 40)
  
  ## 2. Find Clusters
  bcells <- FindClusters(
    object = bcells, 
    reduction.type = "pca", 
    dims.use = use.pcs, 
    resolution = resolution, 
    print.output = FALSE, 
    save.SNN = TRUE)
  
  PrintFindClustersParams(object = bcells)
  
  ## 3. Re-cluster using tSNE
  bcells <- RunTSNE(object = bcells,reduction.use = "pca",dims.use = dims.use,perplexity = 10,do.fast = TRUE)
  
  ## 4. Add 'cluster_loc' category to meta data, new cluster assignments will be in 'ident'
  vec.cell.clus <- bcells@ident
  vec.cell.loc  <- names(vec.cell.clus)
  vec.cell.loc  <- tstrsplit(vec.cell.loc,"_")[[2]]
  
  vec.cell.clus_loc <- vec.cell.loc
  for (i in 1:length(vec.cell.loc)) {
    clus <- as.character(vec.cell.clus[i])
    loc  <- vec.cell.loc[i]
    clus_loc_id <- paste(c(clus,"_",loc),collapse = "")
    vec.cell.clus_loc[i] <- clus_loc_id 
  }
  names(vec.cell.clus_loc) <- names(vec.cell.clus)
  
  # Add Meta Data
  bcells <- AddMetaData(
    object = bcells,
    metadata = vec.cell.clus_loc,
    col.name= "clus_loc")
  
  return(bcells)
}

# Monocle v3
run_monocle_v3 <- function(seurat_obj,DIMS = 100) {
  ### Convert Seurat Object into Monocle Object
  data <- as(as.matrix(seurat_obj@assays$RNA@data), 'sparseMatrix')
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  
  cds <- new_cell_data_set(data,cell_metadata = seurat_obj@meta.data,gene_metadata = fData)
  
  ### Step 1: pre-process
  #Ex: DIMS = 5 CSF B cells, 100 for all B cells
  cds = preprocess_cds(cds, num_dim = DIMS)
  ### Step 2: UMAP reduction
  cds = reduce_dimension(cds, reduction_method = "UMAP")
  ### Step 3: cluster cells
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  ### Step 4: learn Trajectory
  cds <- learn_graph(cds)
  
  return(cds)
}


# Add Monocle Clusters to Seurat Meta Data
add_monocle_results_to_metadata <- function(seurat_obj,monocle_obj){
  louvain_vec <- unlist(monocle_obj@clusters)$UMAP.clusters
  seurat_obj <- AddMetaData(object = seurat_obj,metadata = louvain_vec,col.name = "monocle_louvain_clusters")
  
  return(seurat_obj)
}


# Create a summary table with the cell count breakdown per annotation
create_subtype_summary_df <- function(annot_vec,categ_vec) {
  all_sub_celltypes  <- unique(annot_vec)
  all_clusters       <- as.character(unique(categ_vec))
  
  
  # Create Subtype cell table
  sub.table.df <- data.frame(matrix(nrow = length(all_clusters), ncol = length(all_sub_celltypes)))
  names(sub.table.df) <- all_sub_celltypes
  row.names(sub.table.df) <- all_clusters
  
  subtype.table.counts <- data.frame(table(annot_vec,categ_vec))
  names(subtype.table.counts) <- c("CellType","Clus","Freq")
  
  for (i in 1:nrow(subtype.table.counts)) {
    celltype <- as.character(subtype.table.counts[i,1])
    clus     <- as.character(subtype.table.counts[i,2])
    count    <- as.character(subtype.table.counts[i,3])
    sub.table.df[clus,celltype] <- count
  }
  
  return(sub.table.df)
}

subsetObjByUMAPcoords <- function(seurat_obj, UMAP1=0, UMAP2=0, LESSTHAN_UMAP1 = TRUE,LESSTHAN_UMAP2=TRUE,INVERSE_SUBSET=FALSE) {
  umap_coords <- as.data.frame(Embeddings(object = seurat_obj@reductions[["umap"]]))
  
  # UMAP 1 threshold
  if(LESSTHAN_UMAP1){
    umap_coords <- umap_coords[umap_coords$UMAP_1 <= UMAP1,]
  } else{
    umap_coords <- umap_coords[umap_coords$UMAP_1 >= UMAP1,]
  }
  
  # UMAP 2 threshold
  if(LESSTHAN_UMAP2){
    umap_coords <- umap_coords[umap_coords$UMAP_2 <= UMAP2,]
  }else{
    umap_coords <- umap_coords[umap_coords$UMAP_2 >= UMAP2,]
  }
  
  
  # Get cell barcodes of the 3 cells
  target_cells <- row.names(umap_coords)
  all_cells    <- names(seurat_obj$orig.ident)
  if(INVERSE_SUBSET){
    target_cells <- all_cells[! all_cells %in% target_cells]
  } 
  
  seurat_obj <- seurat_obj[,target_cells]
  
  return(seurat_obj)
}


# DGE ------------------------------------------------------

# Generate table of to DE genes (FindMarkers) in a given Seurat object with an input vector of select cells
#   - Function only works post SingleR & after Meta Data has been added
dge_select_cells <- function(seurat_obj,cells_of_interest, bcells.only = T){
  seurat_obj.interest <- seurat_obj
  if (bcells.only){
    cells_vec <- seurat_obj$subtype.cell.annot[seurat_obj$subtype.cell.annot %in% c("naive B-cells","Memory B-cells","Class-switched memory B-cells","Plasma cells")]
    seurat_obj.interest <- seurat_obj[,names(cells_vec)]
  }
  
  remaining_cells <- colnames(seurat_obj.interest)
  remaining_cells <- remaining_cells[! remaining_cells %in% cells_of_interest]
  
  new_comparison <- character(length = length(colnames(seurat_obj.interest)))
  names(new_comparison) <- colnames(seurat_obj.interest)
  new_comparison[cells_of_interest] <- "1"
  new_comparison[remaining_cells]   <- "0"
  seurat_obj.interest$fuzzy_bcells  <- new_comparison
  
  comp_results<- FindMarkers(seurat_obj.interest,ident.1 = "1",group.by = "fuzzy_bcells")
  return(comp_results)
}

# For a given Seurat object, this code will iterate through all [diagnosis]_[B cell Subtype]_[CSF/PB] comparisons
#  - This function works as long as all meta data has been added to the Seurat Object
#  - returns nested list of FindMarkers results
dge_all_bcell_comparisons <- function(bcells_csfpb,bcell_annot_vec,MIN_CELLS = 3,return.thresh = 0.01,min.pct = 0.1,logfc.threshold = 0.25,test.use = "wilcox"){
  dge_results_list <- list()

  diagnosis_vec <- unique(bcells_csfpb$Diagnostic_Category)
  for (diagnosis in diagnosis_vec) {
    cells_vec <- bcells_csfpb$Diagnostic_Category[bcells_csfpb$Diagnostic_Category %in% diagnosis]
    bcells_csfpb.diag <- bcells_csfpb[,names(cells_vec)]
    
    # Blueprint ENCODE Annotations
    bcell_annot_vec_diag <- bcell_annot_vec[names(cells_vec)]
    bcell.subtype.annot <- unique(bcell_annot_vec_diag)
    for (subtype.annot in bcell.subtype.annot) {
      bcells_vec <- bcell_annot_vec_diag[bcell_annot_vec_diag %in% subtype.annot]
      bcells_csfpb.diag.subtype <- bcells_csfpb.diag[,names(bcells_vec)]
      if(length(unique(bcells_csfpb.diag.subtype$loc)) == 2) {
        csf_count <- length(bcells_csfpb.diag.subtype$loc[bcells_csfpb.diag.subtype$loc == "CSF"])
        pb_count  <- length(bcells_csfpb.diag.subtype$loc[bcells_csfpb.diag.subtype$loc == "PB"])
        if (csf_count < MIN_CELLS || pb_count < MIN_CELLS){next}
        
        # Set all Idents to loc
        Idents(bcells_csfpb.diag.subtype) <- bcells_csfpb.diag.subtype$loc
        comp_label <- paste(c(diagnosis,subtype.annot),collapse = "_")
        cat(paste0(">>> ", comp_label,"\n\n"))
        comp_results_subtype <- FindAllMarkers(bcells_csfpb.diag.subtype,
                                               min.pct = min.pct,logfc.threshold = logfc.threshold,return.thresh = return.thresh,
                                               test.use = test.use)
        
        dge_results_list[[comp_label]] <- comp_results_subtype 
      }
    }
  }
  
  return(dge_results_list)
}

# Compare cells with expanded clonotypes vs. non expanded
# Input: - should be object that is only CSF or only PB
#        - should only have cells with Ig-Seq results
#        - 'imm_thresh' only supports: 1,2,3 (default = 1)
dge_all_igseq_bcell_comparisons <- function(bcells_csfpb,bcell_annot_vec,MIN_CELLS = 3,return.thresh = 0.01,min.pct = 0.1,logfc.threshold = 0.25,test.use = "wilcox",imm_thresh = 1){
  dge_results_list <- list()
  
  diagnosis_vec <- unique(bcells_csfpb$Diagnostic_Category)
  for (diagnosis in diagnosis_vec) {
    cells_vec <- bcells_csfpb$Diagnostic_Category[bcells_csfpb$Diagnostic_Category %in% diagnosis]
    bcells_csfpb.diag <- bcells_csfpb[,names(cells_vec)]
    
    bcell_annot_vec_diag <- bcell_annot_vec[names(cells_vec)]
    bcell.subtype.annot <- unique(bcell_annot_vec_diag)
    for (subtype.annot in bcell.subtype.annot) {
      bcells_vec <- bcell_annot_vec_diag[bcell_annot_vec_diag %in% subtype.annot]
      if (length(bcells_vec) == 0){next}
      bcells_csfpb.diag.subtype <- bcells_csfpb.diag[,names(bcells_vec)]
      
      imm_clonotype_cell_vec <- c()
      if (imm_thresh == 1){
        imm_clonotype_cell_vec <- bcells_csfpb.diag.subtype$expanded_clonotype_status_imm
      } else if (imm_thresh == 2){
        imm_clonotype_cell_vec <- bcells_csfpb.diag.subtype$expanded_clonotype_status_imm2
      } else if (imm_thresh == 3){
        imm_clonotype_cell_vec <- bcells_csfpb.diag.subtype$expanded_clonotype_status_imm3
      }
      
      if(length(unique(imm_clonotype_cell_vec)) == 2) {
        exp_count <- length(imm_clonotype_cell_vec[imm_clonotype_cell_vec == "1"])
        nonexp_count  <- length(imm_clonotype_cell_vec[imm_clonotype_cell_vec == "0"])
        if (exp_count < MIN_CELLS || nonexp_count < MIN_CELLS){next}
        
        # Set all Idents to loc
        Idents(bcells_csfpb.diag.subtype) <- imm_clonotype_cell_vec
        comp_label <- paste(c(diagnosis,subtype.annot),collapse = "_")
        cat(paste0(">>> ", comp_label,"\n\n"))
        comp_results_subtype <- FindAllMarkers(bcells_csfpb.diag.subtype,
                                               min.pct = min.pct,logfc.threshold = logfc.threshold,return.thresh = return.thresh,
                                               test.use = test.use)
        
        dge_results_list[[comp_label]] <- comp_results_subtype 
      }
    }
  }
  
  return(dge_results_list)
}

# Remove Ig genes speicfically V genes from DGE results
filter_out_ig_markers <- function(dge_results){
  genes_to_keep <- dge_results$gene
  genes_mt   <- genes_to_keep[grep("MT-",genes_to_keep)]
  genes_igh  <- genes_to_keep[grep("IGHV",genes_to_keep)]
  genes_igk  <- genes_to_keep[grep("IGKV",genes_to_keep)]
  genes_igl  <- genes_to_keep[grep("IGLV",genes_to_keep)]
  genes_to_rm <- c(genes_mt,genes_igh,genes_igk,genes_igl)

  genes_to_keep <- genes_to_keep[! genes_to_keep %in% genes_to_rm]
  
  dge_results.noig <- dge_results[dge_results$gene %in% genes_to_keep,]
  
  return(dge_results.noig)
}


# QC ------------------------------------------------------

# Use DoubletFinder to annotate which cells have doublets in a given Seurat object
find_doublets <- function(seurat_obj){
  sample_list      <- as.character(unique(seurat_obj$sample.name))
  doublet.rate.vec <- seurat_obj$sample.name

  for (sample in sample_list){
    sample_cell_count <- length(seurat_obj$sample.name[seurat_obj$sample.name == sample])
    cat(paste0(sample,"\t",sample_cell_count,"\n"))
    
    sample_cells      <- names(seurat_obj$sample.name[seurat_obj$sample.name == sample])
    sample_seurat_obj <- seurat_obj[,sample_cells]  # SubsetData

    # Calculate expected doublets
    doublet_rate <- unique(as.numeric(as.character(sample_seurat_obj$DoubletRate)))
    exp_doublets <- round(sample_cell_count * doublet_rate)

    # Run DoubletFinder
    doublet_finder_results <- character(length = ncol(sample_seurat_obj))
    doublet_finder_results[1:length(doublet_finder_results)] <- "Singlet"
    if(exp_doublets != 0){
      ## Figure out how to calculate pK (0.01 used in example), this could effect how many doublets are found
      sample_seurat_obj <- doubletFinder_v3(seu = sample_seurat_obj,pK = doublet_rate,PCs = 1:10,nExp = exp_doublets,sct = TRUE)
      
      sample_output <- sample_seurat_obj@meta.data[,grep("DF.classifications", names(sample_seurat_obj@meta.data))]
      names(sample_output)   <- colnames(sample_seurat_obj)
      doublet_finder_results <- sample_output
    }
    names(doublet_finder_results) <- colnames(sample_seurat_obj)
    doublet.rate.vec[names(doublet.rate.vec) %in% names(doublet_finder_results)] <- doublet_finder_results
  }
  
  ## Add to Meta Data
  seurat_obj$doublet.finder.ann <- doublet.rate.vec
  
  return(seurat_obj)
}

# IGSEQ ANALYSIS ------------------------------------------

# Takes immcantation results and creates a count matrix for making heatmaps
#  - normalizes using Z-scores
#  - supported modes:   "read","umi","cell"
#  - supported usage:   "V_CALL","J_CALL","C_CALL"
create_VJC_usage_matrix <- function(db.sub,splitby = "_", mode = "read", usage = "V_CALL_GENE",patient = "PATIENT",norm=TRUE) {
  no_gene_label    <- "missing"
  db.sub$CELLCOUNT <- replicate(length(db.sub$PATIENT),1)
  
  # Gene gene usage (V,J,C)
  gene_col <- names(db.sub)[names(db.sub) %in% usage]
  gene_usage_vec <- as.data.frame(db.sub)[,gene_col]
  gene_usage_vec[is.na(gene_usage_vec)] <- no_gene_label
  db.sub[,gene_col] <- gene_usage_vec

  count_categ = ""
  if (mode == "read"){
    count_categ = "CONSCOUNT"
  } else if (mode == "umi"){
    count_categ = "UMICOUNT"
  } else if (mode == "cell"){
    count_categ <- "CELLCOUNT"
  }else {
    count_categ = "CONSCOUNT"
  }
  
  # count matrix
  patient_vec_all <- db.sub$PATIENT
  if (patient == "SAMPLE"){patient_vec_all <- db.sub$SAMPLE}
  if (patient == "Paper_ID_Diagnosis_loc") {patient_vec_all <- db.sub$Paper_ID_Diagnosis_loc}
  
  patient_vec <- unique(patient_vec_all)
  gene_vec   <- unique(gene_usage_vec)
  
  count_df <- as.data.frame(matrix(0,nrow = length(gene_vec), ncol = length(patient_vec)))
  names(count_df)     <- patient_vec
  row.names(count_df) <- gene_vec
  
  df <- db.sub[,c(patient,usage,count_categ)]

  # Add up all the read counts
  for (i in 1:nrow(df)) {
    patient <- as.character(df[i,1])
    gene    <- as.character(df[i,2])
    count   <- as.numeric(df[i,3])
    count_df[gene,patient] <- count_df[gene,patient] + count
    
  }
  
  ### NORMALIZE count data frame
  # normalized by column
  count_table <- t(as.matrix(count_df))
  if (norm){
    norm.sum.count_table <- count_table/rowSums(count_table)
    
    
    # Calculate Z-score for each cell
    count_table_zscores <- norm.sum.count_table
    table_row_names <- row.names(count_table_zscores)
    table_col_names <-colnames(count_table_zscores)
    
    for(r in 1:length(row.names(count_table_zscores))){
      row_data <- count_table_zscores[r,]
      std_dev  <- sd(row_data)
      row_mean <- mean(row_data)
      for(c in 1:length(colnames(count_table_zscores))) {
        z_score <- (count_table_zscores[r,c] - row_mean)/std_dev
        count_table_zscores[r,c] <- z_score
      }
    }
    # Inverse so patients are column names again (not pretty)
    count_table <- count_table_zscores
  }
  
  return(count_table)
}

# Exactly like the 'create_VJC_usage_matrix' fuction except it uses clonotype ID as row names instead of patient ID
clonotype_VJC_usage_matrix <- function(db.sub, cell_count_vec,mode = "read", usage = "C_CALL",norm=TRUE) {
  db.sub$CELLCOUNT <- cell_count_vec
  
  # Gene gene usage (V,J,C)
  gene_col <- names(db.sub)[names(db.sub) %in% usage]
  gene_usage_vec <- as.data.frame(db.sub)[,gene_col]
  gene_usage_vec <- tstrsplit(gene_usage_vec, "\\*")[[1]]
  db.sub[,gene_col] <- gene_usage_vec
  
  count_categ = ""
  if (mode == "read"){
    count_categ = "CONSCOUNT"
  } else if (mode == "umi"){
    count_categ = "UMICOUNT"
  } else if (mode == "cell"){
    count_categ <- "CELLCOUNT"
  }else {
    count_categ = "CONSCOUNT"
  }
  
  # count matrix
  clone_vec <- unique(db.sub$CLONE)
  gene_vec    <- unique(gene_usage_vec)
  count_df <- as.data.frame(matrix(0,nrow = length(gene_vec), ncol = length(clone_vec)))
  names(count_df)     <- clone_vec
  row.names(count_df) <- gene_vec
  
  df <- db.sub[,c("CLONE",usage,count_categ)]
  
  # Add up all the read counts
  for (i in 1:nrow(df)) {
    clone_id <- as.character(df[i,1])
    gene     <- as.character(df[i,2])
    count    <- as.numeric(df[i,3])
    count_df[gene,clone_id] <- count_df[gene,clone_id] + count
    
  }
  
  ### NORMALIZE count data frame
  # normalized by column
  count_table <- t(as.matrix(count_df))
  if (norm){
    norm.sum.count_table <- count_table/rowSums(count_table)
    
    
    # Calculate Z-score for each cell
    count_table_zscores <- norm.sum.count_table
    table_row_names <- row.names(count_table_zscores)
    table_col_names <-colnames(count_table_zscores)
    
    for(r in 1:length(row.names(count_table_zscores))){
      row_data <- count_table_zscores[r,]
      std_dev  <- sd(row_data)
      row_mean <- mean(row_data)
      for(c in 1:length(colnames(count_table_zscores))) {
        z_score <- (count_table_zscores[r,c] - row_mean)/std_dev
        count_table_zscores[r,c] <- z_score
      }
    }
    # Inverse so patients are column names again (not pretty)
    count_table <- count_table_zscores
  }
  
  return(count_table)
}


# Used to generate MIXCR output for either the Heavy or Light chain of a given clonotype
runMixcrPerChain <- function(fasta,results_dir,threads = 4,chain = "",MIXCR = "./programs/mixcr-3.0.3/mixcr"){
  name_comp <- chain
  if (nchar(name_comp) > 0){
    name_comp <- paste(c('_',chain),collapse = "")
  }
  thread_str         <- paste(c('-t',threads),collapse = " ")
  aln_vdjca_fh       <- paste(c(results_dir,"/alignments",name_comp,".vdjca"),collapse = "")
  aln_pretty_fh      <- paste(c(results_dir,"/mult_alignment",name_comp,".txt"),collapse = "")
  aln_report_fh      <- paste(c(results_dir,"/alignment_report.txt"),collapse = "")
  assemble_report_fh <- paste(c(results_dir,"/assemble_report.txt"),collapse = "")
  clones_fh          <- paste(c(results_dir,"/clones",name_comp,".clns"),collapse = "")
  clones_txt_fh      <- paste(c(results_dir,"/clones",name_comp,".txt"),collapse = "")
  clonesCDR3_fh      <- paste(c(results_dir,"/clones_cdr3",name_comp,".txt"),collapse = "")
  
  # Run MIXCR
  setwd(results_dir)
  system2(MIXCR,c("align","-f","-s hsa",thread_str,"-OvParameters.geneFeatureToAlign=VTranscript","-r alignment_report.txt",fasta,aln_vdjca_fh))
  system2(MIXCR,c("exportAlignmentsPretty","-f",aln_vdjca_fh,aln_pretty_fh))
  system2(MIXCR,c("assemble","-f","-r",assemble_report_fh,"-OassemblingFeatures=VDJRegion",aln_vdjca_fh,clones_fh))
  system2(MIXCR,c("exportClones","-f",clones_fh,clones_txt_fh))
  system2(MIXCR,c("exportClones","-f","-vHit","-jHit","-dHit","-cHit","-count","-aaFeature","CDR3",clones_fh,clonesCDR3_fh))
  
  # Collect all relavent results
  results.mixcr <- read.delim(clones_txt_fh)
  
  return(results.mixcr)
}


# Create a summary table with the cell count breakdown per annotation
create_VDJ_summary_df <- function(annot_vec,categ_vec) {
  all_sub_celltypes  <- as.character(unique(annot_vec))
  all_clusters       <- as.character(unique(categ_vec))
  
  
  # Create Subtype cell table
  sub.table.df <- data.frame(matrix(nrow = length(all_clusters), ncol = length(all_sub_celltypes)))
  names(sub.table.df) <- all_sub_celltypes
  row.names(sub.table.df) <- all_clusters
  
  subtype.table.counts <- data.frame(table(annot_vec,categ_vec))
  names(subtype.table.counts) <- c("CellType","Clus","Freq")
  
  for (i in 1:nrow(subtype.table.counts)) {
    celltype <- as.character(subtype.table.counts[i,1])
    clus     <- as.character(subtype.table.counts[i,2])
    count    <- as.character(subtype.table.counts[i,3])
    sub.table.df[clus,celltype] <- count
  }
  
  return(sub.table.df)
}


# MANUAL ANNOTATION ---------------------------------------

# Takes a subsetted Seurat object and relabels all non B cells something else
# - this function is for cells that cluster with B cells but have a different annotation
addManualBcellAnnot <- function(seurat_obj, label = "Hybrid B-cells"){
  bcell.subtype.annot <- c("naive B-cells","Memory B-cells","Class-switched memory B-cells","Plasma cells")
  
  manual_annotations <- seurat_obj$subtype.cell.annot
  hybrid_cells_vec <- seurat_obj$subtype.cell.annot[! seurat_obj$subtype.cell.annot %in% bcell.subtype.annot]
  manual_annotations[names(hybrid_cells_vec)] <- label
  
  # Add in meta data for manual annotations
  seurat_obj$manual.cell.annot <- manual_annotations
  
  return(seurat_obj)
}

# PLOTS ---------------------------------------------------
plotPDFsMonocle <- function(monocle_obj, dir_path,sub_dir = "",file_name = "plots.pdf",plot_type="UMAP",subtype="subtype.cell.annot.bencode") {
  if (length(sub_dir) != 0){
    dir_path <- paste(c(dir_path,sub_dir),collapse = '/')
    dir.create(dir_path)
  }
  # List all meta data categories in object
  column_names <- names(pData(monocle_obj))
  
  pdf(file = paste0(dir_path,"/",file_name))
  if (length(column_names[column_names %in% c(subtype)]) != 0){
    p1 <- plot_cells(monocle_obj, color_cells_by = subtype, label_groups_by_cluster=F,label_leaves=F,label_branch_points=F,show_trajectory_graph = F,cell_size = 4)
    print(p1)
  }
  if (length(column_names[column_names %in% c("loc")]) != 0){
    p3 <- plot_cells(monocle_obj, color_cells_by = "loc", label_groups_by_cluster=F,label_leaves=F,label_branch_points=F,show_trajectory_graph = F,cell_size = 2)
    print(p3)
  }
  if (length(column_names[column_names %in% c("seurat_clusters")]) != 0){
    p4 <- plot_cells(monocle_obj, color_cells_by = "seurat_clusters", label_groups_by_cluster=F,label_leaves=T,label_branch_points=T,show_trajectory_graph = F,cell_size = 4)
    print(p4)
  }
  dev.off()
  
}

plotTIFFsMonocle <- function(monocle_obj, dir_path,sub_dir = "",file_name = "plots.tiff",plot_type="UMAP",subtype="subtype.cell.annot.bencode",vec) {
  if (length(sub_dir) != 0){
    dir_path <- paste(c(dir_path,sub_dir),collapse = '/')
    dir.create(dir_path)
  }
  # List all meta data categories in object
  column_names <- names(pData(monocle_obj))
  
  tiff(filename = paste0(dir_path,"/",file_name),height = 700,width = 700)
  if (length(column_names[column_names %in% c(subtype)]) != 0){
    p1 <- plot_cells(monocle_obj, color_cells_by = subtype,label_groups_by_cluster=F,label_leaves=F,label_branch_points=F,label_cell_groups = F,show_trajectory_graph = F,cell_size = 4)
    print(p1)
  }
  dev.off()
  
}



# Generate a volcano plot for a given table of DGE results (output from FindAllMarkers)
volcano_from_dge <- function(dge_results, comp = "CSF", with_labels = TRUE, dir_path = "./",file_name = "volcano_plots.pdf"){
  all_comp_categs <- names(dge_results)
  pdf(file = paste0(dir_path,"/",file_name))
  for (bcell_categ in all_comp_categs) {
    top_genes     <- dge_results[[bcell_categ]]
    top_genes_csf <- top_genes[top_genes$cluster == comp,]
    results_csf   <- mutate(top_genes_csf, sig=ifelse(-log10(top_genes_csf$p_val_adj)>5, "FDR<0.05", "Not Sig"))
    
    # CSF
    p = ggplot(results_csf, aes(avg_logFC, -log10(p_val_adj+10^-20))) +
      geom_point(aes(col=sig)) +
      scale_color_manual(values=c("red", "black")) +
      scale_y_continuous(limits = c(0, 20)) +
      scale_x_continuous(limits = c(-2,2)) +
      ggtitle(bcell_categ) + theme(plot.title = element_text(hjust = 0.5)) + 
      xlab("Log FC") +
      ylab("-log10 p-value") +
      labs(" ")
    
    if (with_labels){
      p <- p+geom_text_repel(data=filter(results_csf, -log10(top_genes_csf$p_val_adj)>5), aes(label=gene), size = 5)
    }
    
    print(p)
  }
  dev.off()
  
}


# Plot optimal number of PCs for re-clustering a Seurat object
PlotOptimalPCsforSeurat <- function(seurat_obj) {
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  
  plot_df <- data.frame(pct = pct, 
                        cumu = cumu, 
                        rank = 1:length(pct))
  # Elbow plot to visualize 
  plot <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + theme_bw()
  
  return(plot)
}



