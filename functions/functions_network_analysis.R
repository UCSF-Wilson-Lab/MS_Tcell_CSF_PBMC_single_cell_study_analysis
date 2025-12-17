# Functions for Cytoscape Network Analysis
# PURPOSE: take Single Cell VDJ output from Immcantation and Convert it
#  - INPUT: mAb/mTCR Table

### SQL command to collect the patient information
library(rjson)     
library(stringr)
library(dplyr)
library(foreach)
library(stringdist)
library(igraph)
library(data.table)
library(Matrix)
# register multicore for foreach
library(doMC)


# Support Functions -------------------

### Input
compress_name <- function(lang_name) {  
  short_name <-""
  for (i in strsplit(lang_name, "_")) { short_name<-substr(i, 0,1) }
  return(paste(short_name, collapse=""))
}

### get_data_bcr.R
# find_identical_CDR3_across_VJ
# basically what it does it to find out for a particular CDR3 AA sequence with a particular B cell subset, 
# whether there is more than one VJ gene pair assigned to it. 
# Theoritically, it should not because the chance is very low due to the diversity of B cell.  
# If yes, keep the most abundant VJ and remove the rest.  Finally we get a filtered dataset from the original one.
# 
find_identical_CDR3_across_VJ<-function(data_in_subset){
  #cdr3_list<-as.character(unique(data_in_subset$CDR3))
  data_final<-foreach(data_in_subset_with_same_cdr3 = split(data_in_subset,data_in_subset$CDR3), .combine='rbind') %do% {
    n_distinct_VJ<-n_distinct(data_in_subset_with_same_cdr3$cluster_criteria)
    if (n_distinct_VJ > 1) {
      max_read_count<-max(data_in_subset_with_same_cdr3$count)
      data_majority<-data_in_subset_with_same_cdr3[data_in_subset_with_same_cdr3$count==max_read_count,]
      
      # the first cluster_criteria has max read_count will be used
      # TODO: the frequency of VJ combination used?
      # 1.
      # if the count of multiple VJ combination are equal to max read cont, 
      # should use the most likely VJ combination
      # 2. 
      # or use the sum of all count?
      if (n_distinct(data_majority$cluster_criteria)>1){
        data_majority<-data_majority[data_majority$cluster_criteria==first(data_majority$cluster_criteria),]
      }
      
    }
    else if (n_distinct_VJ==1) {
      data_majority<-data_in_subset_with_same_cdr3
    }
    
    return (data_majority)
  }
  return(data_final)
}


### edge_list_from_distance.R
# find the clusters for each germline usage and build the edge list out of it,
# then combine all the germline together
cluster_by_distance_tcr <-function(data, distance_method, distance_cutoff){
  CDR3<-paste(data$CDR3_AA, data$CDR3_AA_alpha,sep = ":")
  data_matrix<-stringdistmatrix(CDR3,CDR3,method= distance_method)

  # TODO: distance cutoff, distance method
  
  # two cdr3 is connected if the distance less or equal to 1
  data_matrix<-ifelse(data_matrix<=distance_cutoff,1,0)
  
  ## remove self connected edge
  for (i in 1:nrow(data_matrix)){
    data_matrix[i,i]<-0
  }
  
  rownames(data_matrix)<-data$Node_ID
  colnames(data_matrix)<-data$Node_ID
  
  
  ## format the distance matrix into the adjacency matrix 
  data_adj<-graph.adjacency(data_matrix,mode=c("undirected"))
  
  ## get the edge list from the adjacency matrix
  
  data.edgelist<-as.data.frame(get.edgelist(data_adj))
  
  return(data.edgelist)
}


# Support Functions -------------------

# Collapse mTCR table to 1 line per cell (alpha and Beta)
formatTCRmabTable <- function(mab_table) {
  # Subset only keep cells with Alpha and Beta
  mab_table <- mab_table[mab_table$chain_count == 2,]
  mab_table$vgene_imm <- tstrsplit(mab_table$vgene_imm,"\\*")[[1]]
  mab_table$jgene_imm <- tstrsplit(mab_table$jgene_imm,"\\*")[[1]]
  
  mab_table_beta  <- mab_table[mab_table$chain %in% c("TRB"),]
  row.names(mab_table_beta) <- mab_table_beta$igseq_barcode
  mab_table_alpha <- mab_table[mab_table$chain %in% c("TRA"),]
  row.names(mab_table_alpha) <- mab_table_alpha$igseq_barcode
  
  mab_table_alpha <- mab_table_alpha[,c("vgene_imm","jgene_imm","CDR3_AA")]
  names(mab_table_alpha) <- c("vgene_imm_alpha","jgene_imm_alpha","CDR3_AA_alpha")
  mab_table       <- merge(mab_table_beta,mab_table_alpha, by=0,all=TRUE)
  mab_table <- mab_table[order(mab_table[,2]),]
  row.names(mab_table) <- mab_table[,1]
  mab_table[,1] <- NULL
  
  # Add Cluter Criteria Column (Alpha:Beta)
  cluster_criteria_alpha     <- paste(mab_table$vgene_imm_alpha,mab_table$jgene_imm_alpha,mab_table$CDR3_AA_alpha,sep="_")
  cluster_criteria_beta      <- paste(mab_table$vgene_imm,mab_table$jgene_imm,mab_table$CDR3_AA,sep="_")
  mab_table$cluster_criteria <- paste(cluster_criteria_alpha,cluster_criteria_beta, sep=":")
  
  # Add Node ID for network analysis
  mab_table$Node_ID <-paste(mab_table$sample_igseq,mab_table[,1],sep="_")
  
  return(mab_table)
}



# Clonal Connection Network Functions-------------------

### BCR mAb to GML
createBCRnetworkGML <- function(
  raw_csv_file    = "filtered_contigs_ALLpatients.forSeurat.csv",
  result_dir      = "./Results_Heavy",
  length_cutoff   = as.numeric(8),
  count_cutoff    = as.numeric(1),
  distance_cutoff = as.numeric(2),
  dbname          = "bcr",
  method          = "hamming",
  chain           = "heavy",
  cores           = 12
){
  registerDoMC(cores=12)
  print(Sys.time())
  ptm <- proc.time()
  
  
  ### INPUT/OUTPUT VARIABLES
  bcr_or_tcr <- dbname
  reads<-paste('mixcr', dbname,sep ='_')
  
  config_suffix <-paste(method, distance_cutoff, 'len', length_cutoff, 'count', count_cutoff, sep='_')
  config_suffix <- compress_name(config_suffix)
  
  # Important files for storing output
  data_for_clustering <- file.path(result_dir,paste0('cleaned_data_with_majority_vote_', config_suffix, '.txt'))
  # edge_list_from_distance.R OUTPUT
  overall_cluster_edge_list_file <- file.path(result_dir, paste0("final_edge_list_of_all_", config_suffix , ".txt"))
  # write_all_clusters.R OUTPUT
  all_clusters_gml_file <- file.path(result_dir, paste0("cleaned_clusters_all_", config_suffix, ".gml"))
  all_clusters_csv_file <-file.path(result_dir, paste0("cleaned_clusters_all_", config_suffix ,".csv"))
  
  ### 1. GET DATA (get_data_bcr.R) -----
  # create target dir
  if (!file.exists(result_dir)) {dir.create(result_dir, recursive=TRUE)}
  data<-read.csv(raw_csv_file, header = TRUE,stringsAsFactors=FALSE)
  
  raw_summary_of_sequence_cout <- data %>% group_by(SUBSET) %>% summarise(count = sum(count)) %>% collect()
  raw_summary_of_sequence_cout_path <- file.path(result_dir,'sequence_count_summary_raw.txt')
  write.csv(raw_summary_of_sequence_cout, raw_summary_of_sequence_cout_path,row.names = F)
  
  # clean date
  # filter
  data <- data[data$count>=count_cutoff,]
  # Paired
  data_paired <- data
  # IGH:
  data_igh <- data[str_count(data$vgene,'IGHV')==1 & str_count(data$jgene,'IGHJ') ==1,]
  # IGK:
  data_igk <- data[str_count(data$vgene,'IGKV')==1 & str_count(data$jgene,'IGKJ') ==1,]
  # IGL:
  data_igl <- data[str_count(data$vgene,'IGLV')==1 & str_count(data$jgene,'IGLJ') ==1,]
  
  # Currently only Heavy Chain based clonotypes are supported
  if (chain == "heavy"){
    data <- data_igh
  }else {
    data <- data_igh
  }
  
  # merge row with same vgene, jgene, cdr3, and subset
  # Original way counts were collapsed for sequencing data
  # modification for single cell
  data<- (data %>% group_by(vgene, jgene, CDR3, SUBSET, seurat_converted_cell_barcode) %>% summarise(count = sum(count)) %>% collect())
  
  data$cluster_criteria<-paste(data$vgene,data$jgene,sep="_")
  data$Node_ID<-paste(data$SUBSET,rownames(data),sep="_")
  
  # decide the compartment and subset from barcode
  
  #  Example Subset Format:
  barcode_splited <- strsplit(data$SUBSET,'_')
  patient_id_vec <- sapply(barcode_splited, '[', 1)
  compartment_vec <- sapply(barcode_splited, '[', 3)
  
  # careful, the order of the following three statement matters
  # check whole sample name for CSF 
  data$compartment<-ifelse(substr(compartment_vec,1,3) == 'CSF','CSF','PB')
  data$patient = patient_id_vec
  data$subset<- data$SUBSET
  
  
  data_clean <-foreach(data_in_subset = split(data,data$subset), .combine='rbind') %dopar% {
    find_identical_CDR3_across_VJ(data_in_subset)
  }
  
  clean_summary_of_sequence_cout <- data_clean %>% group_by(SUBSET) %>% summarise(count = sum(count)) %>% collect()
  clean_summary_of_sequence_cout_path <- file.path(result_dir,'sequence_count_summary_clean.txt')
  write.csv(clean_summary_of_sequence_cout, clean_summary_of_sequence_cout_path,row.names = F)
  print(clean_summary_of_sequence_cout)
  write.table(data_clean, data_for_clustering, sep='\t')
  
  
  ### 2. EDGES -----
  edge_list_final <-foreach(data_of_cluster= split(data,data$cluster_criteria), .combine='rbind') %dopar% {
    cluster_by_distance(data_of_cluster,method, distance_cutoff)
  }
  # ? sort list?
  write.table(edge_list_final, overall_cluster_edge_list_file ,sep="\t")
  
  
  ### 3. WRITE CLUSTERS -----
  subsets <- unique(data$subset)
  # do quality control here
  
  #rearrange the column, so first column is node_id
  
  ### ORIGINAL method: for bulk sequencing data
  ### Updated this line for single cell
  vertex<-data[,c("Node_ID","cluster_criteria", "compartment","patient","subset", "vgene","jgene","CDR3","SUBSET","seurat_converted_cell_barcode","count")]
  
  graph_all<- graph.data.frame(edge_list_final, directed=F,vertices=vertex)
  
  # Calculate the maximal (weakly or strongly) connected components of a graph
  vertex$cluster_membership<-as.factor(clusters(graph_all)$membership)
  
  graph_all<- graph.data.frame(edge_list_final, directed=F,vertices=vertex)
  
  write.graph(graph_all, all_clusters_gml_file , format = "gml")
  write.csv(vertex, all_clusters_csv_file )
  
  
  paste(commandArgs(), collapse = " ")
  print(proc.time() - ptm)
}




### TCR mAb table to GML file for Cytoscape
createTCRnetworkGML <- function(
  data_tcr,
  result_dir      = "./Results",
  output_name     = "cleaned_clusters_all",
  count_cutoff    = as.numeric(2), # UMI count cut-off per contig
  distance_cutoff = as.numeric(0),
  method          = "hamming",
  cores           = 12
){
  registerDoMC(cores=cores)
  print(Sys.time())
  ptm <- proc.time()
  
  # Important files for storing output
  overall_cluster_edge_list_file <- file.path(result_dir, paste0("network_edge_list.txt"))
  all_clusters_gml_file <- file.path(result_dir, paste0(output_name, ".gml"))
  all_clusters_csv_file <-file.path(result_dir, paste0(output_name ,".csv"))
  
  ### 1. FORMAT DATA -----
  # create target dir
  if (!file.exists(result_dir)) {dir.create(result_dir, recursive=TRUE)}

  # UMI filter
  data_tcr <- data_tcr[data_tcr$UMIcount>=count_cutoff,]
  
  # Format dataframe so alpha and Beta is on one line
  data_tcr <- formatTCRmabTable(data_tcr)
  
  ### 2. EDGES -----
  edge_list_final <-foreach(data_of_cluster= split(data_tcr,data_tcr$cluster_criteria), .combine='rbind') %dopar% {
    cluster_by_distance_tcr(data_of_cluster,method, distance_cutoff)
  }
  # ? sort list?
  write.table(edge_list_final, overall_cluster_edge_list_file ,sep="\t")
  
  
  ### 3. WRITE CLUSTERS -----

  ### Keep Relavent Columns: Updated this line for single cell
  cols_to_keep <- c("Node_ID","cluster_criteria", "loc","Patient_ID","Paper_ID_Diagnosis_loc", 
                    "MERGED_clonotypeID_imm","vgene_imm","jgene_imm","CDR3_AA","vgene_imm_alpha","jgene_imm_alpha","CDR3_AA_alpha",
                    "subtype.cell.annot.manual","igseq_barcode","seurat_converted_cell_barcode",
                    "clonotype_cell_count_imm","clonotype_CSF_cell_count_imm")
  vertex<-data_tcr[,cols_to_keep]
  vertex$Diagnosis <- tstrsplit(vertex$Paper_ID_Diagnosis_loc,"_")[[2]]
  
  graph_all<- graph.data.frame(edge_list_final, directed=F,vertices=vertex)
  
  write.graph(graph_all, all_clusters_gml_file , format = "gml")
  write.csv(vertex, all_clusters_csv_file )
  
  
  paste(commandArgs(), collapse = " ")
  print(proc.time() - ptm)
}


# GLIPH Network Functions -------------

### TCR mAb table to GML file for Cytoscape
#   - edges are based on GLIPH interaction results
createGLIPHnetworkGML <- function(
  data_tcr,
  result_dir      = "./Results",
  output_name     = "cleaned_clusters_all",
  count_cutoff    = as.numeric(2), # UMI count cut-off per contig
  distance_cutoff = as.numeric(0),
  method          = "hamming",
  gliphCountColumn   = "gliph_csf_interac_count",
  gliphInteracColumn = "gliph_csf_interac",
  cores           = 12
){
  registerDoMC(cores=cores)
  print(Sys.time())
  ptm <- proc.time()
  
  # Important files for storing output
  overall_cluster_edge_list_file <- file.path(result_dir, paste0("network_edge_list.txt"))
  all_clusters_gml_file <- file.path(result_dir, paste0(output_name, ".gml"))
  all_clusters_csv_file <-file.path(result_dir, paste0(output_name ,".csv"))
  
  ### 1. FORMAT DATA -----
  # create target dir
  if (!file.exists(result_dir)) {dir.create(result_dir, recursive=TRUE)}
  
  # UMI filter
  data_tcr <- data_tcr[data_tcr$UMIcount>=count_cutoff,]
  
  # Format dataframe so alpha and Beta is on one line
  data_tcr <- formatTCRmabTable(data_tcr)
  
  ### 2. EDGES -----
 
  ## a. Subset GLIPH results
  interac_col       <- grep(gliphInteracColumn,names(data_tcr))
  interac_count_col <- grep(gliphCountColumn,names(data_tcr))
  data_tcr_gliph    <- data_tcr[data_tcr[,names(data_tcr) %in% gliphCountColumn] > 0,]
  interac_df <- data_tcr_gliph[,c("Node_ID","CDR3_AA","CDR3_AA_alpha",gliphInteracColumn)]
  
  ## b. Create edge list data frame
  edge_list_final <- data.frame(matrix(nrow = 0,ncol = 2))
  for (i in 1:nrow(interac_df)) {
    one_row  <- interac_df[i,]
    interac1 <- one_row$Node_ID
    
    # Find list of all interacting nodes
    interac_cdr3_vec <- one_row[,names(one_row) %in% gliphInteracColumn]
    if(grepl(":",interac_cdr3_vec)){
      interac_cdr3_vec <- as.character(unlist(tstrsplit(interac_cdr3_vec, ":")))
      interac_cdr3_vec <- interac_cdr3_vec[! interac_cdr3_vec %in% c("")]
    }
    
    # list all pair interactions for the row
    for (cdr3 in interac_cdr3_vec) {
      nodes_beta  <- interac_df[interac_df$CDR3_AA %in% cdr3,]
      nodes_alpha <- interac_df[interac_df$CDR3_AA_alpha %in% cdr3,]
      interac2    <- unique(c(nodes_beta$Node_ID,nodes_alpha$Node_ID))
      
      # Skip if interaction is not present
      if(length(interac2) == 0){next}
      
      entry       <- data.frame(matrix(nrow = length(interac2), ncol = 2))
      entry$X1 <- interac1
      entry$X2 <- interac2
      edge_list_final <- rbind(edge_list_final,entry)
    }
  }
  
  
  
  write.table(edge_list_final, overall_cluster_edge_list_file ,sep="\t")
  
  
  ### 3. WRITE CLUSTERS -----
  
  ### Keep Relavent Columns: Updated this line for single cell
  cols_to_keep <- c("Node_ID","cluster_criteria", "loc","Patient_ID","Paper_ID_Diagnosis_loc", 
                    "MERGED_clonotypeID_imm",gliphCountColumn,gliphInteracColumn,"clonotype_cell_count_imm","clonotype_CSF_cell_count_imm",
                    "seurat_converted_cell_barcode","igseq_barcode","vgene_imm","jgene_imm","CDR3_AA","vgene_imm_alpha",
                    "jgene_imm_alpha","CDR3_AA_alpha","subtype.cell.annot.manual")
  vertex<-data_tcr[,cols_to_keep]
  vertex$Diagnosis <- tstrsplit(vertex$Paper_ID_Diagnosis_loc,"_")[[2]]
  
  graph_all<- graph.data.frame(edge_list_final, directed=F,vertices=vertex)
  
  write.graph(graph_all, all_clusters_gml_file , format = "gml")
  write.csv(vertex, all_clusters_csv_file )
  
  
  paste(commandArgs(), collapse = " ")
  print(proc.time() - ptm)
}



### TCR mAb table to GML file for Cytoscape
#   - edges are based on GLIPH interaction results
createGLIPHclonotypeNetworkGML <- function(
  data_tcr,
  result_dir      = "./Results",
  output_name     = "cleaned_clusters_all",
  count_cutoff    = as.numeric(2), # UMI count cut-off per contig
  distance_cutoff = as.numeric(0),
  method          = "hamming",
  gliphCountColumn   = "gliph_csf_interac_count",
  gliphInteracColumn = "gliph_csf_interac",
  cores           = 12
){
  registerDoMC(cores=cores)
  print(Sys.time())
  ptm <- proc.time()
  
  # Important files for storing output
  overall_cluster_edge_list_file <- file.path(result_dir, paste0("network_edge_list.txt"))
  all_clusters_gml_file <- file.path(result_dir, paste0(output_name, ".gml"))
  all_clusters_csv_file <-file.path(result_dir, paste0(output_name ,".csv"))
  
  ### 1. FORMAT DATA -----
  # create target dir
  if (!file.exists(result_dir)) {dir.create(result_dir, recursive=TRUE)}
  
  # UMI filter
  data_tcr <- data_tcr[data_tcr$UMIcount>=count_cutoff,]
  
  # Format dataframe so alpha and Beta is on one line
  data_tcr <- formatTCRmabTable(data_tcr)
  
  # collapse to clonotype resolution rather than cell resolution
  clone_list <- unique(data_tcr$MERGED_clonotypeID_imm)
  # -- Remove columns that are not unique (cell barcodes)
  #    Node ID now becomes Clonotype ID
  
  cols_to_omit <- c()
  
  for (clone in clone_list) {
    clone_df <- data_tcr[data_tcr$MERGED_clonotypeID_imm %in% clone,]
    for (col in names(clone_df)) {
      one_col <- clone_df[,col]
      uni_elements <- length(unique(one_col))
      if(uni_elements > 1){
        cols_to_omit <- c(cols_to_omit,col)
        if (col == "CDR3_AA_alpha"){
          cat(paste0(clone,"\n")) # sanity check 
        }
      }
    }
  }
  cols_to_omit <- unique(cols_to_omit)
  
  # Collapse data frame to one line per clonotype
  data_tcr_fmt <- data_tcr[,! names(data_tcr) %in% cols_to_omit]
  
  data_tcr_fmt <- unique(data_tcr_fmt)
  row.names(data_tcr_fmt) <- data_tcr_fmt$MERGED_clonotypeID_imm
  data_tcr_fmt$Node_ID <- data_tcr_fmt$MERGED_clonotypeID_imm
  
  ### 2. EDGES -----
  
  ## a. Subset GLIPH relavent results
  interac_col       <- grep(gliphInteracColumn,names(data_tcr_fmt))
  interac_count_col <- grep(gliphCountColumn,names(data_tcr_fmt))
  data_tcr_gliph    <- data_tcr_fmt[data_tcr_fmt[,names(data_tcr_fmt) %in% gliphCountColumn] > 0,]
  interac_df <- data_tcr_gliph[,c("Node_ID","CDR3_AA","CDR3_AA_alpha",gliphInteracColumn)]
  
  ## b. Create edge list data frame
  edge_list_final <- data.frame(matrix(nrow = 0,ncol = 2))
  for (i in 1:nrow(interac_df)) {
    one_row  <- interac_df[i,]
    interac1 <- one_row$Node_ID
    
    # Find list of all interacting nodes
    interac_cdr3_vec <- one_row[,names(one_row) %in% gliphInteracColumn]
    if(grepl(":",interac_cdr3_vec)){
      interac_cdr3_vec <- as.character(unlist(tstrsplit(interac_cdr3_vec, ":")))
      interac_cdr3_vec <- interac_cdr3_vec[! interac_cdr3_vec %in% c("")]
    }
    
    # list all pair interactions for the row
    for (cdr3 in interac_cdr3_vec) {
      nodes_beta  <- interac_df[interac_df$CDR3_AA %in% cdr3,]
      nodes_alpha <- interac_df[interac_df$CDR3_AA_alpha %in% cdr3,]
      interac2    <- unique(c(nodes_beta$Node_ID,nodes_alpha$Node_ID))
      
      # Skip if interaction is not present
      if(length(interac2) == 0){next}
      
      entry       <- data.frame(matrix(nrow = length(interac2), ncol = 2))
      entry$X1 <- interac1
      entry$X2 <- interac2
      edge_list_final <- rbind(edge_list_final,entry)
    }
  }
  
  
  
  write.table(edge_list_final, overall_cluster_edge_list_file ,sep="\t")
  
  
  ### 3. WRITE CLUSTERS -----
  
  ### Keep Relavent Columns: Updated this line for single cell
  cols_to_keep <- c("Node_ID","cluster_criteria", "loc","Patient_ID","Paper_ID_Diagnosis_loc", 
                    "MERGED_clonotypeID_imm",gliphCountColumn,gliphInteracColumn,"clonotype_cell_count_imm","clonotype_CSF_cell_count_imm",
                    "vgene_imm","jgene_imm","CDR3_AA","vgene_imm_alpha",
                    "jgene_imm_alpha","CDR3_AA_alpha")
  vertex<-data_tcr_fmt[,cols_to_keep]
  vertex$Diagnosis <- tstrsplit(vertex$Paper_ID_Diagnosis_loc,"_")[[2]]
  
  graph_all<- graph.data.frame(edge_list_final, directed=F,vertices=vertex)
  
  write.graph(graph_all, all_clusters_gml_file , format = "gml")
  write.csv(vertex, all_clusters_csv_file )
  
  
  paste(commandArgs(), collapse = " ")
  print(proc.time() - ptm)
}


# GLIPH 2 Network Function ---------------

createGLIPH2clonotypeNetworkGML <- function(
  data_tcr,
  result_dir      = "./Results",
  output_name     = "cleaned_clusters_all",
  count_cutoff    = as.numeric(2), # UMI count cut-off per contig
  distance_cutoff = as.numeric(0),
  method          = "hamming",
  gliphStatusColumn   = "gliph_status",
  gliphInteracColumn = "gliph_group",
  cores           = 12
){
  registerDoMC(cores=cores)
  print(Sys.time())
  ptm <- proc.time()
  
  # Important files for storing output
  overall_cluster_edge_list_file <- file.path(result_dir, paste0("network_edge_list.txt"))
  all_clusters_gml_file <- file.path(result_dir, paste0(output_name, ".gml"))
  all_clusters_csv_file <-file.path(result_dir, paste0(output_name ,".csv"))
  
  ### 1. FORMAT DATA -----
  # create target dir
  if (!file.exists(result_dir)) {dir.create(result_dir, recursive=TRUE)}
  
  # UMI filter
  data_tcr <- data_tcr[data_tcr$UMIcount>=count_cutoff,]
  
  # Format dataframe so alpha and Beta is on one line
  data_tcr <- formatTCRmabTable(data_tcr)
  clone_list <- unique(data_tcr$MERGED_clonotypeID_imm)
  
  ### ----Add clonotype-Cell type column
  simple_annot <- data_tcr$subtype.cell.annot.manual
  names(simple_annot) <- data_tcr$MERGED_clonotypeID_imm
  for (clone in clone_list) {
    clone_sub <- simple_annot[names(simple_annot) %in% clone]
    uni_elements <- length(unique(clone_sub))
    if(uni_elements > 1){
      cd4_cd8 <- paste(as.character(unique(clone_sub)),collapse = " & ")
      simple_annot[names(simple_annot) %in% clone] <- cd4_cd8
    }
  }
  data_tcr$simple_annot <- simple_annot
  
  
  # collapse to clonotype resolution rather than cell resolution
  # -- Remove columns that are not unique (cell barcodes)
  #    Node ID now becomes Clonotype ID
  
  cols_to_omit <- c()
  
  for (clone in clone_list) {
    clone_df <- data_tcr[data_tcr$MERGED_clonotypeID_imm %in% clone,]
    for (col in names(clone_df)) {
      one_col <- clone_df[,col]
      uni_elements <- length(unique(one_col))
      if(uni_elements > 1){
        cols_to_omit <- c(cols_to_omit,col)
        if (col == "CDR3_AA_alpha"){
          cat(paste0(clone,"\n")) # sanity check 
        }
      }
    }
  }
  cols_to_omit <- unique(cols_to_omit)
  
  # Collapse data frame to one line per clonotype
  data_tcr_fmt <- data_tcr[,! names(data_tcr) %in% cols_to_omit]
  
  data_tcr_fmt <- unique(data_tcr_fmt)
  row.names(data_tcr_fmt) <- data_tcr_fmt$MERGED_clonotypeID_imm
  data_tcr_fmt$Node_ID <- data_tcr_fmt$MERGED_clonotypeID_imm
  
  ### 2. EDGES -----
  
  ## a. Subset GLIPH relavent results
  interac_col        <- grep(gliphInteracColumn,names(data_tcr_fmt))
  interac_status_col <- grep(gliphStatusColumn,names(data_tcr_fmt))
  data_tcr_gliph     <- data_tcr_fmt[data_tcr_fmt[,names(data_tcr_fmt) %in% gliphStatusColumn] %in% "1",]
  interac_df <- data_tcr_gliph[,c("Node_ID","CDR3ab_AA",gliphInteracColumn)]
  
  ## b. Create edge list data frame
  edge_list_final <- data.frame(matrix(nrow = 0,ncol = 2))
  for (i in 1:nrow(interac_df)) {
    one_row  <- interac_df[i,]
    interac1    <- one_row$Node_ID
    interac1_ab <- one_row$CDR3ab_AA
    
    # Find list of all interacting nodes
    interac_grp     <- one_row[,names(one_row) %in% gliphInteracColumn]
    interac_grp_vec <- c()
    if(grepl(":",interac_grp)){
      interac_grp_vec <- as.character(unlist(strsplit(interac_grp, ":")))
    }
    interac_grp_vec <- c(interac_grp_vec,interac_grp)
    
    
    # list all pair interactions for the row
    for (grp in interac_grp_vec) {
      nodes_grp  <- interac_df[interac_df$gliph_group %in% grp,]
      
      # Exclude Interaction 1 from this subset
      nodes_grp   <- nodes_grp[! nodes_grp$CDR3ab_AA %in% interac1_ab,]
      interac2    <- unique(nodes_grp$Node_ID)
      
      # Skip if interaction is not present
      if(length(interac2) == 0){next}
      
      entry       <- data.frame(matrix(nrow = length(interac2), ncol = 2))
      entry$X1 <- interac1
      entry$X2 <- interac2
      edge_list_final <- rbind(edge_list_final,entry)
    }
  }
  
  
  
  write.table(edge_list_final, overall_cluster_edge_list_file ,sep="\t")
  
  
  ### 3. WRITE CLUSTERS -----
  
  ### Keep Columns: Updated this line for single cell
  cols_to_keep <- c("Node_ID","cluster_criteria", "loc","Patient_ID","Paper_ID_Diagnosis_loc","simple_annot", 
                    "MERGED_clonotypeID_imm",gliphStatusColumn,gliphInteracColumn,"clonotype_cell_count_imm","clonotype_CSF_cell_count_imm",
                    "CDR3_AA","CDR3_AA_alpha","clonotype_CSF_freq",
                    "gliph_hla_score","gliph_final_score","gliph_motif_type")
  vertex<-data_tcr_fmt[,cols_to_keep]
  vertex$Diagnosis <- tstrsplit(vertex$Paper_ID_Diagnosis_loc,"_")[[2]]
  
  graph_all<- graph.data.frame(edge_list_final, directed=F,vertices=vertex)
  
  write.graph(graph_all, all_clusters_gml_file , format = "gml")
  write.csv(vertex, all_clusters_csv_file )
  
  
  paste(commandArgs(), collapse = " ")
  print(proc.time() - ptm)
}
