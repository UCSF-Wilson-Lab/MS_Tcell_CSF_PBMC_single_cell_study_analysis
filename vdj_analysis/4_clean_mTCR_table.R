#!/usr/bin/env Rscript

# PURPOSE: filter TCR mAb table
#   - remove cells with incomplete TCR assemblies (particularly cells with incomplete CDR3)
#   - Include MIXCR results for cells it was able to rescue
#   - When including MIXCR results, double check whether clonotype assigment needs to be updated for these cells
library(dplyr)
library(data.table)
library(devtools)
library(ggplot2)
library(stringr)
# contains some MIXCR functions
working_dir <- "./sc_analysis/"
setwd(working_dir)
source("Resources/functions_sc_analysis.R")
source("Resources/functions_network_analysis.R")

## INPUT
current_mab_table <- "./tables/mtcr_table_ALL_TCR.v3.1.csv"
mab_table         <- read.csv(current_mab_table)

incomplete_tcrs_fh <- "./tables/rerun_MIXCR_on_incomplete_TCRs.tsv"
results_incomplete_df <- read.delim(incomplete_tcrs_fh)

## OUTPUT
output_mab_table_fh     <- "./tables/mtcr_table_ALL_TCR.v3.1.CLEAN.csv"
output_mab_table_fh_fmt <- "./tables/fmtVJ_mtcr_table_ALL_TCR.v3.1.csv"
output_dir_clonotype_counts <- "./tables/patient_clonotype_pcts/"
dir.create(output_dir_clonotype_counts,recursive = T)

### 1. Figure out which cells to keep
total_cells  <- length(as.character(unique(results_incomplete_df$igseq_barcode)))

results_na   <- results_incomplete_df[is.na(results_incomplete_df$aaSeqCDR3),]
cells_failed <- as.character(unique(results_na$igseq_barcode))
rescued_df   <- results_incomplete_df[! results_incomplete_df$igseq_barcode %in% cells_failed,]
rescued_df$CDR3_length <- str_count(rescued_df$aaSeqCDR3)
df_sub <- rescued_df[rescued_df$CDR3_length > 21,]

# omit cells with questionable CDR3s (* or _ present)
ambig1      <- rescued_df[grep("\\*",rescued_df$aaSeqCDR3),]
ambig1      <- as.character(unique(ambig1$igseq_barcode))

ambig2      <- rescued_df[grep("\\_",rescued_df$aaSeqCDR3),]
ambig2      <- as.character(unique(ambig2$igseq_barcode))

cells_to_omit <- c(ambig1,ambig2)
rescued_df <- rescued_df[! rescued_df$igseq_barcode %in% cells_to_omit,]

# omit Cells with CDR3 length less than 4 (only 2 cells)
rescued_df    <- rescued_df[rescued_df$CDR3_length > 3,]
rescued_df <- rescued_df %>% mutate_all(as.character)
rescued_cells <- as.character(unique(rescued_df$igseq_barcode))

### 2. Parse mAb table - omit all cells with incomplete TCRs except for cells completed by MIXCR

# Identify all cells with at least one NA for a CDR3 AA
mab_results_na   <- mab_table[is.na(mab_table$CDR3_AA),]
mab_cells_failed <- as.character(unique(mab_results_na$igseq_barcode))
cells_to_omit    <- mab_cells_failed[! mab_cells_failed %in% rescued_cells]


mab_table_clean <- mab_table[! mab_table$igseq_barcode %in% cells_to_omit,]
# Convert all columns to char
mab_table_clean <- mab_table_clean %>% mutate_all(as.character)

# Add in MIXCR results
extra_cells_to_omit <- c()

for (cell in rescued_cells) {
  mab_sub <- mab_table_clean[mab_table_clean$igseq_barcode %in% cell,]
  mab_beta  <- mab_sub[mab_sub$chain %in% "TRB",] 
  mab_alpha <- mab_sub[mab_sub$chain %in% "TRA",]
  
  rescued_sub <- rescued_df[rescued_df$igseq_barcode %in% cell,]
  if(nrow(rescued_sub) < 2){
    extra_cells_to_omit<- c(extra_cells_to_omit,cell)
    next
  }
  
  beta_mixcr  <- rescued_sub[rescued_sub$chain %in% "TRB",]
  alpha_mixcr <- rescued_sub[rescued_sub$chain %in% "TRA",]
  
  # Convert FRW and CDR results to MIXCR results
  
  # Beta
  mab_beta$FRW1 <- beta_mixcr$nSeqFR1
  mab_beta$FRW2 <- beta_mixcr$nSeqFR2
  mab_beta$FRW3 <- beta_mixcr$nSeqFR3
  mab_beta$FRW4 <- beta_mixcr$nSeqFR4
  mab_beta$CDR1 <- beta_mixcr$nSeqCDR1
  mab_beta$CDR2 <- beta_mixcr$nSeqCDR2
  mab_beta$CDR3 <- beta_mixcr$nSeqCDR3
  mab_beta$CDR3_AA <- beta_mixcr$aaSeqCDR3
  mab_sub[mab_sub$chain %in% "TRB",] <- mab_beta
  # Alpha
  mab_alpha$FRW1 <- alpha_mixcr$nSeqFR1
  mab_alpha$FRW2 <- alpha_mixcr$nSeqFR2
  mab_alpha$FRW3 <- alpha_mixcr$nSeqFR3
  mab_alpha$FRW4 <- alpha_mixcr$nSeqFR4
  mab_alpha$CDR1 <- alpha_mixcr$nSeqCDR1
  mab_alpha$CDR2 <- alpha_mixcr$nSeqCDR2
  mab_alpha$CDR3 <- alpha_mixcr$nSeqCDR3
  mab_alpha$CDR3_AA <- alpha_mixcr$aaSeqCDR3
  mab_sub[mab_sub$chain %in% "TRA",] <- mab_alpha
  
  # Finalize cell conversion
  mab_table_clean[mab_table_clean$igseq_barcode %in% cell,] <- mab_sub
}

# Omit any lingering cells that did not have a full MIXCR assembly
mab_table_clean <- mab_table_clean[! mab_table_clean$igseq_barcode %in% extra_cells_to_omit,]

# Add new column to distiguish what alignment method was used
new_col_method        <- rep("IgBLAST",nrow(mab_table_clean))
names(new_col_method) <- mab_table_clean$igseq_barcode

new_col_method[names(new_col_method) %in% rescued_cells] <- "MIXCR"
mab_table_clean$alignment_method <- new_col_method


### 3. Check if MIXCR assembles adjust the clonotype assignment
#   - ONLY applies for this particular version of the TCR table
mab_table_clean_mixcr <- mab_table_clean[mab_table_clean$alignment_method %in% "MIXCR",]
clonotype_list_mixcr  <- unique(mab_table_clean_mixcr$MERGED_clonotypeID_imm)

for (clone in clonotype_list_mixcr) {
  mab_clone     <- mab_table_clean[mab_table_clean$MERGED_clonotypeID_imm %in% clone,]
  mab_clone_fmt <- formatTCRmabTable(mab_clone) # Collapse to 1 row per cell
  
  clonotype_criteria <-  unique(mab_clone_fmt$cluster_criteria)
  
  # MIXCR caused cells to deviate from the current clonotype definition
  # - only one cell that this happens to 
  # - MIXCR create a new clonotype for 1 cell (PB)
  if(length(clonotype_criteria) > 1){
    mab_clone_mixcr <- mab_clone[mab_clone$alignment_method == "MIXCR",]  # 1 cell
    target_cell     <- unique(mab_clone_mixcr$igseq_barcode)
    
    new_cell_count  <- length(target_cell)
    cell_csf_count  <- unique(mab_clone_mixcr$clonotype_CSF_cell_count_imm)
    cell_pb_count   <- unique(mab_clone_mixcr$clonotype_PB_cell_count_imm)
    
    # Change cell counts and expansion status
    mab_clone_mixcr$clonotype_cell_count_imm <- new_cell_count
    if(cell_csf_count > 0){mab_clone_mixcr$clonotype_CSF_cell_count_imm   <- new_cell_count}
    if(cell_pb_count  > 0){mab_clone_mixcr$clonotype_PB_cell_count_imm    <- new_cell_count}
    if(new_cell_count == 1){mab_clone_mixcr$expanded_clonotype_status_imm <- "0"}
    
    # Change Clonotype ID name
    clonotype_id     <- unique(mab_clone_mixcr$MERGED_clonotypeID_imm)
    new_clonotype_id <- paste(c("new",clonotype_id),collapse = "")
    mab_clone_mixcr$MERGED_clonotypeID_imm <- new_clonotype_id
    
    # Insert changes back into mab table
    mab_table_clean[mab_table_clean$igseq_barcode %in% target_cell,] <- mab_clone_mixcr
  }
}


### 4. Go back and modify clonotype cell counts and expansion status
#    - correct the 3 different cell count per clonotype columns since a lot of cells have been omitted
#    - if Total cell count drops to 1, adjust expansion status

# Convert count columns back to numeric
mab_table_clean$clonotype_cell_count_imm     <- as.numeric(mab_table_clean$clonotype_cell_count_imm)
mab_table_clean$clonotype_CSF_cell_count_imm <- as.numeric(mab_table_clean$clonotype_CSF_cell_count_imm)
mab_table_clean$clonotype_PB_cell_count_imm  <- as.numeric(mab_table_clean$clonotype_PB_cell_count_imm)

# Check all expanded clonotypes
mab_table_clean_sub <- mab_table_clean[mab_table_clean$clonotype_cell_count_imm > 1,]
exp_clonotype_list <- unique(mab_table_clean_sub$MERGED_clonotypeID_imm)

for (clone in exp_clonotype_list) {
  mab_clone <- mab_table_clean[mab_table_clean$MERGED_clonotypeID_imm %in% clone,]
  cell_count <- unique(mab_clone$clonotype_cell_count_imm)
  real_cell_count <- length(unique(mab_clone$igseq_barcode))
  
  # >>> Skip if cell count is unchanged
  if (real_cell_count == cell_count){next}
  #cat(paste0(">> ", clone,"\n"))
  
  # Change expansion status if cell count drops to 1
  exp_status <- unique(mab_clone$expanded_clonotype_status_imm)
  if (real_cell_count <= 1){
    exp_status <- "0"
  }
  
  # Determine CSF and PB count
  csf_count <- 0
  pb_count  <- 0
  mab_clone_csf <- mab_clone[mab_clone$loc %in% "CSF",]
  mab_clone_pb  <- mab_clone[mab_clone$loc %in% "PB",]
  
  if(nrow(mab_clone_csf) > 0){
    csf_count <- length(unique(mab_clone_csf$igseq_barcode))
  }
  if (nrow(mab_clone_pb) > 0){
    pb_count <- length(unique(mab_clone_pb$igseq_barcode))
  }
  
  # Correct all clonotype count number columns for these select clonotypes
  mab_clone$expanded_clonotype_status_imm <- exp_status
  mab_clone$clonotype_cell_count_imm      <- real_cell_count
  mab_clone$clonotype_CSF_cell_count_imm  <- csf_count
  mab_clone$clonotype_PB_cell_count_imm   <- pb_count
  mab_table_clean[mab_table_clean$MERGED_clonotypeID_imm %in% clone,] <- mab_clone
}

### 5. spot check Immcantations clonotype assignements and fix them if they are incorrect
mab_table_clean_fmt <- formatTCRmabTable(mab_table_clean)

# Updated Clonotype column
new_clonotype_col        <- mab_table_clean_fmt$MERGED_clonotypeID_imm
names(new_clonotype_col) <- mab_table_clean_fmt$igseq_barcode

target_cells <- c() # cells where the clonotype was changed
all_clones   <- unique(mab_table_clean_fmt$MERGED_clonotypeID_imm)
for (clone in all_clones) {
  mab_table_sub <- mab_table_clean_fmt[mab_table_clean_fmt$MERGED_clonotypeID_imm %in% clone,]
  criteria_list <- unique(mab_table_sub$cluster_criteria)
  
  if (length(criteria_list) > 1){
    # Rename the Clonotype ID for each different set of criteria
    criteria_number <- 1
    for (criteria in criteria_list) {
      criteria_df <- mab_table_sub[mab_table_sub$cluster_criteria %in% criteria,]
      cell_vec    <- unique(criteria_df$igseq_barcode)
      new_id      <- paste(c(clone,criteria_number),collapse = "_")
      
      target_cells <- c(target_cells,cell_vec)
      new_clonotype_col[names(new_clonotype_col) %in% cell_vec] <- new_id
      criteria_number <- criteria_number + 1
    }
    
  }
}

changed_clonotypes <- new_clonotype_col[target_cells]

current_clonotype_col <- mab_table_clean$MERGED_clonotypeID_imm
names(current_clonotype_col) <- mab_table_clean$igseq_barcode
for (cell in names(changed_clonotypes)) {
  new_id <- as.character(changed_clonotypes[cell])
  current_clonotype_col[names(current_clonotype_col) %in% cell] <- new_id
}

mab_table_clean$MERGED_clonotypeID_imm <- current_clonotype_col

# Modify Cell count columns based on new clonotype IDs
# - New clonotype introduced by MIXCR already corrected for cell counts
target_mab_table   <- mab_table_clean[mab_table_clean$igseq_barcode %in% target_cells,]
new_clonotype_list <-  unique(target_mab_table$MERGED_clonotypeID_imm)

# Initialize updated columns and change the counts
new_cell_count_col     <- mab_table_clean$clonotype_cell_count_imm
new_cell_count_col_CSF <- mab_table_clean$clonotype_CSF_cell_count_imm
new_cell_count_col_PB  <- mab_table_clean$clonotype_PB_cell_count_imm
new_expanded_col       <- mab_table_clean$expanded_clonotype_status_imm
names(new_cell_count_col)     <- mab_table_clean$igseq_barcode
names(new_cell_count_col_CSF) <- mab_table_clean$igseq_barcode
names(new_cell_count_col_PB)  <- mab_table_clean$igseq_barcode
names(new_expanded_col)       <- mab_table_clean$igseq_barcode

for (clone in new_clonotype_list) {
  clone_df <- target_mab_table[target_mab_table$MERGED_clonotypeID_imm %in% clone,]
  clone_cell_vec <- unique(clone_df$igseq_barcode)
  loc_vec        <- unique(clone_df$loc)
  expanded_status  <- "1"
  csf_cell_count   <- 0
  pb_cell_count    <- 0
  total_cell_count <- length(unique(clone_df$igseq_barcode))
  if ("CSF" %in% loc_vec){
    csf_df <- clone_df[clone_df$loc %in% "CSF",]
    csf_cell_count <- length(unique(csf_df$igseq_barcode))
  }
  if ("PB" %in% loc_vec){
    pb_df <- clone_df[clone_df$loc %in% "PB",]
    pb_cell_count <- length(unique(pb_df$igseq_barcode))
  }
  if (total_cell_count < 2) {
    expanded_status <- "0"
  }
  
  # Change cell count/expanded status for new manually annotated clonotypes
  new_cell_count_col[names(new_cell_count_col) %in% clone_cell_vec]         <- total_cell_count
  new_cell_count_col_CSF[names(new_cell_count_col_CSF) %in% clone_cell_vec] <- csf_cell_count
  new_cell_count_col_PB[names(new_cell_count_col_PB) %in% clone_cell_vec]   <- pb_cell_count
  new_expanded_col[names(new_expanded_col) %in% clone_cell_vec]             <- expanded_status
}

mab_table_clean$clonotype_cell_count_imm      <- new_cell_count_col
mab_table_clean$clonotype_CSF_cell_count_imm  <- new_cell_count_col_CSF
mab_table_clean$clonotype_PB_cell_count_imm   <- new_cell_count_col_PB
mab_table_clean$expanded_clonotype_status_imm <- new_expanded_col


### 5. Write Final Output
write.csv(mab_table_clean,file = output_mab_table_fh,quote = F,row.names = F)

# Simplify the V and J usage
mab_table_clean$vgene_imm <- tstrsplit(mab_table_clean$vgene_imm,"\\*")[[1]]
mab_table_clean$jgene_imm <- tstrsplit(mab_table_clean$jgene_imm,"\\*")[[1]]

write.csv(mab_table_clean,file = output_mab_table_fh_fmt,quote = F,row.names = F)


### 6. Output per patient tables for deciding expansion thresholds
pt_vec <- unique(mab_table_clean$Patient_ID)

patient_clone_count_list <- list()

for (pt in pt_vec ) {
  mab_pt <- mab_table_clean[mab_table_clean$Patient_ID %in% pt,]
  # Pull cells from all clonotypes observed in this patient
  pt_clones <- unique(mab_pt$MERGED_clonotypeID_imm)
  mab_target_df <- mab_table_clean[mab_table_clean$MERGED_clonotypeID_imm %in% pt_clones,]
  total_cells <- length(unique(mab_target_df$igseq_barcode))
  
  # CSF only
  mab_target_csf  <- mab_target_df[mab_target_df$loc %in% c("CSF"),]
  total_cells_csf <- length(unique(mab_target_csf$igseq_barcode))
  # PB only
  mab_target_pb  <- mab_target_df[mab_target_df$loc %in% c("PB"),]
  total_cells_pb <- length(unique(mab_target_pb$igseq_barcode))
  
  count_df            <- as.data.frame(matrix(0,nrow = length(pt_clones),ncol = 7))
  names(count_df)     <- c("cell_count","CSF_count","PB_count","pct_total","pct_CSF","pct_PB","CSF_PB_ratio")
  row.names(count_df) <- pt_clones
  
  for (clone in pt_clones) {
    mab_clone        <- mab_target_df[mab_target_df$MERGED_clonotypeID_imm %in% clone,]
    clone_cell_count <- length(unique(mab_clone$igseq_barcode))
    
    clone_csf_count <- 0
    clone_pb_count  <- 0
    mab_clone_csf   <- mab_clone[mab_clone$loc %in% c("CSF"),]
    mab_clone_pb    <- mab_clone[mab_clone$loc %in% c("PB"),]
    # CSF
    if(nrow(mab_clone_csf) >0){
      clone_csf_count <- length(unique(mab_clone_csf$igseq_barcode))
    }
    # PB
    if(nrow(mab_clone_pb) >0){
      clone_pb_count <- length(unique(mab_clone_pb$igseq_barcode))
    }
    
    
    # Calulate percentage of total clones in the patient
    clone_pct     <- (clone_cell_count / total_cells) * 100
    clone_pct_csf <- (clone_csf_count / total_cells_csf) * 100
    clone_pct_pb  <- (clone_pb_count / total_cells_pb) * 100
    
    # Calculate CSF:PB ratio
    clone_pct_pb_ratio <- clone_pct_pb
    if(clone_pb_count == 0){
      clone_pct_pb_ratio  <- (1 / total_cells_pb) * 100
    }
    csf_pb_ratio <- clone_pct_csf / clone_pct_pb_ratio
    
    # Add clonotype results to data frame
    results_row <- count_df[row.names(count_df) %in% clone,]
    results_row$pct_total    <- clone_pct
    results_row$pct_CSF      <- clone_pct_csf
    results_row$pct_PB       <- clone_pct_pb
    results_row$cell_count   <- clone_cell_count
    results_row$CSF_count    <- clone_csf_count
    results_row$PB_count     <- clone_pb_count
    results_row$CSF_PB_ratio <- csf_pb_ratio
    
    count_df[row.names(count_df) %in% clone,] <- results_row
  }
  
  patient_clone_count_list[[pt]] <- count_df
}

## Output these tables
#  - per patient tables
#  - also make a combined table with all patients

combined_pt_df        <- as.data.frame(matrix(nrow = 0,ncol = 8))
names(combined_pt_df) <- c("cell_count","CSF_count","PB_count","pct_total","pct_CSF","pct_PB","CSF_PB_ratio","Patient")

for (pt in names(patient_clone_count_list)) {
  pt_df <- patient_clone_count_list[[pt]]
  output_fh <- paste(c(output_dir_clonotype_counts,pt,".csv"),collapse = "")
  
  # Output individual patient files
  write.csv(pt_df,file = output_fh,quote = F,row.names = T)
  
  pt_df$Patient <- pt
  combined_pt_df <- rbind(combined_pt_df,pt_df)
  
}

# Output table with all patients
output_all_fh <- paste(c(output_dir_clonotype_counts,"ALL_PATIENTS_clonotype_counts.csv"),collapse = "")
write.csv(combined_pt_df,file = output_all_fh,quote = F,row.names = T)

cat(paste("\n\n>>> DONE <<<\n\n"))

