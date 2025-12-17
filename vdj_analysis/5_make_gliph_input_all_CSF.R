#!/usr/bin/env Rscript

library(data.table)
library(stringr)
library(dplyr)

# GLIPH and GLIPH2 results directory
working_dir <- "./sc_analysis/"
setwd(working_dir)

results_dir_v2 <- "./gliph_results/all_CSF_GLIPH2/"
dir.create(results_dir_v2,recursive = T)

output_tsv_v2   <- paste(c(results_dir_v2,"TCR_input_table_gliph.allCSF.GLIPH2.tsv"),collapse = "")
output_cdr3b_v2 <- paste(c(results_dir_v2,"ref_CSF.txt"),collapse = "")
output_trb_v2   <- paste(c(results_dir_v2,"ref_CSF_V.txt"),collapse = "")
output_len_v2   <- paste(c(results_dir_v2,"ref_CSF_L.txt"),collapse = "")

# Formatted table - with cell type annotations
mab_table <- read.csv("./tables/mtcr_table_ALL_TCR.v3.1.CLEAN.csv",stringsAsFactors = F)
mab_table <- mab_table[mab_table$chain_count != 1,]
mab_table <- mab_table[! mab_table$subtype.cell.annot.manual %in% "No Annot",]  # Make sure there are only cells that overlap with RNA-Seq

### All CSF
mab_table_sub <- mab_table[mab_table$loc %in% c("CSF"),]
# Add a column for patient:clonotype
mab_table_sub$pt_clonotype <- paste(mab_table_sub$Patient_ID,mab_table_sub$MERGED_clonotypeID_imm,sep = ":")

### 2. Create empty GLIPH data frame
pt_clonotype_list <- unique(mab_table_sub$pt_clonotype)

# Column names: CDR3b		TRBV	TRBJ	CDR3a		TRAV		TRAJ	PatientCounts
input_df <- data.frame(matrix(nrow = length(pt_clonotype_list),ncol = 8))
names(input_df) <- c("CDR3b","TRBV", "TRBJ", "CDR3a", "TRAV", "TRAJ", "Patient","Counts")

row.names(input_df) <- pt_clonotype_list

### 3. Fill Input data frame with Alpha and Beta info from Immcantation Results 
for (clone in pt_clonotype_list) {
  clone_df   <- mab_table_sub[mab_table_sub$pt_clonotype %in% clone,]
  chains     <- unique(clone_df$chain)
  cell_count <- length(unique(clone_df$igseq_barcode))
  
  if(length(chains) != 2){next}
  
  tra_row <- clone_df[clone_df$chain %in% c("TRA"),]
  trb_row <- clone_df[clone_df$chain %in% c("TRB"),]
  
  # Populate input DF
  input_row <- input_df[clone,]
  
  # Beta
  input_row$CDR3b <- unique(trb_row$CDR3_AA)
  vgene           <- unique(tstrsplit(as.character(trb_row$vgene_imm),"\\*")[[1]])
  jgene           <- unique(tstrsplit(as.character(trb_row$jgene_imm),"\\*")[[1]])
  input_row$TRBV  <- vgene
  input_row$TRBJ  <- jgene
  
  # Alpha
  input_row$CDR3a <- unique(tra_row$CDR3_AA)
  vgene           <- unique(tstrsplit(as.character(tra_row$vgene_imm),"\\*")[[1]])
  jgene           <- unique(tstrsplit(as.character(tra_row$jgene_imm),"\\*")[[1]])
  input_row$TRAV  <- vgene
  input_row$TRAJ  <- jgene
  
  # Patient ID ('PatientCounts' column)
  #sample  <- as.character(tstrsplit(cell,"-")[[2]])
  patient   <- unique(clone_df$Patient_ID)
  diagnosis <- unique(clone_df$Paper_ID_Diagnosis_loc)
  diagnosis <- unlist(strsplit(diagnosis,"_"))[2]
  
  input_row$Patient <- paste(c(patient,diagnosis),collapse = ":")
  
  # Cell count
  input_row$Counts <- cell_count
  
  # Add row back to original input_df
  input_df[clone,] <- input_row
}

# Format GLIPH2 input table
# CDR3b  TRBV  TRBJ  CDR3a  subject:condition count
input_df_v2 <- input_df
input_df_v2 <- input_df_v2[,c("CDR3b","TRBV","TRBJ","CDR3a","Patient","Counts")]
names(input_df_v2) <- c("CDR3b","TRBV","TRBJ","CDR3a","subject:condition","count")

### 4. Write TSV for GLIPH v1
input_df$Patient <- tstrsplit(input_df$Patient,":")[[1]]
write.table(input_df,file = output_tsv,row.names = F,quote = F,sep = "\t")


### 5. Write TSV for GLIPH v2
write.table(input_df_v2,file = output_tsv_v2,row.names = F,col.names = F,quote = F,sep = "\t")


### 6. Create CDR3b files for GLIPH2
cdr3b_df <- input_df_v2[,c("CDR3b","TRBV","TRBJ")]
write.table(cdr3b_df,file = output_cdr3b_v2,row.names = F,col.names = F,quote = F,sep = "\t")

### 7. TRBV freq
trb_df        <- as.data.frame(table(cdr3b_df$TRBV))
total_v_usage <- sum(trb_df$Freq)
trb_df$Freq   <- trb_df$Freq /total_v_usage
write.table(trb_df,file = output_trb_v2,row.names = F,col.names = F,quote = F,sep = "\t")

### 8. CDR3b lengths freq
cdr3b_len_df <- as.data.frame(cdr3b_df[,c("CDR3b")])
colnames(cdr3b_len_df) <- c("CDR3b")
cdr3b_len_df$AA_len <- nchar(as.character(cdr3b_len_df$CDR3b))

len_freq_df <- as.data.frame(table(cdr3b_len_df$AA_len))
total       <- sum(len_freq_df$Freq)
len_freq_df$Freq   <- len_freq_df$Freq /total
write.table(len_freq_df,file = output_len_v2,row.names = F,col.names = F,quote = F,sep = "\t")

