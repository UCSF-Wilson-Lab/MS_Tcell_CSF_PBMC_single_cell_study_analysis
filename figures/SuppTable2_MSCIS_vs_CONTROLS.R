library(tidyverse)
library(dplyr)
library(readxl)	

working_dir <- "./sc_analysis/"
setwd(working_dir)								

Supplemental_Table_2_fmtVJ_mab_table_ALL_TCR_v3_1_version_1_xlsb <- read_excel("./tables/SupplementalTable2_fmtVJ_mab_table_ALL_TCR.v3.1.xlsb.xlsx")									

# Keep only CD8 T cells									
CD8_Only_Data <- Supplemental_Table_2_fmtVJ_mab_table_ALL_TCR_v3_1_version_1_xlsb %>%									
  filter(subtype.cell.annot.manual == "CD8 T cells")									

#  Keep only rows where CSF count > 1									
CD8_CSF_Filter <- CD8_Only_Data %>%									
  filter(clonotype_CSF_cell_count_imm > 1) %>% select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc, Paper_ID_Diagnosis_loc)									

#  Keep only rows where PB count > 1									
CD8_PB_Filter <- CD8_Only_Data %>%									
  filter(clonotype_PB_cell_count_imm > 1) %>% select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc, Paper_ID_Diagnosis_loc)									

# Filter for LOC PB and CSF									
CD8_PB_PB_Filter <- CD8_PB_Filter %>%									
  filter(loc == "PB")									

CD8_PB_CSF_Filter <- CD8_PB_Filter %>%									
  filter(loc == "CSF")									

CD8_CSF_PB_Filter <- CD8_CSF_Filter %>%									
  filter(loc == "PB")									

CD8_CSF_CSF_Filter <- CD8_CSF_Filter %>%									
  filter(loc == "CSF")									

#Filtering for MS									
CD8_PB_PB_MS <- CD8_PB_PB_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "RRMS"))									

CD8_PB_CSF_MS <- CD8_PB_CSF_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "RRMS"))									

CD8_CSF_PB_MS <- CD8_CSF_PB_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "RRMS"))									

CD8_CSF_CSF_MS <- CD8_CSF_CSF_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "RRMS"))									

#Filtering for CIS									
CD8_PB_PB_CIS <- CD8_PB_PB_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "CIS"))									

CD8_PB_CSF_CIS <- CD8_PB_CSF_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "CIS"))									

CD8_CSF_PB_CIS <- CD8_CSF_PB_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "CIS"))									

CD8_CSF_CSF_CIS <- CD8_CSF_CSF_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "CIS"))									

#Filtering for HC									
CD8_PB_PB_HC <- CD8_PB_PB_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "HC"))									

CD8_PB_CSF_HC <- CD8_PB_CSF_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "HC"))									

CD8_CSF_PB_HC <- CD8_CSF_PB_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "HC"))									

CD8_CSF_CSF_HC <- CD8_CSF_CSF_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "HC"))									

#Filtering for OND									
CD8_PB_PB_OND <- CD8_PB_PB_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "Neurosarcoidosis|Uveitis"))									

CD8_PB_CSF_OND <- CD8_PB_CSF_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "Neurosarcoidosis|Uveitis"))									

CD8_CSF_PB_OND <- CD8_CSF_PB_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "Neurosarcoidosis|Uveitis"))									

CD8_CSF_CSF_OND <- CD8_CSF_CSF_Filter %>%									
  filter(str_detect(Paper_ID_Diagnosis_loc , "Neurosarcoidosis|Uveitis"))									


#Group by Cell Barcode for the Single Filtered Data									
GID_CD8_CSF_Filter <- CD8_CSF_Filter %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(unique(chain), collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_PB_Filter <- CD8_PB_Filter %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

#Group by Cell Barcode Double Filter									
GID_CD8_PB_PB_Filter <- CD8_PB_PB_Filter %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_PB_CSF_Filter <- CD8_PB_CSF_Filter %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ", "),									
    .groups = "drop"									
  )									

GID_CD8_CSF_PB_Filter <- CD8_CSF_PB_Filter %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(unique(chain), collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_CSF_CSF_Filter <- CD8_CSF_CSF_Filter %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									


#Group by Cell Barcode Double Filter AND MS Subfilter									
GID_CD8_PB_PB_MS <- CD8_PB_PB_MS %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_PB_CSF_MS <- CD8_PB_CSF_MS %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ", "),									
    .groups = "drop"									
  )									

GID_CD8_CSF_PB_MS <- CD8_CSF_PB_MS %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(unique(chain), collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_CSF_CSF_MS <- CD8_CSF_CSF_MS %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

#Group by Cell Barcode Double Filter AND CIS Subfilter									
GID_CD8_PB_PB_CIS <- CD8_PB_PB_CIS %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_PB_CSF_CIS <- CD8_PB_CSF_CIS %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ", "),									
    .groups = "drop"									
  )									

GID_CD8_CSF_PB_CIS <- CD8_CSF_PB_CIS %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(unique(chain), collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_CSF_CSF_CIS <- CD8_CSF_CSF_CIS %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

#Group by Cell Barcode Double Filter AND HC Subfilter									
GID_CD8_PB_PB_HC <- CD8_PB_PB_HC %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_PB_CSF_HC <- CD8_PB_CSF_HC %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ", "),									
    .groups = "drop"									
  )									

GID_CD8_CSF_PB_HC <- CD8_CSF_PB_HC %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(unique(chain), collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_CSF_CSF_HC <- CD8_CSF_CSF_HC %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

#Group by Cell Barcode Double Filter AND OND Subfilter									
GID_CD8_PB_PB_OND <- CD8_PB_PB_OND %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_PB_CSF_OND <- CD8_PB_CSF_OND %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ", "),									
    .groups = "drop"									
  )									

GID_CD8_CSF_PB_OND <- CD8_CSF_PB_OND %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(unique(chain), collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

GID_CD8_CSF_CSF_OND <- CD8_CSF_CSF_OND %>%									
  select(cell_barcode_fmt, vgene_imm, jgene_imm, CDR3_AA, chain, loc) %>%									
  group_by(cell_barcode_fmt) %>%									
  arrange(factor(chain, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    vgene_imm = paste(unique(vgene_imm), collapse = ","),									
    jgene_imm = paste(unique(jgene_imm), collapse = ","),									
    CDR3_AA   = paste(unique(CDR3_AA), collapse = ","),									
    chain     = paste(sort(unique(factor(chain, levels = c("TRA","TRB")))),									
                      collapse = ","),									
    loc       = paste(unique(loc), collapse = ","),									
    .groups = "drop"									
  )									

#Load VDJDG Files									
read.delim("vdjdb_EBV_raw_data.tsv", stringsAsFactors = FALSE)									
read.delim("vdjdb_cmv_raw_data.tsv", stringsAsFactors = FALSE)									

vdjdb_EBV_Original <- read.delim("vdjdb_EBV_raw_data.tsv", stringsAsFactors = FALSE) 									
vdjdb_CMV_Original <- read.delim("vdjdb_cmv_raw_data.tsv", stringsAsFactors = FALSE) 									

GroupedID_VDJDB_EBV <- vdjdb_EBV_Original %>%									
  select(complex.id, Gene, CDR3, V, J, Epitope.species) %>%									
  group_by(complex.id) %>%									
  arrange(factor(Gene, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    Gene            = paste(unique(Gene), collapse = ","),									
    CDR3            = paste(unique(CDR3), collapse = ","),									
    V               = paste(unique(V), collapse = ","),									
    J               = paste(unique(J), collapse = ","),									
    Epitope.species = paste(unique(Epitope.species), collapse = ","),									
    .groups = "drop")									

GroupedID_VDJDB_CMV <- vdjdb_CMV_Original %>%									
  select(complex.id, Gene, CDR3, V, J, Epitope.species) %>%									
  group_by(complex.id) %>%									
  arrange(factor(Gene, levels = c("TRA", "TRB"))) %>%  # Make sure TRA rows are first before summarise									
  summarise(									
    Gene            = paste(unique(Gene), collapse = ","),									
    CDR3            = paste(unique(CDR3), collapse = ","),									
    V               = paste(unique(V), collapse = ","),									
    J               = paste(unique(J), collapse = ","),									
    Epitope.species = paste(unique(Epitope.species), collapse = ","),									
    .groups = "drop"									
  )									

Columns_VDJDB_EBV <- vdjdb_EBV_Original %>%									
  select(complex.id, Gene, CDR3, V, J, Epitope.species)									

Columns_VDJDB_CMV <- vdjdb_CMV_Original %>%									
  select(complex.id, Gene, CDR3, V, J, Epitope.species)									



GroupedID_VDJDB_EBV_AA <- GroupedID_VDJDB_EBV %>%									
  rename(CDR3_AA = CDR3)									

GroupedID_VDJDB_CMV_AA <- GroupedID_VDJDB_CMV %>%									
  rename(CDR3_AA = CDR3)									


SID_VDJDB_EBV_AA <- Columns_VDJDB_EBV %>%									
  rename(CDR3_AA = CDR3)									

SID_VDJDB_CMV_AA <- Columns_VDJDB_CMV %>%									
  rename(CDR3_AA = CDR3)									

#Unique Values only for CDR3_AA 									
unique_GID_VDJDB_CMV_CDR3_AA <- GroupedID_VDJDB_CMV_AA %>%									
  distinct(CDR3_AA, .keep_all = TRUE)									
unique_GID_VDJDB_EBV_CDR3_AA <- GroupedID_VDJDB_EBV_AA %>%									
  distinct(CDR3_AA, .keep_all = TRUE)									

(intersect(GID_CD8_CSF_CSF_Filter$CDR3_AA, unique_GID_VDJDB_EBV_CDR3_AA$CDR3_AA))									
(intersect(GID_CD8_PB_PB_Filter$CDR3_AA, unique_GID_VDJDB_EBV_CDR3_AA$CDR3_AA))									
(intersect(GID_CD8_CSF_CSF_Filter$CDR3_AA, unique_GID_VDJDB_CMV_CDR3_AA$CDR3_AA))									
(intersect(GID_CD8_PB_PB_Filter$CDR3_AA, unique_GID_VDJDB_CMV_CDR3_AA$CDR3_AA))									

Unique_merged_CSF_CSF_EBV <- inner_join(GID_CD8_CSF_CSF_Filter, unique_GID_VDJDB_EBV_CDR3_AA, by = "CDR3_AA")									
Unique_merged_PB_PB_EBV <- inner_join(GID_CD8_PB_PB_Filter, unique_GID_VDJDB_EBV_CDR3_AA, by = "CDR3_AA")									
Unique_merged_CSF_CSF_CMV <- inner_join(GID_CD8_CSF_CSF_Filter, unique_GID_VDJDB_CMV_CDR3_AA, by = "CDR3_AA")									
Unique_merged_PB_PB_CMV <- inner_join(GID_CD8_PB_PB_Filter, unique_GID_VDJDB_CMV_CDR3_AA, by = "CDR3_AA")									

length(intersect(GID_CD8_CSF_CSF_Filter$CDR3_AA, unique_GID_VDJDB_EBV_CDR3_AA$CDR3_AA))									
length(intersect(GID_CD8_PB_PB_Filter$CDR3_AA, unique_GID_VDJDB_EBV_CDR3_AA$CDR3_AA))									
length(intersect(GID_CD8_CSF_CSF_Filter$CDR3_AA, unique_GID_VDJDB_CMV_CDR3_AA$CDR3_AA))									
length(intersect(GID_CD8_PB_PB_Filter$CDR3_AA, unique_GID_VDJDB_CMV_CDR3_AA$CDR3_AA))									

#Merge for Single Filter									
merged_CSF_EBV <- inner_join(GID_CD8_CSF_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
merged_PB_EBV <- inner_join(GID_CD8_PB_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
merged_CSF_CMV <- inner_join(GID_CD8_CSF_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									
merged_PB_CMV <- inner_join(GID_CD8_PB_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									

merged_CSF_CSF_EBV <- inner_join(GID_CD8_CSF_CSF_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									

#Merge for Double Filter									
merged_CSF_CSF_EBV <- inner_join(GID_CD8_CSF_CSF_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
merged_PB_PB_EBV <- inner_join(GID_CD8_PB_PB_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
merged_CSF_CSF_CMV <- inner_join(GID_CD8_CSF_CSF_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									
merged_PB_PB_CMV <- inner_join(GID_CD8_PB_PB_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									

merged_CSF_PB_EBV <- inner_join(GID_CD8_CSF_PB_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
merged_PB_CSF_EBV <- inner_join(GID_CD8_PB_CSF_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
merged_CSF_PB_CMV <- inner_join(GID_CD8_CSF_PB_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									
merged_PB_CSF_CMV <- inner_join(GID_CD8_PB_CSF_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									



#For the merge total data, hits and non hits: Single Filter									
full_merged_CSF_EBV <- full_join(GID_CD8_CSF_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
full_merged_PB_EBV <- full_join(GID_CD8_PB_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
full_merged_CSF_CMV <- full_join(GID_CD8_CSF_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									
full_merged_PB_CMV <- full_join(GID_CD8_PB_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									

#For the merge total data, hits and non hits: Double Filter									
full_merged_CSF_CSF_EBV <- full_join(GID_CD8_CSF_CSF_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
full_merged_PB_PB_EBV <- full_join(GID_CD8_PB_PB_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
full_merged_CSF_CSF_CMV <- full_join(GID_CD8_CSF_CSF_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									
full_merged_PB_PB_CMV <- full_join(GID_CD8_PB_PB_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									

full_merged_CSF_PB_EBV <- full_join(GID_CD8_CSF_PB_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
full_merged_PB_CSF_EBV <- full_join(GID_CD8_PB_CSF_Filter, GroupedID_VDJDB_EBV_AA, by = "CDR3_AA")									
full_merged_CSF_PB_CMV <- full_join(GID_CD8_CSF_PB_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									
full_merged_PB_CSF_CMV <- full_join(GID_CD8_PB_CSF_Filter, GroupedID_VDJDB_CMV_AA, by = "CDR3_AA")									

#Verification to see if there are any hits									
any(GID_CD8_CSF_Filter$CDR3_AA %in% GroupedID_VDJDB_EBV_AA$CDR3_AA)									
any(GID_CD8_PB_Filter$CDR3_AA %in% GroupedID_VDJDB_EBV_AA$CDR3_AA)									
any(GID_CD8_CSF_Filter$CDR3_AA %in% GroupedID_VDJDB_CMV_AA$CDR3_AA)									
any(GID_CD8_PB_Filter$CDR3_AA %in% GroupedID_VDJDB_CMV_AA$CDR3_AA)									

any(GID_CD8_CSF_PB_Filter$CDR3_AA %in% GroupedID_VDJDB_EBV_AA$CDR3_AA)									
any(GID_CD8_CSF_CSF_Filter$CDR3_AA %in% GroupedID_VDJDB_EBV_AA$CDR3_AA)									
any(GID_CD8_PB_PB_Filter$CDR3_AA %in% GroupedID_VDJDB_EBV_AA$CDR3_AA)									
any(GID_CD8_PB_CSF_Filter$CDR3_AA %in% GroupedID_VDJDB_EBV_AA$CDR3_AA)									

any(GID_CD8_CSF_PB_Filter$CDR3_AA %in% GroupedID_VDJDB_CMV_AA$CDR3_AA)									
any(GID_CD8_CSF_CSF_Filter$CDR3_AA %in% GroupedID_VDJDB_CMV_AA$CDR3_AA)									
any(GID_CD8_PB_PB_Filter$CDR3_AA %in% GroupedID_VDJDB_CMV_AA$CDR3_AA)									
any(GID_CD8_PB_CSF_Filter$CDR3_AA %in% GroupedID_VDJDB_CMV_AA$CDR3_AA)									

#merging to see if there are any hits with TRA and TRB Separate									
single_merged_CSF_EBV <- inner_join(CD8_CSF_Filter, SID_VDJDB_EBV_AA, by = "CDR3_AA")									
single_merged_PB_EBV <- inner_join(CD8_PB_Filter, SID_VDJDB_EBV_AA, by = "CDR3_AA")									
