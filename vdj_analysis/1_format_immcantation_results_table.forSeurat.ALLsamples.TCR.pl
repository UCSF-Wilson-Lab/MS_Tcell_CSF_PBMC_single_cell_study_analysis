#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use DataBrowser;

# PURPOSE: Format Immcantation Results for Seurat (cellranger version 3.1)
#  - get clonotype IDs in there (assigned based on immcantation)
#  - create and expanded not expanded column ('Not detected' for cells with no Ig_seq results)
#  - Sample conversion table generated from <create_igseq_results_table_5prime.pl>
# Use all immcantation results instead of a subset
my $immcantation_results_table = "./immcantation_results_v3.1/10X_beta_clone-pass-FILTERED_germ-pass.final.tab";
my $immcantation_alpha_chain_table = "./immcantation_results_v3.1/10X_alpha-parse-FILTERED_germ-pass.final.tab";

my $sample_conv_table = "./metadata/metadata.csv";
my $output_csv_all    = "./tables/FORMATTED_immcantation_results_table.ALL.TCR.v3.1.csv";
my $THRESH            = 1;   # Base Threshold for expanded vs. non-expanded clonotype

# Immcatation alpha beta chain flags
my $TRA_CHAIN = "TRA";
my $TRB_CHAIN = "TRB";


### 1. Store hash of ID conversions
my %sample_conversion;

open(my $IN1, "<", $sample_conv_table) or die "Cannot read sample conversion table\n";
my $header_conv = <$IN1>;
while (my $line = <$IN1>){
  chomp($line);
  my @info = split(",",$line);
  my $id_rnaseq = $info[1];
  my $id_igseq  = $info[2];
  if ($id_rnaseq =~ /coming/ or $id_rnaseq =~ /n\/a/){$id_rnaseq = "";}
  if ($id_igseq =~ /coming/ or $id_igseq =~ /n\/a/)  {$id_igseq  = "";}
  
  if (length($id_rnaseq) > 0 and length($id_igseq) > 0){
    $sample_conversion{$id_igseq} = $id_rnaseq;
  }
  
}
close $IN1;


### 2. Format Immcantation results and store in Hash

### HEADER
#SEQUENCE_ID	SEQUENCE_INPUT	FUNCTIONAL	IN_FRAME	STOP	MUTATED_INVARIANT	INDELS	LOCUS(7)
#V_CALL	D_CALL	J_CALL	SEQUENCE_VDJ	SEQUENCE_IMGT	V_SEQ_START	V_SEQ_LENGTH	V_GERM_START_VDJ	V_GERM_LENGTH_VDJ(16)
#V_GERM_START_IMGT	V_GERM_LENGTH_IMGT	NP1_LENGTH	D_SEQ_START	D_SEQ_LENGTH	D_GERM_START	D_GERM_LENGTH	NP2_LENGTH(23)
#J_SEQ_START	J_SEQ_LENGTH	J_GERM_START	J_GERM_LENGTH	JUNCTION	JUNCTION_LENGTH	GERMLINE_IMGT(31)	
#V_SCORE	V_IDENTITY	V_EVALUE	V_CIGAR	D_SCORED_IDENTITY	D_EVALUE	D_CIGAR	J_SCORE	J_IDENTITY	J_EVALUE(41)	
#J_CIGAR	FWR1_IMGT	FWR2_IMGT	FWR3_IMGT	FWR4_IMGT	CDR1_IMGT	CDR2_IMGT	CDR3_IMGT	CELL	C_CALL(52)	
#CONSCOUNT	UMICOUNT V_CALL_10X	D_CALL_10X	J_CALL_10X	JUNCTION_10X	JUNCTION_10X_AA CLONE	GERMLINE_IMGT_D_MASK(61)	
#GERMLINE_V_CALL	GERMLINE_D_CALL	GERMLINE_J_CALL(64)


# All samples
my %data_hash_all;
my %expanded_all; # cells per clonotype for ALL TCR samples 
my %alpha_beta_clones;  # This hash is purely meant to merge the clonotype IDs of TRB and TRA
my %cell_clonotype;  # Contains heavy chain clonotype assignments per cell

### 1. Process Beta Chain results
open(my $IN2, "<", $immcantation_results_table) or die "ERROR reading input table\n";
my $header = <$IN2>;
while (my $line = <$IN2>){
  chomp($line);
  my @info      = split("\t", $line);
  my $clonotype = $info[60];
  my $cdr3      = $info[59];
  
  my $umi = $info[54];
  
  # Sequences for mTCR creation
  my $cdr1nuc = $info[48];
  my $cdr2nuc = $info[49];
  my $cdr3nuc = $info[50];
  my $frw1 = $info[44];
  my $frw2 = $info[45];
  my $frw3 = $info[46];
  my $frw4 = $info[47];
  # E-values and Identity
  my $veval  = $info[34];
  my $vident = $info[33];
  my $jeval  = $info[42];
  my $jident = $info[41];
  # input contig
  my $input_seq = $info[1];
  my $mab_stuff = "$umi,$veval,$vident,$jeval,$jident,$input_seq,$frw1,$frw2,$frw3,$frw4,$cdr1nuc,$cdr2nuc,$cdr3nuc";

  # Sample
  my $id_str = $info[0];
  my ($sample,$contig) = split(":", $id_str);
  
  my $sample_fmt    = format_sample_name($sample);
  my $sample_rnaseq = "";
  if (defined $sample_conversion{$sample}){$sample_rnaseq = $sample_conversion{$sample};}
  # Cell
  my $cell      = $info[51];
  if ($cell =~ /-/){$cell = (split /-/,$cell)[0];}
  # Seurat converted ID
  my $seurat_id   = "$cell\-$sample_rnaseq";
  my $igseq_id    = "$cell\-$sample";
  my $cell_id_fmt = "$cell\-$sample_fmt";
  
  # Store heavy chain clonotypes
  $cell_clonotype{$igseq_id} = $clonotype;
  
  # VJ CDR3
  my $vcall = $info[8];
  if ($vcall =~ /,/){$vcall = (split /,/, $vcall)[0];}
  my $jcall = $info[10];
  if ($jcall =~ /,/){$jcall = (split /,/, $jcall)[0];}
  my $ccall = $info[52];
  
  my $chain = $TRB_CHAIN;
  my $entry = "$sample,$sample_rnaseq,,$vcall,$jcall,$ccall,$cdr3,$chain,$mab_stuff,$clonotype";
  if ($sample_rnaseq ne ""){
    $entry = "$sample,$sample_rnaseq,$seurat_id,$vcall,$jcall,$ccall,$cdr3,$chain,$mab_stuff,$clonotype";
  }
  
  # All samples
  if ($chain eq $TRB_CHAIN) {
    $alpha_beta_clones{$cell_id_fmt}{$chain}{$clonotype}++;
  }
  
  $data_hash_all{$cell_id_fmt}{$chain}{$entry}++;
  
  
}
close $IN2;


### 2. Process Alpha Chain results

# Light chain HEADER
#SEQUENCE_ID	SEQUENCE_INPUT	FUNCTIONAL	IN_FRAME	STOP	MUTATED_INVARIANT	INDELS	LOCUS(7)	
#V_CALL	D_CALL	J_CALL	SEQUENCE_VDJ	SEQUENCE_IMGT	V_SEQ_START	V_SEQ_LENGTH	V_GERM_START_VDJ	V_GERM_LENGTH_VDJ(16)	
#V_GERM_START_IMGT	V_GERM_LENGTH_IMGT	NP1_LENGTH	D_SEQ_START	D_SEQ_LENGTH	D_GERM_START	D_GERM_LENGTH	NP2_LENGTH(24)	
#J_SEQ_START	J_SEQ_LENGTH	J_GERM_START	J_GERM_LENGTH	JUNCTION	JUNCTION_LENGTH	GERMLINE_IMGT	V_SCORE	V_IDENTITY	V_EVALUE(34)	
#V_CIGAR	D_SCORE	D_IDENTITY	D_EVALUE	D_CIGAR	J_SCORE	J_IDENTITY	J_EVALUE	J_CIGAR	FWR1_IMGT	FWR2_IMGT	FWR3_IMGT	FWR4_IMGT(47)	
#CDR1_IMGT	CDR2_IMGT	CDR3_IMGT	CELL(51)	C_CALL	CONSCOUNT	UMICOUNT	V_CALL_10X	D_CALL_10X	J_CALL_10X	JUNCTION_10X	JUNCTION_10X_AA(59) CLONE

### Alpha HEADER
#SEQUENCE_ID	SEQUENCE_INPUT	FUNCTIONAL	IN_FRAME	STOP	MUTATED_INVARIANT	INDELS	LOCUS(7)
#V_CALL	D_CALL	J_CALL	SEQUENCE_VDJ	SEQUENCE_IMGT	V_SEQ_START	V_SEQ_LENGTH	V_GERM_START_VDJ	V_GERM_LENGTH_VDJ(16)
#V_GERM_START_IMGT	V_GERM_LENGTH_IMGT	NP1_LENGTH	D_SEQ_START	D_SEQ_LENGTH	D_GERM_START	D_GERM_LENGTH	NP2_LENGTH(23)
#J_SEQ_START	J_SEQ_LENGTH	J_GERM_START	J_GERM_LENGTH	JUNCTION	JUNCTION_LENGTH	GERMLINE_IMGT(31)	
#V_SCORE	V_IDENTITY	V_EVALUE	V_CIGAR	D_SCORED_IDENTITY	D_EVALUE	D_CIGAR	J_SCORE	J_IDENTITY	J_EVALUE(41)	
#J_CIGAR	FWR1_IMGT	FWR2_IMGT	FWR3_IMGT	FWR4_IMGT	CDR1_IMGT	CDR2_IMGT	CDR3_IMGT	CELL	C_CALL(52)	
#CONSCOUNT	UMICOUNT V_CALL_10X	D_CALL_10X	J_CALL_10X	JUNCTION_10X	JUNCTION_10X_AA CLONE	GERMLINE_IMGT_D_MASK(61)	
#GERMLINE_V_CALL	GERMLINE_D_CALL	GERMLINE_J_CALL(64)


open(my $IN3, "<",$immcantation_alpha_chain_table) or die "Cannot read lignt chain table\n";
my $header_light = <$IN3>;
while (my $line = <$IN3>){
  chomp($line);
  my @info      = split("\t", $line);
  my $cdr3      = $info[59];
  my $umi       = $info[54];
  
  # Sequences for mTCR creation
  my $cdr1nuc = $info[48];
  my $cdr2nuc = $info[49];
  my $cdr3nuc = $info[50];
  my $frw1 = $info[44];
  my $frw2 = $info[45];
  my $frw3 = $info[46];
  my $frw4 = $info[47];
  # E-values and Identity
  my $veval  = $info[34];
  my $vident = $info[33];
  my $jeval  = $info[42];
  my $jident = $info[41];
  # input contig
  my $input_seq = $info[1];
  my $mab_stuff = "$umi,$veval,$vident,$jeval,$jident,$input_seq,$frw1,$frw2,$frw3,$frw4,$cdr1nuc,$cdr2nuc,$cdr3nuc";
  

  # Sample
  my $id_str = $info[0];
  my ($sample,$contig) = split(":", $id_str);
  
  my $sample_fmt    = format_sample_name($sample);
  my $sample_rnaseq = "";
  if (defined $sample_conversion{$sample}){$sample_rnaseq = $sample_conversion{$sample};}
  # Cell
  my $cell      = $info[51];
  if ($cell =~ /-/){$cell = (split /-/,$cell)[0];}
  # Seurat converted ID
  my $seurat_id   = "$cell\-$sample_rnaseq";
  my $igseq_id    = "$cell\-$sample";
  my $cell_id_fmt = "$cell\-$sample_fmt";
  
  # Get heavy chain clonotype
  my $clonotype = "";
  $clonotype    = $cell_clonotype{$igseq_id} if defined $cell_clonotype{$igseq_id};
  
  # VJ CDR3
  my $vcall = $info[8];
  if ($vcall =~ /,/){$vcall = (split /,/, $vcall)[0];}
  my $jcall = $info[10];
  if ($jcall =~ /,/){$jcall = (split /,/, $jcall)[0];}
  my $ccall = $info[52];
  
  my $chain = $TRA_CHAIN;
  my $entry = "$sample,$sample_rnaseq,,$vcall,$jcall,$ccall,$cdr3,$chain,$mab_stuff,$clonotype";
  if ($sample_rnaseq ne ""){
    $entry = "$sample,$sample_rnaseq,$seurat_id,$vcall,$jcall,$ccall,$cdr3,$chain,$mab_stuff,$clonotype";
  }
  
  if ($chain eq $TRA_CHAIN) {
    $alpha_beta_clones{$cell_id_fmt}{$chain}{$clonotype}++;
  }
  
  $data_hash_all{$cell_id_fmt}{$chain}{$entry}++;
  
  
}
close $IN3;


### 3. Merge Alpha and Beta Clonotypes
#    - exclude cells with 1 chain, only keep cells with paired TRA and TRB
my %merged_clonotypes; # key = cell ID, value = merged TRA/TRB clonotype
my %chain_summary;     # purely for stats to see how many cells are single chain vs paired alpha/beta

foreach my $cell (keys %alpha_beta_clones){
  my $trb_status = 0;
  my $tra_status = 0;
  $trb_status = 1 if defined $alpha_beta_clones{$cell}{$TRB_CHAIN};
  $tra_status = 1 if defined $alpha_beta_clones{$cell}{$TRA_CHAIN};
  
  if ($trb_status == 1 and $tra_status == 1){
    ## PAIRED TRA and TRB per cell
    $chain_summary{"paired"}++;
    
    my $tra_clone = "";
    my $trb_clone = "";
    
    foreach my $clone (keys %{$alpha_beta_clones{$cell}{$TRB_CHAIN}}){
      $trb_clone .= $clone;
    }
    foreach my $clone (keys %{$alpha_beta_clones{$cell}{$TRA_CHAIN}}){
      $tra_clone .= $clone;
    }
    
    # Merge TRA and TRB Clonotypes
    my $combined_clonotype = $tra_clone . $trb_clone;
    if ($tra_clone ne "" and $trb_clone ne ""){
      if ($tra_clone eq $trb_clone){
        $combined_clonotype = $trb_clone;
      } else{
        $combined_clonotype = $tra_clone . ":" . $trb_clone
      }
    }
    
    # >> Store paired cells
    $expanded_all{$combined_clonotype}{$cell}++;
    $merged_clonotypes{$cell} = $combined_clonotype;

  } else{
    ## SINGLE CHAIN
    ## These entries are not kept for further analysis
    foreach my $chain (keys %{$alpha_beta_clones{$cell}}){
      $chain_summary{"single"}{$chain}++;
    }
  }
}


### 3. All Samples: Establish Clonotype status (Expanded vs. Non-Expanded)
my %expanded_status;
my %expanded_status2;
my %expanded_status3;

foreach my $clone (keys %expanded_all){
  my $count     = 0;
  my $count_csf = 0;
  my $count_pb  = 0;
  
  foreach my $cell (keys %{$expanded_all{$clone}}){
    $count++;
    if ($cell =~ /CSF_/){
      $count_csf++;
    } else{
      $count_pb++;
    }
  }
  
  # Use different thresholds to defined clonally expanded
  # thresh = 1
  if ($count > 1){
    $expanded_status{$clone}{"status"}    = 1;
    $expanded_status{$clone}{"count"}     = $count;
    $expanded_status{$clone}{"count_csf"} = $count_csf;
    $expanded_status{$clone}{"count_pb"}  = $count_pb;
  } else{
    $expanded_status{$clone}{"status"}    = 0;
    $expanded_status{$clone}{"count"}     = $count;
    $expanded_status{$clone}{"count_csf"} = $count_csf;
    $expanded_status{$clone}{"count_pb"}  = $count_pb;
  }
  
  # thresh = 2
  if ($count > 2){
    $expanded_status2{$clone}{"status"}    = 1;
    $expanded_status2{$clone}{"count"}     = $count;
    $expanded_status2{$clone}{"count_csf"} = $count_csf;
    $expanded_status2{$clone}{"count_pb"}  = $count_pb;
  } else{
    $expanded_status2{$clone}{"status"}    = 0;
    $expanded_status2{$clone}{"count"}     = $count;
    $expanded_status2{$clone}{"count_csf"} = $count_csf;
    $expanded_status2{$clone}{"count_pb"}  = $count_pb;
  }
  
  # thresh = 3
  if ($count > 3){
    $expanded_status3{$clone}{"status"}    = 1;
    $expanded_status3{$clone}{"count"}     = $count;
    $expanded_status3{$clone}{"count_csf"} = $count_csf;
    $expanded_status3{$clone}{"count_pb"}  = $count_pb;
  } else{
    $expanded_status3{$clone}{"status"}    = 0;
    $expanded_status3{$clone}{"count"}     = $count;
    $expanded_status3{$clone}{"count_csf"} = $count_csf;
    $expanded_status3{$clone}{"count_pb"}  = $count_pb;
  }
}

#browse(\%expanded_status);


### 4. All Samples: Output final table

# Header:
# Patient_ID,seurat_converted_cell_barcode,sample_igseq,sample_gex,vgene,jgene,cgene,clonotypeID,expanded_clonotype_status
open(my $OUT, ">",$output_csv_all) or die "Cannot write output CSV\n";

my $imm_header = "Patient_ID,cell_barcode_fmt,chain_count,sample_igseq,sample_gex,seurat_converted_cell_barcode,vgene_imm,jgene_imm,cgene_imm,CDR3_AA,chain,";
$imm_header   .= "UMIcount,V_Evalue,V_Ident,J_Evalue,J_Ident,input_seq,FRW1,FRW2,FRW3,FRW4,CDR1,CDR2,CDR3,";
$imm_header   .= "clonotypeID_imm,";
$imm_header   .= "clonotype_cell_count_imm,clonotype_CSF_cell_count_imm,clonotype_PB_cell_count_imm,MERGED_clonotypeID_imm,expanded_clonotype_status_imm\n";
print $OUT $imm_header;

foreach my $cell_id (keys %data_hash_all){
  my ($cell,$sample_fmt) = split("-",$cell_id);
  my $patient = (split /_/,$sample_fmt)[0];
  
  # Get Chain count
  my $chain_count = 0;
  foreach my $chain (keys %{$data_hash_all{$cell_id}}){
    foreach my $cell (keys %{$data_hash_all{$cell_id}{$chain}}){
      $chain_count++;
    }
  }
  
  # Get if merged TRA/TRB clonotype
  my $combined_clonotype = "";
  $combined_clonotype    = $merged_clonotypes{$cell_id} if defined $merged_clonotypes{$cell_id};
  
  # Select most reliable heavy and light chain result
  my $best_beta_entry  = "";
  my $heavy_umi_thresh  = 0;
  my $heavy_veval_thresh = 1;
  my $best_alpha_entry  = "";
  my $light_umi_thresh  = 0;
  my $light_veval_thresh = 1;
  foreach my $chain (keys %{$data_hash_all{$cell_id}}){
    foreach my $entry (keys %{$data_hash_all{$cell_id}{$chain}}){
      my $clonotype = (split /,/, $entry)[-1];
      
      my $exp_cell_count = "";
      my $exp_cell_count_csf = "";
      my $exp_cell_count_pb  = "";
      my $exp_status         = "";
      my $exp_status2        = "";
      my $exp_status3        = "";
      $exp_cell_count = $expanded_status{$combined_clonotype}{"count"} if defined $expanded_status{$combined_clonotype}{"count"};
      $exp_cell_count_csf = $expanded_status{$combined_clonotype}{"count_csf"} if defined $expanded_status{$combined_clonotype}{"count_csf"};
      $exp_cell_count_pb  = $expanded_status{$combined_clonotype}{"count_pb"}  if defined $expanded_status{$combined_clonotype}{"count_pb"};
      $exp_status     = $expanded_status{$combined_clonotype}{"status"}  if defined $expanded_status{$combined_clonotype}{"status"};
      $exp_status2    = $expanded_status2{$combined_clonotype}{"status"} if defined $expanded_status2{$combined_clonotype}{"status"};
      $exp_status3    = $expanded_status3{$combined_clonotype}{"status"} if defined $expanded_status3{$combined_clonotype}{"status"};

      my $complete_entry = "$patient,$cell_id,$chain_count,$entry,$exp_cell_count,$exp_cell_count_csf,$exp_cell_count_pb,$combined_clonotype,$exp_status\n";
      
      if ($chain eq $TRB_CHAIN){
        # Beta
        $best_beta_entry   = $complete_entry;
      } else{
        $best_alpha_entry  = $complete_entry;
      }
      
    }
  }
  print $OUT $best_beta_entry if $best_beta_entry ne "";
  print $OUT $best_alpha_entry if $best_alpha_entry ne "";
  
}
close $OUT;

print "\nCREATED:\t$output_csv_all\n";


print "\n\n>>> DONE! <<<\n\n";



#================#
# SUBROUTINES    #
#================#

# Format sample names to facilitate conversion between 10X TCR results and 10X GEX results (5 prime datasets)
sub format_sample_name{
  my ($sample)   = @_;
  my @name_comp  = split("_", $sample);
  my $patient_tp = $name_comp[0];
  if ($patient_tp =~ /[a|d]CSF/) {$patient_tp = (split /CSF/,$patient_tp)[0];}
  if ($patient_tp =~ /[a|d]PBMC/) {$patient_tp = (split /PBMC/,$patient_tp)[0];}

  my $source     = "PB";   # for our data, source is always either PB or CSF
  my $cell_type  = "uns";  # By Default unsorted (uns)
  
  # CSF or PB
  if ($sample =~ /CSF/){$source = "CSF";}
  
  my $formatted_name = "$patient_tp\_$source\_$cell_type";
  
  return $formatted_name;
}


