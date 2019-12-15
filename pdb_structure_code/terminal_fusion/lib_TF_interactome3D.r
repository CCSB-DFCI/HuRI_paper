source("~/R_libs/lib_Interactome3D.r")

faux_switch_partners <- function(data_line)
{
  fields_to_switch1 <- colnames(data_line)[ grep("1$", colnames(data_line)) ]
  fields_to_switch2 <- colnames(data_line)[ grep("2$", colnames(data_line)) ]

  for (i in 1:length(fields_to_switch1))
  {
    field1 <- fields_to_switch1[[i]]
    field2 <- fields_to_switch2[[i]]
    tmp <- data_line[[field1]]
    data_line[[field1]] <- data_line[[field2]]
    data_line[[field2]] <- tmp
  } 
  return (data_line)
}

faux_add_sense_of_pairs_tested <- function(merged_data, tested_info)
{
  new_data <- list()
  for (i in 1:nrow(merged_data))
  {
    this_line <- merged_data[i,]
    P1 <- this_line[["PROT1"]]
    P2 <- this_line[["PROT2"]]
    if (P1 == P2)
    { # return only one pair for homodimers
      # get "co" measure that corresponds to the most complete tail
      if (this_line$co_Cterm_chop1 < this_line$co_Cterm_chop2)
      {
        this_line$co_Cterm_chop2 <- this_line$co_Cterm_chop1
        this_line$co_Cterm_res2 <- this_line$co_Cterm_res1
        this_line$co_d_C2 <- this_line$co_d_C1
      }
      if (this_line$co_Cterm_chop2 < this_line$co_Cterm_chop1)
      {
        this_line$co_Cterm_chop1 <- this_line$co_Cterm_chop2
        this_line$co_Cterm_res1 <- this_line$co_Cterm_res2
      }
      if (this_line$co_Nterm_chop1 < this_line$co_Nterm_chop2)
      {
        this_line$co_Nterm_chop2 <- this_line$co_Nterm_chop1
        this_line$co_Nterm_res2 <- this_line$co_Nterm_res1
      }
      if (this_line$co_Nterm_chop2 < this_line$co_Nterm_chop1)
      {
        this_line$co_Nterm_chop1 <- this_line$co_Nterm_chop2
        this_line$co_Nterm_res1 <- this_line$co_Nterm_res2
      }
      new_data[[i]] <- this_line
    }else
    {
      if ((P1 %in% tested_info$uniprot$tested_non_aa) & (P2 %in% tested_info$uniprot$tested_non_aa))
      {
        new_data[[i]] <- rbind(this_line, faux_switch_partners(this_line))
      }else if (P1 %in% tested_info$uniprot$tested_non_aa) 
      { # pair as it is
        new_data[[i]] <- this_line
      }else if (P2 %in% tested_info$uniprot$tested_non_aa) 
      { # switch info
        new_data[[i]] <- faux_switch_partners(this_line)
      }
    }
  }
  df <- as.data.frame(do.call("rbind", new_data), stringsAsFactors = FALSE)
  return (df)
}

# for each PPI w/ structure, get all pairs that were tested in HuRI (usually each pair A-B has been screened as DB-AD and AD-DB, unless there are autoactivators)
# there are no cases where both proteins are autoactivators
faux_get_tested_pairs_w_structure <- function(intx_bw_tested, tested_info)
{
  tested_pairs <- as.data.frame( do.call("rbind", apply(intx_bw_tested, 1, function(x)
  {
    P1 <- x[["PROT1"]]
    P2 <- x[["PROT2"]]
    if ((P1 %in% tested_info$uniprot$tested_non_aa) & (P2 %in% tested_info$uniprot$tested_non_aa))
    {
      rbind(c(P1, P2), c(P2, P1))
    }else if (P1 %in% tested_info$uniprot$tested_non_aa) 
    {
      c(P1, P2)
    }else if (P2 %in% tested_info$uniprot$tested_non_aa) 
    {
      c(P2, P1)
    }
  })), stringsAsFactors = FALSE)
  colnames(tested_pairs) <- c("DB", "AD")
  tested_pairs$int_ID <- apply(tested_pairs, 1, function(x) { tmp <- sort(c(x[[1]], x[[2]])); sprintf("%s_%s", tmp[[1]], tmp[[2]]) })
  return (tested_pairs)
}


# get PPI models and structures
# interaction_data <- interactome3D_data$interactions$raw_data
faux_get_interactome3D_dimers_bw_huri_tested <- function(tested_info, interaction_data, seq_lengths, B_INCLUDE_HOMODIMERS, B_INCLUDE_PPI_MODELS=FALSE, B_INCLUDE_DDM=FALSE)
{
  PPI_TYPE <- c("Structure")
  if (B_INCLUDE_PPI_MODELS)
  {
    PPI_TYPE <- c(PPI_TYPE, "Model")
  }
  if (B_INCLUDE_DDM)
  {
    PPI_TYPE <- c(PPI_TYPE, "Dom_dom_model")
  }
  intx_bw_tested <- interaction_data[ interaction_data$PROT1 %in% tested_info$uniprot$tested &
                                      interaction_data$PROT2 %in% tested_info$uniprot$tested &
                                      interaction_data$TYPE %in% PPI_TYPE,]
  if (!B_INCLUDE_HOMODIMERS)
  {
    intx_bw_tested <- intx_bw_tested[ intx_bw_tested$PROT1 != intx_bw_tested$PROT2,]
  }
  
  # add protein lengths
  intx_bw_tested$seq_length1 <- seq_lengths[ intx_bw_tested$PROT1 ]
  intx_bw_tested$seq_length2 <- seq_lengths[ intx_bw_tested$PROT2 ]
  
  # calculate number of residues chopped for each tail
  # "co_" identifies complex data; "mo_" will identify monomer data
  intx_bw_tested$co_Nterm_chop1 <- (intx_bw_tested$SEQ_BEGIN1 - 1)
  intx_bw_tested$co_Cterm_chop1 <- (intx_bw_tested$seq_length1 - intx_bw_tested$SEQ_END1)
  intx_bw_tested$co_Nterm_chop2 <- (intx_bw_tested$SEQ_BEGIN2 - 1)
  intx_bw_tested$co_Cterm_chop2 <- (intx_bw_tested$seq_length2 - intx_bw_tested$SEQ_END2)
  
  # remove some useless information (for us)
  fields_to_remove <- c("RANK_MAJOR", "RANK_MINOR", "BIO_UNIT", "MODEL1", "DOMAIN1", "MODEL2", "DOMAIN2", "protein1", "protein2")
  for (this_field in fields_to_remove)
  {
    intx_bw_tested[[ this_field ]] <- NULL
  }
  return (intx_bw_tested)
}

## Adds the residue number of terminal nodes in the structures
# for very few proteins we do not have mappings
faux_add_terminal_residues_struct <- function(intx_bw_tested, interactome3D_version)
{
  residue_mapping_all <- f_read_file(sprintf("~/projects_data/interactome3d/dsysmap_%s.structure_residue.csv.nosync", interactome3D_version), HEADER=TRUE, SKIP=0, SEPARATOR = "|")

  # filter data to make it more efficient
  proteins <- unique(c(intx_bw_tested$PROT1, intx_bw_tested$PROT2))
  pdb_files <- intx_bw_tested$FILENAME
  residue_mapping <- residue_mapping_all[ residue_mapping_all$uniprot_ac %in% proteins & residue_mapping_all$pdb_filename %in% pdb_files, ]
  # generate IDs for residue_mapping
  residue_mapping$res_ID <- sprintf("%s_%s_%s_%s", residue_mapping$uniprot_ac, residue_mapping$pdb_filename, residue_mapping$chain, residue_mapping$res_num)
  # sanity check
  if (sum(duplicated(residue_mapping$res_ID)) > 0) { stop("\nnot unique IDs\n")}
  seq_to_struct <- residue_mapping$struct_res_num
  names(seq_to_struct) <- residue_mapping$res_ID
  
  # get Nterminal residue in structure of prot1 -> always chain A because that's how Roberto generates the complex PDB
  tmp_IDs <- sprintf("%s_%s_%s_%s", intx_bw_tested$PROT1, intx_bw_tested$FILENAME, "A", intx_bw_tested$SEQ_BEGIN1)
  intx_bw_tested$co_Nterm_res1 <- seq_to_struct[ tmp_IDs ]
  # get Cterminal residue in structure of prot1
  tmp_IDs <- sprintf("%s_%s_%s_%s", intx_bw_tested$PROT1, intx_bw_tested$FILENAME, "A", intx_bw_tested$SEQ_END1)
  intx_bw_tested$co_Cterm_res1 <- seq_to_struct[ tmp_IDs ]

  # same for prot2: Nterm
  tmp_IDs <- sprintf("%s_%s_%s_%s", intx_bw_tested$PROT2, intx_bw_tested$FILENAME, "B", intx_bw_tested$SEQ_BEGIN2)
  intx_bw_tested$co_Nterm_res2 <- seq_to_struct[ tmp_IDs ]
  # same for prot2: Cterm
  tmp_IDs <- sprintf("%s_%s_%s_%s", intx_bw_tested$PROT2, intx_bw_tested$FILENAME, "B", intx_bw_tested$SEQ_END2)
  intx_bw_tested$co_Cterm_res2 <- seq_to_struct[ tmp_IDs ]

  ## same for monomers
  pdb_files <- unique(c(intx_bw_tested$mo_Nterm_PDB1, intx_bw_tested$mo_Cterm_PDB1, intx_bw_tested$mo_Nterm_PDB2, intx_bw_tested$mo_Cterm_PDB2))
  residue_mapping <- residue_mapping_all[ residue_mapping_all$uniprot_ac %in% proteins & residue_mapping_all$pdb_filename %in% pdb_files, ]
  # generate IDs for residue_mapping
  residue_mapping$res_ID <- sprintf("%s_%s_%s_%s", residue_mapping$uniprot_ac, residue_mapping$pdb_filename, residue_mapping$chain, residue_mapping$res_num)
  # sanity check
  if (sum(duplicated(residue_mapping$res_ID)) > 0) { stop("\nnot unique IDs\n")}
  seq_to_struct <- residue_mapping$struct_res_num
  names(seq_to_struct) <- residue_mapping$res_ID
  
  # same for monomer 1: Nterm
  tmp_IDs <- sprintf("%s_%s_%s_%s", intx_bw_tested$PROT1, intx_bw_tested$mo_Nterm_PDB1, intx_bw_tested$mo_Nterm_chain1, intx_bw_tested$mo_Nterm_SEQ_BEGIN1)
  intx_bw_tested$mo_Nterm_res1 <- seq_to_struct[ tmp_IDs ]
  # same for monomer 1: Cterm
  tmp_IDs <- sprintf("%s_%s_%s_%s", intx_bw_tested$PROT1, intx_bw_tested$mo_Cterm_PDB1, intx_bw_tested$mo_Cterm_chain1, intx_bw_tested$mo_Cterm_SEQ_END1)
  intx_bw_tested$mo_Cterm_res1 <- seq_to_struct[ tmp_IDs ]
  # same for monomer 2: Nterm
  tmp_IDs <- sprintf("%s_%s_%s_%s", intx_bw_tested$PROT2, intx_bw_tested$mo_Nterm_PDB2, intx_bw_tested$mo_Nterm_chain2, intx_bw_tested$mo_Nterm_SEQ_BEGIN2)
  intx_bw_tested$mo_Nterm_res2 <- seq_to_struct[ tmp_IDs ]
  # same for monomer 2: Cterm
  tmp_IDs <- sprintf("%s_%s_%s_%s", intx_bw_tested$PROT2, intx_bw_tested$mo_Cterm_PDB2, intx_bw_tested$mo_Cterm_chain2, intx_bw_tested$mo_Cterm_SEQ_END2)
  intx_bw_tested$mo_Cterm_res2 <- seq_to_struct[ tmp_IDs ]

  return (intx_bw_tested)
}

#  protein_data <- interactome3D_data$proteins$raw_data
faux_add_best_monomers <- function(intx_bw_tested, protein_data, MIN_COVERAGE, B_INCLUDE_MODELS)
{
  TYPE <- c("Structure")
  if (B_INCLUDE_MODELS) { TYPE <- c(TYPE, "Model") }
  # get best monomers for protein1
  fields_to_keep <- c("FILENAME", "CHAIN", "TYPE", "SEQ_BEGIN", "SEQ_END")

  # PROTEIN1
  best_monomers <- sapply(intx_bw_tested$PROT1, simplify = FALSE,
                          function(this_prot)
                          {
                            monomer_data <- protein_data[ protein_data$UNIPROT_AC == this_prot & 
                                                          protein_data$COVERAGE >= MIN_COVERAGE &
                                                          protein_data$TYPE %in% TYPE,]
                            # select the monomer w/ the most complete N-terminal (lowest SEQ_BEGIN1)
                            # if more than one, select w/ highest SEQ_ID. If still more than one, select w/ highest COVERAGE. If still more than one, pick first
                            if (nrow(monomer_data) > 0)
                            {
                              Nterm_best <- monomer_data[ order(monomer_data$SEQ_BEGIN, -monomer_data$SEQ_IDENT, -monomer_data$COVERAGE)[[1]] , fields_to_keep]
                              Cterm_best <- monomer_data[ order(-monomer_data$SEQ_END, -monomer_data$SEQ_IDENT, -monomer_data$COVERAGE)[[1]] , fields_to_keep]
                              best <- rbind(Nterm_best, Cterm_best)
                              row.names(best) <- c("Nterm", "Cterm")
                              best
                            }
                          })
  intx_bw_tested$mo_Nterm_PDB1 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Nterm", "FILENAME"], NA)})
  intx_bw_tested$mo_Nterm_type1 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Nterm", "TYPE"], NA)})
  intx_bw_tested$mo_Nterm_chain1 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Nterm", "CHAIN"], NA)})
  intx_bw_tested$mo_Nterm_chain1 <- ifelse(intx_bw_tested$mo_Nterm_type1 == "Model", "A", intx_bw_tested$mo_Nterm_chain1)
  intx_bw_tested$mo_Nterm_SEQ_BEGIN1 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Nterm", "SEQ_BEGIN"], NA)})
  intx_bw_tested$mo_Nterm_chop1 <- (intx_bw_tested$mo_Nterm_SEQ_BEGIN1 - 1)
  intx_bw_tested$mo_Cterm_PDB1 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Cterm", "FILENAME"], NA)})
  intx_bw_tested$mo_Cterm_type1 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Cterm", "TYPE"], NA)})
  intx_bw_tested$mo_Cterm_chain1 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Cterm", "CHAIN"], NA)})
  intx_bw_tested$mo_Cterm_chain1 <- ifelse(intx_bw_tested$mo_Cterm_type1 == "Model", "A", intx_bw_tested$mo_Cterm_chain1)
  intx_bw_tested$mo_Cterm_SEQ_END1 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Cterm", "SEQ_END"], NA)})
  intx_bw_tested$mo_Cterm_chop1 <- (intx_bw_tested$seq_length1 - intx_bw_tested$mo_Cterm_SEQ_END1)
  
  ## SAME FOR PROTEIN2
  best_monomers <- sapply(intx_bw_tested$PROT2, simplify = FALSE,
                          function(this_prot)
                          {
                            monomer_data <- protein_data[ protein_data$UNIPROT_AC == this_prot & 
                                                          protein_data$COVERAGE >= MIN_COVERAGE &
                                                          protein_data$TYPE %in% TYPE,]
                            # select the monomer w/ the most complete N-terminal (lowest SEQ_BEGIN1)
                            # if more than one, select w/ highest SEQ_ID. If still more than one, select w/ highest COVERAGE. If still more than one, pick first
                            if (nrow(monomer_data) > 0)
                            {
                              Nterm_best <- monomer_data[ order(monomer_data$SEQ_BEGIN, -monomer_data$SEQ_IDENT, -monomer_data$COVERAGE)[[1]] , fields_to_keep]
                              Cterm_best <- monomer_data[ order(-monomer_data$SEQ_END, -monomer_data$SEQ_IDENT, -monomer_data$COVERAGE)[[1]] , fields_to_keep]
                              best <- rbind(Nterm_best, Cterm_best)
                              row.names(best) <- c("Nterm", "Cterm")
                              best
                            }
                          })
  intx_bw_tested$mo_Nterm_PDB2 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Nterm", "FILENAME"], NA)})
  intx_bw_tested$mo_Nterm_type2 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Nterm", "TYPE"], NA)})
  intx_bw_tested$mo_Nterm_chain2 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Nterm", "CHAIN"], NA)})
  intx_bw_tested$mo_Nterm_chain2 <- ifelse(intx_bw_tested$mo_Nterm_type2 == "Model", "A", intx_bw_tested$mo_Nterm_chain2)
  intx_bw_tested$mo_Nterm_SEQ_BEGIN2 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Nterm", "SEQ_BEGIN"], NA)})
  intx_bw_tested$mo_Nterm_chop2 <- (intx_bw_tested$mo_Nterm_SEQ_BEGIN2 - 1)
  intx_bw_tested$mo_Cterm_PDB2 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Cterm", "FILENAME"], NA)})
  intx_bw_tested$mo_Cterm_type2 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Cterm", "TYPE"], NA)})
  intx_bw_tested$mo_Cterm_chain2 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Cterm", "CHAIN"], NA)})
  intx_bw_tested$mo_Cterm_chain2 <- ifelse(intx_bw_tested$mo_Cterm_type2 == "Model", "A", intx_bw_tested$mo_Cterm_chain2)
  intx_bw_tested$mo_Cterm_SEQ_END2 <- sapply(best_monomers, function(x) { ifelse(!is.null(x), x["Cterm", "SEQ_END"], NA)})
  intx_bw_tested$mo_Cterm_chop2 <- (intx_bw_tested$seq_length2 - intx_bw_tested$mo_Cterm_SEQ_END2)
  
  return (intx_bw_tested)
}

faux_map_iface_residues_onto_monomers <- function(merged_data, iface_res_input_file="output_interface_residues.txt", iface_res_monomer_output_file="iface_residues_on_monomers.txt")
{
  residue_mapping_all <- f_read_file("~/projects_data/interactome3d/dsysmap_2018_04.structure_residue.csv.nosync", HEADER=TRUE, SKIP=0, SEPARATOR = "|")

  ## map interface center residue
  proteins <- unique(c(merged_data$PROT1, merged_data$PROT2))
  pdb_files <- unique(c(merged_data$FILENAME, merged_data$mo_Nterm_PDB1, merged_data$mo_Cterm_PDB1, merged_data$mo_Nterm_PDB2, merged_data$mo_Cterm_PDB2))
  residue_mapping <- residue_mapping_all[ residue_mapping_all$uniprot_ac %in% proteins & residue_mapping_all$pdb_filename %in% pdb_files, ]
  # seq to struct
  residue_mapping$res_ID <- sprintf("%s_%s_%s_%s", residue_mapping$uniprot_ac, residue_mapping$pdb_filename, residue_mapping$chain, residue_mapping$res_num)
  seq_to_struct <- residue_mapping$struct_res_num
  names(seq_to_struct) <- residue_mapping$res_ID
  # struct to seq
  residue_mapping$res_struct_ID <- sprintf("%s_%s_%s_%s", residue_mapping$uniprot_ac, residue_mapping$pdb_filename, residue_mapping$chain, residue_mapping$struct_res_num)
  struct_to_seq <- residue_mapping$res_num
  names(struct_to_seq) <- residue_mapping$res_struct_ID
  
  # sanity check
  if (sum(duplicated(residue_mapping$res_ID)) > 0) { stop("\nnot unique IDs\n")}
  if (sum(duplicated(residue_mapping$res_struct_ID)) > 0) { stop("\nnot unique IDs\n")}
  
  # map interface residue in the complex struct to the sequence (INTERFACE1)
  tmp_IDs <- sprintf("%s_%s_%s_%s", merged_data$PROT1, merged_data$FILENAME, "A", merged_data$iface_res1)
  iface_res_seq1 <- struct_to_seq[ tmp_IDs ]
  # get residue in the monomer structures (N1, C1 use iface_res1)
  tmp_IDs <- sprintf("%s_%s_%s_%s", merged_data$PROT1, merged_data$mo_Nterm_PDB1, merged_data$mo_Nterm_chain1, iface_res_seq1)
  merged_data$mo_Nterm_iface1 <- seq_to_struct[ tmp_IDs ]
  tmp_IDs <- sprintf("%s_%s_%s_%s", merged_data$PROT1, merged_data$mo_Cterm_PDB1, merged_data$mo_Cterm_chain1, iface_res_seq1)
  merged_data$mo_Cterm_iface1 <- seq_to_struct[ tmp_IDs ]

  # map interface residue in the complex struct to the sequence (INTERFACE2)
  tmp_IDs <- sprintf("%s_%s_%s_%s", merged_data$PROT2, merged_data$FILENAME, "B", merged_data$iface_res2)
  iface_res_seq2 <- struct_to_seq[ tmp_IDs ]
  # get residue in the monomer structures (N2, C2 use iface_res2)
  tmp_IDs <- sprintf("%s_%s_%s_%s", merged_data$PROT2, merged_data$mo_Nterm_PDB2, merged_data$mo_Nterm_chain2, iface_res_seq2)
  merged_data$mo_Nterm_iface2 <- seq_to_struct[ tmp_IDs ]
  tmp_IDs <- sprintf("%s_%s_%s_%s", merged_data$PROT2, merged_data$mo_Cterm_PDB2, merged_data$mo_Cterm_chain2, iface_res_seq2)
  merged_data$mo_Cterm_iface2 <- seq_to_struct[ tmp_IDs ]

  ## map all residues at the interface
  iface_data <- f_read_file(iface_res_input_file, SKIP = 0, HEADER=FALSE, SEPARATOR = ",")
  colnames(iface_data) <- c("complex_PDB", "chain1", "iface1", "chain2", "iface2")
  iface_data$complex_PDB <- gsub("\\[|'", "", iface_data$complex_PDB)
  iface_data$chain1 <- gsub("'| ", "", iface_data$chain1)
  iface_data$iface1 <- gsub("'| ", "", iface_data$iface1)
  iface_data$chain2 <- gsub("'| ", "", iface_data$chain2)
  iface_data$iface2 <- gsub("\\]|'| ", "", iface_data$iface2)
  row.names(iface_data) <- iface_data$complex_PDB

  monomer_mapping <- c()
  for (i in 1:nrow(merged_data))
  {
    this_complex <- merged_data$FILENAME[[i]]

    ## map to sequence interface residues of complex subunit 1
    this_P1 <- merged_data$PROT1[[i]]
    this_iface1 <- iface_data$iface1[ iface_data$complex_PDB == this_complex ]
    this_residues1 <- strsplit(this_iface1, ";")[[1]]
    tmp_IDs <- sprintf("%s_%s_%s_%s", this_P1, this_complex, "A", this_residues1)
    iface_res_seq1 <- struct_to_seq[ tmp_IDs ]
    
    # monomer N1: get residue in the monomer structures
    this_mono_N1 <- merged_data$mo_Nterm_PDB1[[i]]
    this_chain_N1 <- merged_data$mo_Nterm_chain1[[i]]
    tmp_IDs <- sprintf("%s_%s_%s_%s", this_P1, this_mono_N1, this_chain_N1, iface_res_seq1)
    this_struct_iface_res_N1 <- unique(seq_to_struct[ tmp_IDs ])

    # monomer C1: get residue in the monomer structures
    this_mono_C1 <- merged_data$mo_Cterm_PDB1[[i]]
    this_chain_C1 <- merged_data$mo_Cterm_chain1[[i]]
    tmp_IDs <- sprintf("%s_%s_%s_%s", this_P1, this_mono_C1, this_chain_C1, iface_res_seq1)
    this_struct_iface_res_C1 <- unique(seq_to_struct[ tmp_IDs ])

    ## map to sequence interface residues of complex subunit 2
    this_P2 <- merged_data$PROT2[[i]]
    this_iface2 <- iface_data$iface2[ iface_data$complex_PDB == this_complex ]
    this_residues2 <- strsplit(this_iface2, ";")[[1]]
    tmp_IDs <- sprintf("%s_%s_%s_%s", this_P2, this_complex, "B", this_residues2)
    iface_res_seq2 <- struct_to_seq[ tmp_IDs ]
    
    # monomer N2: get residue in the monomer structures
    this_mono_N2 <- merged_data$mo_Nterm_PDB2[[i]]
    this_chain_N2 <- merged_data$mo_Nterm_chain2[[i]]
    tmp_IDs <- sprintf("%s_%s_%s_%s", this_P2, this_mono_N2, this_chain_N2, iface_res_seq2)
    this_struct_iface_res_N2 <- unique(seq_to_struct[ tmp_IDs ])
    
    # monomer C2: get residue in the monomer structures
    this_mono_C2 <- merged_data$mo_Cterm_PDB2[[i]]
    this_chain_C2 <- merged_data$mo_Cterm_chain2[[i]]
    tmp_IDs <- sprintf("%s_%s_%s_%s", this_P2, this_mono_C2, this_chain_C2, iface_res_seq2)
    this_struct_iface_res_C2 <- unique(seq_to_struct[ tmp_IDs ])

    monomer_mapping[[ sprintf("%s_N1", this_complex) ]] <- paste(this_struct_iface_res_N1, collapse=";")
    monomer_mapping[[ sprintf("%s_C1", this_complex) ]] <- paste(this_struct_iface_res_C1, collapse=";")
    monomer_mapping[[ sprintf("%s_N2", this_complex) ]] <- paste(this_struct_iface_res_N2, collapse=";")
    monomer_mapping[[ sprintf("%s_C2", this_complex) ]] <- paste(this_struct_iface_res_C2, collapse=";")
  }
  to_print <- cbind(names(monomer_mapping), monomer_mapping)
  f_print_matrix(to_print, OUTPUT_FILE = iface_res_monomer_output_file, colnames = FALSE, rownames = FALSE)
  
  return (merged_data)
}


