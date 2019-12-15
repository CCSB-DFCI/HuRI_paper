f_run_terminal_fusion_interface_intereference_analysis <- function(interactome3D_version    = "2018_04",
                                                                   MONOMER_MIN_COVERAGE     = 50,
                                                                   B_INCLUDE_HOMODIMERS     = TRUE,
                                                                   B_INCLUDE_PPI_MODELS     = TRUE,
                                                                   B_INCLUDE_PPI_DDM        = TRUE,
                                                                   B_INCLUDE_PROTEIN_MODELS = TRUE)
{
  # load huri data
  huri_data <- faux_TF_load_huri_data()
  TMP_OUTPUT_1 <- "part1_output.csv"
  TMP_OUTPUT_2 <- "part2_output.csv"
  TMP_OUTPUT_3 <- "part3_output.csv"
  TMP_OUTPUT_4 <- "part4_output.csv"
  
  # PART1: get structural data. PPI structures/models/DDM + best suited monomers + terminal residue information
  if (!file.exists(TMP_OUTPUT_1))
  {
    # load dimers from interactome3D
    struct_data <- faux_read_structural_data(interactome3D_version)

    intx_bw_tested <- faux_get_interactome3D_dimers_bw_huri_tested(huri_data$tested_info, struct_data$interaction_data, struct_data$seq_length, B_INCLUDE_HOMODIMERS = B_INCLUDE_HOMODIMERS, B_INCLUDE_PPI_MODELS = B_INCLUDE_PPI_MODELS, B_INCLUDE_DDM=B_INCLUDE_PPI_DDM)

    # add info for best monomers (best monomers for Nterm and Cterm may be different)
    intx_bw_tested <- faux_add_best_monomers(intx_bw_tested, struct_data$protein_data, MIN_COVERAGE=MONOMER_MIN_COVERAGE, B_INCLUDE_MODELS=B_INCLUDE_PROTEIN_MODELS)

    # get Nterm and Cterm residues (for complex data and monomers)
    intx_bw_tested <- faux_add_terminal_residues_struct(intx_bw_tested, interactome3D_version)

    f_print_matrix(intx_bw_tested, OUTPUT_FILE=TMP_OUTPUT_1, colnames = TRUE, rownames = FALSE)
  }else
  {
    intx_bw_tested <- f_read_file(TMP_OUTPUT_1, HEADER=TRUE,
                                  colClasses=c(rep("character", 4),            # PROT1/2, TYPE, PDB_ID
                                               "character", rep("numeric", 4), # CHAIN1-SEQ_END1
                                               "character", rep("numeric", 4), # CHAIN2-SEQ_END2
                                               "character",                    # FILENAME
                                               "character", rep("numeric", 2), # int_ID, seq_length1/2
                                               rep("numeric", 4),              # co_chop's
                                               rep("character", 3), rep("numeric", 2), # mo_PDB_N1
                                               rep("character", 3), rep("numeric", 2), # mo_PDB_C1
                                               rep("character", 3), rep("numeric", 2), # mo_PDB_N2
                                               rep("character", 3), rep("numeric", 2), # mo_PDB_C2
                                               rep("character", 8),            # co_term_residues (and +1)
                                               rep("character", 8)             # mo_term_residues (and +1)) )
                                  ) )
  }
  cat("PART 1 done\n")
  
  ## PART2: calculate distances using complex structures (this includes calculating interface residues)
  # calculate distance from terminals to interface using complex subunit structures
  # it also gets interface residues
  if (!file.exists(TMP_OUTPUT_2))
  {
    merged_data <- faux_get_distance_complex_structures(intx_bw_tested)
    f_print_matrix(merged_data, OUTPUT_FILE=TMP_OUTPUT_2, colnames = TRUE, rownames = FALSE)
  }else
  {
    merged_data <- f_read_file(TMP_OUTPUT_2, HEADER=TRUE,
                                             colClasses=c(rep("character", 5),            # FILENAME,PROT1/2, TYPE, PDB_ID
                                                          "character", rep("numeric", 4), # CHAIN1-SEQ_END1
                                                          "character", rep("numeric", 4), # CHAIN2-SEQ_END2
                                                          "character", rep("numeric", 2), # int_ID, seq_length1/2
                                                          rep("numeric", 4),              # co_chop's
                                                          rep("character", 3), rep("numeric", 2), # mo_PDB_N1
                                                          rep("character", 3), rep("numeric", 2), # mo_PDB_C1
                                                          rep("character", 3), rep("numeric", 2), # mo_PDB_N2
                                                          rep("character", 3), rep("numeric", 2), # mo_PDB_C2
                                                          rep("character", 4),            # co_term_residues
                                                          rep("character", 4),            # mo_term_residues
                                                          rep(c("character", "numeric"),2), #iface_res1, iface_dist1 (and 2)
                                                          rep("numeric", 4)              # co_d
                                                     ) )
  }
  cat("PART 2 done\n")
  
  ## PART3: calculate distances using monomer structures (this includes the mapping of interfaces onto monomers)
  # map interface residues to monomer structures
  if (!file.exists(TMP_OUTPUT_3))
  {
    merged_data <- faux_map_iface_residues_onto_monomers(merged_data)

    # calculate distance from terminals to interface using monomer structures
    merged_data <- faux_get_distance_monomer_structures(merged_data)
    f_print_matrix(merged_data, OUTPUT_FILE=TMP_OUTPUT_3, colnames = TRUE, rownames = FALSE)
  }else
  {
    merged_data <- f_read_file(TMP_OUTPUT_3, HEADER=TRUE, colClasses=faux_get_colClasses()[1:114])
  }
  cat("PART 3 done\n")
  
  ## PART4: add directionality (DB or AD), assay information, tail orientation, and get "best" measurements (either complex or monomer data)
  ## For the pairs for which we have PPI structure --> get all tested_pairs in huri
  if (!file.exists(TMP_OUTPUT_4))
  {
    tested_pairs <- faux_add_sense_of_pairs_tested(merged_data, huri_data$tested_info)
    # detectability per assay
    tested_pairs$int_sense <- sprintf("%s_%s", tested_pairs$PROT1, tested_pairs$PROT2)
    tested_pairs$in_v1_assay <- (tested_pairs$int_sense %in% huri_data$found_by_assay$pairs_by_assay[["1"]])
    tested_pairs$in_v2_assay <- (tested_pairs$int_sense %in% huri_data$found_by_assay$pairs_by_assay[["2"]])
    tested_pairs$in_v6_assay <- (tested_pairs$int_sense %in% huri_data$found_by_assay$pairs_by_assay[["6"]])
    tested_pairs$P1_in_huri <- tested_pairs$PROT1 %in% huri_data$huri_positives
    tested_pairs$P2_in_huri <- tested_pairs$PROT2 %in% huri_data$huri_positives

    # get "best" distance --> the one w/ the shortest tail chop
    tested_pairs <- faux_add_best_distance(tested_pairs, B_DB_OR_AD=TRUE, B_N_OR_C=TRUE)
    tested_pairs <- faux_add_best_distance(tested_pairs, B_DB_OR_AD=TRUE, B_N_OR_C=FALSE)
    tested_pairs <- faux_add_best_distance(tested_pairs, B_DB_OR_AD=FALSE, B_N_OR_C=TRUE)
    tested_pairs <- faux_add_best_distance(tested_pairs, B_DB_OR_AD=FALSE, B_N_OR_C=FALSE)
    f_print_matrix(tested_pairs, OUTPUT_FILE=TMP_OUTPUT_4, colnames = TRUE, rownames = FALSE)
  }else
  {
    tested_pairs <- f_read_file(TMP_OUTPUT_4, colClasses="character")
  }
  cat("PART 4 done\n")
  
  ## PART5: analysis
  f_print_matrix(tested_pairs, "results.csv", colnames = TRUE, rownames = FALSE)
  f_complete_analysis(tested_pairs, OUTPUT_FOLDER = "pdfs")
}

faux_read_results <- function(file)
{
  colClasses <- c(rep("character", 5),            # FILENAME....PDB_ID
                  "character", rep("numeric", 4), # CHAIN1-SEQ_END1
                  "character", rep("numeric", 4), # CHAIN2-SEQ_END2
                  "character", rep("numeric", 2), # int_ID, seq_length1/2
                  rep("numeric", 4),              # co_chop's
                  rep("character", 3), rep("numeric", 2), # mo_PDB_N1
                  rep("character", 3), rep("numeric", 2), # mo_PDB_C1
                  rep("character", 3), rep("numeric", 2), # mo_PDB_N2
                  rep("character", 3), rep("numeric", 2), # mo_PDB_C2
                  rep("character", 8),             # co_term_residues (and +1)
                  rep("character", 8),             # mo_term_residues (and +1)
                  "character", "numeric", "character", "numeric", # iface residues
                  rep("numeric", 4),               # co terminal dist
                  rep("character", 4),              # mo iface res
                  rep("numeric", 4),               # mo terminal dist
                  "character", rep("logical", 5),  # int_sense, in_assay1/2/6, P1/2_in_huri
                  rep("numeric",2),  # 4x(best_N1/C1/N2/C2: chop, euc
                  rep("numeric",2),
                  rep("numeric",2),
                  rep("numeric",2) )
  data <- f_read_file(file=file, HEADER=TRUE, colClasses=colClasses, NA_STRINGS = "NA")
  return (data)
}

faux_TF_load_huri_data <- function()
{
  # read HuRI dataset
  huri <- f_load_HuRI(B_HI3 = TRUE) #HURI3_LATEST = TRUE)
  # read interactions found in each of the experimental setting
  found_by_assay <- f_read_interactions_found_per_assay_huri()
  # read which proteins could not be tested fused w/ DBD because of strong auto-activation
  tested_info <- f_read_tested_proteins_huri()
  return (list("huri_positives" = huri$uniprot$in_network,
               "tested_info"    = tested_info,
               "found_by_assay" = found_by_assay))
}

# output info that will be used by python/pymol script to calculate distances from terminals to interface center
faux_get_distance_complex_structures <- function(intx_bw_tested, python_input="input_distances_complex.csv", python_output="output_distances_complex.csv")
{
  fields_to_output <- c("FILENAME", "TYPE", 
                        "co_Nterm_res1", "co_Cterm_res1", "co_Nterm_res2", "co_Cterm_res2") #, # complex data (terminal residues)
  f_print_matrix(intx_bw_tested[,fields_to_output], OUTPUT_FILE = python_input, colnames = TRUE, rownames = FALSE)
  # 1) get closest residue to interface center
  # 2) calculate for all terminal residues their distance (euclidian and surface) to the interface center
  system("python terminal_interface_distance_part1.py > output_distances_complex")
  # read results 
  pymol_results <- f_read_file(python_output, SKIP=0, HEADER=TRUE, colClasses=c(rep("character", 2),    # PDB and type
                                                                                "character", "numeric", # iface res 1
                                                                                "character", "numeric", # iface res 2
                                                                                rep("numeric", 4)))#,      # complex euc dist from terminals to iface center
  pymol_results$PDB_type <- NULL
  fields <- colnames(pymol_results)[ grep("co_d", colnames(pymol_results)) ] 
  for (this_field in fields)
  {
    pymol_results[[ this_field ]][ pymol_results[[ this_field ]] == -1 ] <- NA
  }
  
  merged_data <- merge(intx_bw_tested, pymol_results, by.x = "FILENAME", by.y = "complex_PDB")
  return (merged_data)
}

faux_get_distance_monomer_structures <- function(merged_data, python_input = "input_distances_monomer.csv", python_output = "output_distances_monomer.csv")
{
  # prepare output for python script that calculates distances in monomers
  fields_to_output <- c("FILENAME",
                        "mo_Nterm_PDB1", "mo_Nterm_type1", "mo_Nterm_iface1", "mo_Nterm_res1", # monomer1 Nterm data 
                        "mo_Cterm_PDB1", "mo_Cterm_type1", "mo_Cterm_iface1", "mo_Cterm_res1", # monomer1 Cterm data 
                        "mo_Nterm_PDB2", "mo_Nterm_type2", "mo_Nterm_iface2", "mo_Nterm_res2", # monomer2 Nterm data 
                        "mo_Cterm_PDB2", "mo_Cterm_type2", "mo_Cterm_iface2", "mo_Cterm_res2") # monomer2 Cterm data 
  f_print_matrix(merged_data[,fields_to_output], OUTPUT_FILE = python_input, colnames = TRUE, rownames = FALSE)
  system("python terminal_interface_distance_part2.py > output_distances_monomer")
  # merge results: if dist == -1 --> NA
  pymol_results <- f_read_file(python_output, SKIP=0, HEADER=TRUE, colClasses=c("character", rep("numeric", 4)))
  fields <- colnames(pymol_results)[ grep("mo_d", colnames(pymol_results)) ] 
  for (this_field in fields)
  {
    pymol_results[[ this_field ]][ pymol_results[[ this_field ]] == -1 ] <- NA
  }
  
  merged_data <- merge(merged_data, pymol_results, by.x = "FILENAME", by.y = "complex_PDB")
  
  return (merged_data)
}

faux_read_structural_data <- function(interactome3D_version, species="human")
{
  interactome3D_data <- f_read_Interactome3D_for_species(species, B_REPRESENTATIVE_OR_COMPLETE = TRUE, interactome3D_version = interactome3D_version)
  interaction_data <- interactome3D_data$interactions$raw_data
  protein_data <- interactome3D_data$proteins$raw_data
  return (list("interaction_data" = interaction_data,
               "protein_data"     = protein_data,
               "seq_length"       = interactome3D_data$seq_length))
}

faux_get_colClasses <- function()
{
  colClasses <- c(rep("character", 5),            # FILENAME....PDB_ID
                  "character", rep("numeric", 4), # CHAIN1-SEQ_END1
                  "character", rep("numeric", 4), # CHAIN2-SEQ_END2
                  "character", rep("numeric", 2), # int_ID, seq_length1/2
                  rep("numeric", 4),              # co_chop's
                  rep("character", 3), rep("numeric", 2), # mo_PDB_N1
                  rep("character", 3), rep("numeric", 2), # mo_PDB_C1
                  rep("character", 3), rep("numeric", 2), # mo_PDB_N2
                  rep("character", 3), rep("numeric", 2), # mo_PDB_C2
                  rep("character", 4),             # co_term_residues
                  rep("character", 4),             # mo_term_residues
                  "character", "numeric", "character", "numeric", # iface residues
                  rep("numeric", 4),               # co terminal dist
                  rep("character", 4),              # mo iface res
                  rep("numeric", 4),               # mo terminal dist
                  "character", rep("logical", 5),  # int_sense, in_assay1/2/6, P1/2_in_huri
                  rep("numeric",2),  # 4x(best_N1/C1/N2/C2: chop, euc
                  rep("numeric",2),
                  rep("numeric",2),
                  rep("numeric",2) )
  return (colClasses)
}
