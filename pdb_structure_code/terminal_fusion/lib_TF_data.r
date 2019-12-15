faux_merge_w_monomer_distances_get_monomer_data <- function(merged_data, monomer_dist)
{
  # for each line in merged data, get info for monomer
  merged_monomer <- as.data.frame( t( apply(merged_data, 1, function(x)
  {
    # DB N term
    tmp <- monomer_dist[ monomer_dist$complex_PDB == x[["FILENAME"]] & 
                           monomer_dist$monomer_PDB == x[["DB"]] &
                           monomer_dist$term_ID == "Nterm",]
    int_sense <- x[["int_sense"]]
    monomer_type <- ifelse(nrow(tmp) > 0, tmp$monomer_type[[1]], NA)
    d_DB_Nterm_mono_iface <- ifelse(nrow(tmp) > 0, tmp$d_term_iface[[1]], NA)
    s_DB_Nterm_mono_iface <- ifelse(nrow(tmp) > 0, tmp$s_term_iface[[1]], NA)
    DB_Nterm_res <- ifelse(nrow(tmp) > 0, tmp$term_res[[1]], NA)
    # DB C term
    tmp <- monomer_dist[ monomer_dist$complex_PDB == x[["FILENAME"]] & 
                           monomer_dist$monomer_PDB == x[["DB"]] &
                           monomer_dist$term_ID == "Cterm",]
    d_DB_Cterm_mono_iface <- ifelse(nrow(tmp) > 0, tmp$d_term_iface[[1]], NA)
    s_DB_Cterm_mono_iface <- ifelse(nrow(tmp) > 0, tmp$s_term_iface[[1]], NA)
    DB_Cterm_res <- ifelse(nrow(tmp) > 0, tmp$term_res[[1]], NA)
    # AD N term
    tmp <- monomer_dist[ monomer_dist$complex_PDB == x[["FILENAME"]] & 
                           monomer_dist$monomer_PDB == x[["AD"]] &
                           monomer_dist$term_ID == "Nterm",]
    d_AD_Nterm_mono_iface <- ifelse(nrow(tmp) > 0, tmp$d_term_iface[[1]], NA)
    s_AD_Nterm_mono_iface <- ifelse(nrow(tmp) > 0, tmp$s_term_iface[[1]], NA)
    AD_Nterm_res <- ifelse(nrow(tmp) > 0, tmp$term_res[[1]], NA)
    # AD C term
    tmp <- monomer_dist[ monomer_dist$complex_PDB == x[["FILENAME"]] & 
                           monomer_dist$monomer_PDB == x[["AD"]] &
                           monomer_dist$term_ID == "Cterm",]
    d_AD_Cterm_mono_iface <- ifelse(nrow(tmp) > 0, tmp$d_term_iface[[1]], NA)
    s_AD_Cterm_mono_iface <- ifelse(nrow(tmp) > 0, tmp$s_term_iface[[1]], NA)
    AD_Cterm_res <- ifelse(nrow(tmp) > 0, tmp$term_res[[1]], NA)
    
    c(DB_Nterm_res, d_DB_Nterm_mono_iface, s_DB_Nterm_mono_iface,
      DB_Cterm_res, d_DB_Cterm_mono_iface, s_DB_Cterm_mono_iface,
      AD_Nterm_res, d_AD_Nterm_mono_iface, s_AD_Nterm_mono_iface,
      AD_Cterm_res, d_AD_Cterm_mono_iface, s_AD_Cterm_mono_iface, monomer_type, int_sense)
  }) ), stringsAsFactors=FALSE)
  colnames(merged_monomer) <- c("mono_DB_Nterm_res", "mono_d_DB_Nterm_iface", "mono_s_DB_Nterm_iface",
                                "mono_DB_Cterm_res", "mono_d_DB_Cterm_iface", "mono_s_DB_Cterm_iface",
                                "mono_AD_Nterm_res", "mono_d_AD_Nterm_iface", "mono_s_AD_Nterm_iface",
                                "mono_AD_Cterm_res", "mono_d_AD_Cterm_iface", "mono_s_AD_Cterm_iface", "monomer_type", "int_sense")
  # numeric values for some fields
  for (this_field in c("mono_d_DB_Nterm_iface", "mono_s_DB_Nterm_iface", "mono_d_DB_Cterm_iface", "mono_s_DB_Cterm_iface", "mono_d_AD_Nterm_iface", "mono_s_AD_Nterm_iface", "mono_d_AD_Cterm_iface", "mono_s_AD_Cterm_iface"))
  {
    merged_monomer[[ this_field ]] <- as.numeric(merged_monomer[[ this_field ]])
  }
  return (merged_monomer)
}

faux_merge_w_monomer_distances <- function(merged_data, monomer_dist_file)
{
  ## read distances using monomer data
  monomer_dist <- f_read_file("monomer_distances.csv", HEADER=TRUE, SKIP=0, colClasses=c(rep("character", 6), rep("numeric", 2)))
  
  # get monomer data for each line in merged data
  merged_monomer <- faux_merge_w_monomer_distances_get_monomer_data(merged_data, monomer_dist)
  
  # merge both datasets (merged_data and merged_monomer: they are already sorted)
  merged_data <- merge(merged_data, merged_monomer, by="int_sense")
  
  # add number of chopped residues at the tails
  merged_data$mono_DB_Nterm_chop <- faux_get_residue_number(merged_data$mono_DB_Nterm_res) - 1
  merged_data$mono_DB_Cterm_chop <- merged_data$seq_length1 - faux_get_residue_number(merged_data$mono_DB_Cterm_res)
  merged_data$mono_AD_Nterm_chop <- faux_get_residue_number(merged_data$mono_AD_Nterm_res) - 1
  merged_data$mono_AD_Cterm_chop <- merged_data$seq_length2 - faux_get_residue_number(merged_data$mono_AD_Cterm_res)
  
  # best distance from complex and monomers: that coming from the structure w/ the fewer chopped residues
  merged_data <- faux_add_best_distance(merged_data, B_DB_OR_AD=TRUE, B_N_OR_C=TRUE)
  merged_data <- faux_add_best_distance(merged_data, B_DB_OR_AD=TRUE, B_N_OR_C=FALSE)
  merged_data <- faux_add_best_distance(merged_data, B_DB_OR_AD=FALSE, B_N_OR_C=TRUE)
  merged_data <- faux_add_best_distance(merged_data, B_DB_OR_AD=FALSE, B_N_OR_C=FALSE)
  
  return (merged_data)
}

# add best of monomer and complex distances (we have always data for complex distances)
# also add data for tail orientation
faux_add_best_distance <- function(merged_data, B_DB_OR_AD, B_N_OR_C)
{
#  co_d_C1
  field_template <- "%s_%s_%s%s" # e.g. co_d_C1
  prot_info <- ifelse (B_DB_OR_AD, "1", "2")
  term_info <- ifelse (B_N_OR_C, "N", "C")
  d_complex <- sprintf(field_template, "co", "d", term_info, prot_info)
  complex_tail <- sprintf("co_%sterm_chop%s", term_info, prot_info)
  d_mono <- sprintf(field_template, "mo", "d", term_info, prot_info)
  mono_tail <- sprintf("mo_%sterm_chop%s", term_info, prot_info)
  d_best <- sprintf(field_template, "best", "d", term_info, prot_info)
  best_tail <- sprintf("best_%sterm_chop%s", term_info, prot_info)
  monomer_type <- sprintf("mo_%sterm_type%s", term_info, prot_info)

  # get cases where the complex data has precedence over the monomer data (if same or fewer chopped residues) OR if no data for monomer (check also if monomer is model)
  B_complex_best <- (is.na(merged_data[[ mono_tail ]]) |          # no PDB for monomer (tail chop is NA)
                     is.na(merged_data[[ d_mono ]]) |             # no distance for monomer (probably the interface residue is not in the monomer)
                     merged_data[[ complex_tail ]] <= merged_data[[ mono_tail ]])# |  # fewer chopped residues for complex tails
  merged_data[[ best_tail ]] <- ifelse(B_complex_best, merged_data[[ complex_tail ]], merged_data[[ mono_tail ]])
  merged_data[[ d_best ]] <- ifelse(B_complex_best, merged_data[[ d_complex ]], merged_data[[ d_mono ]])
  return (merged_data)
}

faux_get_residue_number <- function(residue)
{
  to_return <- ifelse(is.na(as.numeric(residue)), as.numeric(substring(residue, 1, nchar(residue)-1)), as.numeric(residue))
  return (to_return)
}

faux_monomers_to_do <- function(protein_info, interactome3D_data)
{
  i_todo <- 1
  to_do <- list()
  for (i_complex in 1:nrow(protein_info))
  {
    this_i <- (i_complex-1)*4
    PDB <- protein_info$PDB_file[[i_complex]]
    P1 <- protein_info$P1[[i_complex]]
    P1_iface_res <- protein_info$A_center_iface[[i_complex]]
    P1_best_monomers <- faux_get_best_monomers(P1, P1_iface_res, interactome3D_data, MIN_COVERAGE = 50, B_INCLUDE_MODELS = TRUE) 
    to_do[[this_i + 1]] <- c(PDB, P1, P1_best_monomers$Nterm[[1]], P1_best_monomers$Nterm[[2]], P1_iface_res, P1_best_monomers$Nterm[[4]], "Nterm")
    to_do[[this_i + 2]] <- c(PDB, P1, P1_best_monomers$Cterm[[1]], P1_best_monomers$Cterm[[2]], P1_iface_res, P1_best_monomers$Cterm[[4]], "Cterm")
    
    P2 <- protein_info$P2[[i_complex]]
    P2_iface_res <- protein_info$B_center_iface[[i_complex]]
    P2_best_monomers <- faux_get_best_monomers(P2, P2_iface_res, interactome3D_data, MIN_COVERAGE = 50, B_INCLUDE_MODELS = TRUE) 
    to_do[[this_i + 3]] <- c(PDB, P2, P2_best_monomers$Nterm[[1]], P2_best_monomers$Nterm[[2]], P2_iface_res, P2_best_monomers$Nterm[[4]], "Nterm")
    to_do[[this_i + 4]] <- c(PDB, P2, P2_best_monomers$Cterm[[1]], P2_best_monomers$Cterm[[2]], P2_iface_res, P2_best_monomers$Cterm[[4]], "Cterm")
  }
  df <- as.data.frame(do.call("rbind", to_do), stringsAsFactors = FALSE)
  colnames(df) <- c("complex_PDB", "monomer", "monomer_PDB", "type", "iface_res", "term_res", "term_type")
  return (df)
}

faux_read_pymol_results <- function(euc_file="euclidian_distances.csv", protein_file="protein_info.csv")
{
  euc_dist <- f_read_file(euc_file, SKIP = 0, HEADER=TRUE, colClasses = c("character", rep("numeric", 14)))
  euc_dist$P1 <- sapply(strsplit(euc_dist$PDB_file, "-"), "[[", 1)
  euc_dist$P2 <- sapply(strsplit(euc_dist$PDB_file, "-"), "[[", 2)
  protein_info <- f_read_file(protein_file, SKIP = 0, HEADER=TRUE, colClasses = c("character", rep( c(rep("numeric", 2), "character", "numeric", "character", "numeric", rep("character", 4)), 2)))
  protein_info$P1 <- sapply(strsplit(protein_info$PDB_file, "-"), "[[", 1)
  protein_info$P2 <- sapply(strsplit(protein_info$PDB_file, "-"), "[[", 2)
  protein_info$seq_length1 <- interactome3D_data$seq_length[ protein_info$P1 ]
  protein_info$seq_length2 <- interactome3D_data$seq_length[ protein_info$P2 ]
  to_return <- list("euc_dist"     = euc_dist,
                    "protein_info" = protein_info)
  return (to_return)
}