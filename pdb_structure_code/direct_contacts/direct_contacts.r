source("lib_aux.r", chdir=TRUE)

## read structural data
species_int3d <- f_read_Interactome3D_interactions_file("interactions.dat")
struct_int3d <- species_int3d$raw_data[ species_int3d$raw_data$TYPE == "Structure",]

## keep the biounit w/ the smaller ID
to_keep <- lapply(split(struct_int3d, struct_int3d$PDB_ID),
                  function(x)
                  {
                    this_biounits <- unique(x[["BIO_UNIT"]])
                    selected_biounit <- ifelse(length(this_biounits) == 1, this_biounits, min(as.numeric(this_biounits)))
                    x[ x[["BIO_UNIT"]] == selected_biounit,]
                  })
struct_int3d <- do.call("rbind", to_keep)
struct_int3d$PDB_BIO <- sprintf("%s_%s", struct_int3d$PDB_ID, struct_int3d$BIO_UNIT)

## consider only one protein per complex chain: gather all info on proteins per complex chain + select the best and remove data from rest
proteins_per_chain1 <- apply(struct_int3d, 1, function(x) { c(x[["PDB_BIO"]], x[["CHAIN1"]], x[["PROT1"]], x[["SEQ_IDENT1"]], x[["COVERAGE1"]]) })
proteins_per_chain2 <- apply(struct_int3d, 1, function(x) { c(x[["PDB_BIO"]], x[["CHAIN2"]], x[["PROT2"]], x[["SEQ_IDENT2"]], x[["COVERAGE2"]]) })
proteins_per_chain <- as.data.frame(t(cbind(proteins_per_chain1, proteins_per_chain2)), stringsAsFactors = FALSE)
colnames(proteins_per_chain) <- c("PDB", "chain", "protein", "seqID", "coverage")
proteins_per_chain$PDB_chain <- sprintf("%s_%s", proteins_per_chain$PDB, proteins_per_chain$chain)
proteins_per_chain$seqID <- as.numeric(proteins_per_chain$seqID)
proteins_per_chain$coverage <- as.numeric(proteins_per_chain$coverage)
proteins_per_chain <- lapply(split(proteins_per_chain, proteins_per_chain$PDB_chain), function(x) { unique(x) })
# select protein w/ highest SEQ ID (first) and coverage (Second)
proteins_per_chain <- lapply(proteins_per_chain, function(x) { x[ order(x[["seqID"]], x[["coverage"]], decreasing = TRUE), ]})
to_remove <- lapply(proteins_per_chain, function(x) { if (nrow(x) > 1) { x[2:nrow(x),] } })
to_remove <- to_remove[ !(sapply(to_remove, is.null)) ]
to_remove <- do.call("rbind", to_remove)
to_remove$to_delete <- sprintf("%s_%s_%s", to_remove$PDB, to_remove$chain, to_remove$protein)

struct_int3d$superID1 <- sprintf("%s_%s_%s", struct_int3d$PDB_BIO, struct_int3d$CHAIN1, struct_int3d$PROT1)
struct_int3d$superID2 <- sprintf("%s_%s_%s", struct_int3d$PDB_BIO, struct_int3d$CHAIN2, struct_int3d$PROT2)
struct_int3d <- struct_int3d[! struct_int3d$superID1 %in% to_remove$to_delete & !struct_int3d$superID2 %in% to_remove$to_delete,]

## get partners per protein
by_prot2 <- split(struct_int3d$PROT1, struct_int3d$PROT2)
by_prot1 <- split(struct_int3d$PROT2, struct_int3d$PROT1)

proteins_w_struct_info <- unique(c(names(by_prot1), names(by_prot2)))
by_prot <- list()
for (this_id in proteins_w_struct_info)
{
  by_prot[[this_id]] <- unique(c(by_prot1[[this_id]], by_prot2[[this_id]]))
}

## split structural data by PDBs
by_PDB <- split(struct_int3d, struct_int3d$PDB_BIO)

## get for each PDB all the contacts for each protein
chains_per_PDB <- list()   # chains that are part of the complex (A1, A2, and B --> 3 different chains: A1, A2, and B)
proteins_per_PDB <- list() # proteins that are part of the complex (A1, A2, and B --> 2 different proteins: A and B)
contacts_per_PDB <- list() # for each complex: for each protein, list of proteins that are in direct contact (no repeats)
real_chains_per_PDB <- list() ### the chain IDs (A,B,C,...) that are part of the complex
for (this_PDB_BIO in names(by_PDB))
{
  this_data <- by_PDB[[this_PDB_BIO]]
  this_data$P1 <- sprintf("%s_%s%s", this_data$PROT1, this_data$CHAIN1, this_data$MODEL1)
  this_data$P2 <- sprintf("%s_%s%s", this_data$PROT2, this_data$CHAIN2, this_data$MODEL2)
  this_chains <- unique(c(this_data$P1, this_data$P2))
  this_proteins <- unique(c(this_data$PROT1, this_data$PROT2))
  chains_per_PDB[[this_PDB_BIO]] <- this_chains
  proteins_per_PDB[[this_PDB_BIO]] <- this_proteins
  real_chains_per_PDB[[this_PDB_BIO]] <- unique(c(this_data$CHAIN1, this_data$CHAIN2))

  # get list of contacts per protein in each PDB
  this_contacts_per_protein <- list()
  for (this_protein in this_proteins)
  {
    this_contacts1 <- this_data$PROT2[ this_data$PROT1 == this_protein]
    this_contacts2 <- this_data$PROT1[ this_data$PROT2 == this_protein]
    this_contacts <- unique(c(this_contacts1, this_contacts2))
    this_contacts_per_protein[[this_protein]] <- this_contacts
  }
  contacts_per_PDB[[this_PDB_BIO]] <- this_contacts_per_protein
}

## we keep only complexes w/ 3+ different proteins
n_uniq_chains_per_PDB <- sapply(real_chains_per_PDB, function(x) { length(unique(x)) })
proteins_per_PDB <- proteins_per_PDB[(sapply(proteins_per_PDB, length)) >= 3 & n_uniq_chains_per_PDB >= 3]

## calculate contacts and no-contacts for all protein pairs within each complex
done <- 0
contacts <- list()
separate <- list()
all_pairs <- list()
for (this_PDB_BIO in names(proteins_per_PDB))
{
  this_PDB_ID <- strsplit(this_PDB_BIO, "_")[[1]]
  this_n_chains <- length(chains_per_PDB[[this_PDB_BIO]])
  this_chains <- paste(chains_per_PDB[[this_PDB_BIO]], collapse=",")

  this_contacts <- contacts_per_PDB[[this_PDB_BIO]]

  ## print list of all possible pairs per complex and whether contact exists or not  
  proteins_in_complex <- proteins_per_PDB[[this_PDB_BIO]]
  pair_info <- list()
  contact_info <- list()
  separate_info <- list()
  for (this_protein in names(this_contacts))
  {
    # get contacts
    this_protein_contacts <- this_contacts[[this_protein]]
    # get no contacts
    this_protein_separate <- proteins_in_complex[ ! proteins_in_complex %in% this_protein_contacts]
    this_contact_pairs <- expand.grid(this_protein, this_protein_contacts)
    this_contact_pairs$ID <- apply(this_contact_pairs, 1, function(x) { tmp <- sort(c(x[[1]], x[[2]])); sprintf("%s_%s", tmp[[1]], tmp[[2]]) })
    this_contact_pairs$in_contact <- TRUE
    this_contact_pairs$support <- this_PDB_BIO
    this_separate_pairs <- expand.grid(this_protein, this_protein_separate)
    this_separate_pairs$ID <- apply(this_separate_pairs, 1, function(x) { tmp <- sort(c(x[[1]], x[[2]])); sprintf("%s_%s", tmp[[1]], tmp[[2]]) })
    if (nrow(this_separate_pairs) > 0)
    {
      this_separate_pairs$in_contact <- FALSE
      this_separate_pairs$support <- this_PDB_BIO
    }
    this_info <- rbind(this_contact_pairs, this_separate_pairs)
    pair_info[[this_protein]] <- this_info
    contact_info[[this_protein]] <- this_contact_pairs
    separate_info[[this_protein]] <- this_separate_pairs
  }
  tmp_contacts <- do.call("rbind", contact_info)
  contacts[[this_PDB_BIO]] <- tmp_contacts[! duplicated(tmp_contacts$ID),]
  tmp_separate <- do.call("rbind", separate_info)
  separate[[this_PDB_BIO]] <- tmp_separate[! duplicated(tmp_separate$ID),]
  tmp_all_pairs <- do.call("rbind", pair_info)
  all_pairs[[this_PDB_BIO]] <- tmp_all_pairs[! duplicated(tmp_all_pairs$ID),]
  
  done <- done+1
  if (done %% 100 == 0)
  {
    cat(sprintf("Done %s out of %s\n", done, length(proteins_per_PDB)))
  }
}

all_contacts <- do.call("rbind", contacts)
all_indirect <- do.call("rbind", separate)

## merge list of direct contacts and get info on which PDB supports it
unique_direct <- all_contacts[!duplicated(all_contacts$ID),(1:4)]
row.names(unique_direct) <- unique_direct$ID
tmp_support <- sapply(split(all_contacts$support, all_contacts$ID), paste, collapse=",")
unique_direct$support <- ""
unique_direct[names(tmp_support),"support"] <- tmp_support

## merge list of no-contacts between proteins of the same complex
unique_indirect <- all_indirect[!duplicated(all_indirect$ID),(1:4)]
row.names(unique_indirect) <- unique_indirect$ID
tmp_support <- sapply(split(all_indirect$support, all_indirect$ID), paste, collapse=",")
unique_indirect$support <- ""
unique_indirect[names(tmp_support),"support"] <- tmp_support

## remove self-pairs info (it is not accurate: self-pairs are evaluated also in complexes w/o homodimers)
DIRECT_FILE <- "direct_pairs.xls"
unique_direct <- unique_direct[ unique_direct$Var1 != unique_direct$Var2,]
cat(file=DIRECT_FILE, append = FALSE, "protein1\tprotein2\tpair_ID\tin_contact\tsupport\n")
cat(file=DIRECT_FILE, append = TRUE, paste(apply(unique_direct, 1, paste, collapse="\t"), collapse="\n"))
cat(file=DIRECT_FILE, append = TRUE, "\n")

## remove indirect pairs that are direct
INDIRECT_FILE <- "indirect_pairs.xls"
unique_indirect <- unique_indirect[ unique_indirect$Var1 != unique_indirect$Var2,]
unique_indirect <- unique_indirect[ !unique_indirect$ID %in% unique_direct$ID,]
cat(file=INDIRECT_FILE, append = FALSE, "protein1\tprotein2\tpair_ID\tin_contact\tsupport\n")
cat(file=INDIRECT_FILE, append = TRUE, paste(apply(unique_indirect, 1, paste, collapse="\t"), collapse="\n"))
cat(file=INDIRECT_FILE, append = TRUE, "\n")
