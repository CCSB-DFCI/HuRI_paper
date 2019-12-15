
f_read_Interactome3D_interactions_file <- function(file)
{
  C_INTERACTIONS_FILE_COLUMNS <- c("character", # PROT1
                                   "character", # PROT2
                                   "numeric",   # RANK_MAJOR
                                   "numeric",   # RANK_MINOR
                                   "character", # TYPE
                                   "character", # PDB_ID
                                   "character",   # BIO_UNIT
                                   "character", # CHAIN1
                                   "character",   # MODEL1
                                   "numeric",   # SEQ_IDENT1
                                   "numeric",   # COVERAGE1
                                   "numeric",   # SEQ_BEGIN1
                                   "numeric",   # SEQ_END1
                                   "character",   # DOMAIN1
                                   "character", # CHAIN2
                                   "character",   # MODEL2
                                   "numeric",   # SEQ_IDENT2
                                   "numeric",   # COVERAGE2
                                   "numeric",   # SEQ_BEGIN2
                                   "numeric",   # SEQ_END2
                                   "character",   # DOMAIN2
                                   "character") # FILENAME

  data <- read.table(file, sep="\t", check.names=FALSE, header=TRUE, comment.char="", skip=0, quote="", colClasses=C_INTERACTIONS_FILE_COLUMNS, strip.white=FALSE, na.strings = NULL)
  data <- unique(data)
  data$protein1 <- apply(data[,c("PROT1", "PROT2")], 1, function(x) { sort(x)[[1]] } )
  data$protein2 <- apply(data[,c("PROT1", "PROT2")], 1, function(x) { sort(x)[[2]] } )
  data$int_ID <- sprintf("%s_%s", data$protein1, data$protein2)
  experimental <- unique(data$int_ID[ data$TYPE == "Structure"])
  model <- unique(data$int_ID[ data$TYPE == "Model"])
  dom_dom_model <- unique(data$int_ID[ data$TYPE == "Dom_dom_model"])
  
  to_return <- list("experimental" = experimental,
                    "model" = model,
                    "dom_dom_model" = dom_dom_model,
                    "raw_data" = data)
  return (to_return)
}