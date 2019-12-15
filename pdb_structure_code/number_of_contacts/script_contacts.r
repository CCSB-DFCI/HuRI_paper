
### read file
f_read_raw_file <- function(file, split_char=NULL)
{
  fc <- file(file)
  lines <- readLines(fc)
  close(fc)
  if (!is.null(split_char))
  {
    mylist <- strsplit(lines, split_char)
  }else
  {
    mylist <- lines
  }
  return (mylist)
}

### get contact raw data 
faux_read_contact_data <- function(contact_data_folder)
{
  contact_data_files <- list.files(contact_data_folder)
  all_data <- list()
  for (this_file in contact_data_files)
  {
    input_file <- sprintf("%s/%s", contact_data_folder, this_file)
    #all_data[[this_file]] <- f_read_file(input_file, colClasses = c("character", "numeric"), SKIP=1, SEPARATOR = ":")[[2]]
    tmp <- f_read_raw_file(input_file, split_char = "\t")
    all_data[[this_file]] <- sapply(tail(tmp,6), function(x) { as.numeric(tail(strsplit(x, " +")[[1]],1)) })
  }
  df <- as.data.frame(do.call("rbind", all_data), stringsAsFactors = FALSE)
  P1 <- sapply(strsplit(row.names(df), "-"), function(x) { x[[1]] })
  P2 <- sapply(strsplit(row.names(df), "-"), function(x) { x[[2]] })
  tmp_df <- as.data.frame(cbind(P1, P2), stringsAsFactors = FALSE)
  tmp_df$int_ID <- apply(tmp_df, 1, function(x) { tmp <- sort(c(x[[1]], x[[2]])); sprintf("%s_%s", tmp[[1]], tmp[[2]]) })
  df <- cbind(tmp_df, df, row.names(df))
  colnames(df) <- c("uniprot1", "uniprot2", "int_ID", "clashes", "disulfide_bridges", "h_bonds", "salt_bridges", "vdw_intx", "total", "PDB_file")
  row.names(df) <- NULL
  return (df)
}

contacts <- faux_read_contact_data(PATH_TO_RESULT_OF_RUN_SCRIPT)

contacts$PDB_file <- sapply(strsplit(as.character(contacts$PDB_file), ".contacts"), "[[", 1)
write.table(contacts, file="number_contacts.csv", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, append=FALSE)

