
f_complete_analysis <- function(tested_pairs, OUTPUT_FOLDER)
{
  dir.create(OUTPUT_FOLDER, recursive = TRUE, showWarnings = FALSE)
  for (B_POSITIVES_IN_HURI_ONLY in c(TRUE, FALSE))
  {
    for (B_COMPLEX_OR_MONOMER_DIST in c(TRUE, FALSE, NA))
    {
      for (MAX_CHOPPED in c(0,1))
      {
          prot_info <- ifelse(is.na(B_COMPLEX_OR_MONOMER_DIST), "best", ifelse(B_COMPLEX_OR_MONOMER_DIST, "co", "mo"))
          bg_info <- ifelse(B_POSITIVES_IN_HURI_ONLY, "positives", "tested")
          OUTPUT_PDF <- sprintf("euc_%s_chop%s.bg_%s.pdf", prot_info, MAX_CHOPPED, bg_info)
          R_object <- sprintf("euc_%s_chop%s.bg_%s.R_object", prot_info, MAX_CHOPPED, bg_info)
          faux_analysis_euclidian_distance(tested_pairs, B_COMPLEX_OR_MONOMER_DIST=B_COMPLEX_OR_MONOMER_DIST, MAX_CHOPPED=MAX_CHOPPED, B_POSITIVES_IN_HURI_ONLY = B_POSITIVES_IN_HURI_ONLY, OUTPUT_FOLDER=OUTPUT_FOLDER, OUTPUT_PDF = OUTPUT_PDF, R_object=R_object)
      }
    }
  }
}

# B_COMPLEX_OR_MONOMER_DIST -> TRUE: complex; FALSE: monomer; NA: combine
faux_analysis_euclidian_distance <- function(tested_pairs, B_COMPLEX_OR_MONOMER_DIST, MAX_CHOPPED, B_POSITIVES_IN_HURI_ONLY, OUTPUT_FOLDER="./", OUTPUT_PDF=NULL, R_object=NULL) # B_TO_IFACE,
{
  faux_analysis_distance_to_iface(tested_pairs, B_COMPLEX_OR_MONOMER_DIST, MAX_CHOPPED, B_POSITIVES_IN_HURI_ONLY, B_RANKSUM_OR_FISHER=TRUE, OUTPUT_FOLDER=OUTPUT_FOLDER, OUTPUT_PDF=sprintf("ranksum_%s", OUTPUT_PDF), R_object=sprintf("ranksum_%s", R_object))
  for (dist_cutoff in c(10,20))
  {
    faux_analysis_distance_to_iface(tested_pairs, B_COMPLEX_OR_MONOMER_DIST, MAX_CHOPPED, B_POSITIVES_IN_HURI_ONLY, B_RANKSUM_OR_FISHER=FALSE, OUTPUT_FOLDER=OUTPUT_FOLDER, FISHER_DISTANCE_CUTOFF=dist_cutoff, OUTPUT_PDF=sprintf("fisher%s_%s", dist_cutoff, OUTPUT_PDF), R_object=sprintf("fisher%s_%s", dist_cutoff, R_object))
  }
}

faux_analysis_distance_to_iface <- function(tested_pairs, 
                                            B_COMPLEX_OR_MONOMER_DIST, MAX_CHOPPED, B_POSITIVES_IN_HURI_ONLY,  
                                            B_RANKSUM_OR_FISHER, FISHER_DISTANCE_CUTOFF,
                                            OUTPUT_FOLDER="./", OUTPUT_PDF=NULL, R_object=NULL)
{
  if (!is.null(OUTPUT_PDF)) { pdf(sprintf("%s/%s", OUTPUT_FOLDER, OUTPUT_PDF), width=10, height=5) }
  par("mfrow"=c(1,4))
  all_res <- list()
  YMAXs <- c(NA) # dummy for fisher test and no YMAX for ranksum
  if (B_RANKSUM_OR_FISHER) { YMAXs <- c(YMAXs, 100) }
  for (this_assay in c(1,2,6))
  {
   for (YMAX in YMAXs)
   {
    for (B_DB_OR_AD in c(TRUE, FALSE))
    {
      for (NTERM_OR_CTERM in c(TRUE, FALSE))
      {
        this_res <- faux_single_analysis(tested_pairs, this_assay, B_DB_OR_AD, NTERM_OR_CTERM, MAX_CHOPPED, B_POSITIVES_IN_HURI_ONLY, B_COMPLEX_OR_MONOMER_DIST, B_RANKSUM_OR_FISHER, FISHER_DISTANCE_CUTOFF, OUTPUT_PDF=NULL, YMAX=YMAX)
        this_ID <- sprintf("%s_%s_%s", as.character(this_assay), as.character(B_DB_OR_AD), as.character(NTERM_OR_CTERM))
        all_res[[ this_ID ]] <- this_res
      }
    }
   }
  }
  if (!is.null(OUTPUT_PDF)) { dev.off() }
  if (!is.null(R_object)) { save(all_res, file=sprintf("%s/%s", OUTPUT_FOLDER, R_object)) }
}


### analysis & plotting function: distance from terminal of one protein (AD or DB) to the interface center
# MAX_CHOPPED: max residues chopped in tail
# assay_ID: 1,2,or 6
# B_DB_OR_AD: analysis on the DB protein or AD protein
faux_single_analysis <- function(tested_pairs, assay_ID, 
                                 B_DB_OR_AD, NTERM_OR_CTERM, MAX_CHOPPED, B_POSITIVES_IN_HURI_ONLY, B_COMPLEX_OR_MONOMER_DIST,
                                 B_RANKSUM_OR_FISHER, FISHER_DISTANCE_CUTOFF=NA,
                                 OUTPUT_PDF=NULL, YMAX=NULL)
{
  assay_field <- sprintf("in_v%s_assay", assay_ID)
  if (B_POSITIVES_IN_HURI_ONLY)
  {
    tested_pairs <- tested_pairs[ tested_pairs$P1_in_huri & tested_pairs$P2_in_huri,]
  }
  # mono_DB_Nterm_chop, complex_DB_Nterm_chop, best_DB_Nterm_chop
  struct_info <- ifelse(is.na(B_COMPLEX_OR_MONOMER_DIST), "best", ifelse(B_COMPLEX_OR_MONOMER_DIST, "co", "mo"))
  prot_info <- ifelse(B_DB_OR_AD, "1", "2")
  term_info <- ifelse(NTERM_OR_CTERM, "N", "C")
  chopped_field <- sprintf("%s_%sterm_chop%s", struct_info, term_info, prot_info)
  dist_info <- "d" 
  distance_field <- sprintf("%s_%s_%s%s", struct_info, dist_info, term_info, prot_info)
  XLAB <- sprintf("Identified in assay %s?", assay_ID)
  all_res <- list()
  if (!is.null(OUTPUT_PDF)) { pdf(OUTPUT_PDF, width=4) }
  # detectability and distance of A_N_terminal to interface
  data_ok <- tested_pairs[ tested_pairs[[ chopped_field ]] <= MAX_CHOPPED & !is.na(tested_pairs[[ chopped_field ]]),]
  detected <- data_ok[ data_ok[[ assay_field ]],]
  missed <- data_ok[ !data_ok[[ assay_field ]],]
  if (B_RANKSUM_OR_FISHER)
  {
    res <- f_calculate_ranksum_test(detected[[ distance_field ]], missed[[ distance_field ]], PLOT = TRUE, ylab=sprintf("[%s_%s%s] Distance to iface", dist_info, term_info, prot_info), xlab=XLAB, xlab1="Yes", xlab2="No", YMIN = 0, YMAX = YMAX)
  }else
  {
    P1 <- sum(detected[[ distance_field ]] < FISHER_DISTANCE_CUTOFF, na.rm=TRUE)
    T1 <- sum(!is.na(detected[[ distance_field ]]))
    P2 <- sum(missed[[ distance_field ]] < FISHER_DISTANCE_CUTOFF, na.rm=TRUE)
    T2 <- sum(!is.na(missed[[ distance_field ]]))
    res <- NA
    if (T1 != 0 & T2 != 0)
    {
      res <- f_calculate_fisher_test_num(P1, T1, P2, T2, bg_includes_group = FALSE, PLOT = TRUE, ylab=sprintf("[%s_%s%s] Fraction at < %s A", dist_info, term_info, prot_info, FISHER_DISTANCE_CUTOFF), xlab=XLAB, xlab1="Yes", xlab2="No", B_SHOW1=FALSE, B_SHOW2=FALSE, B_SHOW3=TRUE)
    }
  }
  if (!is.null(OUTPUT_PDF)) { dev.off() }
  return (res)
}