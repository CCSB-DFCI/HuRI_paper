# CODE: TF (terminal fusion)
# assay version 1: DB CEN N-ter, AD CEN N-ter
# assay version 2: DB CEN N-ter, AD 2u N-ter
# assay version 6: DB CEN N-ter, AD 2u C-ter
source("~/R_libs/lib_common.r")
source("~/R_libs/lib_huri.r")
source("../lib_TF_analysis.r")
source("../lib_TF_interactome3D.r")
source("../lib_TF_data.r")
source("../lib_TF_functions.r")

f_run_terminal_fusion_interface_intereference_analysis(interactome3D_version = "2018_04",
                                                       MONOMER_MIN_COVERAGE = 50,
                                                       B_INCLUDE_HOMODIMERS = TRUE,
                                                       B_INCLUDE_PPI_MODELS = TRUE,
                                                       B_INCLUDE_PPI_DDM = TRUE,
                                                       B_INCLUDE_PROTEIN_MODELS = TRUE,
						                                           B_DSYSMAP_OR_INT3D=FALSE)

