#source("https://bioc.ism.ac.jp/biocLite.R")
#biocLite("Biobase")
library(Biobase)

#load("../external_data/gtex.rdata")
gtex = readRDS("../external_data/gtex.rds")
# read the normalized expression values into a matrix called nmat
nmat = assayData(gtex)[["normalizedMatrix"]]
# get a unique vector of the subtypes
subtypes = unique(pData(gtex)$our_subtypes)
# generate subfolder to save all sample expression files
dir.create(file.path('../external_data/', 'GTEx_sample_data'))

for (subtype in subtypes) {

# get a submatrix with all samples from the subtype
nmat_sub = nmat[,which(pData(gtex)$our_subtypes == subtype)]
write.table(nmat_sub,paste("../external_data/GTEx_sample_data/",subtype,".txt",sep=""),quote=FALSE,sep="\t",col.names=NA)

}
