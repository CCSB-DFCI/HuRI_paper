#install.packages("downloader", repos="http://cran.us.r-project.org")
library(downloader)
tmp = tempfile()
src = "http://networkmedicine.org:3838/gtex_data/gtex_portal_normalized.rds"
download(src,tmp)
obj = readRDS(tmp)
saveRDS(obj,file='../external_data/gtex.rds')
