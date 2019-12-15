####################################################################################################################################
##### Code to reproduce Figure 4a ##################################################################################################
####################################################################################################################################

### R packages
library(RColorBrewer)
library(Homo.sapiens)
library(AnnotationDbi)

### Subcellular localization information: Cell Atlas (Nov 2017)
subcell <- read.table(file = "Subcellular_localization/subcellular_location.tsv",
                      sep = "\t", header = T)
table(subcell$Reliability)
# Approved Supported Uncertain Validated
# 6076      3600       779      1548

########
######## Formatting Cell Atlas data
## Validated, Supported, Approved subcellular location's categories considered as valid per each protein.
subcell <- data.frame(subcell, Subcellular.location = apply(subcell[,c("Validated","Supported","Approved")], 1, function(x)
  paste(x[x != ""],collapse = ";")))

## METACOMPARTMENTs of subcellular locations provided by Cell Atlas (used only to order)
categ <- list(Cytoplasm = c("Cytosol","Cytoplasmic bodies","Rods & Rings"),
              Midbody = c("Midbody","Cytokinetic bridge","Mitotic spindle","Midbody ring"),
              Nucleus = c("Nucleoplasm","Nuclear speckles","Nuclear bodies","Nucleoli","Nucleus",
                          "Nucleoli fibrillar center","Nuclear membrane"),
              Vesicles = c("Golgi apparatus","Vesicles","Lipid droplets","Peroxisomes","Lysosomes",
                           "Endosomes","Aggresome","Endoplasmic reticulum"),
              Membrane = c("Cell Junctions","Plasma membrane","Focal adhesion sites"),
              Mitochondria = c("Mitochondria"),
              Centrosome = c("Centrosome"),
              Cytoskeleton = c("Actin filaments","Microtubules","Microtubule organizing center","Microtubule ends",
                               "Intermediate filaments"))

## Counting proteins per subcellular location of Cell Atlas
p <- sapply(as.character(subcell$Subcellular.location), function(x) unlist(strsplit(x, split = ";", fixed = T)))
freq_subc <- table(unlist(p))
remove_subc <- names(which(table(unlist(p)) < 100))

### ### ###
### READ ALL PROTEIN INTERACTOMES
files <- dir("Interactomes/")
interactomes <- sapply(files, function(x) as.matrix(read.delim(paste("Interactomes/",x,sep =""),
                                                               header = T, sep = "\t")))
names(interactomes) <- unname(sapply(names(interactomes), function(x) unlist(strsplit(x, split = ".", fixed = T))[1]))

### List ALL genes from each INTERACTOME
interactomes_g <- lapply(interactomes, function(x) unique(c(x))) ## per interactome
interactomes_gALL <- unique(unlist(unname(interactomes_g)))

### ### ###
### Subcellular localization information: GO database (Dic 2018)
files <- dir("Subcellular_localization/")[grepl("Proteinlist", dir("Subcellular_localization/"))]

GO <- sapply(files, function(x) {
  p <- as.matrix(read.csv(paste("Subcellular_localization/",x,sep =""), header = FALSE))
  annot <- select(keys = p[,1], x = Homo.sapiens, keytype = "UNIPROT", columns = c("ENSEMBL","SYMBOL"))
  annot <- by(annot$ENSEMBL, annot$UNIPROT, as.character)

  annot <- lapply(as.list(annot), function(x) {
    if(length(x) > 1)
    {
      if(length(x[x %in% interactomes_gALL]) > 0)
        x <- x[x %in% interactomes_gALL]
      else
        x <- x
    }

    return(x)
  })
  p <- c(na.omit(unlist(unname(annot))))

  return(p)
})
names(GO) <- c("Extracellular Region", "Extracellular Vesicle")


### Annotating each gene of interactomes to subcellular location
interactomes_gSubc <- lapply(interactomes_g, function(x) {
  d <- data.frame(x, Subcellular.location = subcell[match(x, subcell$Gene),"Subcellular.location"])
  d$Subcellular.location[is.na(d$Subcellular.location)] <- ""
  colnames(d) <- c("ENSEMBL","Subcellular.location")
  return(d)})

### Creating a binary (1 indicates presence) matrix of protein subcellular location per interactome
interactomes_gSubc <- lapply(interactomes_gSubc, function(d) {
  n <- sapply(as.character(d$Subcellular.location), function(x)
    unlist(strsplit(x, split = ";", fixed = T)))
  names(n) <- d$ENSEMBL

  locations <- unique(unlist(n))

  # Cell Atlas
  n <- lapply(n, function(x) as.integer(locations %in% x))
  n <- do.call(rbind, n)
  colnames(n) <- locations

  # including GO
  n <- cbind(n, do.call(cbind, lapply(GO, function(x) as.integer(rownames(n) %in% x))))

  # including no located annotation
  none <- as.numeric(rowSums(n) == 0)
  n <- cbind(n, None = none)

  return(n)
})

#########################################################################################################
####### CALCULATING ODDS RATIO TO DEFINE BIASES PER INTERACTOME

### Wrapper function to calculate odds ratios
oddsRatioSubc <- function(n, universe, subcell){
  res <- sapply(seq_len(34), function(i) {
    if(i < 33)
      x <- table(interactome = universe %in% rownames(n),
                 location = universe %in% as.character(subcell$Gene)[grepl(
                   colnames(n)[i], as.character(subcell$Subcellular.location), fixed = T)])
    else
      x <- table(interactome = universe %in% rownames(n),
                 location = universe %in% GO[[i-32]])

    p.val <- fisher.test(x)$p.value
    odds <- fisher.test(x)$estimate

    return(c(p.value = p.val, odds = unname(odds)))
  })
  colnames(res) <- colnames(n)[seq_len(34)]

  return(t(res))
}

## CALCULATING ODDS RATIOS
prot_coding <- read.table("Genome/protein_coding_genome.tsv", header = T, sep = "\t")
prot_coding <- as.character(prot_coding$ensembl_gene_id_short)
universe <- intersect(prot_coding, levels(subcell$Gene))

odds <- lapply(interactomes_gSubc, function(y) oddsRatioSubc(y, universe, subcell))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## PLOTTING ODDS RATIOS
pdf("SubcellularLocation_biases.pdf", height = 12, width = 7)
par(mar = c(8,15,7,8))

ord <- c(unlist(unname(categ))[!unlist(unname(categ)) %in% remove_subc], names(GO))
ord_int <- c("HI-III", "BioPlex", "QUBIC", "Lit-BM-17", "CoFrac")

#### Heatmap of ALL interactomes
mx <- log2(do.call(cbind,lapply(odds, function(x) x[rev(ord), "odds"])))[, ord_int]
mx2 <- do.call(cbind,lapply(odds, function(x) x[rev(ord),"p.value"]))[, ord_int]
n <- mx
mx[mx2 > 0.05] <- NA
mx2[] <- ifelse(c(mx2) > 0.05, 1, 0)
n[is.infinite(n)] <- NA
image(t(mx2), axes = F, col = c(adjustcolor("white",0),adjustcolor("grey",0.25)))
image(t(pmax(pmin(mx, 1.0), -1.0)), axes = F, col = c(colorRampPalette(c("blue","white","red"))(32)), add = T,
      zlim = c(-1.0, 1.0))
axis(2, at = seq(0,1,length.out = 22), labels = rev(ord),
     las = 2, lwd = 0)
axis(3, at = seq(0,1,length.out = dim(mx)[2]), labels = colnames(mx),
     las = 2, lwd = 0)

grid(ncol(n),nrow(n),lty = 1)

# LUKE: adding color bar
colfunc <- colorRampPalette(c("blue","white","red"))
legend_image <- as.raster(matrix(colfunc(32), ncol=1))
plot(c(0,2),c(-1,1),type = 'n', axes = F, xlab = '', ylab = '', main = 'Log2 odds ratio')
text(x=1.5, y = seq(-1,1,l=3), labels = seq(-1,1,l=3))
rasterImage(legend_image, 0, 1, 1, -1)


dev.off()


##########################################################################################
##########################################################################################
