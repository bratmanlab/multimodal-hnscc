######################################
## TCGA Analysis w/ Expanded Cohort ##
######################################

# Set directories
scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"

# Load libraries
library(limma)

# # Download data from TCGA ----------
# ## Read manifest from all samples (floor of mouth (OSCC), BRCA, PRAD, LUAD, LUSC, COAD, PAAD)
# manifest <- read.table(paste0(data_dir, "/gdc_manifest_20190710_181822.txt"),
#                        sep="\t", header=T, stringsAsFactors = FALSE)
# 
# ### Subset samples w/ primary tumours only (0A1)
# manifest <- manifest[grepl("-01A-", manifest$filename),]
# 
# ### Identify number of files per tumour type
# manifest_tumourlabel <- unlist(strsplit(manifest$filename, "[.]")) # split by "."
# manifest_tumourlabel <- manifest_tumourlabel[seq(2,length(manifest_tumourlabel),8)] # select only string w/ tumour site
# table(manifest_tumourlabel)
# BRCA: 776, COAD: 294, OSCC: 55, LUAD: 459, LUSC: 365, PAAD: 184, PRAD: 484 #
# 
# ### Index by tumour type and randomly select 55 from each (and all from OSCC)
# TCGA_index <- list()
# tumourlabel <- names(table(manifest_tumourlabel))
# for(i in 1:length(tumourlabel)){
# sublabel <- unlist(strsplit(tumourlabel[i], "_"))[2]
# index <- grepl(sublabel, manifest_tumourlabel)
# index <- which(index == TRUE)
# set.seed(42)
# subindex <- sample(index, 55)
# TCGA_index[[i]] <- subindex
# names(TCGA_index)[i] <- sublabel
# }
# TCGA_index <- unlist(TCGA_index)
# 
# ### Subset manifest file and write
# manifest_sampled <- manifest[TCGA_index,]
# write.table(manifest_sampled,
# file=paste0(data_dir, "/gdc_manifest_20190710_181822_subsampled.txt"),
# sep="\t", quote = FALSE, row.names = FALSE)
# 
# ## Download samples ----------
# ### Performed on Windows CMD
# 
# ## Import files --------
# ### Load libraries
# library(dplyr)
# library(ggplot2)
# library(limma)
# library(reshape2)
# library(minfi)
# 
# library(data.table)
# ### Create functions
# ParseAndPrep <- function(Path){
#   A <- fread(Path, data.table = F)
#   B <- as.character(A[,1])
#   A <- data.matrix(A[,2])
#   rownames(A) <- B
#   return(A)
# }
# 
# Fileprep <- function(Path, pattern){
#   file1 <- list.files(Path, full.names=T)
#   file2 <- c()
#   for(directory in file1){
#     file2 <- c(file2, list.files(directory, pattern=pattern, full.names = TRUE)) #01A denotes primary tumour
#   }
#   return(file2)
# }
# 
# ### Read files and collate
# Combined450K_path <- Fileprep(paste0(data_dir), pattern="01A")
# Combined450K <- lapply(Combined450K_path, ParseAndPrep)
# Combined450K_df <- do.call(cbind, Combined450K)
# 
# #### for colnames
# df_colnames <- c()
# for(path in Combined450K_path){
# split <- unlist(strsplit(path, "/"))[8]
# df_colnames <- c(df_colnames, split)
# }
# 
# colnames(Combined450K_df) <- df_colnames

# Import PBL files from GSE67393 ----------
# NOTE: PBL files were downloaded on GEO portal and deposited into the
#       data/GSE67393 folder
# PBL_path <- list.files(paste0(data_dir, "/GSE67393"),
# full.names = T)
# PBLL_450K <- lapply(PBL_path, ParseAndPrep)

### Pull out number of probeIDs per file
#PBL450K_probeIDlength <- c()
#for(i in 1:length(PBL_450K)){
#PBL450K_probeIDlength <- c(PBL450K_probeIDlength, nrow(PBL_450K[[i]]))
#}

### Remove files with less than 485,577 probes
#index <- PBL450K_probeIDlength == 485577
#PBL_450K <- PBL_450K[index]

### Subset to 55 samples
#set.seed(42)
#PBL_450K <-  PBL_450K[sample(1:length(PBL_450K), 55)]

### Collate to df
#PBL_450K_df <- do.call(cbind, PBL_450K)
#colnames(PBL_450K_df) <- paste("PBL",1:ncol(PBL_450K_df),sep="_")
## Combine w/ PBL --------
#Combined450K_df <- cbind(Combined450K_df, PBL_450K_df)
#Combined450K_df <- Combined450K_df[complete.cases(Combined450K_df),]
#save(Combined450K_df, file=paste0(data_dir, "/Combined450K_df.RData"))

## Identify tumour-specific hyperDMCs -----
load(paste0(data_dir, "/Combined450K_df.RData"))
load(paste0(data_dir, "/All450K_Probes_Info.RData"))

## Set up for differential analysis
Phenotype <- c()
for(label in colnames(Combined450K_df)){
  if(grepl("HNSC", label)) Phenotype <- c(Phenotype, "OSCC")
  if(grepl("BRCA", label)) Phenotype <- c(Phenotype, "BRCA")
  if(grepl("COAD", label)) Phenotype <- c(Phenotype, "COAD")
  if(grepl("LUSC", label)) Phenotype <- c(Phenotype, "LUSC")
  if(grepl("LUAD", label)) Phenotype <- c(Phenotype, "LUAD")
  if(grepl("PRAD", label)) Phenotype <- c(Phenotype, "PRAD")
  if(grepl("PAAD", label)) Phenotype <- c(Phenotype, "PAAD")
  if(grepl("PBL", label)) Phenotype <- c(Phenotype, "PBL")
}

### Perform for each condition
dir.create(paste0(out_dir, "/tcgaHyperDmr"))
setwd(paste0(out_dir, "/tcgaHyperDmr")) # output of DMR analysis

Groups <- levels(factor(Phenotype))
design <- model.matrix(~0 + Phenotype)
colnames(design) <- levels(factor(Phenotype))

source("../../code/MVP_DMR_functions.R")

#### Create funciton
customdiff <- function(Results, cm){
  for(i in 1:numClasses) {
    Tab  <- find.mvp(Combined450K_df,type = "beta", design = design, contrast.matrix = cm, classes = Phenotype, TSBH = F, alpha.TSBH = NULL, coef = i)
    Tab <- Tab$tt
    Results[[i]] <- Tab

    message(paste0(i," done"))

  }

  names(Results) <- make.names(colnames(cm))

  ### Export results

  for (i in 1:length(Results)) {
    CancerClass=unlist(strsplit(names(Results)[i], split="[.]"))[1]
    CancerMedianColumnIndex=which(colnames(Results[[i]])==paste0("median.", CancerClass))

    BloodClass=unlist(strsplit(names(Results)[i], split="[.]"))[4]
    BloodMedianColumnIndex=which(colnames(Results[[i]])==paste0("median.", BloodClass))

    CancerClass_DF=Results[[i]][,c(1,6,CancerMedianColumnIndex,BloodMedianColumnIndex)]


    CancerClass_DF=dplyr::filter(CancerClass_DF, adj.P.Val < 0.01)

    Hyper_Indices=which(CancerClass_DF[,3]-CancerClass_DF[,4] >= 0.25)
    HyperCancer_DF=CancerClass_DF[Hyper_Indices,]

    ## Adding ID to coords mapping here

    HyperProbe_Indices=match(x=as.character(HyperCancer_DF[,1]), table=as.character(All450K_Probes_Info[,4]))
    HyperProbe_Indices=HyperProbe_Indices[!is.na(HyperProbe_Indices)]
    write.table(All450K_Probes_Info[HyperProbe_Indices,], file=paste0(CancerClass,"_vs_",BloodClass,"_HyperDMCs.bed"), col.names=F, row.names=F, sep="\t", quote=F, append=F)

  }
}

#### OSCC-specific
Results_OSCCspecific <- list()
cm <- makeContrasts(OSCC-PBL, OSCC-BRCA, OSCC-COAD, OSCC-LUAD, OSCC-LUSC, OSCC-PAAD, OSCC-PRAD, levels = design)
numClasses=7
customdiff(Results=Results_OSCCspecific, cm=cm)

#### BRCA-specific
Results_Sitespecific <- list()
cm <- makeContrasts(BRCA-PBL, BRCA-OSCC, BRCA-COAD, BRCA-LUAD, BRCA-LUSC, BRCA-PAAD, BRCA-PRAD, levels = design)
numClasses=7
customdiff(Results=Results_Sitespecific, cm=cm)

#### COAD-specific
Results_Sitespecific <- list()
cm <- makeContrasts(COAD-PBL, COAD-OSCC, COAD-BRCA, COAD-LUAD, COAD-LUSC, COAD-PAAD, COAD-PRAD, levels = design)
numClasses=7
customdiff(Results=Results_Sitespecific, cm=cm)

#### LUAD-specific
Results_Sitespecific <- list()
cm <- makeContrasts(LUAD-PBL, LUAD-OSCC, LUAD-BRCA, LUAD-COAD, LUAD-LUSC, LUAD-PAAD, LUAD-PRAD, levels = design)
numClasses=7
customdiff(Results=Results_Sitespecific, cm=cm)

#### LUSC-specific
Results_Sitespecific <- list()
cm <- makeContrasts(LUSC-PBL, LUSC-OSCC, LUSC-BRCA, LUSC-LUAD, LUSC-COAD, LUSC-PAAD, LUSC-PRAD, levels = design)
numClasses=7
customdiff(Results=Results_Sitespecific, cm=cm)

#### PAAD-specific
Results_Sitespecific <- list()
cm <- makeContrasts(PAAD-PBL, PAAD-OSCC, PAAD-BRCA, PAAD-LUAD, PAAD-LUSC, PAAD-COAD, PAAD-PRAD, levels = design)
numClasses=7
customdiff(Results=Results_Sitespecific, cm=cm)

#### PRAD-specific
Results_Sitespecific <- list()
cm <- makeContrasts(PRAD-PBL, PRAD-OSCC, PRAD-BRCA, PRAD-LUAD, PRAD-LUSC, PRAD-PAAD, PRAD-COAD, levels = design)
numClasses=7
customdiff(Results=Results_Sitespecific, cm=cm)

# Load HyperDMCs from h4h cluster analysis -------
### Select site-specific IDS
tumours <- c("BRCA","COAD","LUAD","LUSC","OSCC","PAAD","PRAD")
setwd("../../code")
filedir <- paste0(out_dir, "/tcgaHyperDmr")
sitespecific_vs_all_HyperDMRs <- list()
sitespecific_vs_all_probeIDS <- list()

for(i in 1:length(tumours)){
  sitefiles <- list.files(filedir, pattern=paste0(tumours[i],"_vs"), full.names = T)
  probeID_list <- list()
  windows_list <- list()
  for(n in 1:length(sitefiles)){
    label <- strsplit(strsplit(list.files(filedir, pattern=paste0(tumours[i],"_vs"), full.names = T)[1], split="/")[[1]][7], split="[.]")[[1]][1]
    tmp <- read.table(sitefiles[n], sep="\t", col.names = c("chr","start","end","probeID"))
    probeIDs <- as.character(tmp$probeID)
    windows <- paste(tmp$chr, tmp$start, tmp$end, sep=".")
    probeID_list[[n]] <- probeIDs
    names(probeID_list)[n] <- label
    windows_list[[n]] <- windows
    names(windows_list)[n] <- label
  }
  sitespecific_vs_all_probeIDS[[i]] <- names(table(unlist(probeID_list))[table(unlist(probeID_list)) == 7])
  names(sitespecific_vs_all_probeIDS)[i] <- tumours[i]
  sitespecific_vs_all_HyperDMRs[[i]] <- names(table(unlist(windows_list))[table(unlist(windows_list)) == 7])
  names(sitespecific_vs_all_HyperDMRs)[i] <- tumours[i]
}


#### Evaluate uniqueness
table(table(unlist(sitespecific_vs_all_probeIDS))) # should all be 1 if unique
table(table(unlist(sitespecific_vs_all_HyperDMRs))) # some overlap, subset

#### Combined probes
sitespecific_vs_all_probeIDS_collated <- as.character(unlist(sitespecific_vs_all_probeIDS))
sitespecific_vs_all_HyperDMRs_collated <- as.character(names(table(unlist(sitespecific_vs_all_HyperDMRs))[table(unlist(sitespecific_vs_all_HyperDMRs)) == 1]))

save(sitespecific_vs_all_probeIDS, file=paste0(out_dir, "/sitespecific_vs_all_probeIDs.RData"))
