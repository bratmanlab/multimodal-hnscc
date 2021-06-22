# Identification of informative ctDNA-derived methylation #
###########################################################

# Set directories
scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"
TCGA_dir <- paste0(out_dir, "../TCGAbiolinks")
TCGAdata_dir <- paste0(data_dir, "/TCGA")

# Load libraries
library(TCGAbiolinks)
library(DESeq2)
library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(AnnotationHub)
library(ggbio)
library(Homo.sapiens)
library(glmnet)
data(genesymbol, package = "biovizBase")
source(paste0(scripts, "/ggplot_themeJB.R"))

## Custom functions...
load(paste0(data_dir, "/IDtoGrange.RData"))

## TCGAbiolinks objects...

### NOTE: each object that is loaded in needs to be renamed to avoid removal
### of previous object upon loading of subsequent ones.

### Infinium 450k microarray data of primary tumors
#### NOTE: data was split into two seperate objects due to size constraints
load(paste0(TCGAdata_dir, "/met_hnc_primary_p1.RData"))
met.hnc.primary.p1 <- data
rm(data)

load(paste0(TCGAdata_dir, "/met_hnc_primary_alt.RData"))
met.hnc.primary.p2 <- data
rm(data)

met.hnc.primary.assay <- cbind(assay(met.hnc.primary.p1),
                               assay(met.hnc.primary.p2))

### RNA sequencing data of primary tumors
load(paste0(TCGAdata_dir, "/iso_hnc_primary_p1.RData"))
iso.hnc.primary <- data
rm(data)


## Clinical outcome data provided at 
## https://www.cell.com/cell/fulltext/S0092-8674(18)30229-0
clin.hnc <- fread(file = paste0(TCGAdata_dir, 
                                "/clinical_outcome_hnsc.txt"),
                  sep = "\t", header = TRUE)

## DESeq2 results object...
load(paste0(out_dir, "/DESeq2_cfDNA_DMRs.RData"))

# 2. Download methylation data of adjacent normal tissue by TCGAbiolinks ------

## Create query object
# query.met.hnc.normal <- GDCquery(project = "TCGA-HNSC",
#                                  legacy = TRUE,
#                                  data.category = "DNA methylation",
#                                  platform = "Illumina Human Methylation 450",
#                                  sample.type = c("Solid Tissue Normal"))

## Download and prepare data from query
# GDCdownload(query.met.hnc.normal)
# GDCprepare(query.met.hnc.normal,
#            save = TRUE,
#            summarizedExperiment = TRUE,
#            save.filename = paste0(TCGA_dir, "/met_hnc_normal.RData"))
load(paste0(TCGAdata_dir, "/met_hnc_normal.RData"))
met.hnc.normal <- data
rm(data)

# 3. Select probeIDs significantly elevated in primary tumor tissue -----------

## Extract assay data from loaded RangedSummarizedExperiment objects
met.hnc.primary.assay <- cbind(assay(met.hnc.primary.p1), 
                               assay(met.hnc.primary.p2))
met.hnc.normal.assay <- assay(met.hnc.normal)

## Remove NA values in > 50% of cohort from either 
## primary or normal tissue data
### Determine NAs for each matrix
primary.na_index <- apply(met.hnc.primary.assay, 1, FUN = function(x){
  sum(is.na(x)) > (ncol(met.hnc.primary.assay) / 2)
})

normal.na_index <- apply(met.hnc.normal.assay, 1, FUN = function(x){
  sum(is.na(x)) > (ncol(met.hnc.normal.assay) / 2)
})

### Collate NAs and subset matrices
combined.na_index <- primary.na_index | normal.na_index

met.hnc.primary.assay.nonNA <- met.hnc.primary.assay[!combined.na_index,]
met.hnc.normal.assay.nonNA <- met.hnc.normal.assay[!combined.na_index,]

### Summarize matrices by hyperDMRs
#### Subset probeIDs to regions overlapping with plasma DMRs

##### Create plasma DMR GRange object
hyperDMRs_ROIs <- res_list[["res_padj_hyperDMR"]]
hyperDMRs_ROIs_gr <- IDtoGRange(rownames(hyperDMRs_ROIs))
load(paste0(data_dir, "/probeID_to_300bp_wCpG.RData"))

probeID_to_300bp_gr <- IDtoGRange(as.character(probeID_to_300bp_wCpG$ID))
probeID_to_300bp_gr$probeID <- as.character(probeID_to_300bp_wCpG$probeID)
probeID_to_300bp_gr$ID <- as.character(probeID_to_300bp_wCpG$ID)
probeID_to_300bp_gr <- probeID_to_300bp_gr[probeID_to_300bp_gr$probeID %in%
                                             rownames(met.hnc.primary.assay.nonNA)]

tmp_ind <- 
  match(probeID_to_300bp_gr$probeID, rownames(met.hnc.primary.assay.nonNA))
tmp_ind <- tmp_ind[!is.na(tmp_ind)]


##### Create hm450k assay GRange object and subset
hm450k_assay_gr <- rowRanges(met.hnc.primary.p2)

probeIDs_plasmaDMR_overlap <- findOverlaps(hyperDMRs_ROIs_gr,
                                           hm450k_assay_gr)
probeIDs_plasmaDMR_gr <- 
  hyperDMRs_ROIs_gr[queryHits(probeIDs_plasmaDMR_overlap)]

probeIDs_plasmaDMR_gr$probeIDs <- 
  hm450k_assay_gr[subjectHits(probeIDs_plasmaDMR_overlap)]$probeID

probeIDs_plasmaDMR_gr <- 
  probeIDs_plasmaDMR_gr[probeIDs_plasmaDMR_gr$probeIDs %in%
                          rownames(met.hnc.primary.assay.nonNA)]


#### Collate hm450k assay based on hyperDMRs
met.hnc.DMR.assay <- as.data.frame(matrix(ncol = 
                                            ncol(met.hnc.primary.assay.nonNA),
                                          nrow = 0))
colnames(met.hnc.DMR.assay) <- colnames(met.hnc.primary.assay.nonNA)
hyperDMRs_hm450k_gr <- unique(probeIDs_plasmaDMR_gr[,-1])
for(i in 1:length(hyperDMRs_hm450k_gr)){
  tmp_gr <- hyperDMRs_hm450k_gr[i]
  probeIDs <- subsetByOverlaps(probeIDs_plasmaDMR_gr,
                               tmp_gr)$probeIDs
  tmp_df <- 
    subset(met.hnc.primary.assay.nonNA,
           rownames(met.hnc.primary.assay.nonNA) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.hnc.DMR.assay <- rbind(met.hnc.DMR.assay, tmp_colMeans)
  rownames(met.hnc.DMR.assay)[i] <- paste(seqnames(hyperDMRs_hm450k_gr[i]),
                                          start(hyperDMRs_hm450k_gr[i]),
                                          end(hyperDMRs_hm450k_gr[i]),
                                          sep = ".")
}

met.norm.DMR.assay <- as.data.frame(matrix(ncol = 
                                             ncol(met.hnc.normal.assay.nonNA),
                                           nrow = 0))
colnames(met.norm.DMR.assay) <- colnames(met.hnc.normal.assay.nonNA)
for(i in 1:length(hyperDMRs_hm450k_gr)){
  tmp_gr <- hyperDMRs_hm450k_gr[i]
  probeIDs <- subsetByOverlaps(probeIDs_plasmaDMR_gr,
                               tmp_gr)$probeIDs
  tmp_df <- 
    subset(met.hnc.normal.assay.nonNA,
           rownames(met.hnc.normal.assay.nonNA) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.norm.DMR.assay <- rbind(met.norm.DMR.assay, tmp_colMeans)
  rownames(met.norm.DMR.assay)[i] <- paste(seqnames(hyperDMRs_hm450k_gr[i]),
                                           start(hyperDMRs_hm450k_gr[i]),
                                           end(hyperDMRs_hm450k_gr[i]),
                                           sep = ".")
}

# For each probeID, perform Wilcoxon test
wilcox.pval <- c()
pb <- txtProgressBar(min = 0, max = length(met.hnc.DMR.assay),
                     initial = 0)
for(i in 1:nrow(met.hnc.DMR.assay)){
  results <- wilcox.test(x = as.numeric(met.hnc.DMR.assay[i,]),
                         y = as.numeric(met.norm.DMR.assay[i,]),
                         alternative = "two.sided")
  wilcox.pval <- c(wilcox.pval, results$p.value)
  setTxtProgressBar(pb, i)
}

wilcox.adj.pval <- p.adjust(wilcox.pval, method = "BH")

## Select regions w/ increased mean beta values in primary tumour
## versus adjacent normals (delta beta value > 0.2)
hyperDMC_index <- log2(rowMeans(met.hnc.DMR.assay, na.rm = TRUE) /
                         rowMeans(met.norm.DMR.assay, na.rm = TRUE)) > 1

hyperDMC.sig_index <- wilcox.adj.pval <= 0.05 & hyperDMC_index

met.hnc.TCGA.hyperDMR.assay <- 
  met.hnc.DMR.assay[hyperDMC.sig_index,]

colnames(met.hnc.TCGA.hyperDMR.assay) <- colnames(met.hnc.primary.assay.nonNA)

## Create GRange objects w/ IDs
hyperDMRs_hm450k_gr$ID <- paste(seqnames(hyperDMRs_hm450k_gr),
                                start(hyperDMRs_hm450k_gr),
                                end(hyperDMRs_hm450k_gr),
                                sep = ".")

hyperprobeIDs_plasmaDMR_gr <- 
  hyperDMRs_hm450k_gr[hyperDMRs_hm450k_gr$ID %in%
                        rownames(met.hnc.TCGA.hyperDMR.assay)]

# 4. For each probeID, calculate correlation with gene expression -------------

## Create GRange and assay of RNA expression assay
## Data from https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
## under Genome Annotation Files for Legacy TCGA Data and converted to a 
## GenomicRanges object
load(paste0(TCGAdata_dir, "/TCGA_hg19_June2011_gr.RData"))
RNAseq_gr <- TCGA_hg19_June2011_gr
RNAseq_assay <- iso.hnc.primary
RNAseq_assay <- RNAseq_assay[rownames(RNAseq_assay) %in% 
                               RNAseq_gr$ucsc_id,]
colnames(RNAseq_assay) <- gsub("normalized_count_", "",
                               colnames(RNAseq_assay))

RNAseq_ord <- match(substr(colnames(met.hnc.TCGA.hyperDMR.assay), 1, 12),
                    substr(colnames(RNAseq_assay), 1, 12))

RNAseq_assay <- RNAseq_assay[ , RNAseq_ord]

## For each probeID, evaluate the expression of genes within 2kb
met_vs_exp_cor <- list()
for(i in 1:nrow(met.hnc.TCGA.hyperDMR.assay)){
  # Obtain beta-values of probeID
  matrix_probeID.bval <- met.hnc.TCGA.hyperDMR.assay[i,]
  plasmaDMR_hyperprobeIDs_gr <- IDtoGRange(rownames(matrix_probeID.bval))
  plasmaDMR_hyperprobeIDs_gr$probeID <- 
    paste(seqnames(plasmaDMR_hyperprobeIDs_gr),
          start(plasmaDMR_hyperprobeIDs_gr),
          end(plasmaDMR_hyperprobeIDs_gr),
          sep = ".")
  
  # Find nearest genes by RNA expression data
  probeID_gene.overlap <- 
    suppressWarnings(findOverlaps(plasmaDMR_hyperprobeIDs_gr,
                                  RNAseq_gr,
                                  maxgap = 2000L))
  if(length(probeID_gene.overlap) == 0){
    rho <- NA
    pval <- NA
    tmp_c <- c(rho, pval)
    names(tmp_c) <- c("rho", "p.value")
    tmp_list <- list(tmp_c)
    names(tmp_list) <- "NoOverlap"
    met_vs_exp_cor[[i]] <- tmp_list
    names(met_vs_exp_cor)[i] <- plasmaDMR_hyperprobeIDs_gr$probeID[i]
  } else {
    exp.gene <- RNAseq_gr$ucsc_id[subjectHits(probeID_gene.overlap)]
    matrix_gene.rsem <- subset(RNAseq_assay, 
                               rownames(RNAseq_assay) %in% exp.gene)
    matrix_gene.rsem <- log2(matrix_gene.rsem + 1)
    
    # For each gene, calculate the correlation w/ probeID beta values
    tmp_list <- list()
    for(n in 1:nrow(matrix_gene.rsem)){
      tmp_cor <- cor.test(x = as.numeric(matrix_probeID.bval),
                          y = as.numeric(matrix_gene.rsem[n,]),
                          method = "spearman")
      rho <- c(tmp_cor$estimate)
      pval <- c(tmp_cor$p.value)
      tmp_c <- c(rho, pval)
      names(tmp_c) <- c("rho", "p.value")
      tmp_list[[n]] <- tmp_c
      names(tmp_list)[n] <- rownames(matrix_gene.rsem)[n]
    }
    met_vs_exp_cor[[i]] <- tmp_list
    names(met_vs_exp_cor)[i] <- plasmaDMR_hyperprobeIDs_gr$probeID
  }
}

## Calculate adjusted p-values of correlations
met_vs_exp_pval <- c()
met_vs_exp_rho <- c()
met_vs_exp_probeID <- c()
met_vs_exp_gene <- c()
for(i in 1:length(met_vs_exp_cor)){
  tmp_unlist <- unlist(met_vs_exp_cor[[i]])
  
  grep_pval_index <- grepl("p.value", 
                           names(tmp_unlist))
  met_vs_exp_pval <- c(met_vs_exp_pval, tmp_unlist[grep_pval_index])
  
  grep_rho_index <- grepl("rho",
                          names(tmp_unlist))
  met_vs_exp_rho <- c(met_vs_exp_rho, tmp_unlist[grep_rho_index])
  
  tmp_probeID <- rep(names(met_vs_exp_cor)[i],
                     length(tmp_unlist)/2)
  met_vs_exp_probeID <- c(met_vs_exp_probeID, tmp_probeID)
  
  tmp_gene <- gsub(".rho", "", names(tmp_unlist[grep_rho_index]))
  met_vs_exp_gene <- c(met_vs_exp_gene, tmp_gene)
  
}

## Collate to dataframe
met_vs_exp_df <- data.frame(probeID = met_vs_exp_probeID,
                            ucsc_ID = met_vs_exp_gene,
                            rho = met_vs_exp_rho,
                            pval = met_vs_exp_pval)

## Remove probeIDs w/ no overlap
met_vs_exp_overlap_df <- subset(met_vs_exp_df, 
                                met_vs_exp_df$ucsc_ID != "NoOverlap")

gene_ind <- match(met_vs_exp_overlap_df$ucsc_ID, RNAseq_gr$ucsc_id)
gene_ind <- gene_ind[!is.na(gene_ind)]

met_vs_exp_overlap_df$gene <- 
  RNAseq_gr[gene_ind]$gene

met_vs_exp_overlap_df$refseq <- as.character(RNAseq_gr[gene_ind]$gene)

met_vs_exp_overlap_df$BH <- p.adjust(met_vs_exp_overlap_df$pval, method = "BH")
met_vs_exp_overlap_df$BF <- p.adjust(met_vs_exp_overlap_df$pval, 
                                     method = "bonferroni")

sig_met_vs_exp <- subset(met_vs_exp_overlap_df, 
                         (abs(met_vs_exp_overlap_df$rho) >= 0.3 &
                            met_vs_exp_overlap_df$BF < 0.05))

# 5. Perform Cox regression for all hyperDMC probeIDs -------------------------

## Subset hm450k assay to regions within plasma DMRs and remove NA values

plasmaDMR_hm450k_assay <- met.hnc.TCGA.hyperDMR.assay

## Subset and order TCGA-HNSC clinical data
clin.hnc.sub <- subset(clin.hnc, 
                       clin.hnc$bcr_patient_barcode %in% 
                         substr(colnames(met.hnc.TCGA.hyperDMR.assay), 1, 12))

clin.hnc.sub <- 
  clin.hnc.sub[match(substr(colnames(met.hnc.TCGA.hyperDMR.assay), 1, 12),
                     clin.hnc.sub$bcr_patient_barcode),]

## Summarize to 5-year survival
fiveyr_DSStime <- c()
fiveyr_DSS <- c()
fiveyr_OStime <- c()
fiveyr_OS <- c()
for(i in 1:length(clin.hnc.sub$DSS.time)){
  tmpx <- clin.hnc.sub$DSS.time[i]
  if(is.na(tmpx)){
    tmp_time <- NA
    tmp_stat <- NA
  } else {
    if(tmpx <= 1825){
      tmp_time <- tmpx
      tmp_stat <- clin.hnc.sub$DSS[i]
    } else {
      tmp_time <- 1825
      tmp_stat <- 0
    } 
  }
  fiveyr_DSStime <- c(fiveyr_DSStime, tmp_time)
  fiveyr_DSS <- c(fiveyr_DSS, tmp_stat)
}

for(i in 1:length(clin.hnc.sub$OS.time)){
  tmpx <- clin.hnc.sub$OS.time[i]
  if(is.na(tmpx)){
    tmp_time <- NA
    tmp_stat <- NA
  } else {
    if(tmpx <= 1825){
      tmp_time <- tmpx
      tmp_stat <- clin.hnc.sub$OS[i]
    } else {
      tmp_time <- 1825
      tmp_stat <- 0
    } 
  }
  fiveyr_OStime <- c(fiveyr_OStime, tmp_time)
  fiveyr_OS <- c(fiveyr_OS, tmp_stat)
}

## Create survival object
surv_object <- Surv(time = fiveyr_DSStime, 
                    event = fiveyr_DSS)
surv_df <- as.data.frame(t(plasmaDMR_hm450k_assay))
surv_df$age <- clin.hnc.sub$age_at_initial_pathologic_diagnosis
surv_df$gender <- clin.hnc.sub$gender
surv_df$clinical_stage <- clin.hnc.sub$clinical_stage

surv_OS_object <- Surv(time = fiveyr_OStime,
                       event = fiveyr_OS)

## For each probeID, perform cox analysis and extract the respective
## p-value, hazard ratio, and z-statistic

cox_pval <- c()
cox_HR <- c()
cox_z <- c()

for(i in 1:ncol(surv_df)){
  probeID <- colnames(surv_df)[i]
  
  # create long formula...
  tmp_form <- paste("surv_object ~ ",probeID, 
                    " + age",
                    " + gender",
                    " + clinical_stage",
                    sep = "")
  
  # calculate Cox and extract HR and pvalue
  fit.coxph <- coxph(formula(tmp_form), data = surv_df)
  HR_tmp <- summary(fit.coxph)$coefficients[1,2]
  z_tmp <- summary(fit.coxph)$coefficients[1,4]
  pval_tmp <- summary(fit.coxph)$coefficients[1,5]
  
  cox_HR <- c(cox_HR, HR_tmp)
  cox_pval <- c(cox_pval, pval_tmp)
  cox_z <- c(cox_z, z_tmp)
  
  names(cox_HR)[i] <- probeID
  names(cox_pval)[i] <- probeID
  names(cox_z)[i] <- probeID
}

# 6. Select CpGs w/ significant Cox scores and expression ---------------------

## Select CpGs w/ significant Cox scores
sig_cox_probeIDs <- cox_HR[cox_pval < 0.05]
length(sig_cox_probeIDs) # 33 / 483 DMRs are significant

## Subset significant CpGs w/ significant expression of nearby gene
sig_cox_met_vs_exp <-
  sig_met_vs_exp[sig_met_vs_exp$probeID %in% names(sig_cox_probeIDs),]
nrow(sig_cox_met_vs_exp) 
# 5 / 36 significant probeIDs are also associated with gene expression

# 7. Collate ALL probeIDs within significant regions---------------------------

## Perform KM analysis for DMRs
theme_JBborder <- function(base_size = 9, base_family = "sans", ...) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family) + 
      theme(plot.title = element_text(face = "plain",
                                      size = 9, hjust = 0),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            panel.border = element_rect(colour = "black", fill = NA,
                                        size = 0.13),
            plot.background = element_rect(color = NA),
            axis.title = element_text(size = 9, face = "plain"),
            axis.title.y = element_text(angle=90, vjust = 2,
                                        margin = margin(0,5,0,0)),
            axis.title.x = element_text(vjust = 1,
                                        margin = margin(5,0,0,0)),
            axis.text = element_text(size = 9),
            axis.line = element_line(colour = "black", size = 0.13),
            axis.ticks = element_line(size = 0.13),
            axis.ticks.length = unit(0.15, "cm"),
            panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = margin(0,0,0,0),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face = "plain", size = rel(1)),
            plot.margin = unit(c(1,1,1,1), "mm"),
            strip.background = element_rect(colour = NA, fill = NA),
            strip.text = element_text(face = "plain", vjust = 0.5, hjust = 0.5),
            ...
      ))
}

plasmaDMR_survdf <- surv_df[,colnames(surv_df) %in% names(sig_cox_probeIDs)]

### All regions combined
sig_met_df <- rowMeans(plasmaDMR_survdf, na.rm = TRUE)
sig_met_KM_df <- 
  data.frame(status = ifelse(sig_met_df > 
                               quantile(sig_met_df, 0.5, na.rm = TRUE),
                             "Abv","Blw"))
sig_met_KM_df$status <- factor(sig_met_KM_df$status, 
                               levels = c("Blw", "Abv"))
sig_met_KM_fit <- survfit(surv_object ~ status, data = sig_met_KM_df)
sig_met_OS_KM_fit <- survfit(surv_OS_object ~ status, data = sig_met_KM_df)

sig_met_KMplot <- ggsurvplot(sig_met_KM_fit, pval = TRUE,
                             pval.method = TRUE,
                             ggtheme = theme_JBwBorder(),
                             size = 0.26,
                             palette = palette_JB[c(2,1)], 
                             risk.table.y.text = FALSE,
                             pval.size = (5/14) * 9,
                             legend.title = "Methylation",
                             legend.labs = c("Blw med", "Abv med"),
                             risk.table.fontsize = 4, conf.int = FALSE,
                             xlab = "Time (months)", censor.shape = 124, 
                             censor.size = 2,
                             ylab = "Disease-specific survival",
                             xscale = "d_m", break.time.by=365.25)

TCGA_coxph <- coxph(surv_OS_object ~ status, data = sig_met_KM_df)
cox_TCGA_HR <- summary(TCGA_coxph)$coefficients[2] # HR = 1.66

sig_met_OS_KMplot <- ggsurvplot(sig_met_OS_KM_fit, pval = TRUE,
                                pval.method = TRUE,
                                ggtheme = theme_JBwBorder(),
                                size = 0.26,
                                palette = palette_JB[c(2,1)], 
                                risk.table.y.text = FALSE,
                                pval.size = (5/14) * 9,
                                legend.title = "Methylation",
                                legend.labs = c("Blw med", "Abv med"),
                                risk.table.fontsize = 4, conf.int = FALSE,
                                xlab = "Time (months)", censor.shape = 124, 
                                censor.size = 2,
                                ylab = "Overall survival",
                                xscale = "d_m", break.time.by=365.25)
sig_met_OS_wHR_KMplot <- 
  sig_met_OS_KMplot$plot + 
  annotate("text", x = 190, y = 0.14,
           label = paste0("HR = ", round(cox_TCGA_HR, 2)),
           size = (5/14) * 9)

sig_met_OS_wHR_KMplot

# 9. Prognostic  utility of hyperDMRs in cfDNA cohort -------------------------
## Load matrix w/ RPKM values from cfMeDIP-seq libraries of fragments betwen

## Remove non-essential files to clear up memory
rm(met.hnc.normal.assay)
rm(met.hnc.normal)
rm(met.hnc.primary.assay)
rm(met.hnc.primary.p1)
rm(met.hnc.primary.p2)

## 100 - 150 bp
load(paste0(data_dir, "/cfDNA_Frag100_150_RPKMs.RData"))

## Subset to hyperDMR regions
hyperDMR_IDs <- rownames(res_list$res_padj_hyperDMR)

cfDNA_Frag100_150_RPKMs <- 
  subset(cfDNA_Frag100_150_RPKMs,
         rownames(cfDNA_Frag100_150_RPKMs) %in% hyperDMR_IDs)

## Select CAPP-Seq positive HNSCC samples (baseline [BL]) only
load(paste0(out_dir, "/mclust_wnorm_df.RData"))
ctDNApos_IDs <- substr(rownames(mclust_wnorm_df[mclust_wnorm_df$ctDNA > 0,]),
                       3, 6)
ctDNApos_IDs <- ctDNApos_IDs[-20]

sub_BL_df <- 
  cfDNA_Frag100_150_RPKMs[, grepl(paste0(ctDNApos_IDs, collapse = "|"),
                                  colnames(cfDNA_Frag100_150_RPKMs))]

sub_Norm_df <- 
  cfDNA_Frag100_150_RPKMs[,grepl("Norm", colnames(cfDNA_Frag100_150_RPKMs))]

# For each BL sample, measure fraction of RPKMS from the 5 previously
# identified prognostic region by TCGA analysis
sub_BL_sig_list <- list()
for(i in 1:ncol(sub_BL_df)){
  tmp_df <- sub_BL_df[,i]
  tmp_sub_df <- tmp_df[tmp_df >= 0]
  sub_BL_sig_list[[i]] <- tmp_sub_df
  names(sub_BL_sig_list)[i] <- colnames(sub_BL_df)[i]
}

test <- sapply(sub_BL_sig_list, FUN = function(x){
  sum(x[names(x) %in% sig_cox_met_vs_exp$probeID] / sum(x))
})

names(test) <- substr(names(test), 1, 4)

test <- test[order(names(test))]

## Evaluate relationship between ctDNA abundance and CMS
load(paste0(out_dir, "/mclust_wnorm_df.RData"))
hclust <- mclust_wnorm_df[substr(rownames(mclust_wnorm_df), 3, 6) %in%
                            names(test),]
hclust$cms <- test
hclust$`Combined methylation score` <- ifelse(hclust$cms > quantile(hclust$cms, 0.5),
                    "Abv med", "Blw med")

ttest.hclust <- t.test(meanRPKM ~ `Combined methylation score`, data = hclust)
pval.ttest.hclust <- ttest.hclust$p.value

cappTtest.hclust <- t.test(ctDNA ~ `Combined methylation score`, data = hclust)
pval.cappTtest.hclust <- cappTtest.hclust$p.value

set.seed(775)
hclust_plot <- hclust %>%
  ggplot(aes(x = `Combined methylation score`, y = meanRPKM)) +
  geom_jitter(position = position_jitter(0.4),
              shape = 16, size = 1.5, col = rgb(0,0,0,1/2)) +
  geom_boxplot(outlier.shape = NA,
               show.legend = FALSE, col = "black", fill = NA,
               lwd = 0.13, fatten = 6) + 
  annotate("text", x = 1.5, y = 0.38, 
           label = paste("p = ", round(pval.ttest.hclust, 2))) + 
  theme_JB() + 
  labs(x = "Combined methylation score", y = "Methylation (mean RPKM)") + 
  scale_x_discrete(labels = c("Above median", "Below median"))

set.seed(75)
hclust_capp_plot <- hclust %>%
  ggplot(aes(x = `Combined methylation score`, y = ctDNA)) +
  geom_jitter(position = position_jitter(0.4),
              shape = 16, size = 1.5, col = rgb(0,0,0,1/2)) +
  geom_boxplot(outlier.shape = NA,
               show.legend = FALSE, col = "black", fill = NA,
               lwd = 0.13, fatten = 6) + 
  annotate("text", x = 1.5, y = 5, 
           label = paste("p = ", round(pval.cappTtest.hclust, 2))) + 
  theme_JB() + 
  labs(x = "Combined methylation score", y = "Methylation (mean RPKM)") + 
  scale_x_discrete(labels = c("Above median", "Below median"))

## Create dataframe for kaplan-meier analysis 
KM_cfDNA_df <- data.frame(status = ifelse(test > quantile(test, 0.5),
                                          "Abv med", "Blw med"))
KM_cfDNA_df$status <- factor(KM_cfDNA_df$status,
                             levels = c("Blw med", "Abv med"))

KM_DSS_time <- c(273, 1831, 1669, 280, 1658, 329, 374, 588, 662, 1537,
                 1260, 378, 1397, 358, 1322, 432, 1315, 403, 1247)
KM_DSS <- c(1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0)

KM_surv_obj <- Surv(time = KM_DSS_time, event = KM_DSS)

surv_fit <- survfit(KM_surv_obj ~ status, data = KM_cfDNA_df)

cfDNA_coxph <- coxph(KM_surv_obj ~ status, data = KM_cfDNA_df)
cox_cfDNA_HR <- summary(cfDNA_coxph)$coefficients[2] # HR = 3.78

cfDNA_KM_plot <- ggsurvplot(surv_fit, pval = TRUE,
                            pval.method = TRUE,
                            test.for.trend = FALSE,
                            ggtheme = theme_JBwBorder(),
                            size = 0.26,
                            palette = palette_JB[c(2,1)], 
                            risk.table.y.text = FALSE,
                            pval.size = (5/14) * 9,
                            legend.title = "ctDNA Methylation:",
                            legend.labs = c("Blw med", "Abv med"),
                            risk.table.fontsize = 4, conf.int = FALSE,
                            xlab = "Time (months)", censor.shape = 124, 
                            censor.size = 2,
                            ylab = "Overall survival",
                            xscale = "d_m", break.time.by=365.25)

cfDNA_KM_wHR_plot <- 
  cfDNA_KM_plot$plot + annotate("text", x = 190, y = 0.14,
                                label = paste0("HR = ", round(cox_cfDNA_HR, 2)),
                                size = (5/14) * 9)

# 10. Create plots ------------------------------------------------------------

## Common targets
### Create GRange object based off RNA-seq
load(paste0(TCGAdata_dir, "/TCGA_hg19_June2011_exons_gr.RData"))

### Subset to SHOX2, SEPT9, ONECUT2 and TWIST1
SEPT9_gr <- TCGA_hg19_June2011_exons_gr[grepl("SEPT9",
                                              TCGA_hg19_June2011_exons_gr$gene)]

SHOX2_gr <- TCGA_hg19_June2011_exons_gr[grepl("SHOX2",
                                              TCGA_hg19_June2011_exons_gr$gene)]

ONECUT2_gr <- 
  TCGA_hg19_June2011_exons_gr[grepl("ONECUT2",
                                    TCGA_hg19_June2011_exons_gr$gene)]

TWIST1_gr <-
  TCGA_hg19_June2011_exons_gr[grepl("TWIST1",
                                    TCGA_hg19_June2011_exons_gr$gene)]

### For each gene, highlight hyperDMR and hyperDMC probeIDs

#### SEPT9
SEPT9_DMRs <- subsetByOverlaps(hyperDMRs_ROIs_gr,
                               SEPT9_gr, maxgap = 2000L)

SEPT9_s <- min(start(SEPT9_gr))
SEPT9_r <- (max(end(SEPT9_gr)) - min(start(SEPT9_gr))) / 100
SEPT9_e <- max(end(SEPT9_gr))

SEPT9_ggplot <- SEPT9_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "SEPT9",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(SEPT9_s, 
                                       SEPT9_s + SEPT9_r,
                                       SEPT9_s + SEPT9_r * 2,
                                       SEPT9_s + SEPT9_r, 
                                       SEPT9_s),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(SEPT9_s + SEPT9_r * 1.5, 
                                       SEPT9_s + SEPT9_r * 3, 
                                       SEPT9_s + SEPT9_r * 3 + SEPT9_r,
                                       SEPT9_s + SEPT9_r * 3, 
                                       SEPT9_s + SEPT9_r * 1.5,
                                       SEPT9_s + SEPT9_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(SEPT9_s + SEPT9_r * 1.5 + SEPT9_r * 2, 
                                       SEPT9_s + SEPT9_r * 3 + SEPT9_r * 2, 
                                       SEPT9_s + SEPT9_r * 3 + SEPT9_r * 3,
                                       SEPT9_s + SEPT9_r * 3 + SEPT9_r * 2, 
                                       SEPT9_s + SEPT9_r * 1.5 + SEPT9_r * 2,
                                       SEPT9_s + SEPT9_r * 2.5 + SEPT9_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(SEPT9_s + SEPT9_r * 1.5 + SEPT9_r * 4, 
                                       SEPT9_e, 
                                       SEPT9_e,
                                       SEPT9_s + SEPT9_r * 1.5 + SEPT9_r * 4,
                                       SEPT9_s + SEPT9_r * 2.5 + SEPT9_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(SEPT9_DMRs),
           xmax = end(SEPT9_DMRs), ymin = 0.5, ymax = 12.5, 
           col = NA, fill = "red", alpha = 0.8) + 
  scale_x_sequnit(unit = "Mb") + 
  theme_JB() + 
  theme(axis.line.y = element_blank()) 

#### SHOX2
SHOX2_DMRs <- subsetByOverlaps(hyperDMRs_ROIs_gr,
                               SHOX2_gr, maxgap = 2000L)

SHOX2_s <- min(start(SHOX2_gr))
SHOX2_r <- (max(end(SHOX2_gr)) - min(start(SHOX2_gr))) / 100
SHOX2_e <- max(end(SHOX2_gr))

SHOX2_ggplot <- SHOX2_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "SHOX2",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(SHOX2_e, 
                                       SHOX2_e - SHOX2_r,
                                       SHOX2_e - SHOX2_r * 2,
                                       SHOX2_e - SHOX2_r, 
                                       SHOX2_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(SHOX2_e - SHOX2_r * 1.5, 
                                       SHOX2_e - SHOX2_r * 3, 
                                       SHOX2_e - SHOX2_r * 3 - SHOX2_r,
                                       SHOX2_e - SHOX2_r * 3, 
                                       SHOX2_e - SHOX2_r * 1.5,
                                       SHOX2_e - SHOX2_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(SHOX2_e - SHOX2_r * 1.5 - SHOX2_r * 2, 
                                       SHOX2_e - SHOX2_r * 3 - SHOX2_r * 2, 
                                       SHOX2_e - SHOX2_r * 3 - SHOX2_r * 3,
                                       SHOX2_e - SHOX2_r * 3 - SHOX2_r * 2, 
                                       SHOX2_e - SHOX2_r * 1.5 - SHOX2_r * 2,
                                       SHOX2_e - SHOX2_r * 2.5 - SHOX2_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(SHOX2_e - SHOX2_r * 1.5 - SHOX2_r * 4, 
                                       SHOX2_s, 
                                       SHOX2_s,
                                       SHOX2_e - SHOX2_r * 1.5 - SHOX2_r * 4,
                                       SHOX2_e - SHOX2_r * 2.5 - SHOX2_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(SHOX2_DMRs),
           xmax = end(SHOX2_DMRs), ymin = 0.6, ymax = 3.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  scale_x_sequnit(unit = "Mb") + 
  theme_JB() + 
  theme(axis.line.y = element_blank())

#### ONECUT2
ONECUT2_DMRs <- subsetByOverlaps(hyperDMRs_ROIs_gr,
                                 ONECUT2_gr, maxgap = 2000L)

ONECUT2_s <- min(start(ONECUT2_gr))
ONECUT2_r <- (max(end(ONECUT2_gr)) - min(start(ONECUT2_gr))) / 100
ONECUT2_e <- max(end(ONECUT2_gr))

ONECUT2_ggplot <- ONECUT2_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "ONECUT2",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ONECUT2_s, 
                                       ONECUT2_s + ONECUT2_r,
                                       ONECUT2_s + ONECUT2_r * 2,
                                       ONECUT2_s + ONECUT2_r, 
                                       ONECUT2_s),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ONECUT2_s + ONECUT2_r * 1.5, 
                                       ONECUT2_s + ONECUT2_r * 3, 
                                       ONECUT2_s + ONECUT2_r * 3 + ONECUT2_r,
                                       ONECUT2_s + ONECUT2_r * 3, 
                                       ONECUT2_s + ONECUT2_r * 1.5,
                                       ONECUT2_s + ONECUT2_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ONECUT2_s + ONECUT2_r * 1.5 + ONECUT2_r * 2, 
                                       ONECUT2_s + ONECUT2_r * 3 + ONECUT2_r * 2, 
                                       ONECUT2_s + ONECUT2_r * 3 + ONECUT2_r * 3,
                                       ONECUT2_s + ONECUT2_r * 3 + ONECUT2_r * 2, 
                                       ONECUT2_s + ONECUT2_r * 1.5 + ONECUT2_r * 2,
                                       ONECUT2_s + ONECUT2_r * 2.5 + ONECUT2_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ONECUT2_s + ONECUT2_r * 1.5 + ONECUT2_r * 4, 
                                       ONECUT2_e, 
                                       ONECUT2_e,
                                       ONECUT2_s + ONECUT2_r * 1.5 + ONECUT2_r * 4,
                                       ONECUT2_s + ONECUT2_r * 2.5 + ONECUT2_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(ONECUT2_DMRs),
           xmax = end(ONECUT2_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  scale_x_sequnit(unit = "Mb") +
  theme_JB() + 
  theme(axis.line.y = element_blank())

#### TWIST1
TWIST1_DMRs <- subsetByOverlaps(hyperDMRs_ROIs_gr,
                                TWIST1_gr, maxgap = 2000L)

TWIST1_s <- min(start(TWIST1_gr))
TWIST1_r <- (max(end(TWIST1_gr)) - min(start(TWIST1_gr))) / 100
TWIST1_e <- max(end(TWIST1_gr))

TWIST1_ggplot <- TWIST1_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "TWIST1",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(TWIST1_e, 
                                       TWIST1_e - TWIST1_r,
                                       TWIST1_e - TWIST1_r * 2,
                                       TWIST1_e - TWIST1_r, 
                                       TWIST1_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(TWIST1_e - TWIST1_r * 1.5, 
                                       TWIST1_e - TWIST1_r * 3, 
                                       TWIST1_e - TWIST1_r * 3 - TWIST1_r,
                                       TWIST1_e - TWIST1_r * 3, 
                                       TWIST1_e - TWIST1_r * 1.5,
                                       TWIST1_e - TWIST1_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(TWIST1_e - TWIST1_r * 1.5 - TWIST1_r * 2, 
                                       TWIST1_e - TWIST1_r * 3 - TWIST1_r * 2, 
                                       TWIST1_e - TWIST1_r * 3 - TWIST1_r * 3,
                                       TWIST1_e - TWIST1_r * 3 - TWIST1_r * 2, 
                                       TWIST1_e - TWIST1_r * 1.5 - TWIST1_r * 2,
                                       TWIST1_e - TWIST1_r * 2.5 - TWIST1_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(TWIST1_e - TWIST1_r * 1.5 - TWIST1_r * 4, 
                                       TWIST1_s, 
                                       TWIST1_s,
                                       TWIST1_e - TWIST1_r * 1.5 - TWIST1_r * 4,
                                       TWIST1_e - TWIST1_r * 2.5 - TWIST1_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(TWIST1_DMRs),
           xmax = end(TWIST1_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  scale_x_sequnit(unit = "Mb") +
  theme_JB() + 
  theme(axis.line.y = element_blank())

#### Collate
common_targets_right <- cowplot::plot_grid(SHOX2_ggplot@ggplot, 
                                           ONECUT2_ggplot@ggplot, 
                                           TWIST1_ggplot@ggplot, 
                                           nrow = 3, ncol = 1,
                                           rel_heights = c(1.5, 1, 1))
common_targets_plot <- cowplot::plot_grid(SEPT9_ggplot@ggplot, 
                                          common_targets_right,
                                          nrow = 1, ncol = 2)

## Cox hazard plot
### Number of regions w/ significant cox values
sum(cox_pval < 0.05) # 33 regions

### Plot hazard ratios
#### Extract CpGs w/ significant cox hazard ratios
sig_cox_names <- names(sig_cox_probeIDs)
surv_sig_df <- surv_df[,colnames(surv_df) %in% sig_cox_names]

cox_df <- as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(cox_df) <- c("probeID", "up_lim", "score", "low_lim")
for(i in 1:length(sig_cox_names)){
  tmp_form <- as.formula(paste("surv_object ~ ", sig_cox_names[i]))
  res.cox <- coxph(tmp_form, data = surv_df)
  cox.sum <- summary(res.cox)
  
  tmp_df <- data.frame(probeID = sig_cox_names[i],
                       up_lim = cox.sum$conf.int[4],
                       score = cox.sum$conf.int[1],
                       low_lim = cox.sum$conf.int[3])
  cox_df <- rbind(cox_df, tmp_df)
}

cox_df$probeID <- factor(cox_df$probeID, 
                         levels = cox_df$probeID[order(cox_df$score,
                                                       decreasing = FALSE)])

# Annotated 300-bp windows based on nearby genes
load(paste0(data_dir, "/HNC_vs_Norm_Combined_RPKMs_Output.RData"))

cox_df$gene <- as.character(HNC_vs_Norm_Combined_RPKMs_Output[
  match(cox_df$probeID, HNC_vs_Norm_Combined_RPKMs_Output$ID),
]$gene)

# Manually confirmed transcripts based on genome-browser and corrected
cox_df[cox_df$gene == "AC018470.4",]$gene <- "LINC01305"
cox_df[cox_df$gene == "RP11-548O1.3",]$gene <- "LINC01391"
cox_df[cox_df$gene == "RP11-834C11.12",]$gene <- "AC012531.3"
cox_df[cox_df$gene == "CTD-2368P22.1",]$gene <- "AC008969.1"
cox_df[cox_df$gene == "Charlie8",]$gene <- "STK3"
cox_df[cox_df$gene == "RNU6-1175P",]$gene <- "AP003555.3"
cox_df[cox_df$gene == "GATA2",]$gene <- "GATA2-AS1"

genescore_order <- cox_df %>% 
  group_by(gene) %>% 
  summarise(test = mean(score))

cox_df$gene <- 
  factor(cox_df$gene,
         levels = genescore_order$gene[order(genescore_order$test)])

cox_ggplot <- cox_df %>%
  subset(gene %in% c("OSR1", "GATA2-AS1", "LINC01391",
                     "ZSCAN31", "STK3")) %>%
  ggplot(aes(x = gene, y = score,
             ymin = low_lim, ymax = up_lim)) + 
  geom_pointrange(size = 0.3,
                  shape = 15,
                  alpha = 0.8) + 
  geom_hline(yintercept = 1, lty = 2, col = "lightgrey",
             size = 0.26) +
  scale_y_continuous(trans = scales::log10_trans(),
                     limits = c(0.1, 20)) +
  coord_flip() +
  theme_JB() + 
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9)) +
  labs(x = "Regions (300-bp)", y = "Hazard ratio")

cox_ggplot


## Prognostic gene plots
### STK3 
#### Create data.frame of correaltions
STK3_met_vs_exp_cor <- list()
tmpGRange <- GRanges(seqnames = Rle("chr8"),
                     ranges = IRanges(start = 99465919,
                                      end = 99966918))

STK3_allprobeIDs_gr <- subsetByOverlaps(probeID_to_300bp_gr, tmpGRange)

##### HNC
met.hnc.STK3.assay <- 
  met.hnc.primary.assay.nonNA[rownames(met.hnc.primary.assay.nonNA) %in%
                                STK3_allprobeIDs_gr$probeID,]
met.hnc.STK3.regions.assay <- as.data.frame(matrix(ncol = 
                                                     ncol(met.hnc.STK3.assay),
                                                   nrow = 0))
colnames(met.hnc.STK3.regions.assay) <- 
  colnames(met.hnc.STK3.assay)

STK3_allIDs_gr <- unique(STK3_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(STK3_allIDs_gr)){
  tmp_gr <- STK3_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(STK3_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.hnc.STK3.assay,
           rownames(met.hnc.STK3.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.hnc.STK3.regions.assay <- rbind(met.hnc.STK3.regions.assay, tmp_colMeans)
  rownames(met.hnc.STK3.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                   start(tmp_gr),
                                                   end(tmp_gr),
                                                   sep = ".")
}

##### Norm
met.norm.STK3.assay <- 
  met.hnc.normal.assay.nonNA[rownames(met.hnc.normal.assay.nonNA) %in%
                               STK3_allprobeIDs_gr$probeID,]
met.norm.STK3.regions.assay <- as.data.frame(matrix(ncol = 
                                                      ncol(met.norm.STK3.assay),
                                                    nrow = 0))
colnames(met.norm.STK3.regions.assay) <- 
  colnames(met.norm.STK3.assay)

STK3_allIDs_gr <- unique(STK3_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(STK3_allIDs_gr)){
  tmp_gr <- STK3_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(STK3_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.norm.STK3.assay,
           rownames(met.norm.STK3.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.norm.STK3.regions.assay <- rbind(met.norm.STK3.regions.assay, tmp_colMeans)
  rownames(met.norm.STK3.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                    start(tmp_gr),
                                                    end(tmp_gr),
                                                    sep = ".")
}

STK3_hyperprobeIDs_gr <- 
  subsetByOverlaps(probeID_to_300bp_gr, tmpGRange)
STK3_hyperprobeIDs_gr <- 
  STK3_hyperprobeIDs_gr[STK3_hyperprobeIDs_gr$ID %in%
                          rownames(met.hnc.STK3.regions.assay)]
STK3_hyperprobeIDs_gr <- unique(STK3_hyperprobeIDs_gr[,c(-1)])

for(i in 1:length(STK3_hyperprobeIDs_gr)){
  # Obtain beta-values of probeID
  matrix_probeID.index <- rownames(met.hnc.STK3.regions.assay) %in%
    STK3_hyperprobeIDs_gr$ID[i]
  
  matrix_probeID.bval <- met.hnc.STK3.regions.assay[matrix_probeID.index,]
  
  # Correlate w/ STK3 transcript
  matrix_gene.rsem <- RNAseq_assay[rownames(RNAseq_assay) == "uc003yio.2",]
  matrix_gene.rsem <- log(matrix_gene.rsem + 1)
  
  tmp_cor <- cor.test(x = as.numeric(matrix_probeID.bval), 
                      y = as.numeric(matrix_gene.rsem),
                      method = "spearman")
  
  rho <- c(tmp_cor$estimate)
  pval <- c(tmp_cor$p.value)
  tmp_c <- c(rho, pval)
  names(tmp_c) <- c("rho", "p.value")
  tmp_list <- list(tmp_c)
  names(tmp_list) <- rownames(matrix_gene.rsem)
  
  STK3_met_vs_exp_cor[[i]] <- tmp_list
  names(STK3_met_vs_exp_cor)[i] <- STK3_hyperprobeIDs_gr$ID[i]
  
}

## Calculate adjusted p-values of correlations
STK3_met_vs_exp_pval <- c()
STK3_met_vs_exp_rho <- c()
STK3_met_vs_exp_probeID <- c()
STK3_met_vs_exp_gene <- c()
for(i in 1:length(STK3_met_vs_exp_cor)){
  tmp_unlist <- unlist(STK3_met_vs_exp_cor[[i]])
  
  grep_pval_index <- grepl("p.value", 
                           names(tmp_unlist))
  STK3_met_vs_exp_pval <- c(STK3_met_vs_exp_pval, tmp_unlist[grep_pval_index])
  
  grep_rho_index <- grepl("rho",
                          names(tmp_unlist))
  STK3_met_vs_exp_rho <- c(STK3_met_vs_exp_rho, tmp_unlist[grep_rho_index])
  
  tmp_probeID <- rep(names(STK3_met_vs_exp_cor)[i],
                     length(tmp_unlist)/2)
  STK3_met_vs_exp_probeID <- c(STK3_met_vs_exp_probeID, tmp_probeID)
  
  tmp_gene <- gsub(".rho", "", names(tmp_unlist[grep_rho_index]))
  STK3_met_vs_exp_gene <- c(STK3_met_vs_exp_gene, tmp_gene)
  
}

## Collate to dataframe
STK3_met_vs_exp_df <- data.frame(probeID = STK3_met_vs_exp_probeID,
                                 ucsc_ID = STK3_met_vs_exp_gene,
                                 rho = STK3_met_vs_exp_rho,
                                 pval = STK3_met_vs_exp_pval)

## Get genome position of CpGs
STK3_met_vs_exp_df$pos <- (start(STK3_hyperprobeIDs_gr) + 
                             end(STK3_hyperprobeIDs_gr)) / 2

#### Look at all transcripts
STK3_gr <- TCGA_hg19_June2011_exons_gr[grepl("STK3\\|",
                                             TCGA_hg19_June2011_exons_gr$gene)]

STK3_s <- min(start(STK3_gr))
STK3_r <- (max(end(STK3_gr)) - min(start(STK3_gr))) / 100
STK3_e <- max(end(STK3_gr))

STK3_tmp <- GRanges(seqnames = seqnames(STK3_gr[1]),
                    ranges = IRanges(start = STK3_s, end = STK3_e))
STK3_DMRs <- subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr,
                              hyperDMRs_hm450k_gr[hyperDMRs_hm450k_gr$ID 
                                                  %in% cox_df[cox_df$gene %in% c("STK3", "OSR2 (n = 2)"),]$probeID])

STK3_ggplot <- STK3_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "STK3",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(STK3_e, 
                                       STK3_e - STK3_r,
                                       STK3_e - STK3_r * 2,
                                       STK3_e - STK3_r, 
                                       STK3_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(STK3_e - STK3_r * 1.5, 
                                       STK3_e - STK3_r * 3, 
                                       STK3_e - STK3_r * 3 - STK3_r,
                                       STK3_e - STK3_r * 3, 
                                       STK3_e - STK3_r * 1.5,
                                       STK3_e - STK3_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(STK3_e - STK3_r * 1.5 - STK3_r * 2, 
                                       STK3_e - STK3_r * 3 - STK3_r * 2, 
                                       STK3_e - STK3_r * 3 - STK3_r * 3,
                                       STK3_e - STK3_r * 3 - STK3_r * 2, 
                                       STK3_e - STK3_r * 1.5 - STK3_r * 2,
                                       STK3_e - STK3_r * 2.5 - STK3_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(STK3_e - STK3_r * 1.5 - STK3_r * 4, 
                                       STK3_s, 
                                       STK3_s,
                                       STK3_e - STK3_r * 1.5 - STK3_r * 4,
                                       STK3_e - STK3_r * 2.5 - STK3_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(STK3_DMRs),
           xmax = end(STK3_DMRs), ymin = 0.6, ymax = 3.4, 
           col = "red", fill = "red", alpha = 0.5) + 
  annotate("rect", xmin = 9.945e+07, xmax = 9.997e+07,
           ymin = 1.5, ymax = 2.5, col = "black", size = 0.13,
           fill = NA, lty = 2) + 
  scale_x_sequnit(unit = "Mb") +
  theme_JB() + 
  theme(axis.line.y = element_blank()) + 
  xlim(c(99400000, 100000000))

#### Pair w/ RNA expression for significant transcript
STK3_subggplot <- STK3_gr %>% 
  subset(STK3_gr$pair %in% c("uc003yio.2")) %>%
  ggplot() + 
  geom_alignment(group.selfish = TRUE,
                 main = "STK3",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(STK3_e, 
                                       STK3_e - STK3_r,
                                       STK3_e - STK3_r * 2,
                                       STK3_e - STK3_r, 
                                       STK3_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(STK3_e - STK3_r * 1.5, 
                                       STK3_e - STK3_r * 3, 
                                       STK3_e - STK3_r * 3 - STK3_r,
                                       STK3_e - STK3_r * 3, 
                                       STK3_e - STK3_r * 1.5,
                                       STK3_e - STK3_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(STK3_e - STK3_r * 1.5 - STK3_r * 2, 
                                       STK3_e - STK3_r * 3 - STK3_r * 2, 
                                       STK3_e - STK3_r * 3 - STK3_r * 3,
                                       STK3_e - STK3_r * 3 - STK3_r * 2, 
                                       STK3_e - STK3_r * 1.5 - STK3_r * 2,
                                       STK3_e - STK3_r * 2.5 - STK3_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(STK3_e - STK3_r * 1.5 - STK3_r * 4, 
                                       STK3_s, 
                                       STK3_s,
                                       STK3_e - STK3_r * 1.5 - STK3_r * 4,
                                       STK3_e - STK3_r * 2.5 - STK3_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(STK3_DMRs),
           xmax = end(STK3_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  theme_JB() + 
  theme(axis.line.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")) 


figure_STK3_bot <- STK3_subggplot + 
  scale_x_continuous(expand = c(0,0)) + 
  coord_cartesian(xlim = c(99900000, 99975000),
                  expand = 0)

STK3_hyperProbeIDs <- 
  subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr, STK3_hyperprobeIDs_gr)$ID
STK3_met_vs_exp_df$DMR <- ifelse(STK3_met_vs_exp_df$probeID %in% 
                                   STK3_hyperProbeIDs,
                                 "hyper", "ns")
STK3_assay <- 
  met.hnc.STK3.regions.assay[match(STK3_met_vs_exp_df$probeID, 
                                   rownames(met.hnc.STK3.regions.assay)),]
STK3_upper <- apply(STK3_assay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
STK3_lower <- apply(STK3_assay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
STK3_mean <- apply(STK3_assay, 1, function(x) mean(x, na.rm = TRUE))

STK3_nassay <- 
  met.norm.STK3.regions.assay[match(STK3_met_vs_exp_df$probeID,
                                    rownames(met.norm.STK3.regions.assay)),]
STK3_nupper <- apply(STK3_nassay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
STK3_nlower <- apply(STK3_nassay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
STK3_nmean <- apply(STK3_nassay, 1, function(x) mean(x, na.rm = TRUE))

STK3_met_vs_exp_df$tumour_upper <- STK3_upper
STK3_met_vs_exp_df$tumour_lower <- STK3_lower
STK3_met_vs_exp_df$tumour_mean <- STK3_mean

STK3_met_vs_exp_df$normal_upper <- STK3_nupper
STK3_met_vs_exp_df$normal_lower <- STK3_nlower
STK3_met_vs_exp_df$normal_mean <- STK3_nmean

figure_STK3_top <- STK3_met_vs_exp_df %>%
  ggplot(aes(x = pos, y = rho)) + 
  geom_line(size = 0.26) + 
  geom_point(shape = 15, aes(col = DMR),
             size = 2) +
  geom_hline(yintercept = c(0.3, -0.3),
             lty = 2, color = rgb(0,0,0,1/4),
             size = 0.13) + 
  annotate("rect", xmin = start(STK3_DMRs),
           xmax = end(STK3_DMRs), ymin = -1, ymax = 1, 
           col = NA, fill = "red", alpha = 0.2) +
  labs(y = "Methylation vs. expression correlation") + 
  theme_JB() + 
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.title.y.right = element_text(vjust = 3),
        plot.margin = unit(c(5,5,0,5), "mm")) + 
  scale_y_continuous(limits = c(-1, 1), expand = c(0,0)) + 
  scale_x_continuous(limits = c(99900000, 99975000),
                     expand = c(0,0)) +
  scale_color_manual(values = c("red", rgb(0,0,0,1/4)))

figure_STK3 <- egg::ggarrange(figure_STK3_top, 
                              figure_STK3_bot@ggplot + scale_x_sequnit(unit = "Mb"),
                              nrow = 2, ncol = 1, heights = c(0.9, 0.1))

### ZNF323 
#### Create data.frame of correaltions
ZNF323_met_vs_exp_cor <- list()
ZNF323_allprobeIDs_gr <- 
  probeID_to_300bp_gr[probeID_to_300bp_gr$probeID %in%
                        c("cg22468666", "cg10725301", "cg21750589", "cg18866106", 
                          "cg13961449", "cg03849074",
                          "cg05765371", "cg15743356", "cg21138752", "cg08343881",
                          "cg23084951", "cg06903569", "cg05605299", "cg17968647",
                          "cg13997124", "cg19013262", "cg09996797", "cg07369507",
                          "cg26031541", "cg07206202", "cg04686012")]

##### HNC
met.hnc.ZNF323.assay <- 
  met.hnc.primary.assay.nonNA[rownames(met.hnc.primary.assay.nonNA) %in%
                                ZNF323_allprobeIDs_gr$probeID,]
met.hnc.ZNF323.regions.assay <- as.data.frame(matrix(ncol = 
                                                       ncol(met.hnc.ZNF323.assay),
                                                     nrow = 0))
colnames(met.hnc.ZNF323.regions.assay) <- 
  colnames(met.hnc.ZNF323.assay)

ZNF323_allIDs_gr <- unique(ZNF323_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(ZNF323_allIDs_gr)){
  tmp_gr <- ZNF323_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(ZNF323_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.hnc.ZNF323.assay,
           rownames(met.hnc.ZNF323.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.hnc.ZNF323.regions.assay <- rbind(met.hnc.ZNF323.regions.assay, tmp_colMeans)
  rownames(met.hnc.ZNF323.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                     start(tmp_gr),
                                                     end(tmp_gr),
                                                     sep = ".")
}

##### Norm
met.norm.ZNF323.assay <- 
  met.hnc.normal.assay.nonNA[rownames(met.hnc.normal.assay.nonNA) %in%
                               ZNF323_allprobeIDs_gr$probeID,]
met.norm.ZNF323.regions.assay <- as.data.frame(matrix(ncol = 
                                                        ncol(met.norm.ZNF323.assay),
                                                      nrow = 0))
colnames(met.norm.ZNF323.regions.assay) <- 
  colnames(met.norm.ZNF323.assay)

ZNF323_allIDs_gr <- unique(ZNF323_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(ZNF323_allIDs_gr)){
  tmp_gr <- ZNF323_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(ZNF323_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.norm.ZNF323.assay,
           rownames(met.norm.ZNF323.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.norm.ZNF323.regions.assay <- rbind(met.norm.ZNF323.regions.assay, tmp_colMeans)
  rownames(met.norm.ZNF323.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                      start(tmp_gr),
                                                      end(tmp_gr),
                                                      sep = ".")
}

ZNF323_hyperprobeIDs_gr <- unique(ZNF323_allprobeIDs_gr[,c(-1)])

for(i in 1:length(ZNF323_hyperprobeIDs_gr)){
  # Obtain beta-values of probeID
  matrix_probeID.index <- rownames(met.hnc.ZNF323.regions.assay) %in%
    ZNF323_hyperprobeIDs_gr$ID[i]
  
  matrix_probeID.bval <- met.hnc.ZNF323.regions.assay[matrix_probeID.index,]
  
  # Correlate w/ ZNF323 transcript
  matrix_gene.rsem <- RNAseq_assay[rownames(RNAseq_assay) == "uc003nld.2",]
  matrix_gene.rsem <- log(matrix_gene.rsem + 1)
  
  tmp_cor <- cor.test(x = as.numeric(matrix_probeID.bval), 
                      y = as.numeric(matrix_gene.rsem),
                      method = "spearman")
  
  rho <- c(tmp_cor$estimate)
  pval <- c(tmp_cor$p.value)
  tmp_c <- c(rho, pval)
  names(tmp_c) <- c("rho", "p.value")
  tmp_list <- list(tmp_c)
  names(tmp_list) <- rownames(matrix_gene.rsem)
  
  ZNF323_met_vs_exp_cor[[i]] <- tmp_list
  names(ZNF323_met_vs_exp_cor)[i] <- ZNF323_hyperprobeIDs_gr$ID[i]
  
}

## Calculate adjusted p-values of correlations
ZNF323_met_vs_exp_pval <- c()
ZNF323_met_vs_exp_rho <- c()
ZNF323_met_vs_exp_probeID <- c()
ZNF323_met_vs_exp_gene <- c()
for(i in 1:length(ZNF323_met_vs_exp_cor)){
  tmp_unlist <- unlist(ZNF323_met_vs_exp_cor[[i]])
  
  grep_pval_index <- grepl("p.value", 
                           names(tmp_unlist))
  ZNF323_met_vs_exp_pval <- c(ZNF323_met_vs_exp_pval, tmp_unlist[grep_pval_index])
  
  grep_rho_index <- grepl("rho",
                          names(tmp_unlist))
  ZNF323_met_vs_exp_rho <- c(ZNF323_met_vs_exp_rho, tmp_unlist[grep_rho_index])
  
  tmp_probeID <- rep(names(ZNF323_met_vs_exp_cor)[i],
                     length(tmp_unlist)/2)
  ZNF323_met_vs_exp_probeID <- c(ZNF323_met_vs_exp_probeID, tmp_probeID)
  
  tmp_gene <- gsub(".rho", "", names(tmp_unlist[grep_rho_index]))
  ZNF323_met_vs_exp_gene <- c(ZNF323_met_vs_exp_gene, tmp_gene)
  
}

## Collate to dataframe
ZNF323_met_vs_exp_df <- data.frame(probeID = ZNF323_met_vs_exp_probeID,
                                   ucsc_ID = ZNF323_met_vs_exp_gene,
                                   rho = ZNF323_met_vs_exp_rho,
                                   pval = ZNF323_met_vs_exp_pval)

## Get genome position of CpGs
ZNF323_met_vs_exp_df$pos <- (start(ZNF323_hyperprobeIDs_gr) + 
                               end(ZNF323_hyperprobeIDs_gr)) / 2

#### Look at all transcripts
ZNF323_gr <- TCGA_hg19_June2011_exons_gr[grepl("ZNF323\\|",
                                               TCGA_hg19_June2011_exons_gr$gene)]

ZNF323_s <- min(start(ZNF323_gr))
ZNF323_r <- (max(end(ZNF323_gr)) - min(start(ZNF323_gr))) / 100
ZNF323_e <- max(end(ZNF323_gr))

ZNF323_tmp <- GRanges(seqnames = seqnames(ZNF323_gr[1]),
                      ranges = IRanges(start = ZNF323_s, end = ZNF323_e))
ZNF323_DMRs <- subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr,
                                hyperDMRs_hm450k_gr[hyperDMRs_hm450k_gr$ID 
                                                    %in% cox_df[cox_df$gene %in% c("ZSCAN31"),]$probeID])

ZNF323_ggplot <- ZNF323_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "ZNF323",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ZNF323_e, 
                                       ZNF323_e - ZNF323_r,
                                       ZNF323_e - ZNF323_r * 2,
                                       ZNF323_e - ZNF323_r, 
                                       ZNF323_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ZNF323_e - ZNF323_r * 1.5, 
                                       ZNF323_e - ZNF323_r * 3, 
                                       ZNF323_e - ZNF323_r * 3 - ZNF323_r,
                                       ZNF323_e - ZNF323_r * 3, 
                                       ZNF323_e - ZNF323_r * 1.5,
                                       ZNF323_e - ZNF323_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ZNF323_e - ZNF323_r * 1.5 - ZNF323_r * 2, 
                                       ZNF323_e - ZNF323_r * 3 - ZNF323_r * 2, 
                                       ZNF323_e - ZNF323_r * 3 - ZNF323_r * 3,
                                       ZNF323_e - ZNF323_r * 3 - ZNF323_r * 2, 
                                       ZNF323_e - ZNF323_r * 1.5 - ZNF323_r * 2,
                                       ZNF323_e - ZNF323_r * 2.5 - ZNF323_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ZNF323_e - ZNF323_r * 1.5 - ZNF323_r * 4, 
                                       ZNF323_s, 
                                       ZNF323_s,
                                       ZNF323_e - ZNF323_r * 1.5 - ZNF323_r * 4,
                                       ZNF323_e - ZNF323_r * 2.5 - ZNF323_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(ZNF323_DMRs),
           xmax = end(ZNF323_DMRs), ymin = 0.6, ymax = 3.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  annotate("rect", xmin = 9.945e+07, xmax = 9.997e+07,
           ymin = 1.5, ymax = 2.5, col = "black", size = 0.13,
           fill = NA, lty = 2) + 
  scale_x_sequnit(unit = "Mb") +
  theme_JB() + 
  theme(axis.line.y = element_blank()) + 
  xlim(c(28290000, 28330000))

#### Pair w/ RNA expression for significant transcript
ZNF323_subggplot <- ZNF323_gr %>% 
  subset(ZNF323_gr$pair %in% c("uc003nld.2")) %>%
  ggplot() + 
  geom_alignment(group.selfish = TRUE,
                 main = "ZNF323",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ZNF323_e, 
                                       ZNF323_e - ZNF323_r,
                                       ZNF323_e - ZNF323_r * 2,
                                       ZNF323_e - ZNF323_r, 
                                       ZNF323_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ZNF323_e - ZNF323_r * 1.5, 
                                       ZNF323_e - ZNF323_r * 3, 
                                       ZNF323_e - ZNF323_r * 3 - ZNF323_r,
                                       ZNF323_e - ZNF323_r * 3, 
                                       ZNF323_e - ZNF323_r * 1.5,
                                       ZNF323_e - ZNF323_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ZNF323_e - ZNF323_r * 1.5 - ZNF323_r * 2, 
                                       ZNF323_e - ZNF323_r * 3 - ZNF323_r * 2, 
                                       ZNF323_e - ZNF323_r * 3 - ZNF323_r * 3,
                                       ZNF323_e - ZNF323_r * 3 - ZNF323_r * 2, 
                                       ZNF323_e - ZNF323_r * 1.5 - ZNF323_r * 2,
                                       ZNF323_e - ZNF323_r * 2.5 - ZNF323_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(ZNF323_e - ZNF323_r * 1.5 - ZNF323_r * 4, 
                                       ZNF323_s, 
                                       ZNF323_s,
                                       ZNF323_e - ZNF323_r * 1.5 - ZNF323_r * 4,
                                       ZNF323_e - ZNF323_r * 2.5 - ZNF323_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(ZNF323_DMRs),
           xmax = end(ZNF323_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  theme_JB() + 
  theme(axis.line.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")) 


figure_ZNF323_bot <- ZNF323_subggplot + 
  scale_x_continuous(expand = c(0,0)) + 
  coord_cartesian(xlim = c(28290000, 28330000),
                  expand = 0)

ZNF323_hyperProbeIDs <- 
  subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr, ZNF323_hyperprobeIDs_gr)$ID
ZNF323_met_vs_exp_df$DMR <- ifelse(ZNF323_met_vs_exp_df$probeID %in% 
                                     ZNF323_hyperProbeIDs,
                                   "hyper", "ns")
ZNF323_assay <- 
  met.hnc.ZNF323.regions.assay[match(ZNF323_met_vs_exp_df$probeID, 
                                     rownames(met.hnc.ZNF323.regions.assay)),]
ZNF323_upper <- apply(ZNF323_assay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
ZNF323_lower <- apply(ZNF323_assay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
ZNF323_mean <- apply(ZNF323_assay, 1, function(x) mean(x, na.rm = TRUE))

ZNF323_nassay <- 
  met.norm.ZNF323.regions.assay[match(ZNF323_met_vs_exp_df$probeID,
                                      rownames(met.norm.ZNF323.regions.assay)),]
ZNF323_nupper <- apply(ZNF323_nassay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
ZNF323_nlower <- apply(ZNF323_nassay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
ZNF323_nmean <- apply(ZNF323_nassay, 1, function(x) mean(x, na.rm = TRUE))

ZNF323_met_vs_exp_df$tumour_upper <- ZNF323_upper
ZNF323_met_vs_exp_df$tumour_lower <- ZNF323_lower
ZNF323_met_vs_exp_df$tumour_mean <- ZNF323_mean

ZNF323_met_vs_exp_df$normal_upper <- ZNF323_nupper
ZNF323_met_vs_exp_df$normal_lower <- ZNF323_nlower
ZNF323_met_vs_exp_df$normal_mean <- ZNF323_nmean

figure_ZNF323_top <- ZNF323_met_vs_exp_df %>%
  ggplot(aes(x = pos, y = rho)) + 
  geom_line(size = 0.26) + 
  geom_point(shape = 15, aes(col = DMR),
             size = 2) +
  geom_hline(yintercept = c(0.3, -0.3),
             lty = 2, color = rgb(0,0,0,1/4),
             size = 0.13) + 
  annotate("rect", xmin = start(ZNF323_DMRs),
           xmax = end(ZNF323_DMRs), ymin = -1, ymax = 1, 
           col = NA, fill = "red", alpha = 0.2) +
  labs(y = "Methylation vs. expression correlation") + 
  theme_JB() + 
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.title.y.right = element_text(vjust = 3),
        plot.margin = unit(c(5,5,0,5), "mm")) + 
  scale_x_continuous(limits = c(28290000, 28330000),
                     expand = c(0,0)) +
  scale_color_manual(values = c("red", rgb(0,0,0,1/4)))

figure_ZNF323 <- egg::ggarrange(figure_ZNF323_top, 
                                figure_ZNF323_bot@ggplot + scale_x_sequnit(unit = "Mb"),
                                nrow = 2, ncol = 1, heights = c(0.9, 0.1))

### LINC01391 
#### Create data.frame of correaltions
LINC01391_met_vs_exp_cor <- list()
tmpGRange <- GRanges(seqnames = Rle("chr3"),
                     ranges = IRanges(start = 138652795,
                                      end = 138666762))

LINC01391_allprobeIDs_gr <- subsetByOverlaps(probeID_to_300bp_gr, tmpGRange)

##### HNC
met.hnc.LINC01391.assay <- 
  met.hnc.primary.assay.nonNA[rownames(met.hnc.primary.assay.nonNA) %in%
                                LINC01391_allprobeIDs_gr$probeID,]
met.hnc.LINC01391.regions.assay <- as.data.frame(matrix(ncol = 
                                                          ncol(met.hnc.LINC01391.assay),
                                                        nrow = 0))
colnames(met.hnc.LINC01391.regions.assay) <- 
  colnames(met.hnc.LINC01391.assay)

LINC01391_allIDs_gr <- unique(LINC01391_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(LINC01391_allIDs_gr)){
  tmp_gr <- LINC01391_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(LINC01391_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.hnc.LINC01391.assay,
           rownames(met.hnc.LINC01391.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.hnc.LINC01391.regions.assay <- rbind(met.hnc.LINC01391.regions.assay, tmp_colMeans)
  rownames(met.hnc.LINC01391.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                        start(tmp_gr),
                                                        end(tmp_gr),
                                                        sep = ".")
}

##### Norm
met.norm.LINC01391.assay <- 
  met.hnc.normal.assay.nonNA[rownames(met.hnc.normal.assay.nonNA) %in%
                               LINC01391_allprobeIDs_gr$probeID,]
met.norm.LINC01391.regions.assay <- as.data.frame(matrix(ncol = 
                                                           ncol(met.norm.LINC01391.assay),
                                                         nrow = 0))
colnames(met.norm.LINC01391.regions.assay) <- 
  colnames(met.norm.LINC01391.assay)

LINC01391_allIDs_gr <- unique(LINC01391_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(LINC01391_allIDs_gr)){
  tmp_gr <- LINC01391_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(LINC01391_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.norm.LINC01391.assay,
           rownames(met.norm.LINC01391.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.norm.LINC01391.regions.assay <- rbind(met.norm.LINC01391.regions.assay, tmp_colMeans)
  rownames(met.norm.LINC01391.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                         start(tmp_gr),
                                                         end(tmp_gr),
                                                         sep = ".")
}

LINC01391_hyperprobeIDs_gr <- unique(LINC01391_allprobeIDs_gr[,c(-1)])

for(i in 1:length(LINC01391_hyperprobeIDs_gr)){
  # Obtain beta-values of probeID
  matrix_probeID.index <- rownames(met.hnc.LINC01391.regions.assay) %in%
    LINC01391_hyperprobeIDs_gr$ID[i]
  
  matrix_probeID.bval <- met.hnc.LINC01391.regions.assay[matrix_probeID.index,]
  
  # Correlate w/ LINC01391 transcript
  matrix_gene.rsem <- RNAseq_assay[rownames(RNAseq_assay) == "uc003esv.1",]
  matrix_gene.rsem <- log(matrix_gene.rsem + 1)
  
  tmp_cor <- cor.test(x = as.numeric(matrix_probeID.bval), 
                      y = as.numeric(matrix_gene.rsem),
                      method = "spearman")
  
  rho <- c(tmp_cor$estimate)
  pval <- c(tmp_cor$p.value)
  tmp_c <- c(rho, pval)
  names(tmp_c) <- c("rho", "p.value")
  tmp_list <- list(tmp_c)
  names(tmp_list) <- rownames(matrix_gene.rsem)
  
  LINC01391_met_vs_exp_cor[[i]] <- tmp_list
  names(LINC01391_met_vs_exp_cor)[i] <- LINC01391_hyperprobeIDs_gr$ID[i]
  
}

## Calculate adjusted p-values of correlations
LINC01391_met_vs_exp_pval <- c()
LINC01391_met_vs_exp_rho <- c()
LINC01391_met_vs_exp_probeID <- c()
LINC01391_met_vs_exp_gene <- c()
for(i in 1:length(LINC01391_met_vs_exp_cor)){
  tmp_unlist <- unlist(LINC01391_met_vs_exp_cor[[i]])
  
  grep_pval_index <- grepl("p.value", 
                           names(tmp_unlist))
  LINC01391_met_vs_exp_pval <- c(LINC01391_met_vs_exp_pval, tmp_unlist[grep_pval_index])
  
  grep_rho_index <- grepl("rho",
                          names(tmp_unlist))
  LINC01391_met_vs_exp_rho <- c(LINC01391_met_vs_exp_rho, tmp_unlist[grep_rho_index])
  
  tmp_probeID <- rep(names(LINC01391_met_vs_exp_cor)[i],
                     length(tmp_unlist)/2)
  LINC01391_met_vs_exp_probeID <- c(LINC01391_met_vs_exp_probeID, tmp_probeID)
  
  tmp_gene <- gsub(".rho", "", names(tmp_unlist[grep_rho_index]))
  LINC01391_met_vs_exp_gene <- c(LINC01391_met_vs_exp_gene, tmp_gene)
  
}

## Collate to dataframe
LINC01391_met_vs_exp_df <- data.frame(probeID = LINC01391_met_vs_exp_probeID,
                                      ucsc_ID = LINC01391_met_vs_exp_gene,
                                      rho = LINC01391_met_vs_exp_rho,
                                      pval = LINC01391_met_vs_exp_pval)

## Get genome position of CpGs
LINC01391_met_vs_exp_df$pos <- (start(LINC01391_hyperprobeIDs_gr) + 
                                  end(LINC01391_hyperprobeIDs_gr)) / 2

#### Look at all transcripts
LINC01391_gr <- TCGA_hg19_June2011_exons_gr[grepl("uc003esv.1",
                                                  TCGA_hg19_June2011_exons_gr$pair)]

LINC01391_s <- min(start(LINC01391_gr))
LINC01391_r <- (max(end(LINC01391_gr)) - min(start(LINC01391_gr))) / 100
LINC01391_e <- max(end(LINC01391_gr))

LINC01391_tmp <- GRanges(seqnames = seqnames(LINC01391_gr[1]),
                         ranges = IRanges(start = LINC01391_s, end = LINC01391_e))
LINC01391_DMRs <- subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr,
                                   hyperDMRs_hm450k_gr[hyperDMRs_hm450k_gr$ID 
                                                       %in% cox_df[cox_df$gene %in% c("LINC01391"),]$probeID])

LINC01391_ggplot <- LINC01391_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "LINC01391",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(LINC01391_e, 
                                       LINC01391_e - LINC01391_r,
                                       LINC01391_e - LINC01391_r * 2,
                                       LINC01391_e - LINC01391_r, 
                                       LINC01391_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(LINC01391_e - LINC01391_r * 1.5, 
                                       LINC01391_e - LINC01391_r * 3, 
                                       LINC01391_e - LINC01391_r * 3 - LINC01391_r,
                                       LINC01391_e - LINC01391_r * 3, 
                                       LINC01391_e - LINC01391_r * 1.5,
                                       LINC01391_e - LINC01391_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(LINC01391_e - LINC01391_r * 1.5 - LINC01391_r * 2, 
                                       LINC01391_e - LINC01391_r * 3 - LINC01391_r * 2, 
                                       LINC01391_e - LINC01391_r * 3 - LINC01391_r * 3,
                                       LINC01391_e - LINC01391_r * 3 - LINC01391_r * 2, 
                                       LINC01391_e - LINC01391_r * 1.5 - LINC01391_r * 2,
                                       LINC01391_e - LINC01391_r * 2.5 - LINC01391_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(LINC01391_e - LINC01391_r * 1.5 - LINC01391_r * 4, 
                                       LINC01391_s, 
                                       LINC01391_s,
                                       LINC01391_e - LINC01391_r * 1.5 - LINC01391_r * 4,
                                       LINC01391_e - LINC01391_r * 2.5 - LINC01391_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(LINC01391_DMRs),
           xmax = end(LINC01391_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) +
  scale_x_sequnit(unit = "Mb") +
  theme_JB() + 
  theme(axis.line.y = element_blank()) + 
  xlim(c(138652000, 138664000))

#### Pair w/ RNA expression for significant transcript
LINC01391_subggplot <- LINC01391_gr %>% 
  subset(LINC01391_gr$pair %in% c("uc003esv.1")) %>%
  ggplot() + 
  geom_alignment(group.selfish = TRUE,
                 main = "LINC01391",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(LINC01391_e, 
                                       LINC01391_e - LINC01391_r,
                                       LINC01391_e - LINC01391_r * 2,
                                       LINC01391_e - LINC01391_r, 
                                       LINC01391_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(LINC01391_e - LINC01391_r * 1.5, 
                                       LINC01391_e - LINC01391_r * 3, 
                                       LINC01391_e - LINC01391_r * 3 - LINC01391_r,
                                       LINC01391_e - LINC01391_r * 3, 
                                       LINC01391_e - LINC01391_r * 1.5,
                                       LINC01391_e - LINC01391_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(LINC01391_e - LINC01391_r * 1.5 - LINC01391_r * 2, 
                                       LINC01391_e - LINC01391_r * 3 - LINC01391_r * 2, 
                                       LINC01391_e - LINC01391_r * 3 - LINC01391_r * 3,
                                       LINC01391_e - LINC01391_r * 3 - LINC01391_r * 2, 
                                       LINC01391_e - LINC01391_r * 1.5 - LINC01391_r * 2,
                                       LINC01391_e - LINC01391_r * 2.5 - LINC01391_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(LINC01391_e - LINC01391_r * 1.5 - LINC01391_r * 4, 
                                       LINC01391_s, 
                                       LINC01391_s,
                                       LINC01391_e - LINC01391_r * 1.5 - LINC01391_r * 4,
                                       LINC01391_e - LINC01391_r * 2.5 - LINC01391_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(LINC01391_DMRs),
           xmax = end(LINC01391_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  theme_JB() + 
  theme(axis.line.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")) 


figure_LINC01391_bot <- LINC01391_subggplot + 
  scale_x_continuous(expand = c(0,0)) + 
  coord_cartesian(xlim = c(138652000, 138664000),
                  expand = 0)

LINC01391_hyperProbeIDs <- 
  subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr, LINC01391_hyperprobeIDs_gr)$ID
LINC01391_met_vs_exp_df$DMR <- ifelse(LINC01391_met_vs_exp_df$probeID %in% 
                                        LINC01391_hyperProbeIDs,
                                      "hyper", "ns")
LINC01391_assay <- 
  met.hnc.LINC01391.regions.assay[match(LINC01391_met_vs_exp_df$probeID, 
                                        rownames(met.hnc.LINC01391.regions.assay)),]
LINC01391_upper <- apply(LINC01391_assay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
LINC01391_lower <- apply(LINC01391_assay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
LINC01391_mean <- apply(LINC01391_assay, 1, function(x) mean(x, na.rm = TRUE))

LINC01391_nassay <- 
  met.norm.LINC01391.regions.assay[match(LINC01391_met_vs_exp_df$probeID,
                                         rownames(met.norm.LINC01391.regions.assay)),]
LINC01391_nupper <- apply(LINC01391_nassay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
LINC01391_nlower <- apply(LINC01391_nassay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
LINC01391_nmean <- apply(LINC01391_nassay, 1, function(x) mean(x, na.rm = TRUE))

LINC01391_met_vs_exp_df$tumour_upper <- LINC01391_upper
LINC01391_met_vs_exp_df$tumour_lower <- LINC01391_lower
LINC01391_met_vs_exp_df$tumour_mean <- LINC01391_mean

LINC01391_met_vs_exp_df$normal_upper <- LINC01391_nupper
LINC01391_met_vs_exp_df$normal_lower <- LINC01391_nlower
LINC01391_met_vs_exp_df$normal_mean <- LINC01391_nmean

figure_LINC01391_top <- LINC01391_met_vs_exp_df %>%
  ggplot(aes(x = pos, y = rho)) + 
  geom_line(size = 0.26) + 
  geom_point(shape = 15, aes(col = DMR),
             size = 2) +
  geom_hline(yintercept = c(0.3, -0.3),
             lty = 2, color = rgb(0,0,0,1/4),
             size = 0.13) + 
  annotate("rect", xmin = start(LINC01391_DMRs),
           xmax = end(LINC01391_DMRs), ymin = -1, ymax = 1, 
           col = NA, fill = "red", alpha = 0.2) +
  labs(y = "Methylation vs. expression correlation") + 
  theme_JB() + 
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.title.y.right = element_text(vjust = 3),
        plot.margin = unit(c(5,5,0,5), "mm")) + 
  scale_y_continuous(limits = c(-1, 1), expand = c(0,0)) + 
  scale_x_continuous(limits = c(138652000, 138664000),
                     expand = c(0,0)) +
  scale_color_manual(values = c("red", rgb(0,0,0,1/4)))

figure_LINC01391 <- egg::ggarrange(figure_LINC01391_top, 
                                   figure_LINC01391_bot@ggplot + scale_x_sequnit(unit = "Mb"),
                                   nrow = 2, ncol = 1, heights = c(0.9, 0.1))

### GATA2AS1 
#### Create data.frame of correaltions
GATA2AS1_met_vs_exp_cor <- list()
tmpGRange <- GRanges(seqnames = Rle("chr3"),
                     ranges = IRanges(start = 128208075,
                                      end = 128221190))

GATA2AS1_allprobeIDs_gr <- subsetByOverlaps(probeID_to_300bp_gr, tmpGRange)

##### HNC
met.hnc.GATA2AS1.assay <- 
  met.hnc.primary.assay.nonNA[rownames(met.hnc.primary.assay.nonNA) %in%
                                GATA2AS1_allprobeIDs_gr$probeID,]
met.hnc.GATA2AS1.regions.assay <- as.data.frame(matrix(ncol = 
                                                         ncol(met.hnc.GATA2AS1.assay),
                                                       nrow = 0))
colnames(met.hnc.GATA2AS1.regions.assay) <- 
  colnames(met.hnc.GATA2AS1.assay)

GATA2AS1_allIDs_gr <- unique(GATA2AS1_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(GATA2AS1_allIDs_gr)){
  tmp_gr <- GATA2AS1_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(GATA2AS1_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.hnc.GATA2AS1.assay,
           rownames(met.hnc.GATA2AS1.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.hnc.GATA2AS1.regions.assay <- rbind(met.hnc.GATA2AS1.regions.assay, tmp_colMeans)
  rownames(met.hnc.GATA2AS1.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                       start(tmp_gr),
                                                       end(tmp_gr),
                                                       sep = ".")
}

##### Norm
met.norm.GATA2AS1.assay <- 
  met.hnc.normal.assay.nonNA[rownames(met.hnc.normal.assay.nonNA) %in%
                               GATA2AS1_allprobeIDs_gr$probeID,]
met.norm.GATA2AS1.regions.assay <- as.data.frame(matrix(ncol = 
                                                          ncol(met.norm.GATA2AS1.assay),
                                                        nrow = 0))
colnames(met.norm.GATA2AS1.regions.assay) <- 
  colnames(met.norm.GATA2AS1.assay)

GATA2AS1_allIDs_gr <- unique(GATA2AS1_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(GATA2AS1_allIDs_gr)){
  tmp_gr <- GATA2AS1_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(GATA2AS1_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.norm.GATA2AS1.assay,
           rownames(met.norm.GATA2AS1.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.norm.GATA2AS1.regions.assay <- rbind(met.norm.GATA2AS1.regions.assay, tmp_colMeans)
  rownames(met.norm.GATA2AS1.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                        start(tmp_gr),
                                                        end(tmp_gr),
                                                        sep = ".")
}

GATA2AS1_hyperprobeIDs_gr <- unique(GATA2AS1_allprobeIDs_gr[,c(-1)])

for(i in 1:length(GATA2AS1_hyperprobeIDs_gr)){
  # Obtain beta-values of probeID
  matrix_probeID.index <- rownames(met.hnc.GATA2AS1.regions.assay) %in%
    GATA2AS1_hyperprobeIDs_gr$ID[i]
  
  matrix_probeID.bval <- met.hnc.GATA2AS1.regions.assay[matrix_probeID.index,]
  
  # Correlate w/ GATA2AS1 transcript
  matrix_gene.rsem <- RNAseq_assay[rownames(RNAseq_assay) == "uc003ekp.2",]
  matrix_gene.rsem <- log(matrix_gene.rsem + 1)
  
  tmp_cor <- cor.test(x = as.numeric(matrix_probeID.bval), 
                      y = as.numeric(matrix_gene.rsem),
                      method = "spearman")
  
  rho <- c(tmp_cor$estimate)
  pval <- c(tmp_cor$p.value)
  tmp_c <- c(rho, pval)
  names(tmp_c) <- c("rho", "p.value")
  tmp_list <- list(tmp_c)
  names(tmp_list) <- rownames(matrix_gene.rsem)
  
  GATA2AS1_met_vs_exp_cor[[i]] <- tmp_list
  names(GATA2AS1_met_vs_exp_cor)[i] <- GATA2AS1_hyperprobeIDs_gr$ID[i]
  
}

## Calculate adjusted p-values of correlations
GATA2AS1_met_vs_exp_pval <- c()
GATA2AS1_met_vs_exp_rho <- c()
GATA2AS1_met_vs_exp_probeID <- c()
GATA2AS1_met_vs_exp_gene <- c()
for(i in 1:length(GATA2AS1_met_vs_exp_cor)){
  tmp_unlist <- unlist(GATA2AS1_met_vs_exp_cor[[i]])
  
  grep_pval_index <- grepl("p.value", 
                           names(tmp_unlist))
  GATA2AS1_met_vs_exp_pval <- c(GATA2AS1_met_vs_exp_pval, tmp_unlist[grep_pval_index])
  
  grep_rho_index <- grepl("rho",
                          names(tmp_unlist))
  GATA2AS1_met_vs_exp_rho <- c(GATA2AS1_met_vs_exp_rho, tmp_unlist[grep_rho_index])
  
  tmp_probeID <- rep(names(GATA2AS1_met_vs_exp_cor)[i],
                     length(tmp_unlist)/2)
  GATA2AS1_met_vs_exp_probeID <- c(GATA2AS1_met_vs_exp_probeID, tmp_probeID)
  
  tmp_gene <- gsub(".rho", "", names(tmp_unlist[grep_rho_index]))
  GATA2AS1_met_vs_exp_gene <- c(GATA2AS1_met_vs_exp_gene, tmp_gene)
  
}

## Collate to dataframe
GATA2AS1_met_vs_exp_df <- data.frame(probeID = GATA2AS1_met_vs_exp_probeID,
                                     ucsc_ID = GATA2AS1_met_vs_exp_gene,
                                     rho = GATA2AS1_met_vs_exp_rho,
                                     pval = GATA2AS1_met_vs_exp_pval)

## Get genome position of CpGs
GATA2AS1_met_vs_exp_df$pos <- (start(GATA2AS1_hyperprobeIDs_gr) + 
                                 end(GATA2AS1_hyperprobeIDs_gr)) / 2

#### Look at all transcripts
GATA2AS1_gr <- TCGA_hg19_June2011_exons_gr[grepl("uc003ekp.2",
                                                 TCGA_hg19_June2011_exons_gr$pair)]

GATA2AS1_s <- min(start(GATA2AS1_gr))
GATA2AS1_r <- (max(end(GATA2AS1_gr)) - min(start(GATA2AS1_gr))) / 100
GATA2AS1_e <- max(end(GATA2AS1_gr))

GATA2AS1_tmp <- GRanges(seqnames = seqnames(GATA2AS1_gr[1]),
                        ranges = IRanges(start = GATA2AS1_s, end = GATA2AS1_e))
GATA2AS1_DMRs <- subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr,
                                  hyperDMRs_hm450k_gr[hyperDMRs_hm450k_gr$ID 
                                                      %in% cox_df[cox_df$gene %in% c("GATA2"),]$probeID])

GATA2AS1_ggplot <- GATA2AS1_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "GATA2AS1",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(GATA2AS1_e, 
                                       GATA2AS1_e - GATA2AS1_r,
                                       GATA2AS1_e - GATA2AS1_r * 2,
                                       GATA2AS1_e - GATA2AS1_r, 
                                       GATA2AS1_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(GATA2AS1_e - GATA2AS1_r * 1.5, 
                                       GATA2AS1_e - GATA2AS1_r * 3, 
                                       GATA2AS1_e - GATA2AS1_r * 3 - GATA2AS1_r,
                                       GATA2AS1_e - GATA2AS1_r * 3, 
                                       GATA2AS1_e - GATA2AS1_r * 1.5,
                                       GATA2AS1_e - GATA2AS1_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(GATA2AS1_e - GATA2AS1_r * 1.5 - GATA2AS1_r * 2, 
                                       GATA2AS1_e - GATA2AS1_r * 3 - GATA2AS1_r * 2, 
                                       GATA2AS1_e - GATA2AS1_r * 3 - GATA2AS1_r * 3,
                                       GATA2AS1_e - GATA2AS1_r * 3 - GATA2AS1_r * 2, 
                                       GATA2AS1_e - GATA2AS1_r * 1.5 - GATA2AS1_r * 2,
                                       GATA2AS1_e - GATA2AS1_r * 2.5 - GATA2AS1_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(GATA2AS1_e - GATA2AS1_r * 1.5 - GATA2AS1_r * 4, 
                                       GATA2AS1_s, 
                                       GATA2AS1_s,
                                       GATA2AS1_e - GATA2AS1_r * 1.5 - GATA2AS1_r * 4,
                                       GATA2AS1_e - GATA2AS1_r * 2.5 - GATA2AS1_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(GATA2AS1_DMRs),
           xmax = end(GATA2AS1_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) +
  scale_x_sequnit(unit = "Mb") +
  theme_JB() + 
  theme(axis.line.y = element_blank()) + 
  xlim(c(128208075, 128221190))

#### Pair w/ RNA expression for significant transcript
GATA2AS1_subggplot <- GATA2AS1_gr %>% 
  subset(GATA2AS1_gr$pair %in% c("uc003ekp.2")) %>%
  ggplot() + 
  geom_alignment(group.selfish = TRUE,
                 main = "GATA2-AS1",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(GATA2AS1_e, 
                                       GATA2AS1_e - GATA2AS1_r,
                                       GATA2AS1_e - GATA2AS1_r * 2,
                                       GATA2AS1_e - GATA2AS1_r, 
                                       GATA2AS1_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(GATA2AS1_e - GATA2AS1_r * 1.5, 
                                       GATA2AS1_e - GATA2AS1_r * 3, 
                                       GATA2AS1_e - GATA2AS1_r * 3 - GATA2AS1_r,
                                       GATA2AS1_e - GATA2AS1_r * 3, 
                                       GATA2AS1_e - GATA2AS1_r * 1.5,
                                       GATA2AS1_e - GATA2AS1_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(GATA2AS1_e - GATA2AS1_r * 1.5 - GATA2AS1_r * 2, 
                                       GATA2AS1_e - GATA2AS1_r * 3 - GATA2AS1_r * 2, 
                                       GATA2AS1_e - GATA2AS1_r * 3 - GATA2AS1_r * 3,
                                       GATA2AS1_e - GATA2AS1_r * 3 - GATA2AS1_r * 2, 
                                       GATA2AS1_e - GATA2AS1_r * 1.5 - GATA2AS1_r * 2,
                                       GATA2AS1_e - GATA2AS1_r * 2.5 - GATA2AS1_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(GATA2AS1_e - GATA2AS1_r * 1.5 - GATA2AS1_r * 4, 
                                       GATA2AS1_s, 
                                       GATA2AS1_s,
                                       GATA2AS1_e - GATA2AS1_r * 1.5 - GATA2AS1_r * 4,
                                       GATA2AS1_e - GATA2AS1_r * 2.5 - GATA2AS1_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(GATA2AS1_DMRs),
           xmax = end(GATA2AS1_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  theme_JB() + 
  theme(axis.line.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")) 


figure_GATA2AS1_bot <- GATA2AS1_subggplot + 
  scale_x_continuous(expand = c(0,0)) + 
  coord_cartesian(xlim = c(128208075, 128221190),
                  expand = 0)

GATA2AS1_hyperProbeIDs <- 
  subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr, GATA2AS1_hyperprobeIDs_gr)$ID
GATA2AS1_met_vs_exp_df$DMR <- ifelse(GATA2AS1_met_vs_exp_df$probeID %in% 
                                       GATA2AS1_hyperProbeIDs,
                                     "hyper", "ns")
GATA2AS1_assay <- 
  met.hnc.GATA2AS1.regions.assay[match(GATA2AS1_met_vs_exp_df$probeID, 
                                       rownames(met.hnc.GATA2AS1.regions.assay)),]
GATA2AS1_upper <- apply(GATA2AS1_assay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
GATA2AS1_lower <- apply(GATA2AS1_assay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
GATA2AS1_mean <- apply(GATA2AS1_assay, 1, function(x) mean(x, na.rm = TRUE))

GATA2AS1_nassay <- 
  met.norm.GATA2AS1.regions.assay[match(GATA2AS1_met_vs_exp_df$probeID,
                                        rownames(met.norm.GATA2AS1.regions.assay)),]
GATA2AS1_nupper <- apply(GATA2AS1_nassay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
GATA2AS1_nlower <- apply(GATA2AS1_nassay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
GATA2AS1_nmean <- apply(GATA2AS1_nassay, 1, function(x) mean(x, na.rm = TRUE))

GATA2AS1_met_vs_exp_df$tumour_upper <- GATA2AS1_upper
GATA2AS1_met_vs_exp_df$tumour_lower <- GATA2AS1_lower
GATA2AS1_met_vs_exp_df$tumour_mean <- GATA2AS1_mean

GATA2AS1_met_vs_exp_df$normal_upper <- GATA2AS1_nupper
GATA2AS1_met_vs_exp_df$normal_lower <- GATA2AS1_nlower
GATA2AS1_met_vs_exp_df$normal_mean <- GATA2AS1_nmean

figure_GATA2AS1_top <- GATA2AS1_met_vs_exp_df %>%
  ggplot(aes(x = pos, y = rho)) + 
  geom_line(size = 0.26) + 
  geom_point(shape = 15, aes(col = DMR),
             size = 2) +
  geom_hline(yintercept = c(0.3, -0.3),
             lty = 2, color = rgb(0,0,0,1/4),
             size = 0.13) + 
  annotate("rect", xmin = start(GATA2AS1_DMRs),
           xmax = end(GATA2AS1_DMRs), ymin = -1, ymax = 1, 
           col = NA, fill = "red", alpha = 0.2) +
  labs(y = "Methylation vs. expression correlation") + 
  theme_JB() + 
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.title.y.right = element_text(vjust = 3),
        plot.margin = unit(c(5,5,0,5), "mm")) + 
  scale_y_continuous(limits = c(-1, 1), expand = c(0,0)) + 
  scale_x_continuous(limits = c(128208075, 128221190),
                     expand = c(0,0)) +
  scale_color_manual(values = c("red", rgb(0,0,0,1/4)))

figure_GATA2AS1 <- egg::ggarrange(figure_GATA2AS1_top, 
                                  figure_GATA2AS1_bot@ggplot + scale_x_sequnit(unit = "Mb"),
                                  nrow = 2, ncol = 1, heights = c(0.9, 0.1))

### OSR1 
#### Create data.frame of correaltions
OSR1_met_vs_exp_cor <- list()
tmpGRange <- GRanges(seqnames = Rle("chr2"),
                     ranges = IRanges(start = 19551000,
                                      end = 19559000))

OSR1_allprobeIDs_gr <- subsetByOverlaps(probeID_to_300bp_gr, tmpGRange)

##### HNC
met.hnc.OSR1.assay <- 
  met.hnc.primary.assay.nonNA[rownames(met.hnc.primary.assay.nonNA) %in%
                                OSR1_allprobeIDs_gr$probeID,]
met.hnc.OSR1.regions.assay <- as.data.frame(matrix(ncol = 
                                                     ncol(met.hnc.OSR1.assay),
                                                   nrow = 0))
colnames(met.hnc.OSR1.regions.assay) <- 
  colnames(met.hnc.OSR1.assay)

OSR1_allIDs_gr <- unique(OSR1_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(OSR1_allIDs_gr)){
  tmp_gr <- OSR1_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(OSR1_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.hnc.OSR1.assay,
           rownames(met.hnc.OSR1.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.hnc.OSR1.regions.assay <- rbind(met.hnc.OSR1.regions.assay, tmp_colMeans)
  rownames(met.hnc.OSR1.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                   start(tmp_gr),
                                                   end(tmp_gr),
                                                   sep = ".")
}

##### Norm
met.norm.OSR1.assay <- 
  met.hnc.normal.assay.nonNA[rownames(met.hnc.normal.assay.nonNA) %in%
                               OSR1_allprobeIDs_gr$probeID,]
met.norm.OSR1.regions.assay <- as.data.frame(matrix(ncol = 
                                                      ncol(met.norm.OSR1.assay),
                                                    nrow = 0))
colnames(met.norm.OSR1.regions.assay) <- 
  colnames(met.norm.OSR1.assay)

OSR1_allIDs_gr <- unique(OSR1_allprobeIDs_gr[,c(-1,-2)])

for(i in 1:length(OSR1_allIDs_gr)){
  tmp_gr <- OSR1_allIDs_gr[i]
  probeIDs <- subsetByOverlaps(OSR1_allprobeIDs_gr,
                               tmp_gr)$probeID
  tmp_df <- 
    subset(met.norm.OSR1.assay,
           rownames(met.norm.OSR1.assay) %in% probeIDs)
  tmp_colMeans <- colMeans(tmp_df, na.rm = TRUE)
  met.norm.OSR1.regions.assay <- rbind(met.norm.OSR1.regions.assay, tmp_colMeans)
  rownames(met.norm.OSR1.regions.assay)[i] <- paste(seqnames(tmp_gr),
                                                    start(tmp_gr),
                                                    end(tmp_gr),
                                                    sep = ".")
}

OSR1_hyperprobeIDs_gr <- unique(OSR1_allprobeIDs_gr[,c(-1)])

for(i in 1:length(OSR1_hyperprobeIDs_gr)){
  # Obtain beta-values of probeID
  matrix_probeID.index <- rownames(met.hnc.OSR1.regions.assay) %in%
    OSR1_hyperprobeIDs_gr$ID[i]
  
  matrix_probeID.bval <- met.hnc.OSR1.regions.assay[matrix_probeID.index,]
  
  # Correlate w/ OSR1 transcript
  matrix_gene.rsem <- RNAseq_assay[rownames(RNAseq_assay) == "uc002rdc.2",]
  matrix_gene.rsem <- log(matrix_gene.rsem + 1)
  
  tmp_cor <- cor.test(x = as.numeric(matrix_probeID.bval), 
                      y = as.numeric(matrix_gene.rsem),
                      method = "spearman")
  
  rho <- c(tmp_cor$estimate)
  pval <- c(tmp_cor$p.value)
  tmp_c <- c(rho, pval)
  names(tmp_c) <- c("rho", "p.value")
  tmp_list <- list(tmp_c)
  names(tmp_list) <- rownames(matrix_gene.rsem)
  
  OSR1_met_vs_exp_cor[[i]] <- tmp_list
  names(OSR1_met_vs_exp_cor)[i] <- OSR1_hyperprobeIDs_gr$ID[i]
  
}

## Calculate adjusted p-values of correlations
OSR1_met_vs_exp_pval <- c()
OSR1_met_vs_exp_rho <- c()
OSR1_met_vs_exp_probeID <- c()
OSR1_met_vs_exp_gene <- c()
for(i in 1:length(OSR1_met_vs_exp_cor)){
  tmp_unlist <- unlist(OSR1_met_vs_exp_cor[[i]])
  
  grep_pval_index <- grepl("p.value", 
                           names(tmp_unlist))
  OSR1_met_vs_exp_pval <- c(OSR1_met_vs_exp_pval, tmp_unlist[grep_pval_index])
  
  grep_rho_index <- grepl("rho",
                          names(tmp_unlist))
  OSR1_met_vs_exp_rho <- c(OSR1_met_vs_exp_rho, tmp_unlist[grep_rho_index])
  
  tmp_probeID <- rep(names(OSR1_met_vs_exp_cor)[i],
                     length(tmp_unlist)/2)
  OSR1_met_vs_exp_probeID <- c(OSR1_met_vs_exp_probeID, tmp_probeID)
  
  tmp_gene <- gsub(".rho", "", names(tmp_unlist[grep_rho_index]))
  OSR1_met_vs_exp_gene <- c(OSR1_met_vs_exp_gene, tmp_gene)
  
}

## Collate to dataframe
OSR1_met_vs_exp_df <- data.frame(probeID = OSR1_met_vs_exp_probeID,
                                 ucsc_ID = OSR1_met_vs_exp_gene,
                                 rho = OSR1_met_vs_exp_rho,
                                 pval = OSR1_met_vs_exp_pval)

## Get genome position of CpGs
OSR1_met_vs_exp_df$pos <- (start(OSR1_hyperprobeIDs_gr) + 
                             end(OSR1_hyperprobeIDs_gr)) / 2

#### Look at all transcripts
OSR1_gr <- TCGA_hg19_June2011_exons_gr[grepl("uc002rdc.2",
                                             TCGA_hg19_June2011_exons_gr$pair)]

OSR1_s <- min(start(OSR1_gr))
OSR1_r <- (max(end(OSR1_gr)) - min(start(OSR1_gr))) / 100
OSR1_e <- max(end(OSR1_gr))

OSR1_tmp <- GRanges(seqnames = seqnames(OSR1_gr[1]),
                    ranges = IRanges(start = OSR1_s, end = OSR1_e))
OSR1_DMRs <- subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr,
                              hyperDMRs_hm450k_gr[hyperDMRs_hm450k_gr$ID 
                                                  %in% cox_df[cox_df$gene %in% c("OSR1"),]$probeID])

OSR1_ggplot <- OSR1_gr %>% ggplot() + 
  geom_alignment(group.selfish = FALSE,
                 stat = "stepping",
                 aes(group = pair),
                 main = "OSR1",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(OSR1_e, 
                                       OSR1_e - OSR1_r,
                                       OSR1_e - OSR1_r * 2,
                                       OSR1_e - OSR1_r, 
                                       OSR1_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(OSR1_e - OSR1_r * 1.5, 
                                       OSR1_e - OSR1_r * 3, 
                                       OSR1_e - OSR1_r * 3 - OSR1_r,
                                       OSR1_e - OSR1_r * 3, 
                                       OSR1_e - OSR1_r * 1.5,
                                       OSR1_e - OSR1_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(OSR1_e - OSR1_r * 1.5 - OSR1_r * 2, 
                                       OSR1_e - OSR1_r * 3 - OSR1_r * 2, 
                                       OSR1_e - OSR1_r * 3 - OSR1_r * 3,
                                       OSR1_e - OSR1_r * 3 - OSR1_r * 2, 
                                       OSR1_e - OSR1_r * 1.5 - OSR1_r * 2,
                                       OSR1_e - OSR1_r * 2.5 - OSR1_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(OSR1_e - OSR1_r * 1.5 - OSR1_r * 4, 
                                       OSR1_s, 
                                       OSR1_s,
                                       OSR1_e - OSR1_r * 1.5 - OSR1_r * 4,
                                       OSR1_e - OSR1_r * 2.5 - OSR1_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(OSR1_DMRs),
           xmax = end(OSR1_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) +
  scale_x_sequnit(unit = "Mb") +
  theme_JB() + 
  theme(axis.line.y = element_blank()) + 
  xlim(c(19551000, 19559000))

#### Pair w/ RNA expression for significant transcript
OSR1_subggplot <- OSR1_gr %>% 
  subset(OSR1_gr$pair %in% c("uc002rdc.2")) %>%
  ggplot() + 
  geom_alignment(group.selfish = TRUE,
                 main = "OSR1",
                 xlab = "Genomic location",
                 gap.geom = "segment",
                 col = "steelblue",
                 fill = "steelblue") + 
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(OSR1_e, 
                                       OSR1_e - OSR1_r,
                                       OSR1_e - OSR1_r * 2,
                                       OSR1_e - OSR1_r, 
                                       OSR1_e),
                                 y = c(0, 0, 0.2, 0.4, 0.4)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(OSR1_e - OSR1_r * 1.5, 
                                       OSR1_e - OSR1_r * 3, 
                                       OSR1_e - OSR1_r * 3 - OSR1_r,
                                       OSR1_e - OSR1_r * 3, 
                                       OSR1_e - OSR1_r * 1.5,
                                       OSR1_e - OSR1_r * 2.5),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(OSR1_e - OSR1_r * 1.5 - OSR1_r * 2, 
                                       OSR1_e - OSR1_r * 3 - OSR1_r * 2, 
                                       OSR1_e - OSR1_r * 3 - OSR1_r * 3,
                                       OSR1_e - OSR1_r * 3 - OSR1_r * 2, 
                                       OSR1_e - OSR1_r * 1.5 - OSR1_r * 2,
                                       OSR1_e - OSR1_r * 2.5 - OSR1_r * 2),
                                 y = c(0, 0, 0.2, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  geom_polygon(aes(x = x, y = y),
               inherit.aes = FALSE, 
               data = data.frame(x = c(OSR1_e - OSR1_r * 1.5 - OSR1_r * 4, 
                                       OSR1_s, 
                                       OSR1_s,
                                       OSR1_e - OSR1_r * 1.5 - OSR1_r * 4,
                                       OSR1_e - OSR1_r * 2.5 - OSR1_r * 4),
                                 y = c(0, 0, 0.4, 0.4, 0.2)),
               fill = "darkblue", col = NA) +
  annotate("rect", xmin = start(OSR1_DMRs),
           xmax = end(OSR1_DMRs), ymin = 0.6, ymax = 1.4, 
           col = NA, fill = "red", alpha = 0.5) + 
  theme_JB() + 
  theme(axis.line.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")) 


figure_OSR1_bot <- OSR1_subggplot + 
  scale_x_continuous(expand = c(0,0)) + 
  coord_cartesian(xlim = c(19551000, 19559000),
                  expand = 0)

OSR1_hyperProbeIDs <- 
  subsetByOverlaps(hyperprobeIDs_plasmaDMR_gr, OSR1_hyperprobeIDs_gr)$ID
OSR1_met_vs_exp_df$DMR <- ifelse(OSR1_met_vs_exp_df$probeID %in% 
                                   OSR1_hyperProbeIDs,
                                 "hyper", "ns")
OSR1_assay <- 
  met.hnc.OSR1.regions.assay[match(OSR1_met_vs_exp_df$probeID, 
                                   rownames(met.hnc.OSR1.regions.assay)),]
OSR1_upper <- apply(OSR1_assay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
OSR1_lower <- apply(OSR1_assay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
OSR1_mean <- apply(OSR1_assay, 1, function(x) mean(x, na.rm = TRUE))

OSR1_nassay <- 
  met.norm.OSR1.regions.assay[match(OSR1_met_vs_exp_df$probeID,
                                    rownames(met.norm.OSR1.regions.assay)),]
OSR1_nupper <- apply(OSR1_nassay, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
OSR1_nlower <- apply(OSR1_nassay, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
OSR1_nmean <- apply(OSR1_nassay, 1, function(x) mean(x, na.rm = TRUE))

OSR1_met_vs_exp_df$tumour_upper <- OSR1_upper
OSR1_met_vs_exp_df$tumour_lower <- OSR1_lower
OSR1_met_vs_exp_df$tumour_mean <- OSR1_mean

OSR1_met_vs_exp_df$normal_upper <- OSR1_nupper
OSR1_met_vs_exp_df$normal_lower <- OSR1_nlower
OSR1_met_vs_exp_df$normal_mean <- OSR1_nmean

figure_OSR1_top <- OSR1_met_vs_exp_df %>%
  ggplot(aes(x = pos, y = rho)) + 
  geom_line(size = 0.26) + 
  geom_point(shape = 15, aes(col = DMR),
             size = 2) +
  geom_hline(yintercept = c(0.3, -0.3),
             lty = 2, color = rgb(0,0,0,1/4),
             size = 0.13) + 
  annotate("rect", xmin = start(OSR1_DMRs),
           xmax = end(OSR1_DMRs), ymin = -1, ymax = 1, 
           col = NA, fill = "red", alpha = 0.2) +
  labs(y = "Methylation vs. expression correlation") + 
  theme_JB() + 
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.title.y.right = element_text(vjust = 3),
        plot.margin = unit(c(5,5,0,5), "mm")) + 
  scale_y_continuous(limits = c(-1, 1), expand = c(0,0)) + 
  scale_x_continuous(limits = c(19551000, 19559000),
                     expand = c(0,0)) +
  scale_color_manual(values = c("red", rgb(0,0,0,1/4)))

figure_OSR1 <- egg::ggarrange(figure_OSR1_top, 
                              figure_OSR1_bot@ggplot + scale_x_sequnit(unit = "Mb"),
                              nrow = 2, ncol = 1, heights = c(0.9, 0.1))

