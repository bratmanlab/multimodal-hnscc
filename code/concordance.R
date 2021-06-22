# cfMeDIP-seq and CAPP-Seq concordance analysis source code #
#############################################################

# Set directories
scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"
fig_dir <- "../figures"

# Load libraries
library(survminer)
library(survival)
library(ggdendro)
library(DESeq2)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(egg)
library(pROC)
library(cutpointr)
library(tidyverse)
source(paste0(scripts, "/ggplot_themeJB.R"))

# Create CAPP-seq vs. cfMeDIP-seq abundance data.frame
## Convert matrix counts to RPKM values of cfMeDIP-seq samples
load(paste0(data_dir, "/cfDNA_Frag100_150_Counts.RData"))

## Save total reads for division
RPKM_denom <- colSums(cfDNA_Frag100_150_Counts) / 1000000

## Subset to PBLdepleted windows
load(paste0(out_dir, "/PBLdepleted_IDs.RData"))

cfDNA_CombinedRPKMs_PBLdepleted <- 
  cfDNA_Frag100_150_Counts[rownames(cfDNA_Frag100_150_Counts) %in% 
                                    PBLdepleted_IDs,]

cfDNA_CombinedRPKMs_PBLdepleted <- t(t(cfDNA_CombinedRPKMs_PBLdepleted / 0.3) /
                                       RPKM_denom)

rm(cfDNA_Frag100_150_Counts)

## Subset to hyperDMRs
load(paste0(out_dir, "/DESeq2_cfDNA_DMRs.RData"))

cfDNA_CombinedRPKMs_hyperDMRs <- 
  cfDNA_CombinedRPKMs_PBLdepleted[rownames(cfDNA_CombinedRPKMs_PBLdepleted) %in%
                                    rownames(res_list$res_padj_hyperDMR),]

mat <- cfDNA_CombinedRPKMs_hyperDMRs

## Create data.frame with phenotypic information of samples
pheno.v <- data.frame(condition = ifelse(grepl("Norm", colnames(mat)),
                                         "Norm", "HNC"))
rownames(pheno.v) <- paste0("HN", substr(colnames(mat), 1, 4))

### Add ctDNA abundance information
load(paste0(out_dir, "/cfDNA_PBLfiltered_SNVs_meanMAF.RData"))
meanMAF <- cfDNA_PBLfiltered_SNVs_meanMAF$meanMAF
names(meanMAF) <- cfDNA_PBLfiltered_SNVs_meanMAF$upn

donorMAF <- rep(0, sum(pheno.v$condition == "Norm"))
names(donorMAF) <- rownames(subset(pheno.v, pheno.v$condition == "Norm"))

meanMAF <- c(meanMAF, donorMAF)

pheno.v$ctDNA <- meanMAF[match(rownames(pheno.v), names(meanMAF))]

## Compare metrics of each cluster between HNSCC cfDNA samples
mclust_df <- pheno.v[pheno.v$condition == "HNC",]
HNC_mat <- mat[,!grepl("Norm", colnames(mat))]
mclust_df$meanRPKM <- colMeans(HNC_mat)

### Create data.frame containing healthy donor RPKMs
mclust_wnorm_df <- pheno.v

cfDNA_hyperDMR_colMeans <- colMeans(cfDNA_CombinedRPKMs_hyperDMRs)
names(cfDNA_hyperDMR_colMeans) <- 
  paste0("HN", substr(colnames(cfDNA_CombinedRPKMs_hyperDMRs), 1, 4))

mclust_wnorm_df$meanRPKM <- 
  cfDNA_hyperDMR_colMeans[match(rownames(mclust_wnorm_df), 
                                names(cfDNA_hyperDMR_colMeans))]

max_norm_RPKM <-
  max(mclust_wnorm_df[mclust_wnorm_df$condition == "Norm",]$meanRPKM)

mvm_cor <- 
  cor.test(x = mclust_wnorm_df[mclust_wnorm_df$condition == "HNC",]$meanRPKM,
           y = mclust_wnorm_df[mclust_wnorm_df$condition == "HNC",]$ctDNA)

MvR_cor <- round(mvm_cor$estimate, 2)
MvR_pval <- round(mvm_cor$p.value, 9)

### Create ggplot
hyperDMR_corplot <- mclust_wnorm_df %>%
  ggplot(aes(x = meanRPKM, y = ctDNA, colour = condition)) + 
  geom_smooth(data = subset(mclust_wnorm_df, mclust_wnorm_df$condition == "HNC"),
              aes(x = meanRPKM, y = ctDNA),
              method = lm, fill = "lightgrey", fullrange = TRUE,
              show.legend = FALSE, col = "red", size = 0.13) + 
  geom_point(size = 1.5,
             stroke = 1,
             aes(shape = condition),
             alpha = 0.5) +
  annotate(geom = "text",
           x = 0.03, y = 4.5,
           label = paste0("R = ",MvR_cor,"\np = ",
                          MvR_pval),
           hjust = 0,
           size = (5/14) * 9) + 
  geom_vline(xintercept = max_norm_RPKM, 
             lty = 2, col = "blue", size = 0.13) +
  coord_cartesian(ylim = c(-0.1, 5), xlim = c(0, 1.2),
                  expand = FALSE) +
  labs(y = "Mutant allele frequency (%)",
       x = "Methylation (mean RPKM)") + 
  theme_JB() + 
  scale_colour_JB("Condition:") + 
  scale_shape_manual("Condition:", values = c(16, 4))

# Load table w/ patient clinical information and methylation clusters
clinical_relapse <- 
  read.csv(paste0(data_dir, 
                  "/HNC_prognosis_info_w_relapse_updated_nov182019.csv"))

# Name column w/ sample identifiers
colnames(clinical_relapse)[1] <- "HN_ID"

clinical_relapse$hyperDMR_meanRPKM <- mclust_df$meanRPKM

ctDNA_positive <- clinical_relapse[clinical_relapse$ctDNA > 0,]
ctDNA_negative <- clinical_relapse[clinical_relapse$ctDNA == 0,]

KM_df <- as.data.frame(matrix(nrow = 30, ncol=0))
rownames(KM_df) <- clinical_relapse$HN_ID

futime <- c()
fustat <- c()

# Evaluate association of ctDNA detection with clinical stage
clin2 <- mclust_wnorm_df[!grepl("Norm", mclust_wnorm_df$condition),]
clin2$stage <- clinical_relapse$Stage
clin2$detection <- ifelse(clin2$ctDNA >0 & clin2$meanRPKM > max_norm_RPKM,
                          "positve", "negative")
clin2$stage2 <- ifelse(clin2$stage %in% c("I", "II"),
                       "early", "advanced")

# Fisher's exact test
clin2Data <- table(clin2$detection, clin2$stage2)
clin2Test <- fisher.test(clin2Data)


## Construct kaplan-meier based on overall survival
for(i in 1:nrow(clinical_relapse)){
  if(clinical_relapse$Vital.Status[i] == "Dead"){
    tmp <- c()
    tmp <- c(tmp, as.Date(clinical_relapse$Last.FU[i],
                          format = "%m-%d-%Y") - 
               as.Date(clinical_relapse$Date.Dx[i],
                       format = "%m-%d-%Y"))
    
    futime <- c(futime, min(tmp))
    fustat <- c(fustat, 1)
  } else {
    futime <- c(futime, as.Date(clinical_relapse$Last.FU[i], 
                                format = "%m-%d-%Y") -
                  as.Date(clinical_relapse$Date.Dx[i], 
                          format = "%m-%d-%Y"))
    fustat <- c(fustat, 0)
  }
}

KM_df$futime <- futime
KM_df$fustat <- fustat

## Stratify patients based on ctDNA detection by CAPP-Seq
KM_df$ctDNA <- ifelse(clinical_relapse$ctDNA > 0, "positive", "negative")
KM_df$ctDNA <- factor(KM_df$ctDNA, levels = c("positive", "negative"))

## Stratify patients based on ctDNA methylation by cfMeDIP-seq
meth <- 
  rownames(mclust_wnorm_df[mclust_wnorm_df$meanRPKM >
                             max(mclust_wnorm_df[mclust_wnorm_df$condition == 
                                                   "Norm",]$meanRPKM),])

KM_df$meth_pos <- ifelse(rownames(KM_df) %in% substr(meth, 3, 6),
                         "positive", "negative")

## Stratify patients based on detection of both metrics
KM_df$both_pos <- ifelse(KM_df$ctDNA == "positive" & KM_df$meth_pos == "positive",
                         "positive", "negative")

## Stratify patients by early/advanced stage (I+II vs. III+IVa)
KM_df$stage <- clin2$stage2

## Fit kaplan meier based on each
surv_object <- Surv(time = KM_df$futime, event = KM_df$fustat)
### ctDNA
ctDNA_fit <- survfit(surv_object ~ ctDNA, data = KM_df)
### methylation 
meth_pos_fit <- survfit(surv_object ~ meth_pos, data = KM_df)
### both
both_pos_fit <- survfit(surv_object ~ both_pos, data = KM_df)
### stage
stage_fit <- survfit(surv_object ~ stage, data = KM_df)

## Get hazard ratios of each
### ctDNA
#### Rearrange factors for HR analysis
KM_df$ctDNA <- factor(KM_df$ctDNA, levels = c("negative", "positive"))
ctDNA_coxph <- coxph(surv_object ~ ctDNA, data = KM_df)
ctDNA_HR <- summary(ctDNA_coxph)$coefficient[2]

### methylation
#### Rearrange factors for meth_pos analysis
KM_df$meth_pos <- factor(KM_df$meth_pos, levels = c("negative", "positive"))
meth_pos_coxph <- coxph(surv_object ~ meth_pos, data = KM_df)
meth_pos_HR <- summary(meth_pos_coxph)$coefficient[2]

### both
#### Rearrange factors for both_pos analysis
KM_df$both_pos <- factor(KM_df$both_pos, levels = c("negative", "positive"))
both_pos_coxph <- coxph(surv_object ~ both_pos, data = KM_df)
both_pos_HR <- summary(both_pos_coxph)$coefficient[2]

### stage
#### Rearrange factors for both_pos analysis
KM_df$stage <- factor(KM_df$stage, levels = c("early", "advanced"))
stage_coxph <- coxph(surv_object ~ stage, data = KM_df)
stage_HR <- summary(stage_coxph)$coefficient[2]

# Create ggplots
## ggplot theme for KM curves
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

## ctDNA
ggsurv_BL_ctDNA <- 
  ggsurvplot(ctDNA_fit, data = KM_df, pval = TRUE, risk.table = TRUE,
             tables.height = 0.3, ggtheme = theme_JBborder(), size = 0.26,
             palette = palette_JB[1:2], risk.table.y.text = FALSE,
             pval.size = (5/14) * 9, pval.method = TRUE,
             legend.title = "CAPP-Seq", legend.labs = c("positive", "negative"),
             risk.table.fontsize = 4, conf.int = FALSE, 
             risk.table.col = "strata",
             xlab = "Time (months)", censor.shape = 124, censor.size = 2,
             ylab = "Overall survival",
             xscale = "d_m", break.time.by = 365.25) 

ggsurv_BL_ctDNA <- 
  ggarrange(ggsurv_BL_ctDNA$plot +
              annotate("text", x = 180, y = 0.13,
                       label = paste0("HR = ", round(ctDNA_HR, 2)),
                       size = (5/14) * 9),
            ncol = 1, nrow = 1,
            draw = FALSE)

## methylation clusters
ggsurv_BL_meth_pos <- 
  ggsurvplot(meth_pos_fit, data = KM_df, pval = TRUE, risk.table = TRUE,
             tables.height = 0.3, ggtheme = theme_JBborder(), size = 0.26,
             palette = palette_JB[2:1], risk.table.y.text = FALSE,
             pval.size = (5/14) * 9, pval.method = TRUE,
             legend.title = "cfMeDIP-seq", 
             legend.labs = c("negative", "positive"),
             risk.table.fontsize = 4, conf.int = FALSE, 
             risk.table.col = "strata",
             xlab = "Time (months)", censor.shape = 124, censor.size = 2,
             ylab = "Overall survival",
             xscale = "d_m", break.time.by = 365.25) 

ggsurv_BL_meth_pos <- 
  ggarrange(ggsurv_BL_meth_pos$plot +
              annotate("text", x = 180, y = 0.13,
                       label = paste0("HR = ", round(meth_pos_HR, 2)),
                       size = (5/14) * 9),
            ncol = 1, nrow = 1,
            draw = FALSE)

## both
ggsurv_BL_both_pos <- 
  ggsurvplot(both_pos_fit, data = KM_df, pval = TRUE, risk.table = TRUE,
             tables.height = 0.3, ggtheme = theme_JBborder(), size = 0.26,
             palette = palette_JB[2:1], risk.table.y.text = FALSE,
             pval.size = (5/14) * 9, pval.method = TRUE,
             legend.title = "CAPP-Seq + cfMeDIP-seq", 
             legend.labs = c("negative", "positive"),
             risk.table.fontsize = 4, conf.int = FALSE, 
             risk.table.col = "strata",
             xlab = "Time (months)", censor.shape = 124, censor.size = 2,
             ylab = "Overall survival",
             xscale = "d_m", break.time.by = 365.25) 

ggsurv_BL_both_pos <- 
  ggarrange(ggsurv_BL_both_pos$plot +
              annotate("text", x = 180, y = 0.13,
                       label = paste0("HR = ", round(both_pos_HR, 2)),
                       size = (5/14) * 9),
            ncol = 1, nrow = 1,
            draw = FALSE)

## stage
ggsurv_BL_stage <- 
  ggsurvplot(stage_fit, data = KM_df, pval = TRUE, risk.table = TRUE,
             tables.height = 0.3, ggtheme = theme_JBborder(), size = 0.26,
             palette = palette_JB[1:2], risk.table.y.text = FALSE,
             pval.size = (5/14) * 9, pval.method = TRUE,
             legend.title = "Clinical stage", 
             legend.labs = c("advanced (III/IVa)", "early (I/II)"),
             risk.table.fontsize = 4, conf.int = FALSE, 
             risk.table.col = "strata",
             xlab = "Time (months)", censor.shape = 124, censor.size = 2,
             ylab = "Overall survival",
             xscale = "d_m", break.time.by = 365.25) 

ggsurv_BL_stage <- 
  ggarrange(ggsurv_BL_stage$plot +
              annotate("text", x = 180, y = 0.13,
                       label = paste0("HR = ", round(stage_HR, 2)),
                       size = (5/14) * 9),
            ncol = 1, nrow = 1,
            draw = FALSE)

ggsave(plot = as_ggplot(ggsurv_BL_stage), height = 3.54, width = 7,
              dpi = 300, filename = paste0("baseline_stage_OS_KMplot.pdf"),
              useDingbats = FALSE)

# Perform Area Under the Curve (AUC) analysis
## Custom function for AUROC and eventual cross-validation 
AUC_calc <- function(MAFvsDMR_obj, ctDNA_only=FALSE){
  
  # ctDNA_only = TRUE
  if(isTRUE(ctDNA_only)){
    ctDNApos_index <- (MAFvsDMR_obj$condition == "HNC" & MAFvsDMR_obj$ctDNA > 0) | MAFvsDMR_obj$condition == "Norm"
    MAFvsDMR_obj <- MAFvsDMR_obj[ctDNApos_index,]
  }
  
  # Create ROC_medip functin
  ROC_medip <- function(MAFvsDMR_obj){
    
    # Order samples from highest to lowest RPKM
    MAFvsDMR_obj_ordered <- MAFvsDMR_obj[order(MAFvsDMR_obj$meanRPKM, decreasing = TRUE),]
    
    # Calculate category and predition
    ## Category
    category <- c()
    for(i in 1:nrow(MAFvsDMR_obj_ordered)){
      if(MAFvsDMR_obj_ordered$condition[i] == "HNC"){
        category <- c(category, 1)
      } else {
        category <- c(category, 0)
      }
    }
    ## Prediction
    prediction <- rev(seq_along(category))
    
    # Calculate ROC
    roc_obj <- pROC::roc(category, prediction)
    
    # Return object
    return(roc_obj)
  }
  
  # Calculate ROC with function
  rocObj <- ROC_medip(MAFvsDMR_obj = MAFvsDMR_obj)
  
  # Return object
  return(rocObj)
  
}

## Calculate AUC for ctDNA positive samples only
hyperDMR_AUC_ctDNApos <- AUC_calc(MAFvsDMR_obj = mclust_wnorm_df, 
                                  ctDNA_only = TRUE)

## Calculate AUC for all HNSCC samples
hyperDMR_AUC_all <- AUC_calc(MAFvsDMR_obj = mclust_wnorm_df,
                             ctDNA_only = FALSE)

## Plot results
### Combine AUC analysis of ctDNA only and all conditions
roc.list <- list(hyperDMR_AUC_ctDNApos, hyperDMR_AUC_all)

AUROC_plot <- ggroc(roc.list, size = 0.52) + 
  theme_JB() + 
  scale_colour_JB("HNSCC samples:", labels = c("SNV only", "All")) + 
  ylab("Sensitivity") + 
  xlab("Specificity") + 
  annotate(geom = "text", x = 0.1, y = 0.03,
           label = paste0("AUC = ",round(hyperDMR_AUC_all$auc, 3)),
           col = palette_JB[2],
           size = (5/14) * 9) + 
  annotate(geom = "text", x = 0.1, y = 0.1,
           label = paste0("AUC = ", round(hyperDMR_AUC_ctDNApos$auc, 3)),
           col = palette_JB[1], 
           size = (5/14) * 9) + 
  annotate(geom = "segment", x = 1, xend = 0, y = 0, yend = 1,
           colour = "red", size = 0.13, lty = 2)

# Concordance of DMR calling and comparison to CAPP-Seq metrics

## Create 10 training and test cohorts by bootstrap
### For each bootstrap, select HNSCC patients with detectable cTDNA
SNVplusNorm <- mclust_wnorm_df[(mclust_wnorm_df$condition == "HNC" & 
                                  mclust_wnorm_df$ctDNA > 0) | 
                                 mclust_wnorm_df$condition == "Norm",]
SNVplusNorm <- SNVplusNorm[order(SNVplusNorm$ctDNA),]
rownames(SNVplusNorm) <- substr(rownames(SNVplusNorm), 3, 6)

### Assign labels (donors = 0, increasing ctDNA abundance = 1 - 5)
SNVplusNorm$rank <- c(rep(0, 20), rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5))

### Create function for selection of samples for training and testing sets 
### using a 60/40% split respectively
subsample_60_40_DMR <- function(){
  # Select 12 normals
  tmp0 <- rownames(SNVplusNorm[SNVplusNorm$rank == 0,])[sample(1:20, 12)]
  # Select 3 samples at rank 1
  tmp1 <- rownames(SNVplusNorm[SNVplusNorm$rank == 1,])[sample(1:5, 3)]
  # Select 3 samples at rank 2
  tmp2 <- rownames(SNVplusNorm[SNVplusNorm$rank == 2,])[sample(1:5, 3)]
  # Select 3 samples at rank 3
  tmp3 <- rownames(SNVplusNorm[SNVplusNorm$rank == 3,])[sample(1:5, 3)]
  # Select 3 samples at rank 4
  tmp4 <- rownames(SNVplusNorm[SNVplusNorm$rank == 4,])[sample(1:5, 3)]
  
  # Add to training samples
  training_samples <- c(tmp0, tmp1, tmp2, tmp3, tmp4)
  
  # Place remaining in validation samples
  validation_samples <- rownames(SNVplusNorm[!(rownames(SNVplusNorm) %in% 
                                                 training_samples),])
  
  output_list <- list(training_samples, validation_samples)
  names(output_list) <- c("training", "validation")
  return(output_list)
}

### Join each bootstrap as list
set.seed(42)
bootstrap_60_40 <- replicate(50, subsample_60_40_DMR())

### Load cfDNA counts from HNSCC and healthy donor samples

load(paste0(data_dir, "/cfDNA_CombinedCounts_AllWindows.RData"))

# Subset matrix on PBLdepleted regions previously generated in
# PBL_depleted_windows.R script
load(paste0(out_dir, "/PBLdepleted_IDs.RData"))
cfDNA_CombinedCounts_PBLdepleted <-
  cfDNA_CombinedCounts_AllWindows[rownames(cfDNA_CombinedCounts_AllWindows) %in%
                                    PBLdepleted_IDs,]

rm(cfDNA_CombinedCounts_AllWindows)

## Create function for replicates
AUC_cross_calc <- function(ID){
  
  # Create nested functions
  ## DMR function:
  DESeq2_DMR_figs <- function(ID, grepl_list = FALSE, ...){
    
    # Subset matrix by IDs and filtered samples ------
    
    # grepl_list = FALSE
    if(isFALSE(grepl_list)){
      grepl_list <- paste0(c(bootstrap_60_40[[1,1]], bootstrap_60_40[[2,1]]), 
                           collapse = "|")
    } else {
      grepl_list <- paste0(grepl_list, collapse = "|")
    }
    
    HNC_vs_Norm_CombinedCounts_IDsubset <- 
      cfDNA_CombinedCounts_PBLdepleted[
        rownames(cfDNA_CombinedCounts_PBLdepleted) %in% ID, 
        grepl(grepl_list, colnames(cfDNA_CombinedCounts_PBLdepleted))]
    
    # Perform DESeq2 analysis ------
    ## Create colData object
    coldata <- 
      data.frame(condition = 
                   ifelse(grepl("Norm", 
                                colnames(HNC_vs_Norm_CombinedCounts_IDsubset)),
                          "Norm", "HNC"))
    rownames(coldata) <- colnames(HNC_vs_Norm_CombinedCounts_IDsubset)
    
    ## Perform DMR analysis
    dds <- DESeqDataSetFromMatrix(countData = HNC_vs_Norm_CombinedCounts_IDsubset,
                                  colData = coldata,
                                  design = ~ condition)
    dds <- dds[rowSums(counts(dds)) >= 10,]
    dds <- DESeq(dds)
    
    res <- results(dds, contrast = c("condition", "HNC", "Norm"))
    
    ## Select significant hyper/hypomethylated regions
    res_padj <- subset(res, res$padj < 0.1)
    res_padj_hyperDMR <- subset(res_padj, res_padj$log2FoldChange > 0)
    res_padj_hypoDMR <- subset(res_padj, res_padj$log2FoldChange < 0)
    
    
    # Return output as list -----
    
    res_list <- list(res, res_padj_hyperDMR, res_padj_hypoDMR)
    names(res_list) <- c("res", "res_padj_hyperDMR", "res_padj_hypoDMR")
    
    return(res_list)
    
  }
  
  ## MAFvsDMR function:
  MAFvsDMR <- function(res_list, grepl_list = FALSE){
    
    dmr_index <- rownames(cfDNA_CombinedRPKMs_PBLdepleted) %in%
      rownames(res_list[['res_padj_hyperDMR']])
    
    hyperDMR <- colMeans(cfDNA_CombinedRPKMs_PBLdepleted[dmr_index,])
    
    MeanMAF_vs_DMR_tmp <- mclust_wnorm_df
    MeanMAF_vs_DMR_tmp$meanRPKM <- hyperDMR
    
    if(isFALSE(grepl_list) == FALSE){
      grepl_list <- grepl(paste0(grepl_list, collapse = "|"), 
                          rownames(MeanMAF_vs_DMR_tmp))
      MeanMAF_vs_DMR_tmp <- MeanMAF_vs_DMR_tmp[grepl_list, ]
    }
    
    return(MeanMAF_vs_DMR_tmp)
  }
  
  # Perform DESeq2 analysis for each replicate
  AUC_crossval_list <- list()
  
  for(i in 1:ncol(bootstrap_60_40)){
    
    cat("performing DESeq2 analysis on replicate", i, "...\n")
    AUC_crossval_list[[i]] <- DESeq2_DMR_figs(ID, grepl_list = bootstrap_60_40[[1,i]])
    
  }
  
  names(AUC_crossval_list) <- paste0("replicate_", 1:ncol(bootstrap_60_40))
  
  # Remove cases with 0 DMRs
  hyperDMR_ind <- sapply(AUC_crossval_list, FUN = function(x){
    nrow(x$res_padj_hyperDMR) > 1
  })
  
  # Create MAFvsDMR object in validation cohort of each replicate
  MAFvsDMR_list <- list()
  
  for(i in 1:length(AUC_crossval_list[hyperDMR_ind])){
    
    cat("performing MAFvsDMR function on replicate", i, "...\n")
    MAFvsDMR_list[[i]] <- MAFvsDMR(AUC_crossval_list[hyperDMR_ind][[i]], 
                                   grepl_list = bootstrap_60_40[,hyperDMR_ind][[2,i]])
    
  }
  
  names(MAFvsDMR_list) <- paste0("replicate_", 1:ncol(bootstrap_60_40[,hyperDMR_ind]))
  
  # Generate roc curve for each replicate
  ROC_medip <- function(MAFvsDMR_obj){
    
    # Order samples from highest to lowest RPKM
    MAFvsDMR_obj_ordered <- MAFvsDMR_obj[order(MAFvsDMR_obj$meanRPKM, decreasing = TRUE),]
    
    # Calculate category and predition
    ## Category
    category <- c()
    for(i in 1:nrow(MAFvsDMR_obj_ordered)){
      if(MAFvsDMR_obj_ordered$condition[i] == "HNC"){
        category <- c(category, 1)
      } else {
        category <- c(category, 0)
      }
    }
    ## Prediction
    prediction <- rev(seq_along(category))
    
    # Calculate ROC
    roc_obj <- pROC::roc(category, prediction)
    
    # Return object
    return(roc_obj)
  }
  
  # Calculate ROC with function
  rocObj_replicates <- lapply(MAFvsDMR_list, FUN=function(x) ROC_medip(x))
  
  # Calculate AUC for each replicate
  auc_replicates <- sapply(MAFvsDMR_list, FUN=function(x) pROC::auc(AUC_calc(x)))
  
  # Create return list
  output_list <- list(AUC_crossval_list, MAFvsDMR_list, rocObj_replicates, 
                      auc_replicates)
  names(output_list) <- c("res", "MAFvsDMR", "roc", "auc")
  
  # Return object
  return(output_list)
  
}

AUC_crossval <- AUC_cross_calc(PBLdepleted_IDs)

save(AUC_crossval, file = paste0(out_dir, "/AUC_crossval.RData"))

## Calculate median AUC value
medauc <- median(AUC_crossval$auc)

## Plot results
crossval_df <- data.frame(AUC = AUC_crossval$auc)

set.seed(42)
crossval_plot <- crossval_df %>%
  ggplot(aes(x = "", y = AUC)) +
  geom_jitter(position = position_jitter(0.35),
              shape = 16, size = 1.5, col = "black",
              alpha = 0.2) +
  geom_boxplot(col = "black", fill = NA,
               lwd = 0.13, fatten = 6,
               outlier.shape = NA) +
  annotate("text", x = 0.95, y = 0,
           label = paste0("Median = ",
                          round(medauc, 3)),
           size = (5/14) * 9) +
  theme_JB() +
  theme(axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank()) +
  scale_fill_JB() +
  theme(legend.position = "none") +
  xlab("") +
  coord_cartesian(ylim = c(0,1))

## Create average AUROC plot
### Calculate median ROC
require(purrr)
mean_roc <- function(data, cutoffs = seq(from = -5, to = 5, by = 0.5)) {
  map_df(cutoffs, function(cp) {
    out <- cutpointr(data = data, x = PredictionValues, class = RealClass,
                     subgroup = Sample, method = oc_manual, cutpoint = cp,
                     pos_class = "pathogenic", direction = ">=")
    data.frame(cutoff = cp, 
               sensitivity = mean(out$sensitivity),
               specificity = mean(out$specificity))
  })
}

predictions_50_samples <- map2_dfr(AUC_crossval$MAFvsDMR, names(AUC_crossval$MAFvsDMR), function(x,y){
  data.frame(Sample = rep(y, nrow(x)),
             PredictionValues = x[,3],
             RealClass = x[,1])
})

mean_roc2 <- function(data, cutoffs = seq(from = 0, to = 4, by = 0.01)) {
  map_df(cutoffs, function(cp) {
    out <- cutpointr(data = data, x = PredictionValues, class = RealClass,
                     subgroup = Sample, method = oc_manual, cutpoint = cp,
                     pos_class = "HNC", direction = ">=")
    data.frame(cutoff = cp, 
               sensitivity = mean(out$sensitivity),
               specificity = mean(out$specificity))
  })
}

mr2 <- mean_roc2(predictions_50_samples)

meanRocPlot <- cutpointr(data = predictions_50_samples, 
                         x = PredictionValues, class = RealClass, subgroup = Sample,
                         pos_class = "HNC", direction = ">=") %>% 
  plot_roc(display_cutpoint = FALSE, size = 0.52) + theme_JB() + 
  theme(legend.position="none", plot.title = element_blank(),
        plot.subtitle = element_blank()) + 
  scale_color_manual(values = rep(rgb(0,0,0,1/20), nrow(predictions_50_samples))) + 
  geom_line(data = mr2, mapping = aes(x = 1 - specificity, y = sensitivity), 
            color = "darkblue", size = 1)

ExtendedFigure8 <- egg::ggarrange(plots = list(meanRocPlot + theme(panel.background = element_blank()), 
                                  crossval_plot + theme(panel.background = element_blank()),
                                  as_ggplot(ggsurv_BL_stage) + theme(panel.background = element_blank())),
                                  ncol = 2, widths = c(0.7, 0.25))

