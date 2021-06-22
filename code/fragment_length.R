# cfMeDIP-seq and CAPP-Seq fragment length analysis #
#####################################################

# Set directories
scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"

# Load libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Repitools)
source(paste0(scripts, "/ggplot_themeJB.R"))

# Create function for extracting fragment length info from hyperDMRs
fun_fraglen <- function(res_list, file_dir, MAFvsDMR_obj){
  # List files
  sample_dir <- list.files(file_dir)
  
  # For each file - load - subset by hyperDMRs - remove
  all_hyperDMR_frag <- list()
  
  for(i in 1:length(sample_dir)){
    # Create path to load
    sample_full <- paste0(file_dir, "/", sample_dir[i])
    R_obj <- gsub(".RData", "", sample_dir[i])
    
    # Load
    load(sample_full)
    
    # Create tmp file and remove original
    tmp <- get(R_obj)
    rm(list=R_obj)
    
    # Subset tmp file by hyperDMRs of interest
    tmp_nametoID <- sub(":", ".", names(tmp))
    tmp_nametoID <- sub("-", ".", tmp_nametoID)
    hyperDMR_ID <- rownames(res_list[['res_padj_hyperDMR']])
    hyperDMR_index <- tmp_nametoID %in% hyperDMR_ID
    
    tmp_hyperDMR <- tmp[hyperDMR_index]
    
    # Add to list
    
    all_hyperDMR_frag[[R_obj]] <- tmp_hyperDMR
    
  }
  
  # Divide into HNC and norm samples
  hnc_index <- grepl("BL", names(all_hyperDMR_frag))
  norm_index <- grepl("Norm", names(all_hyperDMR_frag))
  
  norm_list  <- all_hyperDMR_frag[norm_index]
  hnc_list <- all_hyperDMR_frag[hnc_index]
  
  # Get fragment size distribution across all windows for each sample
  ## Perform on grouped norms
  n <- 1
  grouped_norm_list <- norm_list[[n]]
  while(n <= length(norm_list)){
    grouped_norm_list <- Map(c, grouped_norm_list, norm_list[[n]])
    n <- n + 1
  }
  
  norm_allfrag <- as.numeric(unlist(grouped_norm_list))
  
  ## Perform on each hnc
  hnc_allfrag <- lapply(hnc_list, function(x) as.numeric(unlist(x)))
  
  # Rank hnc_allfrag based on hyperDMR signal
  MAFvsDMR_grepl <- paste0(substr(names(hnc_allfrag), 8, 11), collapse = "|")
  MAFvsDMR_index <- grepl(MAFvsDMR_grepl, rownames(MAFvsDMR_obj))
  
  MAFvsDMR_obj <- MAFvsDMR_obj[MAFvsDMR_index,]
  MAFvsDMR_obj <- MAFvsDMR_obj[order(MAFvsDMR_obj$meanRPKM, decreasing = TRUE),]
  
  hnc_dmr_rank <- c()
  for(i in 1:nrow(MAFvsDMR_obj)){
    ident <- substr(rownames(MAFvsDMR_obj)[i], 3, 6)
    whichtmp <- which(substr(names(hnc_allfrag), 8, 11) %in% ident)
    hnc_dmr_rank <- c(hnc_dmr_rank, whichtmp)
  }
  
  hnc_allfrag <- hnc_allfrag[hnc_dmr_rank]
  
  # Create MedianFragment size vs. RPKM
  medfrag_vs_hyperDMR <- data.frame(medfrag = sapply(hnc_allfrag, median),
                                    hyperDMR = MAFvsDMR_obj$meanRPKM)
  
  # Collate to list
  hyperDMR_frag_output <- list(norm_allfrag, hnc_allfrag, medfrag_vs_hyperDMR)
  names(hyperDMR_frag_output) <- c("norm_allfrag", "hnc_allfrag", "medfrag_vs_hyprDMR")
  
  # Return object
  return(hyperDMR_frag_output)
  
}

# Create data.frame with ctDNA abundance as calculated by CAPP-Seq as well as
# methylation (mean RPKM) as calculated by cfMeDIP-seq
## Convert matrix counts to RPKM values of cfMeDIP-seq samples
load(paste0(data_dir, "/cfDNA_CombinedCounts_AllWindows.RData"))

## Save total reads for division
RPKM_denom <- colSums(cfDNA_CombinedCounts_AllWindows) / 1000000

## Subset to PBLdepleted windows
load(paste0(out_dir, "/PBLdepleted_IDs.RData"))

cfDNA_CombinedRPKMs_PBLdepleted <- 
  cfDNA_CombinedCounts_AllWindows[rownames(cfDNA_CombinedCounts_AllWindows) %in% 
                                    PBLdepleted_IDs,]

cfDNA_CombinedRPKMs_PBLdepleted <- t(t(cfDNA_CombinedRPKMs_PBLdepleted / 0.3) /
                                       RPKM_denom)

rm(cfDNA_CombinedCounts_AllWindows)

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

### Save ctDNA vs RPKM object
save(mclust_wnorm_df, file = paste0(out_dir, "/mclust_wnorm_df.RData"))

# Run function on directory containing output files containing fragment length
# info across all 300-bp windows from cfMeDIP-seq libraries
load(paste0(out_dir, "/DESeq2_cfDNA_DMRs.RData"))

## NOTE: fragment length data was generated via bash/SAMtools and deposited into
## the local directory ~/data/insert_size
frag_dir <- paste0(data_dir, "/insert_size")

hyperDMR_frag <- fun_fraglen(res_list = res_list, 
                             file_dir = frag_dir,
                             MAFvsDMR_obj = mclust_wnorm_df)

# Subset below 220 for 1st peak only
hyperDMR_frag_1stPeak <- lapply(hyperDMR_frag[1:2], FUN = function(x){
  lapply(x, FUN = function(x){
    x[x < 220]
  })
})

# Plot fragment length distribution across samples
hyperDMR_medfrag_c <- c(list(hyperDMR_frag_1stPeak[['norm_allfrag']]), 
                        hyperDMR_frag_1stPeak[['hnc_allfrag']][30:1])

names(hyperDMR_medfrag_c)[1] <- "Norm"
names(hyperDMR_medfrag_c)[2:31] <- 
  paste0("HN", substr(names(hyperDMR_medfrag_c)[2:31], 8, 11))

hyperDMR_medfrag <- melt(hyperDMR_medfrag_c)

hyperDMR_medfrag$L1 <- as.factor(hyperDMR_medfrag$L1)
hyperDMR_medfrag$L1 <- factor(hyperDMR_medfrag$L1, levels = names(hyperDMR_medfrag_c))

hyperDMR_medfrag$condition <- 
  as.factor(ifelse(grepl("Norm", hyperDMR_medfrag$L1), "Grouped norm", "HNC"))
hyperDMR_medfrag$condition <- 
  factor(hyperDMR_medfrag$condition, levels = c("HNC", "Grouped norm"))

## ggplot
hyperDMR_medfrag_plot <- hyperDMR_medfrag %>%
  ggplot(aes(x = L1, y = value, col = condition)) + 
  geom_jitter(position = position_jitter(0.3),
              shape = 16, size = 0.25) + 
  geom_boxplot(outlier.shape = NA, coef = 0,
               col = "black", fill = NA,
               lwd = 0.13, fatten = 6) + 
  geom_hline(yintercept = median(hyperDMR_medfrag[hyperDMR_medfrag$L1 == "Norm",]$value),
             colour = "blue", linetype = "longdash", size = 0.13) + 
  theme_JB() + 
  scale_color_manual(values = c(alpha(palette_JB[1], 0.05), alpha(palette_JB[2], 0.01))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5,5,15,5), "mm")) +
  ylim(c(100, 220)) + 
  labs(y = "Fragment length", x = "Samples") + 
  annotation_custom(legendGrob(labels = c("HNSCC", "Healthy donors"),
                               nrow = 1, ncol = 2,
                               pch = 16, gp = gpar(col = c(palette_JB[1],palette_JB[2]),
                                                   cex = 0.75)),
                    xmin = 0, xmax = 34, ymin = 65, ymax = 65) +
  coord_cartesian(xlim = c(0,34), clip = "off")

# Repeat for fragments overlapping with identified SNVs
mut_insertsize <- read.csv(paste0(data_dir, 
                                  "/fragment_length_wildtype_mutant_reads.csv"), 
                           stringsAsFactors = TRUE)

# Generate ID of Identified Mutaions
load(paste0(out_dir, "/cfDNA_DepthFiltered_coding_SNVs.RData"))
cfDNA_DepthFiltered_coding_SNVs$ID <- paste(cfDNA_DepthFiltered_coding_SNVs$upn, 
                                   cfDNA_DepthFiltered_coding_SNVs$chr,
                                   cfDNA_DepthFiltered_coding_SNVs$pos,
                                   cfDNA_DepthFiltered_coding_SNVs$ref,
                                   cfDNA_DepthFiltered_coding_SNVs$alt, sep=".")

# Select muations present in insert-size data 
cfDNA_DepthFiltered_coding_SNVs_wInset <- mut_insertsize[mut_insertsize$ID %in% 
                                                  cfDNA_DepthFiltered_coding_SNVs$ID,]
cfDNA_DepthFiltered_coding_SNVs_wInset$upn <- factor(cfDNA_DepthFiltered_coding_SNVs_wInset$upn)

## Create list version as well
### Mutant
cfDNA_DepthFiltered_coding_SNVs_wMutInset_list <- list()
for(i in 1:length(levels(cfDNA_DepthFiltered_coding_SNVs_wInset$upn))){
  tmp <- cfDNA_DepthFiltered_coding_SNVs_wInset[
    levels(cfDNA_DepthFiltered_coding_SNVs_wInset$upn)[i] == 
      cfDNA_DepthFiltered_coding_SNVs_wInset$upn,]
  
  cfDNA_DepthFiltered_coding_SNVs_wMutInset_list[[i]] <- 
    tmp$medianFragmentLengthMutant
  names(cfDNA_DepthFiltered_coding_SNVs_wMutInset_list)[i] <- 
    levels(cfDNA_DepthFiltered_coding_SNVs_wInset$upn)[i]
}

### Wildtype
cfDNA_DepthFiltered_coding_SNVs_wWTInset_list <- list()
for(i in 1:length(levels(cfDNA_DepthFiltered_coding_SNVs_wInset$upn))){
  tmp <- cfDNA_DepthFiltered_coding_SNVs_wInset[
    levels(cfDNA_DepthFiltered_coding_SNVs_wInset$upn)[i] == 
      cfDNA_DepthFiltered_coding_SNVs_wInset$upn,]
  
  cfDNA_DepthFiltered_coding_SNVs_wWTInset_list[[i]] <- 
    tmp$medianFragmentLengthWildtype
  
  names(cfDNA_DepthFiltered_coding_SNVs_wWTInset_list)[i] <- 
    levels(cfDNA_DepthFiltered_coding_SNVs_wInset$upn)[i]
}

## Create data.frame to plot results
MAF_index <- c()
for(i in 1:length(cfDNA_DepthFiltered_coding_SNVs_wMutInset_list)){
  tmp <- which(substr(names(cfDNA_DepthFiltered_coding_SNVs_wMutInset_list)[i], 
                      3, 6) == 
                 substr(rownames(mclust_wnorm_df[
                   order(mclust_wnorm_df$ctDNA, 
                         decreasing = TRUE),
                   ]), 3, 6))
  
  MAF_index <- c(MAF_index, tmp)
}

WT_list <- melt(cfDNA_DepthFiltered_coding_SNVs_wWTInset_list)
WT_list$condition <- rep("WT", nrow(WT_list))

Mut_list <- melt(cfDNA_DepthFiltered_coding_SNVs_wMutInset_list)
Mut_list$condition <- rep("Mut", nrow(Mut_list))

SNV_medfrag <- rbind(WT_list, Mut_list)

SNV_medfrag$L1 <- as.factor(SNV_medfrag$L1)
SNV_medfrag$L1 <- factor(SNV_medfrag$L1, 
                         levels = levels(SNV_medfrag$L1)[MAF_index])

set.seed(42) # for jitter

## Create ggplot
SNV_medfrag_plot <- SNV_medfrag %>%
  ggplot(aes(x = L1, y = value, fill = condition, color = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2),
              size = 1.5, aes(col = condition)) + 
  geom_boxplot(position = "identity",
               outlier.shape = NA,
               coef = 0, fill = NA, 
               lwd = 0.13, fatten = 6) + 
  scale_fill_manual("Allele:", labels = c("SNV", "Wild-type"),
                    values = alpha(palette_JB)) + 
  scale_color_manual("Allele:", labels = c("SNV", "Wild-type"),
                     values = alpha(palette_JB)) + 
  theme_JB() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  labs(y = "Median fragment length", x = "Samples")

# Evaluate Median fragment length vs. MAF by CAPP-Seq
SNV_medfrag_median <- SNV_medfrag %>%
  subset(condition == "Mut") %>%
  group_by(L1) %>%
  summarize(medfrag = median(value))

MAF_tmp <- mclust_wnorm_df[match(SNV_medfrag_median$L1, 
                                 rownames(mclust_wnorm_df)),]$ctDNA

SNV_medfrag_median$MAF <- MAF_tmp

plot_pval <- cor.test(x = SNV_medfrag_median$MAF,
                      y = SNV_medfrag_median$medfrag)$p.value
plot_pval <- round(plot_pval, 2)

MAF_vs_SNVmedfrag <- SNV_medfrag_median %>%
  ggplot(aes(x = MAF, y = medfrag)) +
  geom_smooth(method = "lm", color = palette_JB[1], size = 0.26,
              fill = "lightgrey") + 
  geom_point(size = 1.5, col = palette_JB[1]) +
  labs(x = "Mutant allele fraction * 100 (%)",
       y = "Median fragment length - all SNVs") + 
  theme_JB() + 
  annotate("text", x = 4.5, y = 113,
           label = paste0("p = ", plot_pval),
           size = (5/14) * 9)

# Evaluate fragment length vs. RPKM of hyper-DMRs
meanRPKM_vs_medfrag <- mclust_wnorm_df[mclust_wnorm_df$ctDNA > 0,]
med_ind <- match(rownames(meanRPKM_vs_medfrag), names(hyperDMR_medfrag_c))
med_c <- sapply(hyperDMR_medfrag_c[med_ind], median)
meanRPKM_vs_medfrag$medfrag <- med_c

plot_cor <- cor.test(x = meanRPKM_vs_medfrag$meanRPKM,
                     y = meanRPKM_vs_medfrag$medfrag)

## Plot
meanRPKM_vs_medfrag_plot <- meanRPKM_vs_medfrag %>% 
  ggplot(aes(x = meanRPKM, y = medfrag)) +
  geom_smooth(method = lm, fill = "lightgrey", col = "red",
              size = 0.13, fullrange = TRUE) + 
  geom_point(shape = 16, colour = palette_JB[1]) + 
  annotate("text", x = 0.01, y = 147,
           label = paste0("R = ",round(plot_cor$estimate, 2),
                          "\np = ",round(plot_cor$p.value, 3)),
           hjust = 0, size = (5/14) * 9) + 
  labs(y = "Median fragment length",
       x = "Methylation (mean RPKM)") + 
  theme_JB()

# Compare fragment length between CAPP-Seq and cfMeDIP-seq profiles
SNV_medfrag_test <- SNV_medfrag
SNV_medfrag_test$value <- as.numeric(SNV_medfrag_test$value)
SNV_medfrag_test$L1 <- as.character(SNV_medfrag_test$L1)

CAPP_vs_MeDIP_medfrag <- 
  SNV_medfrag_test[SNV_medfrag_test$condition == "Mut",] %>%
  group_by(L1) %>%
  summarize(mut = median(value))

median_hyperDMR_frag <- sapply(hyperDMR_frag_1stPeak$hnc_allfrag, median)
names(median_hyperDMR_frag) <- 
  paste0("HN", substr(names(median_hyperDMR_frag), 8, 11))
median_hyperDMR_frag <- median_hyperDMR_frag[names(median_hyperDMR_frag) %in%
                                               CAPP_vs_MeDIP_medfrag$L1]

CAPP_vs_MeDIP_medfrag$meth <- 
  median_hyperDMR_frag[order(names(median_hyperDMR_frag))]

## Subset to patients above median methylation
rpkm_ordered <- meanRPKM_vs_medfrag[order(meanRPKM_vs_medfrag$meanRPKM),]

CAPP_vs_MeDIP_filtered <- CAPP_vs_MeDIP_medfrag[CAPP_vs_MeDIP_medfrag$L1 %in%
                                 rownames(rpkm_ordered)[11:20],]
filtered_cor <- cor.test(x = CAPP_vs_MeDIP_filtered$meth,
                      y = CAPP_vs_MeDIP_filtered$mut)

## Create ggplot
CAPP_vs_MeDIP_medfrag_plot <- CAPP_vs_MeDIP_filtered %>%
  ggplot(aes(x = meth, y = mut)) + 
  geom_smooth(method = lm, se = TRUE,
              col = "red", size = 0.13, fill = "lightgrey",
              fullrange = TRUE) + 
  geom_point(shape = 16, colour = palette_JB[1],
             size = 1.5) + 
  annotate("text", x = 142, y = 168,
           label = paste0("R = ", round(filtered_cor$estimate, 2),
                          "\np = ", round(filtered_cor$p.value, 4)),
           hjust = 0, size = (5/14) * 9) + 
  theme_JB() + 
  labs(y = "SNV median fragment length",
       x = "HyperDMR median fragment length") +
  scale_x_continuous(breaks = seq(130, 170, by = 5))


# Subset analysis to short (100-150) and large (151-220) fragments
## Load output file of subsetted BAM files processed by MEDIPS R package
load(paste0(out_dir, "/PBLdepleted_IDs.RData"))
load(paste0(out_dir, "/DESeq2_cfDNA_DMRs.RData"))

# NOTE: cfMeDIP-seq profiles were subsetted between 100-150 and 151-220 bp
## via pysam
load(paste0(data_dir, "/cfDNA_Frag100_150_Counts.RData"))
cfDNA_Frag100_150_PBLdepleted_Counts <- 
  cfDNA_Frag100_150_Counts[rownames(cfDNA_Frag100_150_Counts) %in%
                             PBLdepleted_IDs,]
rm(cfDNA_Frag100_150_Counts)

load(paste0(data_dir, "/cfDNA_Frag151_220_Counts.RData"))
cfDNA_Frag151_220_PBLdepleted_Counts <- 
  cfDNA_Frag151_220_Counts[rownames(cfDNA_Frag151_220_Counts) %in%
                             PBLdepleted_IDs,]
rm(cfDNA_Frag151_220_Counts)

load(paste0(data_dir, "/cfDNA_Frag100_220_Counts.RData"))
cfDNA_Frag100_220_PBLdepleted_Counts <- 
  cfDNA_Frag100_220_Counts[rownames(cfDNA_Frag100_220_Counts) %in%
                             PBLdepleted_IDs,]
rm(cfDNA_Frag100_220_Counts)

## Generate a null distribution of expected counts at random windows w/
## identical number and CpG content as hyperDMRs
### Create GenomicRanges object for genome-wide bins and calculate
### CpG density
gw300bp <- genomeBlocks(seqlengths(Hsapiens), width = 300)
gw300bp$CpG <- cpgDensityCalc(gw300bp, Hsapiens)
gw300bp$ID <- paste(seqnames(gw300bp),
                    start(gw300bp), 
                    end(gw300bp),
                    sep = ".")
gw300bp_PBLdepleted <- gw300bp[gw300bp$ID %in% PBLdepleted_IDs]

hyperDMR_CpG <- 
  gw300bp_PBLdepleted[gw300bp_PBLdepleted$ID %in% 
                        rownames(res_list$res_padj_hyperDMR)]$CpG

### Calculate CpG density distribution of hyperDMRs
CpG_table <- table(hyperDMR_CpG)

### Create function to take random PBLdepleted windows of with identical 
### numbers and CpG distributions as hyperDMRs
null_exp <- function(gr, CpG_table, counts, n_replicates){
  tmp_df <- as.data.frame(matrix(nrow = 0, ncol = ncol(counts)))
  for(r in 1:n_replicates){
    tmp_ind <- c()
    for(i in 1:length(CpG_table)){
      tmp_samp <- sample(which(gr$CpG == names(CpG_table)[i]),
                         CpG_table[i])
      tmp_ID <- gr[tmp_samp]$ID
      tmp_ind <- c(tmp_ind, tmp_ID)
    }
    tmp_obs <- counts[rownames(counts) %in% tmp_ind,]
    tmp_obs_colMeans <- colMeans(tmp_obs)
    tmp_df <- rbind(tmp_df, tmp_obs_colMeans)
  }
  colnames(tmp_df) <- colnames(counts)
  return(tmp_df)
}

### Estimate distribution of expected counts across the null distiribution
#### Short
set.seed(42)
exp_short_df <- null_exp(gr = gw300bp_PBLdepleted,
                         CpG_table = CpG_table,
                         counts = cfDNA_Frag100_150_PBLdepleted_Counts,
                         n_replicates = 30)

##### Determine observed value by taking mean counts across hyperDMRs
obs_short <- colMeans(cfDNA_Frag100_150_PBLdepleted_Counts[
  rownames(cfDNA_Frag100_150_PBLdepleted_Counts) %in%
    rownames(res_list$res_padj_hyperDMR),
  ]) 

##### Calculate significance of enrichment
enrc_short_sig <- c()
for(i in 1:ncol(exp_short_df)){
  exp_distr <- exp_short_df[,i]
  obs_tmp <- obs_short[i]
  
  b <- sum(exp_distr > obs_tmp)
  pval <- (b + 1) / length(exp_distr)
  
  enrc_short_sig <- c(enrc_short_sig, pval)
  names(enrc_short_sig)[i] <- substr(colnames(exp_short_df)[i], 1, 6)
}

##### Calculate percent enrichment
enrc_short_percent <- obs_short / colMeans(exp_short_df)

#### Long
set.seed(42)
exp_long_df <- null_exp(gr = gw300bp_PBLdepleted,
                         CpG_table = CpG_table,
                         counts = cfDNA_Frag151_220_PBLdepleted_Counts,
                         n_replicates = 30)

##### Determine observed value by taking mean counts across hyperDMRs
obs_long <- colMeans(cfDNA_Frag151_220_PBLdepleted_Counts[
  rownames(cfDNA_Frag151_220_PBLdepleted_Counts) %in%
    rownames(res_list$res_padj_hyperDMR),
  ]) 

##### Calculate significance of enrichment
enrc_long_sig <- c()
for(i in 1:ncol(exp_long_df)){
  exp_distr <- exp_long_df[,i]
  obs_tmp <- obs_long[i]
  
  b <- sum(exp_distr > obs_tmp)
  pval <- (b + 1) / length(exp_distr)
  
  enrc_long_sig <- c(enrc_long_sig, pval)
  names(enrc_long_sig)[i] <- substr(colnames(exp_long_df)[i], 1, 6)
}

##### Calculate percent enrichment
enrc_long_percent <- obs_long / colMeans(exp_long_df)

#### Combined
set.seed(42)
exp_combined_df <- null_exp(gr = gw300bp_PBLdepleted,
                        CpG_table = CpG_table,
                        counts = cfDNA_Frag100_220_PBLdepleted_Counts,
                        n_replicates = 30)

##### Determine observed value by taking mean counts across hyperDMRs
obs_combined <- colMeans(cfDNA_Frag100_220_PBLdepleted_Counts[
  rownames(cfDNA_Frag100_220_PBLdepleted_Counts) %in%
    rownames(res_list$res_padj_hyperDMR),
  ]) 

##### Calculate significance of enrichment
enrc_combined_sig <- c()
for(i in 1:ncol(exp_combined_df)){
  exp_distr <- exp_combined_df[,i]
  obs_tmp <- obs_combined[i]
  
  b <- sum(exp_distr > obs_tmp)
  pval <- (b + 1) / length(exp_distr)
  
  enrc_combined_sig <- c(enrc_combined_sig, pval)
  names(enrc_combined_sig)[i] <- substr(colnames(exp_combined_df)[i], 1, 6)
}

##### Calculate percent enrichment
enrc_combined_percent <- obs_combined / colMeans(exp_combined_df)

#### Calculate foldchange in enrichment between short and long/combined
short_vs_long_enrc_FC <- enrc_short_percent / enrc_long_percent
short_vs_long_enrc_FC <- short_vs_long_enrc_FC[order(short_vs_long_enrc_FC,
                                                     decreasing = TRUE)]

short_vs_combined_enrc_FC <- enrc_short_percent / enrc_combined_percent
short_vs_combined_enrc_FC <- short_vs_combined_enrc_FC[
  order(short_vs_combined_enrc_FC,
        decreasing = TRUE)
]

### Get range of enrichment for patients with detectable ctDNA
ctDNApos_HN <- rownames(mclust_wnorm_df[mclust_wnorm_df$ctDNA > 0,])
ctDNApos_HN <- substr(ctDNApos_HN, 3, 6)

short_enrc_ctDNApos <- short_vs_combined_enrc_FC[
  grepl(paste0(ctDNApos_HN, collapse = "|"), names(short_vs_combined_enrc_FC))
]

## Create ggplot of output
short_vs_combined_df <- 
  data.frame(condition = ifelse(grepl("Norm",names(short_vs_combined_enrc_FC)),
                                "Norm", "HNSCC"),
             enr_FC = (short_vs_combined_enrc_FC - 1) * 100)

df_length <- 1:nrow(short_vs_combined_df)

short_enrc_plot <- short_vs_combined_df %>%
  ggplot(aes(x = 1:50,
             y = enr_FC,
             fill = condition)) + 
  geom_bar(stat = "identity", position = "dodge", 
           size = 0.13, width = 0.75,
           col = NA) + 
  labs(y = "100-150/100-220 bp hyperDMR enrichment (%)",
       x = "samples") + 
  lims(y = c(0,2)) + 
  theme_JB() + 
  scale_y_continuous(expand = c(0,0), limits = c(-25, 100)) + 
  scale_x_discrete(expand = c(0,1)) +
  scale_fill_manual(values = palette_JB)

