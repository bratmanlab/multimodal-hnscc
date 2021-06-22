# Identification of PBL-depleted windows source code #
######################################################

# Set directories
scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"

# Load libraries
library(GenomicRanges)
library(Repitools)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggridges)
library(RColorBrewer)
library(stats)
library(DESeq2)
library(ggplot2)
source(paste0(scripts, "/ggplot_themeJB.R"))

# Load RData object containing MeDIP-seq PBL RPKMs summarized across 300-bp 
# windows. Columns represent samples collected from HNSCC patients at 
# diagnosis as well as healthy donors. Rows represent location of 300-bp
# window denoted as chromosome.start.end
load(paste0(data_dir, "/PBLs_CombinedRPKMs_Size300_NonXY.RData"))

# Subset data.frame to samples from healthy donors only
Norm_PBL_NonXY_RPKMs <- PBLs_CombinedRPKMs_Size300_NonXY[
  , grepl("Norm", colnames(PBLs_CombinedRPKMs_Size300_NonXY))
  ]
rm(PBLs_CombinedRPKMs_Size300_NonXY)

# Summarize median RPKM of each 300-bin across healthy donor samples
Norm_PBL_NonXY_rowMedianRPKMs <- Biobase::rowMedians(Norm_PBL_NonXY_RPKMs)
names(Norm_PBL_NonXY_rowMedianRPKMs) <- rownames(Norm_PBL_NonXY_RPKMs)

# Load RData object containing MeDIP-seq FaDu RPKMs summarized across 300-bp
# windows.
load(paste0(data_dir, "/FaDu_NonXY_RPKMs.RData"))

# Calculate CpG density for each 300-bp bin
## Function to convert row IDs of previous data.frames
IDtoGrange <- function(ID){
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  Names <- unlist(strsplit(ID, "[.]"))
  Chr <- Names[seq(1,length(Names), by=3)]
  Start <- as.numeric(Names[seq(2,length(Names), by=3)])
  End <- as.numeric(Names[seq(3,length(Names), by=3)])
  
  GRanges(seqnames=Rle(Chr),
          ranges=IRanges(start=Start, end=End))
}

NonXY_300bp_gr <- IDtoGrange(names(Norm_PBL_NonXY_rowMedianRPKMs))
NonXY_300bp_gr$CpG <- cpgDensityCalc(NonXY_300bp_gr, Hsapiens)
NonXY_300bp_gr$ID <- paste(seqnames(NonXY_300bp_gr),
                           start(NonXY_300bp_gr),
                           end(NonXY_300bp_gr),
                           sep = ".")

# Create data.frame with healthy donor and FaDu RPKMs and CpG density
gw300bp_IDswCpGs_NonXY <- data.frame(ID = NonXY_300bp_gr$ID,
                                     CpG_count = NonXY_300bp_gr$CpG)
gw300bp_IDswCpGs_NonXY$FaDu_RPKMs <- as.numeric(FaDu_NonXY_RPKMs)
gw300bp_IDswCpGs_NonXY$Norm_rowMedianRPKMs <- 
  as.numeric(Norm_PBL_NonXY_rowMedianRPKMs)

# Summarize RPKMs >= x CpGs for cfDNA and FaDu
CpG_list <- unique(gw300bp_IDswCpGs_NonXY$CpG_count)
RPKM_get <- c()
for(i in CpG_list){
  df_tmp <- gw300bp_IDswCpGs_NonXY[gw300bp_IDswCpGs_NonXY$CpG_count >= i,]
  tmp <- median(df_tmp$Norm_rowMedianRPKMs)
  RPKM_get <- c(RPKM_get, tmp)
  
}

FaDu_RPKM_get <- c()
for(i in CpG_list){
  df_tmp <- gw300bp_IDswCpGs_NonXY[gw300bp_IDswCpGs_NonXY$CpG_count >= i,]
  tmp <- median(df_tmp$FaDu_RPKMs)
  FaDu_RPKM_get <- c(FaDu_RPKM_get, tmp)
  
}

names(RPKM_get) <- CpG_list
RPKM_get <- RPKM_get[order(as.numeric(names(RPKM_get)))]

names(FaDu_RPKM_get) <- CpG_list
FaDu_RPKM_get <- FaDu_RPKM_get[order(as.numeric(names(FaDu_RPKM_get)))]

# Collate >= x CpGs for cfDNA and FaDu for ggplot
MedianRPKM <- data.frame(condition = c(rep("FaDu", length(FaDu_RPKM_get)), 
                                       rep("Healthy donor", length(RPKM_get))),
                         MedianRPKM = c(FaDu_RPKM_get, RPKM_get),
                         nCpG = as.numeric(c(names(FaDu_RPKM_get), names(RPKM_get))))

MedianRPKM <- MedianRPKM[MedianRPKM$nCpG <= 40,]

# Create ggplot
FaDu_vs_Norm_nCpGMedianRPKM <- MedianRPKM %>%
  ggplot(aes(x = nCpG, y = MedianRPKM, color = condition, group = condition)) + 
  geom_line(size = 0.13) + 
  geom_point(size = 1.5) + 
  theme_JB() + 
  theme(panel.background = element_blank()) + 
  scale_color_manual(values = c(palette_JB[1], palette_JB[2])) + 
  scale_x_continuous(breaks = seq(0, length(unique(as.character(MedianRPKM$nCpG))), by = 4)) + 
  labs(y = "Median methylation\n(RPKM)", x = "# of CpGs >= x")

# Create table w/ number of 300-bp windows per >= n CpGs
gw300bp_table <- table(gw300bp_IDswCpGs_NonXY$CpG_count)
getnCpGs_n <- sum(gw300bp_table) / 1000000
names(getnCpGs_n) <- 0
for(i in 1:40){
  tmp <- sum(gw300bp_table[-c(1:i)]) / 1000000
  getnCpGs_n <- c(getnCpGs_n, tmp)
  names(getnCpGs_n)[i + 1] <- i
}
gw300bp_df <- data.frame(getnCpGs = as.numeric(names(getnCpGs_n)),
                         n = as.numeric(getnCpGs_n))
gw300bp_df <- gw300bp_df[1:41,]

nWindows_vs_nCpG <- gw300bp_df %>%
  ggplot(aes(x = getnCpGs, y = n))  + 
  geom_line(size = 0.13) + 
  geom_point(size = 1.5) + 
  theme_JB() + 
  theme(panel.background = element_blank(), 
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  scale_color_manual(values = "black") + 
  scale_x_continuous(breaks = seq(0, 40, by = 4)) + 
  labs(y = "No. 300-bp windows (M)", x = ">= n CpGs")

FaDu_vs_Norm_nCpGMedianRPKM <- egg::ggarrange(nWindows_vs_nCpG, FaDu_vs_Norm_nCpGMedianRPKM,
                       nrow = 2, ncol = 1, draw = FALSE)



# Select 300-bp windows w/ >= 8 CpGs
CpGgt8_IDs <- 
  gw300bp_IDswCpGs_NonXY[gw300bp_IDswCpGs_NonXY$CpG_count >= 8, ]$ID
## Save for later use
save(CpGgt8_IDs, file = paste0(out_dir, "/CpGgt8_IDs.RData"))

# Reload RData object with MeDIP-seq profiles of HNSCC and healthy donor PBL
# RPKMs
rm(Norm_PBL_NonXY_RPKMs) # remove due to memory limitations

load(paste0(data_dir, "/PBLs_CombinedRPKMs_Size300_NonXY.RData"))

PBLs_CombinedRPKMs_Size300_CpGgt8 <- 
  PBLs_CombinedRPKMs_Size300_NonXY[
    rownames(PBLs_CombinedRPKMs_Size300_NonXY) %in% CpGgt8_IDs,
    ]
rm(PBLs_CombinedRPKMs_Size300_NonXY)

# Load RData object with cfMeDIP-seq profiles of HNSCC and healthy donor cfDNA
# RPKMs and subset to 300-bp bins with >= 8 CpGs
load(paste0(data_dir, "/cfDNA_CombinedRPKMs_Size300_NonXY.RData"))

cfDNA_CombinedRPKMs_Size300_CpGgt8 <-
  cfDNA_CombinedRPKMs_Size300_NonXY[
    rownames(cfDNA_CombinedRPKMs_Size300_NonXY) %in% CpGgt8_IDs,
  ]

rm(cfDNA_CombinedRPKMs_Size300_NonXY)

# Calculate correlation of cfDNA w/ FaDu, unmatched PBLs, and matched PBLs
## custom function to incoporate Cook's Distance cutoff between data.frames
cooks_cor <- function(cfDNA_mat, PBL_mat, cutoff){
  
  tmp_mat <- matrix(nrow = ncol(cfDNA_mat), ncol = ncol(PBL_mat))
  colnames(tmp_mat) <- colnames(cfDNA_mat)
  rownames(tmp_mat) <- colnames(PBL_mat)
  
  # for each sample, filter regions above cooks cutoff
  for(i in 1:ncol(cfDNA_mat)){
    tmp_cor_list <- c()
    for(n in 1:ncol(cfDNA_mat)){
      df <- data.frame(cfDNA=cfDNA_mat[,i], PBL=PBL_mat[,n])
      tmp_lm <- lm(PBL ~ cfDNA, data=df)
      tmp_cooks <- cooks.distance(tmp_lm)
      
      ## remove windows with cooks cutoff > 3 * mean
      filtered_df <- df[rownames(df) %in% names(tmp_cooks[tmp_cooks < cutoff]),]
      
      ## calculate cor
      filtered_cor <- cor(filtered_df$cfDNA, filtered_df$PBL)
      
      ## add to list
      tmp_cor_list <- c(tmp_cor_list, filtered_cor)
    }
    
    tmp_mat[,i] <- tmp_cor_list
    
  }
  
  return(tmp_mat)
  
}

## cfDNA vs. all PBLs
## NOTE: order of samples for both cfDNA and PBL data.frames must be in the same
cooks_cutoff <- 1 # maximum Cook's Distance threshold
CpGgt8_cfDNA_PBL_cooks_cor <- 
  cooks_cor(cfDNA_mat =  cfDNA_CombinedRPKMs_Size300_CpGgt8,
            PBL_mat = PBLs_CombinedRPKMs_Size300_CpGgt8,
            cutoff = cooks_cutoff)

colnames(CpGgt8_cfDNA_PBL_cooks_cor) <- paste0("pb_",1:50)
rownames(CpGgt8_cfDNA_PBL_cooks_cor) <- paste0("cf_",1:50)

CpGgt8_cfDNA_PBL_cooks_cor_matched <- diag(CpGgt8_cfDNA_PBL_cooks_cor) 
CpGgt8_cfDNA_PBL_cooks_cor_unmatched <- c(CpGgt8_cfDNA_PBL_cooks_cor)
CpGgt8_cfDNA_PBL_cooks_cor_unmatched <- CpGgt8_cfDNA_PBL_cooks_cor_unmatched[
  !(CpGgt8_cfDNA_PBL_cooks_cor_unmatched %in% 
      CpGgt8_cfDNA_PBL_cooks_cor_matched)
  ]

## cfDNA vs. FaDu
FaDu_CpGgt8 <- FaDu_NonXY_RPKMs[names(FaDu_NonXY_RPKMs) %in% CpGgt8_IDs]

tmp_mat <- c()

### for each sample, filter regions above cook's cutoff
for(i in 1:ncol(cfDNA_CombinedRPKMs_Size300_CpGgt8)){
  
  df <- data.frame(cfDNA = cfDNA_CombinedRPKMs_Size300_CpGgt8[,i], 
                   PBL = FaDu_CpGgt8)
  tmp_lm <- lm(PBL ~ cfDNA, data = df)
  tmp_cooks <- cooks.distance(tmp_lm)
  
  ## remove windows with cooks cutoff > 3 * mean
  filtered_df <- df[rownames(df) %in% names(tmp_cooks[tmp_cooks < 
                                                        cooks_cutoff]),]
  
  ## calculate cor
  filtered_cor <- cor(filtered_df$cfDNA, filtered_df$PBL)
  
  
  tmp_mat <- c(tmp_mat, filtered_cor)
}

## For unmatched PBLs, calculate median correlation of each sample
perSampleCorUnmatched <- c()
for(i in 1:50){
  tmp_mat2 <- CpGgt8_cfDNA_PBL_cooks_cor[i,]
  tmp_mat2 <- tmp_mat2[-i] # remove matched observation
  tmp_median <- median(tmp_mat2)
  perSampleCorUnmatched <- c(perSampleCorUnmatched, tmp_median)
  names(perSampleCorUnmatched)[i] <- paste0("cfDNA_",i)
}

## Collate correlations into data.frame
density_df <- 
  data.frame(source = c(rep("FaDu", length(tmp_mat)),
                        rep("Unmatched\nPBLs", 
                            length(perSampleCorUnmatched)),
                        rep("Matched\nPBLs",
                            length(CpGgt8_cfDNA_PBL_cooks_cor_matched))),
             cor = c(tmp_mat, 
                     perSampleCorUnmatched, 
                     CpGgt8_cfDNA_PBL_cooks_cor_matched))

density_df$source <- factor(density_df$source, 
                            levels = c("FaDu", "Unmatched\nPBLs", 
                                       "Matched\nPBLs"))

## Store median correlations for each comparison
cfDNAvsMatched <- median(CpGgt8_cfDNA_PBL_cooks_cor_matched)
cfDNAvsUnmatched <- median(perSampleCorUnmatched)
cfDNAvsFaDu <- median(tmp_mat)

## Create ggplot
set.seed(42)
density_plot <- density_df %>%
  subset(source != "FaDu") %>%
  ggplot(aes(x = cor, y = source, fill = source, col = source)) + 
  geom_density_ridges(size = 0.13, jittered_points = TRUE,
                      aes(point_shape = 16, 
                          point_alpha = c(0.5, 0.5, 0.5)[source])) +
  theme_JB() + 
  theme(axis.title.y = element_blank(), legend.position = "none",
        axis.text = element_text(size = 8), 
        axis.title.x = element_text(size = 8),
        plot.margin = margin(c(5,5,5,5))) + 
  scale_color_manual("vs. cfDNA", values = brewer.pal(3, "Set1")[2:3]) + 
  scale_fill_manual("vs. cfDNA", values = alpha(brewer.pal(3, "Set1")[2:3], 0.2)) + 
  lims(x = c(0.5,1)) + 
  coord_cartesian(expand = TRUE) + 
  labs(x = "Pearson's correlation") 

# Load cfMeDIP-seq libraries of HNSCC patients and healthy donors with 
# MEDEStrand converted absolute methylation levels
load(paste0(data_dir, "/PBL_CombinedBinMethyl_Size300_AllWindows.RData"))
load(paste0(data_dir, "/cfDNA_CombinedBinMethyl_Size300_AllWindows.RData"))

# Subset absolute methylation to regions w/ >= 8 CpGs
PBL_CombinedBinMethyl_Size300_CpGgt8 <- 
  PBL_CombinedBinMethyl_Size300_AllWindows[
    rownames(PBL_CombinedBinMethyl_Size300_AllWindows) %in% CpGgt8_IDs,]
cfDNA_CombinedBinMethyl_Size300_CpGgt8 <-
  cfDNA_CombinedBinMethyl_Size300_AllWindows[
    rownames(cfDNA_CombinedBinMethyl_Size300_AllWindows) %in% CpGgt8_IDs,]

rm(PBL_CombinedBinMethyl_Size300_AllWindows)
rm(cfDNA_CombinedBinMethyl_Size300_AllWindows)

# Select regions w/ median absolute methylation across healthy donors less
# than 0.1 in PBLs
Norm_PBL_BinMethyl_Size300_CpGgt8 <- PBL_CombinedBinMethyl_Size300_CpGgt8[
  , grepl("Norm", colnames(PBL_CombinedBinMethyl_Size300_CpGgt8))
]
Norm_PBLdepleted_Size300_CpGgt8 <- Norm_PBL_BinMethyl_Size300_CpGgt8[
  Biobase::rowMedians(Norm_PBL_BinMethyl_Size300_CpGgt8) <= 0.1,
]
PBLdepleted_IDs <- as.character(rownames(Norm_PBLdepleted_Size300_CpGgt8))

# Subset median absolute methylation across HNSCC PBLs based on previously
# identified "PBLdepleted" windows
HNC_PBL_BinMethyl_Size300_CpGgt8 <- PBL_CombinedBinMethyl_Size300_CpGgt8[ 
  , !grepl("Norm", colnames(PBL_CombinedBinMethyl_Size300_CpGgt8))
]
HNC_PBLdepleted_Size300_CpGgt8 <- HNC_PBL_BinMethyl_Size300_CpGgt8[
  Biobase::rowMedians(HNC_PBL_BinMethyl_Size300_CpGgt8) <= 0.1,
]

# Collate to data.frame for ggplot
HNC_and_Norm_binMethyl_PBLdepletion_data <- 
  as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(HNC_and_Norm_binMethyl_PBLdepletion_data) <- 
  c("sample", "condition", "depletion", "binMethyl")

for(i in 1:ncol(Norm_PBL_BinMethyl_Size300_CpGgt8)){
  
  # sample ID
  sample_ID <- as.character(i)
  
  # before depletion
  bf_depl <- as.numeric(Norm_PBL_BinMethyl_Size300_CpGgt8[,i])
  
  # after depletion
  af_depl <- as.numeric(Norm_PBLdepleted_Size300_CpGgt8[,i])
  
  # create variable vectors
  sample_ID_c <- rep(sample_ID, length(bf_depl) + length(af_depl))
  depl_c <- c(rep("before", length(bf_depl)), rep("after", length(af_depl)))
  condition_c <- rep("norm", length(sample_ID_c))
  
  # create data.frame for rbind
  tmp_df <- data.frame(sample = sample_ID_c, condition = condition_c, 
                       depletion = depl_c, binMethyl = c(bf_depl, af_depl))
  
  # add to data.frame
  HNC_and_Norm_binMethyl_PBLdepletion_data <- 
    rbind(HNC_and_Norm_binMethyl_PBLdepletion_data, tmp_df)
  
}

for(i in 1:ncol(HNC_PBL_BinMethyl_Size300_CpGgt8)){
  
  # sample ID
  sample_ID <- as.character(i)
  
  # before depletion
  bf_depl <- as.numeric(HNC_PBL_BinMethyl_Size300_CpGgt8[,i])
  
  # after depletion
  af_depl <- as.numeric(HNC_PBLdepleted_Size300_CpGgt8[,i])
  
  # create variable vectors
  sample_ID_c <- rep(sample_ID, length(bf_depl) + length(af_depl))
  depl_c <- c(rep("before", length(bf_depl)), rep("after", length(af_depl)))
  condition_c <- rep("HNC", length(sample_ID_c))
  
  # create data.frame for rbind
  tmp_df <- data.frame(sample = sample_ID_c, condition = condition_c, 
                       depletion = depl_c, binMethyl = c(bf_depl, af_depl))
  
  # add to data.frame
  HNC_and_Norm_binMethyl_PBLdepletion_data <- 
    rbind(HNC_and_Norm_binMethyl_PBLdepletion_data, tmp_df)
  
}

HNC_and_Norm_binMethyl_PBLdepletion_data$depletion <- 
  as.factor(HNC_and_Norm_binMethyl_PBLdepletion_data$depletion)
HNC_and_Norm_binMethyl_PBLdepletion_data$depletion <- 
  relevel(HNC_and_Norm_binMethyl_PBLdepletion_data$depletion, "before")

# Create ggplot
HNC_and_Norm_binMethyl_PBLdepletion_data$condition <- 
  as.character(HNC_and_Norm_binMethyl_PBLdepletion_data$condition)

HNC_and_Norm_binMethyl_PBLdepletion_data$condition <- 
  ifelse(HNC_and_Norm_binMethyl_PBLdepletion_data$condition == "norm", 
         "Healthy donor", "HNSCC")

HNC_and_Norm_binMethyl_PBLdepletion_data$condition <- 
  as.factor(HNC_and_Norm_binMethyl_PBLdepletion_data$condition)

HNC_and_Norm_binMethyl_PBLdepletion_data$sample <- 
  factor(HNC_and_Norm_binMethyl_PBLdepletion_data$sample, 
         levels = c(paste(1:30)))

set.seed(42)
HNC_and_Norm_binMethyl_PBLdepletion_plot <- 
  HNC_and_Norm_binMethyl_PBLdepletion_data[sample(1:nrow(HNC_and_Norm_binMethyl_PBLdepletion_data), 100000),] %>%
  ggplot(aes(x = sample, y = binMethyl)) + 
  geom_point(pch = 16, position = position_jitterdodge(), 
             alpha = 0.05, aes(col = depletion),
             size = 0.4) + 
  geom_boxplot(aes(fill = depletion), col = "black", outlier.shape = NA,
               lwd = 0.13, fatten = 6, alpha = 0) + 
  theme_JB() + 
  scale_colour_manual(values = c(palette_JB[2], palette_JB[1])) + 
  labs(y = "MeDEStrand\nAbsolute methylation") + 
  facet_grid(cols = vars(condition), scales = "free", space = "free_x") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Perform differential methylation analysis on HNSCC vs. healthy donor PBLs
## Load matrix containing counts from HNSCC and healthy donor PBLs
load(paste0(data_dir, "/PBL_CombinedCounts_AllWindows.RData"))

## Subset to PBLdepleted regions
PBL_CombinedCounts_PBLdepleted <- PBL_CombinedCounts_AllWindows[
  rownames(PBL_CombinedCounts_AllWindows) %in% PBLdepleted_IDs,
]

## Perform differential methylation analysis
coldata <- 
  data.frame(condition = 
               ifelse(grepl("Norm", colnames(PBL_CombinedCounts_PBLdepleted)),
                                         "Norm", "HNC"))
rownames(coldata) <- colnames(PBL_CombinedCounts_PBLdepleted)

dds <- DESeqDataSetFromMatrix(countData = PBL_CombinedCounts_PBLdepleted,
                              colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 20,]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "HNC", "Norm"), cooksCutoff = 1,
               independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH")

res_padj <- subset(res, res$padj < 0.05)
res_padj_hyperDMR <- subset(res_padj, res_padj$log2FoldChange > 0)
res_padj_hypoDMR <- subset(res_padj, res_padj$log2FoldChange < 0)

PBL_res_df <- 
  data.frame(x = res$log2FoldChange,
             y = -log10(res$pvalue),
             label = ifelse(is.na(res$padj), "ns",
                            ifelse(res$padj <= 0.05 &
                                     res$log2FoldChange < 0,
                                   "hypo", ifelse(res$padj <= 0.05 & 
                                                    res$log2FoldChange > 0,
                                                  "hyper", "ns"))))

PBL_res_plot <- PBL_res_df %>%
  ggplot(aes(x = x, y = y, col = label)) + 
  geom_point(size = 1.5) + 
  theme_JB() + 
  scale_color_manual(values = c("blue", rgb(0,0,0,1/50))) + 
  xlim(-5, 5) + 
  ylim(0, 10) + 
  labs(x = "log2FC(HNSCC/healthy donor)",
       y = "-log10(pvalue)")

# Save RData object of PBL-depleted windows
save(PBLdepleted_IDs, file = paste0(out_dir, "/PBLdepleted_IDs.RData"))
