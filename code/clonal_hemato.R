# Analysis of clonal hematopoeisis source code #
################################################

# Set directories
scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"

# Load libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)

# Load table generated from CAPP-Seq and iDES analysis:
# $ upn - patient identifier
# $ gene - gene symbol
# $ var.class - class of identified variant
# $ chr - chromosome
# $ pos - position
# $ ref - reference allele
# $ alt - alternative allele
# $ DcsSc.n.read.PBL - number of duplex consensus strand (w/ singleton
#                      correction) reads w/ alternative allele in matched PBLs
# $ allUniqDcs.vaf.PBL - variant allele fraction of alternative allele across
#                        total reads in matched PBLs
# $ allUniqDcs.n.read.PBL - number of total reads w/ alterantive allele in 
#                           matched PBLs
# $ allUniqDcs.depths.PBL - sequencing depth of total reads covering 
#                           alterantive allele in matched PBLs
# $ DcsSc.n.read.Baseline - same as DcsSc.n.read.PBL but in plasma sample
# $ allUniqDcs.vaf.Baseline - same as allUniqDcs.vaf.PBL but in plasma sample
# $ allUniqDcs.n.read.Baseline - same as allUniqDcs.n.read.PBL but in 
#                                plasma sample
# $ allUniqDcs.depths.Baseline - same as allUniqDcs.depths.PBL but in 
#                                plasma sample

vafAllUniqDcs <- read.table(file = paste0(data_dir, "/vafAllUniqDcs.txt"), 
                            header=T, stringsAsFactors = F)

# Create unique IDs for each SNV
# @ i.e. sample_ID.chromosome.basepair_position.reference_allele.mutant_allele
vafAllUniqDcs$ID <- paste(vafAllUniqDcs$upn, vafAllUniqDcs$chr, 
                          vafAllUniqDcs$pos, vafAllUniqDcs$ref, 
                          vafAllUniqDcs$alt, 
                          sep=".")

# Select SNVs present in matched cfDNA only, remove SNVs at extremes of
# read coverage 
cfDNA_SNVs <-  vafAllUniqDcs[vafAllUniqDcs$DcsSc.n.read.PBL == 0 &
                               vafAllUniqDcs$DcsSc.n.read.Baseline >= 3,]
q15_cfDNA <- quantile(vafAllUniqDcs$allUniqDcs.depths.Baseline, 0.15)
q85_cfDNA <- quantile(vafAllUniqDcs$allUniqDcs.depths.Baseline, 0.85)
q10_PBL <- quantile(vafAllUniqDcs$allUniqDcs.depths.PBL, 0.10)
q85_PBL <- quantile(vafAllUniqDcs$allUniqDcs.depths.PBL, 0.85)
cfDNA_index <- cfDNA_SNVs$allUniqDcs.depths.Baseline >= q15_cfDNA &
  cfDNA_SNVs$allUniqDcs.depths.Baseline <= q85_cfDNA
PBL_index <- cfDNA_SNVs$allUniqDcs.depths.PBL >= q10_PBL &
  cfDNA_SNVs$allUniqDcs.depths.PBL <= q85_PBL
cfDNA_DepthFiltered_SNVs <- cfDNA_SNVs[cfDNA_index & PBL_index,]

# Remove CNVs and non-coding variants
SNV_index <- !(cfDNA_DepthFiltered_SNVs$var.class %in% 
                 c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
                   "In_Frame_Ins", "Intron", "3'UTR", "IGR", "RNA"))
cfDNA_DepthFiltered_coding_SNVs <- cfDNA_DepthFiltered_SNVs[SNV_index,]

# Remove outlier genes w/ reduced read depth
cfDNA_DepthFiltered_coding_SNVs <- cfDNA_DepthFiltered_coding_SNVs[
  !(cfDNA_DepthFiltered_coding_SNVs$gene %in% c("FAM47A", "ZNF449")),]

# Repeat w/o PBL filter
all_cfDNA_index <- vafAllUniqDcs$allUniqDcs.depths.Baseline >= q15_cfDNA &
  vafAllUniqDcs$allUniqDcs.depths.Baseline <= q85_cfDNA

all_PBL_index <- vafAllUniqDcs$allUniqDcs.depths.PBL >= q10_PBL &
  vafAllUniqDcs$allUniqDcs.depths.PBL <= q85_PBL

all_vaf <- vafAllUniqDcs[all_cfDNA_index & all_PBL_index,]


# Remove CNVs and non-coding variants
all_SNV_index <- !(all_vaf$var.class %in% 
                     c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
                       "In_Frame_Ins", "Intron", "3'UTR", "IGR", "RNA"))
all_vaf_coding_SNVs <- all_vaf[all_SNV_index,]

# Remove outlier genes w/ reduced read depth
all_vaf_coding_SNVs <- all_vaf_coding_SNVs[
  !(all_vaf_coding_SNVs$gene %in% c("FAM47A", "ZNF449")),]


# Create ggplot
before_PBLfilter <- 
  all_vaf_coding_SNVs[all_vaf_coding_SNVs$DcsSc.n.read.Baseline >= 3 &
                        all_vaf_coding_SNVs$allUniqDcs.vaf.Baseline <= 0.1,]

before_PBLfilter <- 
  data.frame(upn = rep(before_PBLfilter$upn, 2),
             DCs = c(before_PBLfilter$DcsSc.n.read.Baseline,
                     before_PBLfilter$DcsSc.n.read.PBL),
             vaf = c(before_PBLfilter$allUniqDcs.vaf.Baseline,
                     before_PBLfilter$allUniqDcs.vaf.PBL),
             source = c(rep("cfDNA", length(before_PBLfilter$upn)),
                                          rep("PBL", 
                                              length(before_PBLfilter$upn))))

after_PBLfilter <- 
  all_vaf_coding_SNVs[all_vaf_coding_SNVs$DcsSc.n.read.PBL == 0 &
                        all_vaf_coding_SNVs$DcsSc.n.read.Baseline >= 3,]

after_PBLfilter <- 
  data.frame(upn = rep(after_PBLfilter$upn, 2),
             DCs = c(after_PBLfilter$DcsSc.n.read.Baseline,
                                      after_PBLfilter$DcsSc.n.read.PBL),
             vaf = c(after_PBLfilter$allUniqDcs.vaf.Baseline,
                                      after_PBLfilter$allUniqDcs.vaf.PBL),
             source = c(rep("cfDNA", length(after_PBLfilter$upn)),
                        rep("PBL", length(after_PBLfilter$upn))))

before_after_PBLfilter <- rbind(before_PBLfilter, after_PBLfilter)
before_after_PBLfilter$filter <- 
  c(rep("Before\nPBL filter", nrow(before_PBLfilter)),
    rep("After\nPBL filter", nrow(after_PBLfilter)))

before_after_PBLfilter$filter <- factor(before_after_PBLfilter$filter,
                                        levels = c("Before\nPBL filter", 
                                                   "After\nPBL filter"))

HN_IDs_ordered <- all_vaf_coding_SNVs %>% 
  subset(DcsSc.n.read.Baseline >= 3 & allUniqDcs.vaf.Baseline <= 0.1) %>%
  group_by(upn) %>%
  summarize(meanMAF = mean(allUniqDcs.vaf.Baseline))

HN_IDs_ordered <- HN_IDs_ordered[order(HN_IDs_ordered$meanMAF,
                                       decreasing = TRUE),]$upn

before_after_PBLfilter$upn <- factor(before_after_PBLfilter$upn,
                                     levels = c(HN_IDs_ordered))

cfDNA_vs_PBLs_SNVs <- all_vaf_coding_SNVs %>%
  ggplot(aes(x = allUniqDcs.vaf.Baseline, y = allUniqDcs.vaf.PBL)) + 
  geom_abline(size = 0.26, col = "red") + 
  geom_point(size = 1.5, col = "black", alpha = 0.2, shape = 16) + 
  theme_JB() + 
  labs(x = "Mutant allele fraction - cfDNA", 
       y = "Mutant allele fraction - PBLs") +
  theme(panel.background = element_blank())

cfDNA_vs_PBLs_SNVs_filtered <- 
  all_vaf_coding_SNVs[all_vaf_coding_SNVs$DcsSc.n.read.Baseline >= 3 &
                        all_vaf_coding_SNVs$allUniqDcs.vaf.Baseline <= 0.1,] %>%
  ggplot(aes(x = allUniqDcs.vaf.Baseline * 100, y = allUniqDcs.vaf.PBL * 100)) + 
  geom_smooth(method = "lm", 
              data = all_vaf_coding_SNVs[
                (all_vaf_coding_SNVs$DcsSc.n.read.Baseline >= 3 &
                   all_vaf_coding_SNVs$allUniqDcs.vaf.Baseline <= 0.1) &
                  all_vaf_coding_SNVs$allUniqDcs.vaf.PBL > 0,],
              aes(x = allUniqDcs.vaf.Baseline * 100, 
                  y = allUniqDcs.vaf.PBL * 100),
              inherit.aes = FALSE,
              size = 0.26,
              col = rgb(1,0,0,1/4),
              fullrange = TRUE, se = FALSE) + 
  geom_point(size = 1.5, col = "red", alpha = 0.5, shape = 16) + 
  theme_JB() + 
  labs(x = "Mutant allele fraction - cfDNA", 
       y = "Mutant allele fraction - PBLs") + 
  theme(panel.background = element_blank()) + 
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 10),
                     breaks = c(0, 0.1, 1, 10)) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 10),
                     breaks = c(0, 0.1, 1, 10)) + 
  annotate("rect", xmin = 0.085, xmax = 11.5, ymin = -0.03, ymax = 0.03,
           fill = NA, col = "red", lty = 2, lwd = 0.13) + 
  annotate("text", x = 5, y = 2, label = "R = 0.94", size = (4/15) * 9)

set.seed(43)

before_after_PBLfilter_plot <- before_after_PBLfilter %>% 
  ggplot(aes(x = upn, y = vaf * 100, fill = source)) + 
  geom_point(position = position_jitterdodge(0.4),
             size = 1.5, aes(col = source)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.13, fatten = 3) + 
  facet_grid(rows = vars(filter), scales = "free_y",
             space = "free_y") + 
  scale_fill_manual(values = c(rgb(0,0,0,0), rgb(0,0,0,0))) + 
  scale_color_manual(values = c(palette_JB[1], palette_JB[2])) + 
  theme_JB() + 
  theme(axis.text.x.bottom = element_blank(),
        panel.grid.major.x = element_line(size = 0.13, color = rgb(0,0,0,1/10)),
        plot.margin = unit(c(1,1,1,1), "mm"),
        axis.title.x = element_text(vjust = 1,
                                    margin = margin(5,0,0,0)),
        legend.key.size = unit(0.4, "cm")) + 
  labs(x = "HNSCC samples", y = "Mutant allele fraction (%)") + 
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 10),
                     breaks = c(0, 0.1, 1, 10))


vaf_bygenes <- 
  all_vaf_coding_SNVs[all_vaf_coding_SNVs$DcsSc.n.read.Baseline >= 3 &
                        all_vaf_coding_SNVs$allUniqDcs.vaf.Baseline <= 0.1,]
vaf_bygenes_CHIP <- vaf_bygenes[vaf_bygenes$allUniqDcs.n.read.PBL > 0,]

## Create oncoprint
DepthFiltered_SNP_variants <- vaf_bygenes_CHIP

HN_list <- unique(vaf_bygenes_CHIP$upn)
gene_list <- unique(DepthFiltered_SNP_variants$gene)
mat <- matrix(nrow=length(HN_list), ncol=length(gene_list))
colnames(mat) <- gene_list
rownames(mat) <- HN_list

for(hn in 1:nrow(mat)){
  hn_tmp <- rownames(mat)[hn]
  muts_tmp <- 
    unique(DepthFiltered_SNP_variants[DepthFiltered_SNP_variants$upn == 
                                        hn_tmp,1:3])
  if(nrow(muts_tmp) >= 1){
    for(genes in 1:ncol(mat)){
      if(colnames(mat)[genes] %in% muts_tmp$gene){
        dir <- which(muts_tmp$gene %in% colnames(mat)[genes])
        vars <- muts_tmp[dir,]$var.class[1]
        mat[hn,genes] <- vars
      } else {
        mat[hn,genes] <- NA
      }
    }
  } else {
    mat[hn,] <- rep(NA,ncol(mat))
  }
}

get_type_fun <- function(x){
  strsplit(x,";")[[1]]
}

col <- c(Missense_Mutation = "red", Nonsense_Mutation = "blue", 
         Splice_Site = "green", Silent = "black")

mat_melted <- melt(t(mat))

FullOncoprint_data <- mat_melted
colnames(FullOncoprint_data) <- c("cols", "rows", "value")

## Create components of oncoprint...
gg_hm = FullOncoprint_data %>%
  ggplot(aes(x = rows, y = cols, fill = value)) + 
  geom_tile(colour = "lightgrey", size = 0.13) +
  theme_JB() + 
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        legend.position = "bottom") + 
  scale_fill_JB(name = "SNV:", 
                labels = c("Missense", "Nonsense", 
                           "Silent", "None")) + 
  scale_x_discrete(limits = rev(levels(FullOncoprint_data$rows)))

gg_rows = FullOncoprint_data %>%
  group_by(rows) %>%
  summarize(value = sum(!is.na(value))) %>%
  ggplot(aes(x = rows, y = value)) + 
  geom_bar(stat = "identity", position = "dodge",
           fill = "white", col = "black",
           size = 0.13) + 
  scale_y_continuous(breaks = seq(0,12,2)) + 
  labs(y = "# SNVs") + 
  theme_JB() + 
  theme(axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin =  unit(c(0,5,-10,0), "mm")) + 
  scale_x_discrete(limits = rev(levels(FullOncoprint_data$rows)))

gg_cols = FullOncoprint_data %>%
  group_by(cols) %>%
  summarize(value = sum(!is.na(value))) %>%
  ggplot(aes(x = cols, y = value)) + 
  geom_bar(stat = "identity", position = "dodge",
           fill = "white", col = "black", size = 0.13) + 
  scale_y_continuous(position = "right",
                     breaks = seq(0,12,2)) + 
  labs(y = "# samples") + 
  coord_flip() + 
  theme_JB() + 
  theme(axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,-5), "mm"))

gg_empty = FullOncoprint_data %>%
  ggplot(aes(x = cols, y = value)) + 
  geom_blank() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0,5,-10,-5), "mm")) + 
  scale_x_discrete(limits = rev(levels(FullOncoprint_data$rows)))

a <- egg::ggarrange(gg_rows, gg_empty, gg_hm, gg_cols,
                    nrow = 2, ncol = 2, widths = c(3,1), heights = c(0.5,3),
                    draw = FALSE)

CHIP_oncoprint_plot <- as_ggplot(a)

# Save object for subsequent scripts
cfDNA_PBLfiltered_SNVs <- after_PBLfilter[after_PBLfilter$source == 
                                            "cfDNA",]

save(cfDNA_PBLfiltered_SNVs, file = 
       paste0(out_dir, "/cfDNA_PBLfiltered_SNVs.RData"))
save(cfDNA_DepthFiltered_coding_SNVs,
     file = paste0(out_dir, "/cfDNA_DepthFiltered_coding_SNVs.RData"))
