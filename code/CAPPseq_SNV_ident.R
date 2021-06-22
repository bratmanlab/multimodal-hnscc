# Tumor-naive detection of pre-treatment ctDNA by CAPP-Seq source code #
########################################################################

# Set directories
scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"

# Load libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)

# Load RData object w/ SNVs after PBL filtering
# source(paste0(scripts, "/clonal_hemato.R"))
load(paste0(out_dir, "/cfDNA_PBLfiltered_SNVs.RData"))
#load(paste0(out_dir, "/all_vaf_coding_SNVs.RData"))

# Calculate mean MAF of identified SNVs for each patient
cfDNA_PBLfiltered_SNVs_meanMAF <- cfDNA_PBLfiltered_SNVs %>%
  group_by(upn) %>%
  summarise(meanMAF = mean(vaf) * 100)

cfDNA_PBLfiltered_SNVs_meanMAF <- as.data.frame(cfDNA_PBLfiltered_SNVs_meanMAF)
cfDNA_PBLfiltered_SNVs_meanMAF$upn <- 
  as.character(cfDNA_PBLfiltered_SNVs_meanMAF$upn)

# Add patients w/ non-dectable SNVs by CAPP-Seq
# NOTE: table was previously generated in workflow.R, uncomment code if this
#       script is being ran independently
# dsDNA_logs <- read.table(file = paste0(data_dir, "/hnc_DNAperPlasma.txt"),
#                          sep = "\t", header = TRUE)
# dsDNA_logs$Timepoint2 <- factor(dsDNA_logs$Timepoint2, 
#                                 levels = c("Diagnosis", "Post-surgery",
#                                            "Mid-radiotherapy",
#                                            "Post-treatment-1", 
#                                            "Post-treatment-2", "Normal"))

HN_IDs <- subset(dsDNA_logs, dsDNA_logs$Timepoint2 == "Diagnosis")$HNC

for(i in 1:length(HN_IDs)){
  if(sum(grepl(HN_IDs[i], cfDNA_PBLfiltered_SNVs_meanMAF$upn)) == 0){
    tmp_df <- data.frame(upn = paste0("HN", HN_IDs[i]), meanMAF = 0)
    cfDNA_PBLfiltered_SNVs_meanMAF <- 
      rbind(cfDNA_PBLfiltered_SNVs_meanMAF, tmp_df)
  }
}

# Reorganize data.frame for ggplot
cfDNA_PBLfiltered_SNVs_meanMAF$upn <- 
  factor(cfDNA_PBLfiltered_SNVs_meanMAF$upn,
         levels = cfDNA_PBLfiltered_SNVs_meanMAF$upn[
           order(cfDNA_PBLfiltered_SNVs_meanMAF$meanMAF, 
                 decreasing = TRUE)
           ]
         )

# Create OncoPrint data.frame
HN_list <- cfDNA_PBLfiltered_SNVs_meanMAF$upn
## all_vaf_coding_SNVs previously generated from clonal_hemato.R script

gene_list <- unique(cfDNA_DepthFiltered_coding_SNVs$gene) 
mat <- matrix(nrow=length(HN_list), ncol=length(gene_list))
colnames(mat) <- gene_list
rownames(mat) <- HN_list

for(hn in 1:nrow(mat)){
  hn_tmp <- rownames(mat)[hn]
  muts_tmp <- 
    unique(cfDNA_DepthFiltered_coding_SNVs[cfDNA_DepthFiltered_coding_SNVs$upn == 
                                        hn_tmp,1:3])
  if(nrow(muts_tmp) >= 1){
    for(genes in 1:ncol(mat)){
      if(colnames(mat)[genes] %in% muts_tmp$gene){
        dir <- which(muts_tmp$gene %in% colnames(mat)[genes])
        vars <- muts_tmp[dir,]$var.class
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

sample_index <- tapply(1:nrow(mat_melted), factor(mat_melted$Var2), 
                       function(x) sum(!is.na(mat_melted[x,]$value)))
sample_index <- sample_index[order(sample_index, decreasing = TRUE)]
sample_index <- match(names(sample_index), levels(mat_melted$Var2))

mat_melted$Var2 <- factor(mat_melted$Var2, 
                          levels(mat_melted$Var2)[sample_index])

best_gene_index <- tapply(1:nrow(mat_melted), factor(mat_melted$Var2), 
                          function(x) sum(!is.na(mat_melted[x,]$value) &
                                            mat_melted[x,]$Var1 == "TP53"))
best_gene_index <- best_gene_index[order(best_gene_index, decreasing = TRUE)]
best_gene_index <- match(names(best_gene_index), levels(mat_melted$Var2))

mat_melted$Var2 <- factor(mat_melted$Var2, levels(mat_melted$Var2)[best_gene_index])

gene_index <- tapply(1:nrow(mat_melted), factor(mat_melted$Var1), function(x) sum(!is.na(mat_melted[x,]$value)))
gene_index <- gene_index[order(gene_index, decreasing = TRUE)]
gene_index <- match(names(gene_index), levels(mat_melted$Var1))

mat_melted$Var1 <- factor(mat_melted$Var1, 
                          levels(mat_melted$Var1)[rev(gene_index)])

test <- sapply(table(mat_melted[!is.na(mat_melted$value),]$Var1), 
               function(x) seq(1,x,by=1))
test2 <- sapply(table(as.character(mat_melted[!is.na(mat_melted$value),]$Var2)), 
                function(x) seq(1,x,by=1))

gene_percent <- 
  round(sapply(tapply(1:nrow(mat_melted[!is.na(mat_melted$value),]), 
                      factor(mat_melted[!is.na(mat_melted$value),]$Var1),
                      function(x) sapply(x, length)), sum) / 32 * 100)
gene_percent_df <- data.frame(label = gene_percent, gene = names(gene_percent))

colnames(mat_melted) <- c("cols", "rows", "value")

FullOncoprint_data <- mat_melted

# Create ggplot
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
                           "Silent", "Splice site", "None"))

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
        plot.margin =  unit(c(0,5,-10,0), "mm"))

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
        plot.margin = unit(c(0,5,-10,-5), "mm"))

a <- egg::ggarrange(gg_rows, gg_empty, gg_hm, gg_cols,
                    nrow = 2, ncol = 2, widths = c(3,1), heights = c(0.5,3),
                    draw = FALSE)

FullOncoprint_plot <- as_ggplot(a)

## Generate reduced OncoPrint containing genes of interest
## Small oncoprint
SmallOncoprint_data <- FullOncoprint_data[grepl("TP53|PIK3CA|FAT1|TP63|NOTCH1|MYC|CASP8|NFE2L2", 
                                                FullOncoprint_data$cols),]
# create "other" row
OtherOncoprint <- 
  FullOncoprint_data[!grepl("TP53|PIK3CA|FAT1|TP63|NOTCH1|MYC|CASP8|NFE2L2", 
                            FullOncoprint_data$cols),]
OtherOncoprint$cols <- "Other"
OtherOncoprint$value <- ifelse(is.na(OtherOncoprint$value), NA, "Other")

SmallOncoprint_data <- rbind(SmallOncoprint_data, OtherOncoprint)

SmallOncoprint_data$cols <- factor(SmallOncoprint_data$cols,
                                   levels = rev(c("TP53", "PIK3CA", "FAT1", 
                                                  "TP63", "NOTCH1", "MYC", 
                                                  "CASP8", "NFE2L2", "Other")))
gg_hm = SmallOncoprint_data %>%
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
  scale_fill_manual(name = "SNV",
                    values = c("#ef3b2c","#386cb0","grey90","#fdb462", "#7fc97f", "#a6cee3","#fb9a99","#984ea3","#ffff33"),
                    labels = c("Missense", "Nonsense", 
                               "Other", "Silent", "Splice site", "None"))

gg_rows = SmallOncoprint_data %>%
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
        plot.margin =  unit(c(0,5,-10,0), "mm"))

gg_cols <- SmallOncoprint_data %>%
  subset(cols != "Other") %>%
  group_by(cols) %>%
  summarize(value = sum(!is.na(value)))

gg_cols <- rbind(gg_cols, data.frame(cols = "Other",
                                     value = 0))

gg_cols = gg_cols %>%
  ggplot(aes(x = cols, y = value)) + 
  geom_bar(stat = "identity", position = "dodge",
           fill = "white", col = "black", size = 0.13) + 
  scale_y_continuous(position = "right",
                     breaks = seq(0,12,4)) + 
  labs(y = "# samples") + 
  coord_flip() + 
  theme_JB() + 
  theme(axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,-5), "mm"))

gg_empty = SmallOncoprint_data %>%
  ggplot(aes(x = cols, y = value)) + 
  geom_blank() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0,5,-10,-5), "mm"))

a <- egg::ggarrange(gg_rows, 
                    gg_empty, gg_hm, gg_cols,
                    nrow = 2, ncol = 2, widths = c(9,2), heights = c(2,3),
                    draw = FALSE)

SmallOncoprint_plot <- as_ggplot(a)

# Create ggplot
MAFperPatient_plot <- cfDNA_PBLfiltered_SNVs_meanMAF %>%
  ggplot(aes(x = upn, y = meanMAF)) + 
  geom_bar(stat = "identity", position = "dodge",
           col = palette_JB[1], fill = palette_JB[1], size = 0.13,
           width = 0.75) + 
  geom_hline(yintercept = median(cfDNA_PBLfiltered_SNVs_meanMAF$meanMAF), 
             size = 0.2, lty = 2, col = "grey30") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 5)) + 
  scale_x_discrete(expand = c(0,1)) +
  theme_JB() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_text(margin = margin(15,0,0,0))) + 
  labs(y = "Mean mutant allele fraciton (%)",
       x = "Samples")

# Save objects
save(cfDNA_PBLfiltered_SNVs_meanMAF, 
     file = paste0(out_dir, "/cfDNA_PBLfiltered_SNVs_meanMAF.RData"))
