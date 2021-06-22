# Longitduinal analysis #
#########################

# Set directories
scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"

# 1. Load libraries -----------------------------------------------------------
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(data.table)
library(egg)
library(DT)
source(paste0(scripts, "/ggplot_themeJB.R"))

# 2. Generate table w/ clinical info and mean RPKM across hyperDMRs -----------

## Load and subset clinical and RPKM tables
Clinical_info <- read.table(paste0(data_dir, "/Clinical_info.txt"),
                            sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE)

load(paste0(data_dir,
            "/cfDNA_BLwPostTx_Frag100_150_CombinedRPKMs_hyperDMRs.RData"))

## Norm RPKMs...
cfDNA_Norm_RPKMs_PBLdepleted <- 
  cfDNA_BLwPostTx_Frag100_150_CombinedRPKMs_hyperDMRs[
    ,grepl("Norm", colnames(cfDNA_BLwPostTx_Frag100_150_CombinedRPKMs_hyperDMRs))
    ]

## Baseline RPKMs...
cfDNA_Baseline_RPKMs_PBLdepleted <- 
  cfDNA_BLwPostTx_Frag100_150_CombinedRPKMs_hyperDMRs[
    ,grepl("BL", colnames(cfDNA_BLwPostTx_Frag100_150_CombinedRPKMs_hyperDMRs))
    ]

## Post-treatment RPKMs...
cfDNA_PostTx_RPKMs_PBLdepleted <- 
  cfDNA_BLwPostTx_Frag100_150_CombinedRPKMs_hyperDMRs[
    ,!grepl("Norm|BL", colnames(cfDNA_BLwPostTx_Frag100_150_CombinedRPKMs_hyperDMRs))
    ]

PostTx_index <- unique(substr(colnames(cfDNA_PostTx_RPKMs_PBLdepleted), 1, 4))
Baseline_index <- substr(colnames(cfDNA_Baseline_RPKMs_PBLdepleted), 1, 4)

cfDNA_PostTx_RPKMs_PBLdepleted <- 
  cfDNA_PostTx_RPKMs_PBLdepleted[ 
    ,grepl(paste0(Baseline_index, collapse = "|"), 
           colnames(cfDNA_PostTx_RPKMs_PBLdepleted))]

cfDNA_Baseline_RPKMs_PBLdepleted <- 
  cfDNA_Baseline_RPKMs_PBLdepleted[ 
    ,grepl(paste0(PostTx_index, collapse = "|"), 
           colnames(cfDNA_Baseline_RPKMs_PBLdepleted))]

CombinedTx_hyperDMR_RPKMs <- cbind(cfDNA_Baseline_RPKMs_PBLdepleted,
                                   cfDNA_PostTx_RPKMs_PBLdepleted,
                                   cfDNA_Norm_RPKMs_PBLdepleted)

CombinedTx_hyperDMR_MeanRPKMs <- colMeans(CombinedTx_hyperDMR_RPKMs)

cfDNA_Norm_RPKMs_PBLdepleted <- 
  CombinedTx_hyperDMR_MeanRPKMs[grepl("Norm", 
                                      names(CombinedTx_hyperDMR_MeanRPKMs))]

Norm_meanRPKM <- mean(cfDNA_Norm_RPKMs_PBLdepleted)
Norm_maxRPKM <- max(cfDNA_Norm_RPKMs_PBLdepleted)
Norm_minRPKM <- min(cfDNA_Norm_RPKMs_PBLdepleted)

hyperDMR_colMeans <- c()
for(i in 1:nrow(Clinical_info)){
  grepl_index <- paste(Clinical_info$HN[i], Clinical_info$Timepoint[i],
                       sep = "_")
  tmp <- 
    CombinedTx_hyperDMR_MeanRPKMs[grepl(grepl_index,
                                        names(CombinedTx_hyperDMR_MeanRPKMs))]
  hyperDMR_colMeans <- c(hyperDMR_colMeans, tmp)
}

Clinical_info$meanRPKM <- hyperDMR_colMeans

## Estimate ctDNA percentage based on comparison w/ CAPP-Seq data -------------
load(paste0(out_dir, "/mclust_wnorm_df.RData"))

MAF_ind <- mclust_wnorm_df$ctDNA[
  match(substr(colnames(cfDNA_Baseline_RPKMs_PBLdepleted), 1, 4),
        substr(rownames(mclust_wnorm_df), 3, 6))]

cor_plot <- 
  cor.test(x = colMeans(cfDNA_Baseline_RPKMs_PBLdepleted), 
           y = MAF_ind)

MvR_cor <- cor_plot$estimate
MvR_pval <- cor_plot$p.value

lm_df <- data.frame(ctDNA = MAF_ind,
                    hyperDMR = colMeans(cfDNA_Baseline_RPKMs_PBLdepleted))

linearMod <- lm(ctDNA ~ hyperDMR, 
                data = lm_df)
print(linearMod)
# formula = y = 3.2394x - 0.1233

RPKMtoCtDNA <- function(x){
  3.2394 * x - 0.1233
}

Clinical_info$MeDIP_ctDNA <- RPKMtoCtDNA(Clinical_info$meanRPKM)

## Determine lower limit of detection based off healthy donors...
MeDIP_ctDNA_LLOD <- RPKMtoCtDNA(Norm_maxRPKM)
MeDIP_ctDNA_LLOD # approximately 0.2%

## Plot
mclust_wnorm_100_150_df <- 
  data.frame(meanRPKM = c(colMeans(cfDNA_Baseline_RPKMs_PBLdepleted),
                          cfDNA_Norm_RPKMs_PBLdepleted),
             ctDNA = c(MAF_ind, rep(0, length(cfDNA_Norm_RPKMs_PBLdepleted))),
             condition = c(rep("HNC", ncol(cfDNA_Baseline_RPKMs_PBLdepleted)),
                           rep("Norm", length(cfDNA_Norm_RPKMs_PBLdepleted))))

hyperDMR_shortfrag_corplot <- 
  mclust_wnorm_100_150_df %>%
  ggplot(aes(x = meanRPKM, y = ctDNA, colour = condition)) + 
  geom_smooth(data = subset(mclust_wnorm_100_150_df, 
                            mclust_wnorm_100_150_df$condition == "HNC"),
              aes(x = meanRPKM, y = ctDNA),
              method = lm, fill = "lightgrey", fullrange = TRUE,
              show.legend = FALSE, col = "red", size = 0.13) + 
  geom_point(size = 1.5,
             stroke = 1,
             aes(shape = condition),
             alpha = 0.5) +
  annotate(geom = "text",
           x = 0.03, y = 4.5,
           label = paste0("R = ",round(MvR_cor, 2),
                          "\np = ", round(MvR_pval, 9)),
           hjust = 0,
           size = (5/14) * 9) + 
  geom_hline(yintercept = MeDIP_ctDNA_LLOD, 
             lty = 2, col = "blue", size = 0.13) +
  coord_cartesian(ylim = c(-0.1, 5), expand = FALSE) +
  xlim(c(0.01, 1.25)) + 
  labs(y = "Mutant allele frequency (%)",
       x = "Methylation (mean RPKM)") + 
  theme_JB() + 
  scale_colour_JB("Condition:") + 
  scale_shape_manual("Condition:", values = c(16, 4))

# 3. Create list w/ ctDNA converted values at each timepoint ------------------

HN_IDs <- 
  mclust_wnorm_100_150_df %>%
  subset(mclust_wnorm_100_150_df$condition == "HNC") %>%
  subset(mclust_wnorm_100_150_df$ctDNA > 0 &
           mclust_wnorm_100_150_df$meanRPKM > Norm_maxRPKM) %>%
  rownames() %>%
  substr(1, 4)

HN_list <- c()
for(i in HN_IDs){
  tmp_df <- Clinical_info[grepl(i, Clinical_info$HN),]
  start <- as.Date(tmp_df[tmp_df$Timepoint == "BL",]$Visit_date)
  start_RPKM <- tmp_df[tmp_df$Timepoint == "BL",]$MeDIP_ctDNA
  times <- as.Date(tmp_df[tmp_df$Timepoint != "BL",]$Visit_date)
  times_RPKM <- tmp_df[tmp_df$Timepoint != "BL",]$MeDIP_ctDNA
  days_postBL <- times - start
  days_BL <- 1
  if(tmp_df$Date_of_Surgery[1] != ""){
    days_surgery <- as.Date(tmp_df$Date_of_Surgery[1]) - start 
  } else {
    days_surgery <- NA
  }
  if(tmp_df$Date_of_Failure[1] != ""){
    days_failure <- as.Date(tmp_df$Date_of_Failure[1]) - start 
  } else {
    days_failure <- NA
  }
  if(tmp_df$RT_start[1] != ""){
    days_RTstart <- as.Date(tmp_df$RT_start[1]) - start
    days_RTend <- as.Date(tmp_df$RT_end[1]) - start
  } else {
    days_RTstart <- NA
    days_RTend <- NA
  }
  days_lastFU <- as.Date(tmp_df$Last_FU[1]) - start
  tmp_df2 <- data.frame(time = c(days_BL, days_postBL),
                        MeDIP_ctDNA = c(start_RPKM, times_RPKM))
  dates <- c(as.numeric(days_surgery), as.numeric(days_RTstart), 
             as.numeric(days_RTend), as.numeric(days_failure),
             as.numeric(days_lastFU))
  names(dates) <- c("surgery", "RTstart", "RTend", "relapse", "lastFU")
  tmp_list <- list(tmp_df2, dates)
  HN_list[[paste(i)]] <- tmp_list
}


# 4. Create longitudinal plots for each sample --------------------------------
long_df <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(long_df) <- c("HN", "time", "ctDNA", "condition")

test_HN <- HN_list$`2187`
test_HN_tmp <- data.frame(HN = rep("2187", nrow(test_HN[[1]]) + 
                                     length(test_HN[[2]])),
                          time = c(test_HN[[1]]$time,
                                   as.numeric(test_HN[[2]])),
                          ctDNA = c(test_HN[[1]]$MeDIP_ctDNA,
                                    rep(NA, length(test_HN[[2]]))),
                          condition = c(rep("cfDNA", nrow(test_HN[[1]])), 
                                        names(test_HN[[2]])))

for(i in 1:length(HN_list)){
  test_HN <- HN_list[[i]]
  HN <- names(HN_list)[i]
  
  test_HN_tmp <- data.frame(HN = rep(HN, nrow(test_HN[[1]]) + 
                                       length(test_HN[[2]])),
                            time = c(test_HN[[1]]$time,
                                     as.numeric(test_HN[[2]])),
                            ctDNA = c(test_HN[[1]]$MeDIP_ctDNA,
                                      rep(NA, length(test_HN[[2]]))),
                            condition = c(rep("cfDNA", nrow(test_HN[[1]])), 
                                          names(test_HN[[2]])))
  long_df <- rbind(long_df, test_HN_tmp)
}

plot_list <- list()
for(i in 1:length(names(HN_list))){
  test_HN_tmp <- subset(long_df, as.character(long_df$HN) ==
                          names(HN_list)[i])
  
  cfDNA_sub <- subset(test_HN_tmp,
                      test_HN_tmp$condition == "cfDNA")
  cfDNA_sub$stat <- ifelse(cfDNA_sub$ctDNA > MeDIP_ctDNA_LLOD,
                           "pos", "neg")
  cfDNA_sub$stat <- factor(cfDNA_sub$stat, levels = c("pos", "neg"))
  ctDNA2 <- c()
  for(n in 1:nrow(cfDNA_sub)){
    if(cfDNA_sub$ctDNA[n] > MeDIP_ctDNA_LLOD){
      ctDNA2 <- c(ctDNA2, cfDNA_sub$ctDNA[n])
    } else {
      ctDNA2 <- c(ctDNA2, 0)
    }
  }
  cfDNA_sub$ctDNA2 <- ctDNA2
  
  if(length(unique(cfDNA_sub$stat)) == 2){
    col_stat <- c("black", "white")
  } else {
    if(sum(unique(cfDNA_sub$stat) == "pos") == 1){
      col_stat <- c("black", "white")
    } else {
    col_stat <- c("white", "black")
    }
  }
  
  plot_list[[i]] <- ggplot(data = test_HN_tmp) + 
    geom_rect(data = subset(test_HN_tmp,
                            test_HN_tmp$condition %in% 
                              c("surgery", "relapse", "lastFU")),
              aes(xmin = time, xmax = time + (time[3] * 0.0125), 
                  ymin = -Inf, ymax = Inf),
              fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/2), "black")) + 
    geom_rect(data = subset(test_HN_tmp,
                            test_HN_tmp$condition %in%
                              c("RTstart", "RTend")),
              aes(xmin = time[1], xmax = time[2],
                  ymin = -Inf, ymax = Inf),
              fill = rgb(0,1,0,1/4)) + 
    geom_line(data = cfDNA_sub,
              aes(x = time, y = ctDNA2),
              size = 0.26) + 
    geom_point(data = cfDNA_sub,
               aes(x = time, y = ctDNA2, fill = stat), 
               size = 1,
               shape = 21,
               col = "black",
               inherit.aes = TRUE) + 
    geom_hline(yintercept = MeDIP_ctDNA_LLOD,
               col = "blue", lty = 2, size = 0.13) + 
    theme_JBwBorder() + 
    scale_x_continuous(expand = c(0, 20)) + 
    scale_y_continuous(trans = scales::pseudo_log_trans(0.05, 10),
                       breaks = c(0, seq(0.2, 1, by = 0.1),
                                  seq(2, 10, length.out = 9)),
                       labels = c("N.D", "0.2", rep("",7), "1", rep("", 8),
                                  "10"),
                       expand = c(0,0)) + 
    coord_cartesian(ylim = c(0, 10), clip = "off") + 
    scale_fill_manual(values = col_stat) + 
    labs(y = "ctDNA (%)", 
         x = "Time (days)") + 
    theme(axis.line.x.bottom = element_blank(),
          legend.position = "none") + 
    ggtitle(paste(names(HN_list)[i]))
}

names(plot_list) <- names(HN_list)

AllPlots <- ggarrange(plots = plot_list, nrow = 6, ncol = 3,
                      newpage = TRUE)

# 6. Define kinetics ----------------------------------------------------------

Clinical_info_lastFU <- 
  as.data.frame(Clinical_info %>% group_by(HN) %>% top_n(1, Visit_date))
Clinical_info_lastFU <- 
  subset(Clinical_info_lastFU, Clinical_info_lastFU$HN %in% HN_IDs)
Clinical_info_lastFU <-
  rbind(Clinical_info[Clinical_info$Timepoint == "BL" & Clinical_info$HN %in% HN_IDs,],
        Clinical_info_lastFU)
Clinical_info_lastFU$Timepoint <- 
  ifelse(Clinical_info_lastFU$Timepoint == "BL",
         "Pre-treatment", "Mid-/post-treatment")
Clinical_info_lastFU$Timepoint <- 
  factor(Clinical_info_lastFU$Timepoint,
         levels = c("Pre-treatment", "Mid-/post-treatment"))

preTreat <- Clinical_info_lastFU[Clinical_info_lastFU$Timepoint == "Pre-treatment",]$MeDIP_ctDNA
postTreat <- Clinical_info_lastFU[Clinical_info_lastFU$Timepoint == "Mid-/post-treatment",]$MeDIP_ctDNA

kinetic <- ((postTreat/preTreat) - 1) * 100
names(kinetic) <- unique(Clinical_info_lastFU$HN)
kinetic

annt.kinetic <- c()
for(i in 1:length(kinetic)){
  tmp <- kinetic[i]
  tmp.id <- names(tmp)
  tmp.df <- Clinical_info_lastFU[Clinical_info_lastFU$HN == tmp.id,]
  
  if(tmp.df[tmp.df$Timepoint == "Mid-/post-treatment",]$MeDIP_ctDNA < MeDIP_ctDNA_LLOD){
    tmp.annt <- "Complete clearance"
  } else {
    tmp.annt <- ifelse(tmp <= -90,
                       "Partial clearance",
                       "No clearance")
  }
  annt.kinetic <- c(annt.kinetic, tmp.annt)
  }

names(annt.kinetic) <- names(kinetic)

# 5. Kinetic examples ---------------------------------------------------------

complete_clearance <- plot_list$`2291` + ggtitle("Complete clearance")
partial_clearance <- plot_list$`2423` + ggtitle("Partial clearance")
stable_disease <- plot_list$`2298` + ggtitle("No clearance")
long_legend <- 
  cowplot::get_legend(plot_list$`2464` + theme(legend.position = "bottom"))

kinetics_ex <- cowplot::plot_grid(complete_clearance, partial_clearance, 
                                stable_disease,
                                ncol = 3)
kinetics_ex_wlegend <- 
  cowplot::plot_grid(long_legend, kinetics_ex, nrow = 2, rel_heights = c(1,9))

# 9. BL vs. last-FU -----------------------------------------------------------

Clinical_info_lastFU$delta <-
  annt.kinetic
Clinical_info_lastFU$delta <- ifelse(Clinical_info_lastFU$delta == "No clearance",
                                     "No clearance", "Partial/complete clearance")

MeDIP_ctDNA_mod <- c()
for(i in 1:nrow(Clinical_info_lastFU)){
  tmp <- Clinical_info_lastFU$MeDIP_ctDNA[i]
  if(tmp < MeDIP_ctDNA_LLOD){
    MeDIP_ctDNA_mod <- c(MeDIP_ctDNA_mod, 0)
  } else {
    MeDIP_ctDNA_mod <- c(MeDIP_ctDNA_mod, tmp)
  }
}

Clinical_info_lastFU$MeDIP_ctDNA <- MeDIP_ctDNA_mod

BL_lastFU_ggplot <- Clinical_info_lastFU %>%
  ggplot(aes(x = Timepoint, y = MeDIP_ctDNA, group = HN)) + 
  geom_line(size = 0.26, aes(col = delta)) + 
  geom_hline(yintercept = MeDIP_ctDNA_LLOD, lty = 2, col = "blue", size = 0.39) + 
  scale_x_discrete(expand = c(0,0.1)) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(0.05, 10),
                     breaks = c(0, seq(0.2, 1, length.out = 9),
                                seq(2, 10, length.out = 9)),
                     labels = c("N.D", "0.2",
                                rep("",7), "1", rep("", 8),
                                "10"),
                     limit = c(0, 10)) + 
  scale_color_manual(name = "Kinetic:", values = c(palette_JB[c(1)], rgb(0,0,0,1/10))) + 
  labs(y = "Methylated ctDNA (%)") + 
  theme_JB() + 
  theme(axis.title.x = element_blank())

BL_lastFU_ggplot

# 10. RFS analysis ------------------------------------------------------------

## Create data.frame for Kaplan-meier
KM_df <- as.data.frame(matrix(nrow = 18, ncol=0))
tmp_mat <- HN_list

futime <- c()
fustat <- c()
for(i in 1:length(tmp_mat)){
  if(is.na(tmp_mat[[i]][[2]]['relapse'])){
    futime <- c(futime, tmp_mat[[i]][[2]]['lastFU'])
    fustat <- c(fustat, 0)
  } else {
    futime <- c(futime, tmp_mat[[i]][[2]]['relapse'])
    fustat <- c(fustat, 1)
  }
}

relapsePos_HNs <- 
  Clinical_info_lastFU[Clinical_info_lastFU$delta == "No clearance",]$HN

KM_df$futime <- futime
KM_df$fustat <- fustat
rownames(KM_df) <- names(tmp_mat)
KM_df$rx <- ifelse(rownames(KM_df) %in%
                     relapsePos_HNs,
                   "A", "B")

surv_object <- Surv(time = KM_df$futime, event = KM_df$fustat)
fit1 <- survfit(surv_object ~ rx, data = KM_df)
KM_df$rx <- factor(KM_df$rx, levels = c("B", "A"))
rx_coxph <- coxph(surv_object ~ rx, data = KM_df)
rx_HR <- summary(rx_coxph)$coefficient[2]

ggsurv_postTx <- ggsurvplot(fit1, data = KM_df, pval = TRUE, risk.table = TRUE,
                            tables.height = 0.3, ggtheme = theme_JBwBorder(), 
                            size = 0.26,palette = palette_JB[1:2], 
                            risk.table.y.text = FALSE,
                            pval.size = (5/14) * 9,
                            legend.title = "Kinetics", 
                            legend.labs = c("No clearance", "clearance"),
                            risk.table.fontsize = 4, conf.int = FALSE, 
                            risk.table.col = "strata",
                            xlab = "Time (days)", censor.shape = 124, 
                            censor.size = 2,
                            ylab = "Recurrence-free survival",
                            pval.coord = c(100, 0.05),
                            xscale = "d_m", break.time.by = 365.25) 
ggsurv_postTx <- egg::ggarrange(ggsurv_postTx$plot, 
                                ncol = 1, nrow = 1,
                                draw = FALSE)

# Compare day to recurrence across patients based on clearance / no clearance
tru_pos <- HN_list[names(HN_list) %in% rownames(KM_df[KM_df$fustat == 1 & 
                                                      KM_df$rx == "A",])]

fal_neg <- HN_list[names(HN_list) %in% rownames(KM_df[KM_df$fustat == 1 & 
                                                        KM_df$rx == "B",])]

## Get range of last collection to relapse
tru_pos_lc2rel <- c()
for(i in 1:length(tru_pos)){
  tmp <- tru_pos[[i]]
  lc <- tmp[[1]]$time[which(tmp[[1]]$time == max(tmp[[1]]$time))]
  rel <- tmp[[2]][4]
  lc2rel <- rel - lc
  tru_pos_lc2rel <- c(tru_pos_lc2rel, lc2rel)
  names(tru_pos_lc2rel)[i] <- names(tru_pos)[i]
}

### Days to Months
tru_pos_lc2rel <- tru_pos_lc2rel / 365.25 * 12

fal_neg_lc2rel <- c()
for(i in 1:length(fal_neg)){
  tmp <- fal_neg[[i]]
  lc <- tmp[[1]]$time[which(tmp[[1]]$time == max(tmp[[1]]$time))]
  rel <- tmp[[2]][4]
  lc2rel <- rel - lc
  fal_neg_lc2rel <- c(fal_neg_lc2rel, lc2rel)
  names(fal_neg_lc2rel)[i] <- names(fal_neg)[i]
}

fal_neg_lc2rel <- fal_neg_lc2rel / 365.25 * 12

t.test(x = tru_pos_lc2rel, y = fal_neg_lc2rel[-1])
