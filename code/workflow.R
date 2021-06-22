###############################################################################
##
## Title: Workflow
##
## Author: Justin Matthew Burgener
##
###############################################################################
##
## Description: Summary code that sources individual scripts used for 
##              subsequent generation of data. Subsequent code is sourced
##              according to order of data presented in manuscript
##
###############################################################################

# Load libraries and set directories ------------------------------------------

scripts <- "../code"
data_dir <- "../data"
out_dir <- "../results"
fig_dir <- "../figures"

# Libraries

library(ggplot2)
library(dplyr)
library(cowplot)
source(paste0(scripts, "/ggplot_themeJB.R"))

# 1. Isolation cfDNA yields from plasma samples -------------------------------

# load table w/ sample information for HNSCC and healthy donor plasma:
# $ HNC - identifier for each patient and healthy donor
# $ Timepoint - time of collection for HNC patients according to study 
#               definitions, or label of healthy donors (i.e. "Normal")
# $ dsDNApermlPlasma - total yield of cfDNA normalized to volume of plasma
# $ Timepoint2 - updated labels
# $ Final - Y/N on whether sample is final collection for HNC patient

dsDNA_logs <- read.table(file = paste0(data_dir, "/hnc_DNAperPlasma.txt"),
                         sep = "\t", header = TRUE)
dsDNA_logs$Timepoint2 <- factor(dsDNA_logs$Timepoint2, 
                                levels = c("Diagnosis", "Post-surgery",
                                           "Mid-radiotherapy",
                                           "Post-treatment-1", 
                                           "Post-treatment-2", "Normal"))

# generate plot

set.seed(42)
dsDNA_logs_plot <- dsDNA_logs %>%
  ggplot(aes(x = Timepoint2, y = dsDNApermLPlasma)) + 
  geom_jitter(position = position_jitter(),
              aes(col = Timepoint2), size = 1.5) + 
  geom_boxplot(fill = NA, fatten = 6, lwd = 0.13, outlier.shape = NA) + 
  theme_JB() + 
  scale_color_manual(values = c(palette_JB[1:6][-2], palette_JB[2])) + 
  theme(legend.position = "none", 
        panel.grid.major.y = element_line(size = 0.13,
                                          colour = rgb(0,0,0,1/20)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) + 
  labs(y = "ng cfDNA per mL plasma") +
  scale_y_continuous(breaks = seq(0, 70, by = 10))

# 2. Analysis of clonal hematopoiesis -----------------------------------------
current_ls <- ls()
current_ls <- c(current_ls, "current_ls")

source(paste0(scripts, "/clonal_hemato.R"))

current_ls <- c(current_ls, "cfDNA_vs_PBLs_SNVs_filtered", "CHIP_oncoprint_plot",
                "before_after_PBLfilter_plot", "before_after_PBLfilter_plot",
                "cfDNA_DepthFiltered_coding_SNVs")
rm(list = grep(paste0(current_ls, collapse = "|"), ls(),
               value = TRUE, invert = TRUE))
   
before_after_PBLfilter_plot <- before_after_PBLfilter_plot +
  theme(plot.margin = unit(c(5,5,5,5), "mm")) 

# 3. CAPP-Seq SNV analysis ----------------------------------------------------

source(paste0(scripts, "/CAPPseq_SNV_ident.R"))

current_ls <- c(current_ls, "FullOncoprint_plot", "SmallOncoprint_plot",
                "MAFperPatient_plot")
rm(list = grep(paste0(current_ls, collapse = "|"), ls(),
               value = TRUE, invert = TRUE))

# 4. Identification of PBL-depleted regions ------------------------------------

source(paste0(scripts, "/PBL_depleted_windows.R"))

current_ls <- c(current_ls, "FaDu_vs_Norm_nCpGMedianRPKM", "density_plot",
                "HNC_and_Norm_binMethyl_PBLdepletion_plot",
                "PBL_res_plot")

rm(list = grep(paste0(current_ls, collapse = "|"), ls(),
               value = TRUE, invert = TRUE))

# 5. cfMeDIP-seq differential methylation analysis ----------------------------

source(paste0(scripts, "/TCGA_site_specific_DMR.R"))
source(paste0(scripts, "/cfMeDIP_DMR.R"))

current_ls <- c(current_ls, "cfDNA_DMRs_plot", "hyperDMR_reducedDMRs_barplot",
                "CGI_DMRandEnrichment_plot", "CGI_DMRandhypo_enrichment_plot",
                "TCGA_DMRandhyper_enrichment_plot")

rm(list = grep(paste0(current_ls, collapse = "|"), ls(),
               value = TRUE, invert = TRUE))

# 6. cfMeDIP-seq and CAPP-Seq fragment length analysis ------------------------

source(paste0(scripts, "/fragment_length.R"))

current_ls <- c(current_ls, "SNV_medfrag_plot", "hyperDMR_medfrag_plot",
                "CAPP_vs_MeDIP_medfrag_plot", "short_enrc_plot",
                "MAF_vs_SNVmedfrag", "meanRPKM_vs_medfrag_plot")

rm(list = grep(paste0(current_ls, collapse = "|"), ls(),
               value = TRUE, invert = TRUE))

# 7. cfMeDIP-seq and CAPP-Seq concordance and prognosis analysis --------------

source(paste0(scripts, "/concordance.R"))

source(paste0(scripts, "/TCGA_analysis.R"))

current_ls <- c(current_ls, "hyperDMR_corplot", "ggsurv_BL_ctDNA", 
                "ggsurv_BL_meth_pos", "ggsurv_BL_both_pos", 
                "common_targets_plot", "cox_ggplot", "figure_STK3",
                "figure_OSR1", "figure_GATA2AS1", "figure_LINC01391",
                "figure_ZNF323", "sig_met_OS_wHR_KMplot", 
                "cfDNA_KM_wHR_plot", "ggsurv_BL_stage", "ExtendedFigure8",
                "hclust_capp_plot", "hclust_plot")

rm(list = grep(paste0(current_ls, collapse = "|"), ls(),
               value = TRUE, invert = TRUE))

# 8. Longitudinal analysis ----------------------------------------------------

source(paste0(scripts, "/longitudinal_analysis.R"))

current_ls <- c(current_ls, "hyperDMR_shortfrag_corplot", "AllPlots",
                "kinetics_ex_wlegend", "BL_lastFU_ggplot", 
                "ggsurv_postTx")

rm(list = grep(paste0(current_ls, collapse = "|"), ls(),
               value = TRUE, invert = TRUE))

# 10. Save and export plots for figure collation ---------------------------------------
## Save
save(dsDNA_logs_plot, cfDNA_vs_PBLs_SNVs_filtered, CHIP_oncoprint_plot,
     before_after_PBLfilter_plot, FullOncoprint_plot, SmallOncoprint_plot,
     MAFperPatient_plot, FaDu_vs_Norm_nCpGMedianRPKM, density_plot,
     HNC_and_Norm_binMethyl_PBLdepletion_plot, PBL_res_plot, cfDNA_DMRs_plot,
     hyperDMR_reducedDMRs_barplot, CGI_DMRandEnrichment_plot,
     CGI_DMRandhypo_enrichment_plot,
     hyperDMR_corplot, SNV_medfrag_plot, hyperDMR_medfrag_plot,
     CAPP_vs_MeDIP_medfrag_plot, short_enrc_plot,
     MAF_vs_SNVmedfrag, meanRPKM_vs_medfrag_plot, 
     ggsurv_BL_meth_pos, ggsurv_BL_ctDNA, ggsurv_BL_both_pos,
     common_targets_plot,
     cox_ggplot, figure_STK3, figure_OSR1, figure_GATA2AS1, figure_LINC01391,
     figure_ZNF323, sig_met_OS_wHR_KMplot, cfDNA_KM_wHR_plot,
     hyperDMR_shortfrag_corplot, AllPlots, kinetics_ex_wlegend,
     BL_lastFU_ggplot, ggsurv_postTx, ExtendedFigure8, hclust_capp_plot,
     hclust_plot, TCGA_DMRandhyper_enrichment_plot,
     file = paste0(fig_dir, "/all_manuscript_figures.R"))

## Load for future analysis
# load(paste0(fig_dir, "/all_manuscript_figures.R"))

# Export individual figures

## Main Figures

# Figure 2b
ggsave(plot = cfDNA_vs_PBLs_SNVs_filtered,
       height = 2.8, width = 2.8, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_2b.pdf"),
       useDingbats=FALSE)

# Figure 2c
ggsave(plot = CHIP_oncoprint_plot, height = 2.8, width = 2.4, dpi = 300,
       filename=paste0(fig_dir, "/figure_2c.pdf"),
       useDingbats=FALSE)

# Figure 2d
ggsave(plot = before_after_PBLfilter_plot +
         theme(plot.margin = unit(c(0,0,0,0), "mm")),
       height = 3.54, width = 4.8, dpi = 300,
       filename=paste0(fig_dir,
                  "/figure_2d.pdf"),
       useDingbats=FALSE)

# Figure 2e
ggsave(plot = SmallOncoprint_plot, height = 2.8, width = 3.2, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_2e.pdf"),
       useDingbats=FALSE)

# Figure 2f
ggsave(plot = MAFperPatient_plot, height = 3.54, width = 3.54, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_2f.pdf"),
       useDingbats=FALSE)

# Figure 3b*
ggsave(plot = density_plot, height = 2.8, width = 2.8, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_3b.pdf"),
       useDingbats=FALSE)

# Figure 3c
ggsave(plot = HNC_and_Norm_binMethyl_PBLdepletion_plot +
         theme(axis.text = element_text(size = 8)),
       height = 2.8, width = 5.6, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_3c.pdf"),
       useDingbats=FALSE)

# Figure 3d
ggsave(plot = cfDNA_DMRs_plot + theme(axis.text = element_text(size = 8),
                                      panel.background = element_blank()),
       height = 2.8, width = 2.8, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_3d.pdf"),
       useDingbats=FALSE)

# Figure 3e
ggsave(plot = as_ggplot(CGI_DMRandEnrichment_plot),
       height = 2.8, width = 1.4, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_3e.pdf"),
       useDingbats=FALSE)

# Figure 3f
ggsave(plot = as_ggplot(TCGA_DMRandhyper_enrichment_plot),
       height = 2.8, width = 1.4, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_3f.pdf"),
       useDingbats=FALSE)

# Figure 4a
ggsave(plot = SNV_medfrag_plot,
       height = 2.8, width = 4, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_4a.pdf"),
       useDingbats=FALSE)

# Figure 4b
ggsave(plot = hyperDMR_medfrag_plot,
       height = 2.8, width = 4, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_4b.pdf"),
       useDingbats=FALSE)

# Figure 4c
ggsave(plot = CAPP_vs_MeDIP_medfrag_plot,
       height = 3.2, width = 3.2, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_4c.pdf"),
       useDingbats=FALSE)

# Figure 4d
ggsave(plot = short_enrc_plot,
       height = 3.2, width = 4.8, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_4d.pdf"),
       useDingbats=FALSE)

# Figure 5a
ggsave(plot = hyperDMR_corplot,
       height = 3.54, width = 3.54, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_5a.pdf"),
       useDingbats=FALSE)

# Figure 5b
ggsave(plot = as_ggplot(ggsurv_BL_both_pos),
       height = 3.54, width = 3.54, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_5b.pdf"),
       useDingbats=FALSE)


# Figure 5d
ggsave(plot = cox_ggplot,
       height = 3.54, width = 2.4, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_5d.pdf"),
       useDingbats=FALSE)

# Figure 5e
ggsave(plot = sig_met_OS_wHR_KMplot,
       height = 3.2, width = 3.2, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_5e.pdf"),
       useDingbats=FALSE)

# Figure 5f
ggsave(plot = cfDNA_KM_wHR_plot,
       height = 3.2, width = 3.2, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_5f.pdf"),
       useDingbats=FALSE)

# Figure 6a
ggsave(plot = kinetics_ex_wlegend,
       height = 2.8, width = 8, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_6a.pdf"),
       useDingbats=FALSE)

# Figure 6b
ggsave(plot = BL_lastFU_ggplot,
       height = 3.54, width = 3.54, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_6b.pdf"),
       useDingbats=FALSE)

# Figure 6c
ggsave(plot = as_ggplot(ggsurv_postTx),
       height = 3.54, width = 3.54, dpi = 300,
       filename=paste0(fig_dir,
                       "/figure_6c.pdf"),
       useDingbats=FALSE)

### Supplementary Figures

# Supplementary Figure 1b
ggsave(plot = dsDNA_logs_plot, height = 3.54, width = 8, dpi = 300,
       filename=paste0(fig_dir, "/sup_figure_1b.pdf"),
       useDingbats=FALSE)

# Supplementary Figure 2 not provided

# Supplementary Figure 3
ggsave(plot = FullOncoprint_plot, height = 7.8, width = 7.2, dpi = 300,
       filename=paste0(fig_dir,
                       "/sup_figure_3.pdf"),
       useDingbats=FALSE)

# Supplementary Figure 4A
ggsave(plot = as_ggplot(FaDu_vs_Norm_nCpGMedianRPKM), 
       height = 3.54, width = 3.54,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_4a.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 4B
ggsave(plot = PBL_res_plot, height = 3.54, width = 3.54,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_4b.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 5A
ggsave(plot = hyperDMR_reducedDMRs_barplot, height = 3.54, width = 3.54,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_5a.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 5B
ggsave(plot = as_ggplot(CGI_DMRandhypo_enrichment_plot), 
       height = 3.54, width = 3.54,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_5b.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 6 (not provided)

# Supplementary Figure 7A
ggsave(plot = MAF_vs_SNVmedfrag, height = 3.54, width = 3.54,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_7a.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 7B*
ggsave(plot = meanRPKM_vs_medfrag_plot, height = 3.54, width = 3.54,
       dpi = 300, filename = paste0(fig_dir,
                                  "/sup_figure_7b.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 8
ggsave(plot = as_ggplot(ExtendedFigure8),
       height = 7, width = 8,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_fig_8.pdf"),
       useDingbats = FALSE)


# Supplementary Figure 9a
ggsave(plot = common_targets_plot, height = 3.54, width = 3.54,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_9a.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 9d
ggsave(plot = as_ggplot(figure_STK3), height = 2.8, width = 2.8,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_9d.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 9e
ggsave(plot = as_ggplot(figure_OSR1), height = 2.8, width = 2.8,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_9e.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 9b
ggsave(plot = as_ggplot(figure_GATA2AS1), height = 2.8, width = 2.8,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_9b.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 9f
ggsave(plot = as_ggplot(figure_LINC01391), height = 2.8, width = 2.8,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_9f.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 9c
ggsave(plot = as_ggplot(figure_ZNF323), height = 2.8, width = 2.8,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_9c.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 9g
ggsave(plot = hclust_plot, height = 2.8, width = 2.8,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_9g.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 9h
ggsave(plot = hclust_capp_plot, height = 2.8, width = 2.8,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_9h.pdf"),
       useDingbats = FALSE)

# Supplementary Figure 10
ggsave(plot = as_ggplot(AllPlots), height = 10, width = 8.5,
       dpi = 300, filename = paste0(fig_dir,
                                    "/sup_figure_10.pdf"),
       useDingbats = FALSE)

# Session Info
sessionInfo()
