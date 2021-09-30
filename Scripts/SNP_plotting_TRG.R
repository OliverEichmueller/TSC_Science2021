## plotting SNPs of TSC2 TRG
##------ Wed Sep 29 08:31:34 2021 ------##
## Oliver Eichmueller

# load in libraries ------------------------------------------------------------
library(tidyr)
library(dplyr)
library(janitor)
library(reshape2)
library(ggplot2)
library(ggpubr)

# Set OutPut Directory ---------------------------------------------------------
OutDir <- 'path_to/OutDir'

# load in SNP tables ------------------------------------------------------------
control6  <- read.delim("path_to_SNPs/6/control.alleleCount")
tumor6    <- read.delim("path_to_SNPs/6/tumor.alleleCount")
control7  <- read.delim("path_to_SNPs/7/control.alleleCount")
tumor7    <- read.delim("path_to_SNPs/7/tumor.alleleCount")
control8  <- read.delim("path_to_SNPs/8/control.alleleCount")
tumor8    <- read.delim("path_to_SNPs/8/tumor.alleleCount")
control11 <- read.delim("path_to_SNPs/11/control.alleleCount")
tumor11   <- read.delim("path_to_SNPs/11/tumor.alleleCount")
control12 <- read.delim("path_to_SNPs/12/control.alleleCount")
tumor12   <- read.delim("path_to_SNPs/12/tumor.alleleCount")
control13 <- read.delim("path_to_SNPs/13/control.alleleCount")
tumor13   <- read.delim("path_to_SNPs/13/tumor.alleleCount")

# load in SNP reference tables -------------------------------------------------
lovd_snps_tsc2 <-  
  read.delim("path_to_SNPs/lovd/lovd.nl.TSC2.hg38.bed", header = F)
lovd_snps_tsc1 <-  
  read.delim("path_to_SNPs/lovd/lovd.nl.TSC1.hg38.bed", header = F)

lovd_snps_tsc2_filt <- lovd_snps_tsc2 %>%  mutate(database = "lovd")
lovd_snps_tsc1_filt <- lovd_snps_tsc1 %>% mutate(database = "lovd")

lovd_snps_combined <- bind_rows(lovd_snps_tsc1_filt, lovd_snps_tsc2_filt)

# define function to check single position -------------------------------------
check_snp_pos <- function(control_data, tumor_data, position){
  control_data[,-1] %>% melt(id.vars = c(1,6)) %>% mutate(genotype = "control") %>%
    mutate(rel.C = value/ Good_depth, POS_base = paste(POS, variable, sep = "_")) %>% filter(rel.C >0) %>%
    full_join(tumor_data[,-1] %>% melt(id.vars = c(1,6)) %>% mutate(genotype = "tumor") %>% mutate(rel.T = value/ Good_depth, POS_base = paste(POS, variable, sep = "_")) %>% filter(rel.T >0), 
              by = "POS_base") %>%
    left_join(lovd_snps_combined %>% select(V2, database), by = c("POS.x" = "V2")) %>% replace_na(replace = list(database = "clinivar")) %>% filter(POS.x == position)
  
}

# define function to filter and plot SNPs --------------------------------------
filter_plot_snps_tsc1_2 <- function(control_data, tumor_data, dp.filter = 10){
  pl1 = 
    control_data[,-1] %>% melt(id.vars = c(1,6)) %>% 
    mutate(genotype = "control") %>%
    mutate(rel.C = value/ Good_depth, 
           POS_base = paste(POS, variable, sep = "_")) %>% 
    filter(rel.C >0) %>%
    full_join(tumor_data[,-1] %>% melt(id.vars = c(1,6)) %>% 
                mutate(genotype = "tumor") %>% 
                mutate(rel.T = value/ Good_depth, 
                       POS_base = paste(POS, variable, sep = "_")) %>% 
                filter(rel.T >0), 
              by = "POS_base") %>%
    full_join(lovd_snps_combined %>% 
                select(V2, database), by = c("POS.x" = "V2")) %>% 
    replace_na(replace = list(database = "clinivar")) %>% 
    mutate( database = case_when(POS.x == 2055441 ~ "Patient", 
                                 TRUE ~ database))%>%
    mutate(ratio = rel.T-rel.C, gene = case_when(POS.x < 5e7 ~ "TSC2", 
                                                 TRUE ~ "TSC1")) %>% 
    filter(Good_depth.x>dp.filter & Good_depth.y>dp.filter) %>%
    ggplot(aes(x = POS.x, y = ratio, color = database)) + 
    geom_point() + scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = .25, alpha = .5, color = "grey") + 
    geom_hline(yintercept = -.25, alpha = .5, color = "grey")+
    geom_hline(yintercept = .5, alpha = .75, color = "grey") + 
    geom_hline(yintercept = -.5, alpha = .75, color = "grey")+
    facet_grid(cols = vars(gene), scales = "free_x")+
    ylim(-1,1) + theme_pubclean()
  
  return(pl1)
}

# define function to filter SNPs -----------------------------------------------
filter_snps_tsc1_2 <- function(control_data, tumor_data, 
                               dp.filter = 10, tumor = "Tumor6"){
  pl1 = control_data[,-1] %>% melt(id.vars = c(1,6)) %>% 
    mutate(genotype = "control") %>%
    mutate(rel.C = value/ Good_depth, 
           POS_base = paste(POS, variable, sep = "_")) %>% 
    filter(rel.C >0) %>%
    full_join(tumor_data[,-1] %>% melt(id.vars = c(1,6)) %>% 
                mutate(genotype = "tumor") %>% 
                mutate(rel.T = value/ Good_depth, 
                       POS_base = paste(POS, variable, sep = "_")) %>% 
                filter(rel.T >0), 
              by = "POS_base") %>%
    full_join(lovd_snps_combined %>% 
                select(V2, database), by = c("POS.x" = "V2")) %>% 
    replace_na(replace = list(database = "clinivar")) %>% 
    mutate( database = case_when(POS.x == 2055441 ~ "Patient", 
                                 TRUE ~ database))%>%
    mutate(ratio = rel.T-rel.C, 
           gene = case_when(POS.x < 5e7 ~ "TSC2", 
                            TRUE ~ "TSC1"), tumor = tumor) %>% 
    filter(Good_depth.x>dp.filter & Good_depth.y>dp.filter)
  
  return(pl1)
}

# filter SNPs ------------------------------------------------------------------
tumor6_filter <- filter_snps_tsc1_2(control6, tumor6, tumor = "Tumor6")
tumor7_filter <- filter_snps_tsc1_2(control7, tumor7, tumor = "Tumor7")
tumor8_filter <- filter_snps_tsc1_2(control8, tumor8, tumor = "Tumor8")
tumor11_filter <- filter_snps_tsc1_2(control11, tumor11, tumor = "Tumor11")
tumor12_filter <- filter_snps_tsc1_2(control12, tumor12, tumor = "Tumor12")
tumor13_filter <- filter_snps_tsc1_2(control13, tumor13, tumor = "Tumor13")

# plot SNPs --------------------------------------------------------------------
pl1 <- 
  bind_rows(tumor6_filter, tumor11_filter, tumor12_filter, tumor13_filter) %>% 
  mutate(tumor = 
           factor(tumor, 
                  levels = c("Tumor12", "Tumor6", "Tumor11", "Tumor13"))) %>%
  ggplot(aes(x = POS.x, y = ratio, color = database)) + 
  geom_hline(yintercept = .25, alpha = .5, color = "grey") + 
  geom_hline(yintercept = -.25, alpha = .5, color = "grey")+
  geom_hline(yintercept = .5, alpha = .75, color = "grey") + 
  geom_hline(yintercept = -.5, alpha = .75, color = "grey")+
  geom_hline(yintercept = 1, alpha = 1, color = "black") + 
  geom_hline(yintercept = -1, alpha = 1, color = "black")+
  geom_point(size = .8, alpha = .5) + scale_color_brewer(palette = "Set1") +
  facet_grid(cols = vars(gene), rows = vars(tumor), scales = "free_x")+
  ylim(-1,1) + theme_pubclean()

# save plot --------------------------------------------------------------------
pdf(file = paste0(OutDir, "/snp_freq_TRG.pdf"), width = 3, height = 5)
pl1
dev.off()