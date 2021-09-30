## Analysis of interneuron maturation
##------ Tue Sep 28 18:43:43 2021 ------##
## Oliver Eichmueller
library(zoo)
library(dplyr)
library(monocle3)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(grid)

# wrapper for heatmap plotting ----------

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

# set Output Directory ---------------------------------------------------------
OutDir <- 'path_to/OutDir'

# plot statistics for all organoid integration ---------------------------------

stat_all_orgs <- readRDS(file = "path_to/statistics_AllOrgs.Rds")

# filter for only d110 organoids
stat_all_orgs_filter <- stat_all_orgs %>% 
  filter(media_genotype != "Old.Tumor", clusters_low %in% c(6,8,9,16)) %>%
  mutate(clusters_low = factor(clusters_low),
         media = factor(as.vector(media)))

# create media annotation
media_label <- c("High Nutrient d110", "Low Nutrient d110")
names(media_label) <- c("HN", "LN")


pl1 <- ggplot(stat_all_orgs_filter
              , aes(x = clusters_low, y = relative, fill = media_genotype)) +
  ggnewscale::new_scale_fill()  +
  geom_col( aes(x = clusters_low, y = relative, fill = media_genotype), 
            position = "dodge", width = 0.8) +
  geom_text(data = stat_all_orgs_filter[
    stat_all_orgs_filter$media_genotype %in% c("HN.TSC2", "LN.TSC2"),]
            , aes(x = clusters_low, y = relative +2, label = round(relative)), 
    color = "black", nudge_x = .2) +
  geom_text(data = stat_all_orgs_filter[
    stat_all_orgs_filter$media_genotype %in% c("HN.Ctrl", "LN.Ctrl"),]
            , aes(x = clusters_low, y = relative +2, label = round(relative)), 
    color = "black", nudge_x = -.2) +
  scale_fill_manual("Dataset", 
                    values = c(rgb(102, 153, 102, maxColorValue = 255), 
                               rgb(51, 102, 153, maxColorValue = 255), 
                               rgb(102, 153, 102, maxColorValue = 255), 
                               rgb(51, 102, 153, maxColorValue = 255), 
                               rgb(51, 102, 153, maxColorValue = 255))) +
  theme_light() +
  theme(strip.text.y = element_text(face = "bold", size = 10, color = "black")
        , strip.background.y = element_rect(fill = 
                                              c(alpha("lightgrey", alpha = .5), 
                                                alpha("blue", alpha = .5), 
                                                alpha("green", alpha = .5)))) +
  xlab("Clusters") + ylab("Percentage of Dataset")+
  scale_y_continuous(limits = c(0,15))+
  facet_grid(rows = vars(media), labeller = labeller(media = media_label)) +
  ggtitle("Percentage of Interneuron Clusters of whole dataset")

ggsave(paste0(OutDir, '/Perc_IN_AllClustering.pdf'),
       device = "pdf", plot = pl1)



# load dataset of all organoids with pseudotime --------------------------------

cds.integration <- readRDS(file = 'path_to/tsc_paper_integration_all_pseudotime.Rds')

# choose immature to mature IN
clip_cds_IN <- monocle3::choose_graph_segments(cds.integration)
clip_cds_IN <- clip_cds[,clip_cds_IN@colData$barcode]

plot_cells(clip_cds_IN, color_cells_by = "media_genotype", show_trajectory_graph = F
           , label_groups_by_cluster = F, labels_per_group = 0)

# downsample barcodes based on smallest d110 dataset
barcodes_downsampled <- cds.integration@colData %>% as.data.frame() %>%
  filter(media_genotype != "Old.Tumor") %>%
  dplyr::group_by(media_genotype) %>%
  dplyr::sample_n(4816) %>%
  magrittr::use_series(barcode)

# subset IN lineage dataset for downsampled barcodes
barcodes_downsampled_found <- intersect(barcodes_downsampled, 
                                        colnames(clip_cds_IN))

clip_cds_IN_downsampled <- clip_cds_IN[,barcodes_downsampled_found]
clip_cds_IN_downsampled@colData$orig.ident <- 
  as.vector(clip_cds_IN_downsampled@colData$orig.ident)

# re-calculate UMAP
clip_cds_IN_downsampled <- reduce_dimension(clip_cds_IN_downsampled, 
                                            max_components = 2)

# re-cluster IN lineage dataset
clip_cds_IN_downsampled <- cluster_cells(clip_cds_IN_downsampled, k = 8)
clip_cds_IN_downsampled@colData$clusters_new <- clip_cds_IN_downsampled@clusters$UMAP$clusters

# plot UMAPs
plot_cells(clip_cds_IN_downsampled, color_cells_by = "cluster", cell_size = 2,
           alpha = .5, show_trajectory_graph = FALSE, label_cell_groups = FALSE)
ggsave(paste0(OutDir, '/UMAP_IN_SubClustering.png'),
       device = "png")

plot_cells(clip_cds_IN_downsampled, color_cells_by = "media_genotype", cell_size = 2,
           alpha = .5, show_trajectory_graph = FALSE, label_cell_groups = FALSE)+   
  scale_color_manual(values = pals::brewer.set1(4))
ggsave(paste0(OutDir, 'UMAP_IN__media_genotype_SubClustering.png'),
       device = "png")

# learn graph
clip_cds_IN_downsampled <- 
  learn_graph(clip_cds_IN_downsampled, use_partition = F,
              learn_graph_control = list(minimal_branch_len = 12,
                                         rann.k = NULL,
                                         orthogonal_proj_tip = FALSE,
                                         geodesic_distance_ratio = 1/3,
                                         euclidean_distance_ratio = 1))

# order cells in pseudotime
clip_cds_IN_downsampled <- order_cells(clip_cds_IN_downsampled)

# UMAP of pseudotime
plot_cells(clip_cds_IN_downsampled, color_cells_by = "pseudotime",
           cell_size = 2, show_trajectory_graph = T)
ggsave(paste0(OutDir, 'UMAP_IN_Pseudotime_SubClustering.png'),
       device = "png")

# subset only tumor to tuber IN
clip_cds_IN_downsampled_subset <- choose_cells(clip_cds_IN_downsampled)

# graph test on tumor to tuber IN
graph_test_newIN_subset <- graph_test(clip_cds_IN_downsampled_subset)

write.csv(graph_test_newIN_subset, file = paste0(OutDir, 'graph_test_pseudotime_IN.csv'))

# filter gois
goi <- graph_test_newIN_subset %>% 
  filter(status %in% "OK", q_value <1e-20|q_value==0) %>%
  magrittr::use_series(gene_short_name)

# Order cells from LN to HN cluster
clip_cds_IN_downsampled_subset_neworder <- 
  order_cells(clip_cds_IN_downsampled_subset)

# UMAPs of re-ordered subset
plot_cells(clip_cds_IN_downsampled_subset_neworder, 
           color_cells_by = "pseudotime",
           show_trajectory_graph = F, cell_size = 2) 
ggsave(paste0(OutDir, 'UMAP_INsubset_pseudotime_SubClustering.png'),
       device = "png")

plot_cells(clip_cds_IN_downsampled_subset_neworder, color_cells_by = "cluster",
           show_trajectory_graph = F, cell_size = 2) 
ggsave(paste0(OutDir, 'UMAP_INsubset_cluster_SubClustering.png'),
       device = "png")

# bin along pseudotime
pseudotime_bin <- pseudotime(clip_cds_IN_downsampled_subset_neworder) %>% 
  data.frame(pseudotime_bin = .) %>%
  mutate(barcode = row.names(.), pseudotime_bin = floor(pseudotime_bin)) %>%
  select(barcode, pseudotime_bin)

# add pseudotime bin to coldata
clip_cds_IN_downsampled_subset_neworder@colData <- 
  clip_cds_IN_downsampled_subset_neworder@colData  %>% as.data.frame() %>%
  left_join(pseudotime_bin, by = "barcode") %>%
  DataFrame(row.names = .$barcode)

# convert pseudotime bin to factor
clip_cds_IN_downsampled_subset_neworder@colData <- 
  clip_cds_IN_downsampled_subset_neworder@colData %>% as.data.frame() %>%
  select(barcode, Size_Factor, orig.ident, clusters_low, media_genotype, 
         clusters_new, pseudotime_bin) %>%
  mutate(pseudotime_bin = factor(pseudotime_bin)) %>%
  DataFrame(row.names = .$barcode)

# perform sliding average and plot from tumor to tuber IN
# set window and step
window <- 2
step <- 1

# aggregate gois of graph test per pseudotime bin
dat_IN_subset <- aggregate_gene_expression(clip_cds_IN_downsampled_subset_neworder[goi,]
                                           , cell_group_df = pseudotime_bin
                                           , scale_agg_values = T
                                           , max_agg_value = 1)


dat_IN_subset_backup <- dat_IN_subset

dat_IN_subset <- as.matrix(dat_IN_subset_backup)

# order using sliding average
dat_IN_subset_ordered <- 
  dat_IN_subset[order(apply(t(rollapply(t(dat_IN_subset), width=window, by=step, FUN=mean)), 1, which.max)), ]

# generate annotation for groups
groupings <- 
  read.delim('path_to/Overrepresentation.clusterProfiler3.18.1.graph_test.splitByConsecutiveMaxc.gene2set.tsv')

groupings$gene_id <- stringr::str_replace(groupings$gene_id, pattern = "-", "\\.")
groupings_sub <- groupings[1:360,]
row.names(groupings_sub) <- groupings_sub$gene_id
groupings_sub <- groupings_sub %>% select(-gene_id)
groupings_color <- groupings_sub
groupings_color <- list(set = c(`1` = pals::brewer.set2(3)[1],
                                `2` = pals::brewer.set2(3)[2],
                                `3` = pals::brewer.set2(3)[3]))

# select genes for annotation
goi_plot <- c("GRIA1", "GRIA2", "DPP10", "GABRA2", "GRIP1", "ARX", "STMN1", "SNAP25", "GRIN2B",
              "GABARAPL2", "RPL22", "RPS8", "RPS23", "LAPTM4B")

# plot heatmap
ph1 <- pheatmap(dat_IN_subset_ordered
                , cluster_rows = F, cluster_cols = F
                , color = rev(pals::brewer.rdbu(100))
                , annotation_row = groupings_sub
                , annotation_colors = groupings_color
                , scale = "row", show_rownames = T, fontsize_row = 15, angle_col = 0, fontsize_col = 20)

# annotated selected genes
pl1 <- add.flag(ph1, goi_plot, repel.degree = .1)

# save plot
ggsave(paste0(OutDir, 'Heatmap_pseudotime_ordered.pdf'),
       device = "pdf", plot = pl1)
ggsave(paste0(OutDir, 'Heatmap_pseudotime_ordered.png'),
       device = "png", plot = pl1)


# calculate statistics and perform plotting
clip_cds_IN_downsampled_stat <- clip_cds_IN_downsampled@colData %>%
  as.data.frame() %>%
  dplyr::group_by(media_genotype) %>%
  dplyr:: count(clusters_new) %>%
  mutate(relative = round(n/4816*100, 2),
         media = stringr::str_extract(media_genotype, "LN|HN"))


pl1 <- ggplot(clip_cds_IN_downsampled_stat
              , aes(x = clusters_new, y = relative, fill = media_genotype)) +
  ggnewscale::new_scale_fill()  +
  geom_col( aes(x = clusters_new, y = relative, fill = media_genotype), 
            position = "dodge", width = 0.8) +
  geom_text(data = clip_cds_IN_downsampled_stat[
    clip_cds_IN_downsampled_stat$media_genotype %in% c("HN.TSC2", "LN.TSC2"),]
            , aes(x = clusters_new, y = relative +2, label = relative), 
    color = "black", nudge_x = .2) +
  geom_text(data = clip_cds_IN_downsampled_stat[
    clip_cds_IN_downsampled_stat$media_genotype %in% c("HN.Ctrl", "LN.Ctrl"),]
            , aes(x = clusters_new, y = relative +2, label = relative), 
    color = "black", nudge_x = -.2) +
  scale_fill_manual("Dataset", 
                    values = c(rgb(102, 153, 102, maxColorValue = 255), 
                               rgb(51, 102, 153, maxColorValue = 255), 
                               rgb(102, 153, 102, maxColorValue = 255), 
                               rgb(51, 102, 153, maxColorValue = 255), 
                               rgb(51, 102, 153, maxColorValue = 255))) +
  theme_light() +
  theme(strip.text.y = element_text(face = "bold", size = 10, color = "black")
        , strip.background.y = 
          element_rect(fill = c(alpha("lightgrey", alpha = .5), 
                                alpha("blue", alpha = .5), 
                                alpha("green", alpha = .5)))) +
  xlab("Clusters") + ylab("Percentage of Dataset")+
  scale_y_continuous(limits = c(0,15))+
  facet_grid(rows = vars(media), labeller = labeller(media = media_label)) +
  ggtitle("Percentage of Interneuron Subclusters of whole dataset")


ggsave(paste0(OutDir, 'Perc_IN_SubClustering.pdf'),
       device = "pdf", plot = pl1)

write.csv(clip_cds_IN_downsampled_stat, 
          paste0(OutDir, 'clip_cds_IN_downsampled_stat.csv'))

# plot UMAPs for individual genes
pl1 <- plot_cells(clip_cds_IN_downsampled, genes = "LAPTM4B",
                  show_trajectory_graph = F, cell_size = 2, norm_method = "size_only") +
  facet_wrap(~media_genotype)
ggsave(paste0(OutDir, 'UMAP_LAPTM4B_IN_downsampled.png'),
       device = "png", plot = pl1)

pl1 <- plot_cells(clip_cds_IN_downsampled, genes = "GABARAPL2",
                  show_trajectory_graph = F, cell_size = 2, norm_method = "size_only") +
  facet_wrap(~media_genotype)
ggsave(paste0(OutDir, 'UMAP_GABARAPL2_IN_downsampled.png'),
       device = "png", plot = pl1)

pl1 <- plot_cells(clip_cds_IN_downsampled, genes = "RPL22",
                  show_trajectory_graph = F, cell_size = 2, norm_method = "size_only") +
  facet_wrap(~media_genotype)
ggsave(paste0(OutDir, 'UMAP_RPL22_IN_downsampled.png'),
       device = "png", plot = pl1)

pl1 <- plot_cells(clip_cds_IN_downsampled, genes = "RPS10",
                  show_trajectory_graph = F, cell_size = 2, norm_method = "size_only") +
  facet_wrap(~media_genotype)
ggsave(paste0(OutDir, 'UMAP_RPS10_IN_downsampled.png'),
       device = "png", plot = pl1)

# save files and objects
saveRDS(clip_cds_IN, file = paste0(OutDir, "/CLIP_cds_IN.Rds"))
saveRDS(clip_cds_IN_downsampled, 
        file = paste0(OutDir, "CLIP_cds_in_downsampled.Rds"))
saveRDS(clip_cds_IN_downsampled_subset, 
        file = paste0(OutDir, "CLIP_cds_in_downsampled_subset.Rds"))
saveRDS(clip_cds_IN_downsampled_subset_neworder, 
        file = paste0(OutDir, "CLIP_cds_in_downsampled_subset_neworder.Rds"))
write.csv(dat_IN_subset_ordered %>% as.data.frame()%>% 
            mutate(gene = row.names(.)), 
          paste0(OutDir, "dat_IN_subset_ordered.csv"))
