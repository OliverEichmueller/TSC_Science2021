## Pseudotime analysis and plotting of organoid and fetal datasets
##------ Tue Sep 28 16:34:12 2021 ------##
## Oliver Eichmueller

# load required libraries ------------------------------------------------------
try (source ('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R'), silent= F)
require(MarkdownReports)# devtools::install_github(repo = "vertesy/MarkdownReports")
library(Seurat)
library(monocle3)
library(reticulate)
library(pals)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(zoo)

# wrapper for heatmap plotting -------------------------------------------------

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

# read in monocle object -------------------------------------------------------
# Integration of organoid data (d110 H- and L-medium, d220 tumors)
cds.integration <- readRDS("path_to/AllOrgs_plusDiss_monocle.Rds")

cds.integration@colData$media_genotype <- 
  factor(cds.integration@colData$media_genotype, 
         levels = rev(levels(cds.integration@colData$media_genotype)))

# learn trajectory graph
cds.integration <- learn_graph(cds = cds.integration, use_partition = T, 
                               verbose = T, close_loop = F, 
                               learn_graph_control = 
                                 list(minimal_branch_len = 12, 
                                      orthogonal_proj_tip = F, 
                                      rann.k = 200))
# subset for ventral cells only
cds.integration_ventral <- 
  cds.integration[,cds.integration@colData$clusters_low %in% 
                    c(5,6,7,8,9,11,13,15,16,17)]

plot_cells(cds.integration_ventral,
           color_cells_by = "clusters_low",
           show_trajectory_graph = T,
           group_label_size = 10, alpha = .3)

# order cells in pseudotime
cds.integration_ventral <- order_cells(cds.integration_ventral)

plot_cells(cds.integration_ventral,
           color_cells_by = "pseudotime",
           show_trajectory_graph = T)

# subset for OPC lineage of ventral cells
# use barcodes of subsetting to extract cells from original object
cds_ventral_OPC <- monocle3::choose_graph_segments(cds.integration_ventral)
cds_ventral_neurog <- monocle3::choose_graph_segments(cds.integration_ventral)

cds_ventral_subset_OPC <- 
  cds.integration_ventral[,colnames(cds_ventral_OPC)]
cds_ventral_subset_neurog <- 
  cds.integration_ventral[,colnames(cds_ventral_neurog)]

# perform graph testing on subsetted objects
graph_test_OPC <- graph_test(cds_ventral_subset_OPC)

graph_test_ventral_neurog <- graph_test(cds_ventral_subset_neurog)

# subset graph test on monans_I > .1
graph_test_ventral_neurog_subset <- graph_test_ventral_neurog %>% arrange(desc(morans_I)) %>%
  filter(morans_I > .1) %>%
  select(gene_short_name)

graph_test_OPC_subset <- graph_test_OPC %>% arrange(desc(morans_I)) %>%
  filter(morans_I > .1) %>%
  select(gene_short_name)

# bin along pseudotime, for ventral pseudotime val. are floored and halfed
binned_pseudotime_ventral <- 
  data.frame(floor(pseudotime(cds_ventral_subset_neurog)/2)) %>% 
  rename(bin=1) %>% mutate(barcode = row.names(.)) %>%
  select(barcode,bin)

# bin along pseudotime, for OPC pseudotime val. are floored and multipl. w/ 5
binned_pseudotime_OPC <- 
  data.frame(floor(pseudotime(cds_ventral_subset_OPC)*5)) %>% 
  rename(bin=1) %>% mutate(barcode = row.names(.)) %>%
  select(barcode,bin)

# plot OPC pseudotime on organoids ---------------------------------------------
# set window and step for sliding average
window <- 2
step <- 1

# aggregate expression using bins and genes from graph testing
dat_OPC <- aggregate_gene_expression(cds_ventral_subset_OPC[graph_test_OPC_subset$gene_short_name,]
                                     , cell_group_df = binned_pseudotime_OPC
                                     , scale_agg_values = T
                                     , max_agg_value = 1)

dat_OPC_backup <- dat_OPC

dat_OPC <- as.matrix(dat_OPC_backup)

# order genes using a sliding average
dat_OPC_ordered <- dat_OPC[order(apply(t(rollapply(t(dat_OPC), 
                                                   width=window, by=step, 
                                                   FUN=mean)), 1, which.max)), ]

# select genes of interest for plotting
goi <- c("EDNRB", "PTGDS", "PDGFRA", "PCDH15", "OLIG1", "OLIG2", "EGFR",
         "HOPX", "TNC", "VIM", "CD9")

# plot pseudotime heatmap
ph1 <- pheatmap(dat_OPC_ordered
                , cluster_rows = F, cluster_cols = F
                , color = rev(pals::brewer.rdbu(100))
                , scale = "row", show_rownames = T, 
                fontsize_row = 15, angle_col = 0, fontsize_col = 20)

# annotate only selected genes
add.flag(ph1, goi, repel.degree = .1)

# save heatmap
png(paste0(OutDir, "Pseudotime_Heatmap_OPC_Org.png"),
    width = 750, height = 500)
add.flag(ph1, goi, repel.degree = .1)
dev.off()

# subset for control and TSC mutant cells 
ctrl_cells_OPC <- cds_ventral_subset_OPC@colData$barcode[
  cds_ventral_subset_OPC@colData$media_genotype %in% c("LN.Ctrl", "HN.Ctrl")]
TSC_cells_OPC <- cds_ventral_subset_OPC@colData$barcode[
  !(cds_ventral_subset_OPC@colData$media_genotype %in% c("LN.Ctrl", "HN.Ctrl"))]

# plot genes along pseudotime for mutant
pl1 <- 
  plot_genes_in_pseudotime(cds_ventral_subset_OPC[
    c("EDNRB", "PTGDS", "PDGFRA","OLIG2", "EGFR"),TSC_cells_OPC]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("EDNRB", "PTGDS", "EGFR", 
                                                  "PDGFRA", "OLIG2")) + 
  scale_color_viridis_c(alpha = .4)

png(paste0(OutDir, "Pseudotime_genes_OPC_MUT_CLIP_Org.png"), width = 500, 
    height = 750)
pl1
dev.off()

# plot genes along pseudotime for control
pl1 <- 
  plot_genes_in_pseudotime(cds_ventral_subset_OPC[
    c("EDNRB", "PTGDS", "PDGFRA","OLIG2", "EGFR"),ctrl_cells_OPC]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("EDNRB", "PTGDS", "EGFR", 
                                                  "PDGFRA", "OLIG2")) + 
  scale_color_viridis_c(alpha = .4)

png(paste0(OutDir, "Pseudotime_genes_OPC_CTRL_CLIP_Org.png"), width = 500, 
    height = 750)
pl1
dev.off()

# save with alpha = 1 for scale bar
pl1 <- 
  plot_genes_in_pseudotime(cds_ventral_subset_OPC[
    c("EDNRB", "PTGDS", "PDGFRA","OLIG2", "EGFR"),ctrl_cells_OPC]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("EDNRB", "PTGDS", "EGFR", 
                                                  "PDGFRA", "OLIG2")) + 
  scale_color_viridis_c(alpha = .4)

png(paste0(OutDir, "Pseudotime_genes_OPC_scale_CLIP_Org.png"), width = 500, 
    height = 750)
pl1
dev.off()

# save table
write.csv(dat_OPC_ordered, file = paste0(OutDir, 
                                         "Pseudotime_Heatmap_OPC_Org.csv"))

# plot CLIP Neurogenic lineage pseudotime on organoids -------------------------
# set window and step for sliding average
window <- 3
step <- 1

# aggregate expression using bins and genes from graph testing
dat_neurog <- 
  aggregate_gene_expression(cds_ventral_subset_neurog
                            [graph_test_ventral_neurog_subset$gene_short_name,], 
                            cell_group_df = binned_pseudotime_ventral, 
                            scale_agg_values = F, max_agg_value = 1)

dat_neurog_backup <- dat_neurog

dat_neurog <- as.matrix(dat_neurog_backup)

# order genes using a sliding average
dat_neurog_ordered <- 
  dat_neurog[order(apply(t(rollapply(t(dat_neurog), width=window, by=step, 
                                     FUN=mean)), 1, which.max)), ]
# select genes of interest for plotting
goi <- c("EDNRB", "PTGDS", "SCGN", "EGFR", "DLX1", "DLX5", "SP8"
         , "DLX6.AS1", "ASCL1", "HOPX", "TNC", "VIM", "GFAP", "S100B")

# plot pseudotime heatmap
ph1 <- pheatmap(dat_neurog_ordered
                , cluster_rows = F, cluster_cols = F
                , color = rev(pals::brewer.rdbu(100))
                , scale = "row", show_rownames = T, fontsize_row = 15
                , angle_col = 0, fontsize_col = 20)

# annotate only selected genes
add.flag(ph1, goi, repel.degree = .1)

# save heatmap
png(paste0(OutDir, "Pseudotime_Heatmap_neurogenic_CLIP_Org.png"),width = 750, 
    height = 500)
add.flag(ph1, goi, repel.degree = .1)
dev.off()

# select control and TSC mutant cells of CLIP neurogenic lineage
ctrl_cells_neurog <- cds_ventral_subset_neurog@colData$barcode[
  cds_ventral_subset_neurog@colData$media_genotype %in% c("LN.Ctrl", "HN.Ctrl")]
TSC_cells_neurog <- cds_ventral_subset_neurog@colData$barcode[
  !(cds_ventral_subset_neurog@colData$media_genotype %in% 
      c("LN.Ctrl", "HN.Ctrl"))]

# plot gois along pseudotime for TSC mutant cells
pl1 <- plot_genes_in_pseudotime(cds_ventral_subset_neurog[
  c("GFAP", "HOPX", "SCGN", "PTGDS", "EDNRB", "EGFR", "DLX2", "DLX6.AS1")
                                                          ,TSC_cells_neurog]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("GFAP", "HOPX","EDNRB", 
                                                  "PTGDS", "EGFR", "DLX2", 
                                                  "DLX6.AS1", "SCGN")) + 
  scale_color_viridis_c(alpha = .4)

png(paste0(OutDir, "Pseudotime_genes_neurogenic_MUT_CLIP_Org.png"), width = 500, 
    height = 1000)
pl1
dev.off()

# plot gois along pseudotime for control cells
pl1 <- plot_genes_in_pseudotime(cds_ventral_subset_neurog[
  c("GFAP", "HOPX", "SCGN", "PTGDS", "EDNRB", "EGFR", "DLX2", "DLX6.AS1")
                                                          ,ctrl_cells_neurog]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("GFAP", "HOPX","EDNRB", 
                                                  "PTGDS", "EGFR", "DLX2", 
                                                  "DLX6.AS1", "SCGN")) + 
  scale_color_viridis_c(alpha = .4)

png(paste0(OutDir, "Pseudotime_genes_neurogenic_CTRL_CLIP_Org.png"), width = 500
    , height = 1000)
pl1
dev.off()

# save with alpha = 1 for scale bar
pl1 <- plot_genes_in_pseudotime(cds_ventral_subset_neurog[
  c("GFAP", "HOPX", "SCGN", "PTGDS", "EDNRB", "EGFR", "DLX2", "DLX6.AS1")
                                                          ,ctrl_cells_neurog]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("GFAP", "HOPX","EDNRB", 
                                                  "PTGDS", "EGFR", "DLX2", 
                                                  "DLX6.AS1", "SCGN")) + 
  scale_color_viridis_c(alpha = 1)
png(paste0(OutDir, "Pseudotime_genes_neurogenic_scale_Org_210405.png"), width = 500, height = 1000)
pl1
dev.off()

# save table
write.csv(dat_neurog_ordered, file = 
            paste0(OutDir, "Pseudotime_Heatmap_neurogenic_CLIP_Org_210405.csv"))

# save objects
saveRDS(cds.integration, file = 
          paste0(OutDir, "tsc_paper_integration_all_pseudotime.Rds"))
saveRDS(cds.integration_ventral, file = 
          paste0(OutDir, "tsc_paper_integration_ventral_pseudotime.Rds"))
saveRDS(cds_ventral_subset_OPC, file = 
          paste0(OutDir, "tsc_paper_integration_ventral_OPC_pseudotime.Rds"))
saveRDS(cds_ventral_subset_neurog, file = 
          paste0(OutDir, "tsc_paper_integration_ventral_neurog_pseudotime.Rds"))

# save graph test results
write.csv(graph_test_OPC, file = 
            paste0(OutDir, "tsc_paper_integration_OPC_pseudotime_graphtest.Rds"))
write.csv(graph_test_ventral_neurog, file = 
            paste0(OutDir, "tsc_paper_integration_CGE_pseudotime_graphtest.Rds"))


# plot umap with overview of OPC vs Neurogenic trajectory ----------------------
# select UMAP coordinates of OPC lineage with pseudotime values
OPC_UMAP_coord <- reducedDims(cds_ventral_subset_OPC)$UMAP %>% 
  data.frame() %>%
  rename(UMAP_1 = 1, UMAP_2 = 2) %>%
  mutate(barcode = row.names(.)) %>%
  left_join(pseudotime(cds_ventral_subset_OPC) %>% data.frame() %>% 
              rename(pseudotime = 1) %>% 
              mutate(barcode = row.names(.)),  by = "barcode") %>%
  select(barcode, UMAP_1, UMAP_2, pseudotime)

# select UMAP coordinates of CLIP neurogenic lienage with pseudotime values
neurog_UMAP_coord <-  reducedDims(cds_ventral_subset_neurog)$UMAP %>% 
  data.frame() %>%
  rename(UMAP_1 = 1, UMAP_2 = 2) %>%
  mutate(barcode = row.names(.)) %>%
  left_join(pseudotime(cds_ventral_subset_neurog) %>% data.frame() %>% 
              rename(pseudotime = 1) %>% 
              mutate(barcode = row.names(.)),  by = "barcode") %>%
  select(barcode, UMAP_1, UMAP_2, pseudotime)

# UMAP coordinates of all ventral cells
All_UMAP_coord <- reducedDims(cds.integration_ventral)$UMAP%>% 
  data.frame() %>%
  rename(UMAP_1 = 1, UMAP_2 = 2)

# plot UMAP color coded for OPC and neurogenic pseudotime
pl1 <- ggplot() +
  geom_point(data = All_UMAP_coord, aes(x = UMAP_1, y = UMAP_2), 
             color = "light grey", alpha = .2) +
  geom_point(data = neurog_UMAP_coord, aes(x = UMAP_1, y = UMAP_2, 
                                           color = pseudotime/2), alpha = .2) +
  scale_color_gradientn("Neurogenic", colors = 
                          pals::kovesi.diverging_linear_bjr_30_55_c53(100)) + 
  ggnewscale::new_scale_color() +
  geom_point(data = OPC_UMAP_coord, aes(x = UMAP_1, y = UMAP_2, 
                                        color = pseudotime*5), alpha = .2) +
  scale_color_gradient("OPC", low = "#002AD7", high = "dark green") + 
  ggnewscale::new_scale_color() +
  theme(panel.grid = element_blank(), axis.line = element_line(size = 1), 
        panel.background = element_blank()
        , axis.text = element_text(size = 10, face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"))

# save plot
png(paste0(OutDir, "Pseudotime_UMAP_Org.png"), width = 750, height = 500)
pl1
dev.off()

# Reproduce on fetal integration -----------------------------------------------
# set Output Directory 
OutDir <- 'path_to/OutDir'
dir.create(OutDir)
setwd(OutDir)

# load monocle object fetal integration
cds.integration_fetal <- readRDS("path_to/AllOrg_Diss_Fetal_monocle.Rds")

# subset for ventral lineage only
ventral_withfetal <- cds.integration_fetal@colData$barcode[
  cds.integration_fetal@colData$clusters_low %in% c(4,5,7,10,13,15,17,18)]
cds.integration_fetal_ventral <- cds.integration_fetal[,ventral_withfetal]

# recreate one single partition
vec <- rep(1, length(ventral_withfetal))
names(vec) <- ventral_withfetal
vec <- factor(vec)

cds.integration_fetal_ventral@clusters$UMAP$partitions <- vec

# subset for fetal cells only
cds_withfetal_ventral_fetal_only <- 
  cds.integration_fetal_ventral[,cds.integration_fetal_ventral@colData$barcode[
    cds.integration_fetal_ventral@colData$orig.ident %in% "ORC.reanalysis"]]

# subset partitions otherwise learn_graph fails
cds_withfetal_ventral_fetal_only@clusters$UMAP$partitions  <- 
  cds_withfetal_ventral_fetal_only@clusters$UMAP$partitions[
    cds_withfetal_ventral@colData$barcode[
      cds_withfetal_ventral@colData$orig.ident %in% "ORC.reanalysis"]]


# learn trajectory graph on fetal only
cds_withfetal_ventral_fetal_only <- 
  learn_graph(cds = cds_withfetal_ventral_fetal_only, use_partition = F, 
              verbose = T, close_loop = F, learn_graph_control = 
                list(minimal_branch_len = 5, orthogonal_proj_tip = F, 
                     rann.k = NULL))

# create new cell type annotation in fetal only 
# focusing on OPC, oRG and oRG/Astrocyte
cds_withfetal_ventral_fetal_only@colData$Subtype_2 <- 
  cds_withfetal_ventral_fetal_only@colData$Subtype
cds_withfetal_ventral_fetal_only@colData$Subtype_2[
  !(cds_withfetal_ventral_fetal_only@colData$Subtype %in% 
      c("OPC", "oRG", "oRG/Astrocyte"))] = "other fetal"

plot_cells(cds_withfetal_ventral_fetal_only, color_cells_by = "Subtype_2"
           , label_cell_groups = F, show_trajectory_graph = F) + 
  scale_color_manual(values = as.vector(pals::alphabet(25)))

# create new cell type annotation in organoid integration 
# focusing on OPC, oRG and oRG/Astrocyte
cds.integration_fetal_ventral@colData$Subtype_2 <- cds.integration_fetal_ventral@colData$Subtype
cds.integration_fetal_ventral@colData$Subtype_2[
  !(cds.integration_fetal_ventral@colData$Subtype %in% 
      c("OPC", "oRG", "oRG/Astrocyte", NA))] = "other fetal"
cds.integration_fetal_ventral@colData$Subtype_2[
  cds.integration_fetal_ventral@colData$Subtype_2 %in% c(NA)] = "Organoid"

plot_cells(cds.integration_fetal_ventral, color_cells_by = "Subtype_2"
           , label_cell_groups = F, show_trajectory_graph = F) 


# Subset lineages and create pseudotime plots for fetal ------------------------

# subset for OPC lineage
cds_fetal_OPC <- 
  monocle3::choose_graph_segments(cds_withfetal_ventral_fetal_only)
cds_fetal_subset_OPC <- cds_withfetal_ventral_fetal_only[,colnames(cds_fetal_OPC)]

# subset for neurogenic lineage
cds_fetal_neurog <- 
  monocle3::choose_graph_segments(cds_withfetal_ventral_fetal_only)
cds_fetal_subset_neurog <- cds_withfetal_ventral_fetal_only[,colnames(cds_fetal_neurog)]

# perform graph testing on OPC clustering
graph_test_fetal_OPC <- graph_test(cds_fetal_subset_OPC)

# subset graph testing of OPC
graph_test_fetal_OPC_subset <- 
  graph_test_fetal_OPC %>% 
  arrange(desc(morans_I)) %>%
  filter(morans_I > .1) %>%
  select(gene_short_name)

# perform graph testing on neurogenic cells
graph_test_fetal_neurog <- graph_test(cds_fetal_subset_neurog)

# subset graph testing of neurogenic cells
graph_test_fetal_neurog_subset <- 
  graph_test_fetal_neurog %>% 
  arrange(desc(morans_I)) %>%
  filter(morans_I > .1) %>%
  select(gene_short_name)

# bin pseudotime of OPC
binned_pseudotime_OPC_fetal <- 
  data.frame(floor(pseudotime(cds_fetal_subset_OPC)*3)) %>% 
  rename(bin=1) %>% mutate(barcode = row.names(.)) %>%
  select(barcode,bin)

# bin pseudotime of neurogenic
binned_pseudotime_ventral_fetal <- 
  data.frame(floor(pseudotime(cds_fetal_subset_neurog)/2)) %>% 
  rename(bin=1) %>% mutate(barcode = row.names(.)) %>%
  select(barcode,bin)

# plot OPC pseudotime on fetal cells -------------------------------------------
# set window and step for sliding average
window <- 2
step <- 1

# aggregate expression based on graph testing and binned pseudotime
dat_OPC_fetal <- 
  aggregate_gene_expression(cds_fetal_subset_OPC[
    graph_test_fetal_OPC_subset$gene_short_name,], 
    cell_group_df = binned_pseudotime_OPC_fetal, 
    scale_agg_values = T, max_agg_value = 1)

dat_OPC_fetal_backup <- dat_OPC_fetal

dat_OPC_fetal <- as.matrix(dat_OPC_fetal_backup)

# order genes using a sliding average
dat_OPC_fetal_ordered <- 
  dat_OPC_fetal[order(apply(t(rollapply(t(dat_OPC_fetal), 
                                        width=window, by=step, FUN=mean)), 1, 
                            which.max)), ]

# select genes of interest for plotting
goi <- c("EDNRB", "PTGDS", "PDGFRA", "PCDH15", "OLIG1", "OLIG2", "EGFR",
         "HOPX", "TNC", "VIM", "CD9")

# plot heatmap
ph1 <- pheatmap(dat_OPC_fetal_ordered
                , cluster_rows = F, cluster_cols = F
                , color = rev(pals::brewer.rdbu(100))
                , scale = "row", show_rownames = T, fontsize_row = 15,
                angle_col = 0, fontsize_col = 20)

# annotate only selected genes
add.flag(ph1, goi, repel.degree = .1)

# save heatmap
png(paste0(OutDir, "Pseudotime_Heatmap_OPC_fetal.png"),width = 750, 
    height = 500)
add.flag(ph1, goi, repel.degree = .1)
dev.off()

# plot gois along OPC pseudotime for fetal cells
pl1 <- plot_genes_in_pseudotime(cds_fetal_subset_OPC[
  c("EDNRB", "PTGDS", "PDGFRA","OLIG2", "EGFR"),]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("EDNRB", "PTGDS", 
                                                  "EGFR", "PDGFRA", "OLIG2")) + 
  scale_color_viridis_c(alpha = .4)

png(paste0(OutDir, "Pseudotime_genes_OPC_fetal.png"), width = 500, height = 750)
pl1
dev.off()

# plot neurogenic pseudotime on fetal cells ------------------------------------
# set window and step for sliding average
window <- 3
step <- 1

# aggregate expression of graph test genes per bin
dat_neurog_fetal <- 
  aggregate_gene_expression(cds_fetal_subset_neurog[
    graph_test_fetal_neurog_subset$gene_short_name,], 
    cell_group_df = binned_pseudotime_ventral_fetal, 
    scale_agg_values = F, max_agg_value = 1)

dat_neurog_fetal_backup <- dat_neurog_fetal

dat_neurog_fetal <- as.matrix(dat_neurog_fetal_backup)

# order genes using a sliding average
dat_neurog_fetal_ordered <- 
  dat_neurog_fetal[order(apply(t(rollapply(t(dat_neurog_fetal), width=window,
                                           by=step, FUN=mean)), 1, which.max)),]

# select genes of interest for plotting
goi <-  c("EDNRB", "PTGDS", "SCGN", "EGFR", "DLX1", "DLX5", "SP8"
          , "DLX6.AS1", "ASCL1", "HOPX", "TNC", "VIM", "GFAP", "S100B")

# plot heatmap
ph1 <- pheatmap(dat_neurog_fetal_ordered
                , cluster_rows = F, cluster_cols = F
                , color = rev(pals::brewer.rdbu(100))
                , scale = "row", show_rownames = T, fontsize_row = 15, 
                angle_col = 0, fontsize_col = 20)

# annotate selected gois
add.flag(ph1, goi, repel.degree = .1)

# save plot
png(paste0(OutDir, "Pseudotime_Heatmap_neurogenic_fetal.png"),width = 750, 
    height = 500)
add.flag(ph1, goi, repel.degree = .1)
dev.off()

# plot gois along neurogenicpseudotime for fetal cells
pl1 <- 
  plot_genes_in_pseudotime(cds_fetal_subset_neurog[
    c("GFAP", "HOPX", "SCGN", "PTGDS", "EDNRB", "EGFR", "DLX2", "DLX6.AS1"),]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("GFAP", "HOPX","EDNRB", 
                                                  "PTGDS", "EGFR", "DLX2", 
                                                  "DLX6.AS1", "SCGN")) + 
  scale_color_viridis_c(alpha = .4)

png(paste0(OutDir, "Pseudotime_genes_neurogenic_fetal.png"), width = 500,
    height = 1000)
pl1
dev.off()

# plot with alpha = 1 for color scale
pl1 <- 
  plot_genes_in_pseudotime(cds_fetal_subset_neurog[
    c("GFAP", "HOPX", "SCGN", "PTGDS", "EDNRB", "EGFR", "DLX2", "DLX6.AS1"),]
                                , cell_size = 2, vertical_jitter = .1
                                , color_cells_by = "pseudotime"
                                , panel_order = c("GFAP", "HOPX","EDNRB", 
                                                  "PTGDS", "EGFR", "DLX2", 
                                                  "DLX6.AS1", "SCGN")) + 
  scale_color_viridis_c(alpha = 1)

png(paste0(OutDir, "Pseudotime_genes_neurogenic_fetal_scale.png"), width = 500,
    height = 1000)
pl1
dev.off()

# save tables
write.csv(dat_neurog_fetal_ordered
          , file = paste0(OutDir, "neurog_lineage_fetal.csv"))
write.csv(dat_OPC_fetal_ordered
          , file = paste0(OutDir, "OPC_lineage_fetal.csv"))
write.csv(graph_test_fetal_neurog
          , file = paste0(OutDir, "graph_test_fetal_neurog_lineage_fetal.csv"))
write.csv(graph_test_fetal_OPC
          , file = paste0(OutDir, "graph_test_fetal_OPC_lineage_fetal.csv"))

# plot UMAP color coded for OPC and neurogenic pseudotime
# extract UMAP coordinates and pseudotime values for fetal OPC lineage
OPC_UMAP_coord_fetal <- reducedDims(cds_fetal_subset_OPC)$UMAP %>% 
  data.frame() %>%
  rename(UMAP_1 = 1, UMAP_2 = 2) %>%
  mutate(barcode = row.names(.)) %>%
  left_join(pseudotime(cds_fetal_subset_OPC) %>% data.frame() %>% 
              rename(pseudotime = 1) %>% 
              mutate(barcode = row.names(.)),  by = "barcode") %>%
  select(barcode, UMAP_1, UMAP_2, pseudotime)

# extract UMAP coordinates and pseudotime values for fetal neurogenic lineage
neurog_UMAP_coord_fetal <-  reducedDims(cds_fetal_subset_neurog)$UMAP %>% 
  data.frame() %>%
  rename(UMAP_1 = 1, UMAP_2 = 2) %>%
  mutate(barcode = row.names(.)) %>%
  left_join(pseudotime(cds_fetal_subset_neurog) %>% data.frame() %>% 
              rename(pseudotime = 1) %>% 
              mutate(barcode = row.names(.)),  by = "barcode") %>%
  select(barcode, UMAP_1, UMAP_2, pseudotime)


# all UMAP coordinates with organoid cells
All_UMAP_coord_fetal <- reducedDims(cds.integration_fetal_ventral)$UMAP%>% 
  data.frame() %>%
  rename(UMAP_1 = 1, UMAP_2 = 2)

# All umap coordinates without organoid cells
All_UMAP_coord_fetalonly <- reducedDims(cds_withfetal_ventral_fetal_only)$UMAP%>% 
  data.frame() %>%
  rename(UMAP_1 = 1, UMAP_2 = 2)

# plot color coded UMAP
pl1 <- ggplot() +
  geom_point(data = All_UMAP_coord_fetal, aes(x = UMAP_1, y = UMAP_2), 
             color = "light grey", alpha = .2) +
  geom_point(data = All_UMAP_coord_fetalonly, aes(x = UMAP_1, y = UMAP_2), 
             color = "dark grey", alpha = 1) +
  geom_point(data = OPC_UMAP_coord_fetal, aes(x = UMAP_1, y = UMAP_2, 
                                              color = pseudotime*5), alpha = 1) +
  scale_color_gradient("OPC", low = "#002AD7", high = "dark green") + 
  ggnewscale::new_scale_color() +
  geom_point(data = neurog_UMAP_coord_fetal, 
             aes(x = UMAP_1, y = UMAP_2, color = pseudotime/2), alpha = 1) +
  scale_color_gradientn("Neurogenic", colors = 
                          pals::kovesi.diverging_linear_bjr_30_55_c53(100)) + 
  ggnewscale::new_scale_color() +
  theme(panel.grid = element_blank(), axis.line = element_line(size = 1), 
        panel.background = element_blank()
        , axis.text = element_text(size = 10, face = "bold"), axis.title = 
          element_text(size = 15, face = "bold"))

# save UMAP plot
png(paste0(OutDir, "Pseudotime_UMAP_fetal.png"), width = 750, height = 500)
pl1
dev.off()

# save UMAP of Subtypes
png(paste0(OutDir, "Subtype_UMAP_fetal_210405.png"), width = 750, height = 500)
plot_cells(cds.integration_fetal_ventral, color_cells_by = "Subtype_2",
           show_trajectory_graph = F, cell_size = 1.5
           , label_cell_groups = F) + 
  scale_color_manual("Subtypes", values = c("dark green", "dark blue", 
                                            "dark red", "light grey", "orange"))
dev.off()

# save objects
saveRDS(cds_withfetal_ventral_fetal_only, 
        paste0(OutDir, "cds_ventral_fetal_only.Rds"))
saveRDS(cds_fetal_subset_neurog, 
        paste0(OutDir, "cds_ventral_fetal_only_neurog.Rds"))
saveRDS(cds_fetal_subset_OPC, paste0(OutDir, "cds_ventral_fetal_only_OPC.Rds"))


