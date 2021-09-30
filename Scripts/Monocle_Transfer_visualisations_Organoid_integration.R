## Transfer to monocle and visualisation of organoid integration
##------ Tue Sep 28 15:14:28 2021 ------##
## Oliver Eichmueller

# load required libraries ------------------------------------------------------
library(Seurat)
library(monocle3)
library(reticulate)
library(pals)
library(ggplot2)
library(ggnewscale)
library(dplyr)

# set Output Directory ---------------------------------------------------------
OutDir <- 'path_to/OutDir'

# read in seurat object --------------------------------------------------------
# Integration of organoid data (d110 H- and L-medium, d220 tumors)

sc.obj.integration <- readRDS("path_to/LN_HN_d220_integration.Rds")

# transfer RNA assay from seurat to monocle ====================================


# extract genes
genes <- as.data.frame(rownames(sc.obj.integration@assays$RNA), 
                       row.names = rownames(sc.obj.integration@assays$RNA))
colnames(genes) <- "gene_short_name"

# extract cells
cells <- as.data.frame(
  sc.obj.integration@assays[["RNA"]]@counts@Dimnames[[2]], 
  row.names = sc.obj.integration@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cells) <- "barcode"

# extract expression matrix
expression_matrix <- sc.obj.integration@assays[["RNA"]]@counts
expression_matrix <- expression_matrix[rownames(sc.obj.integration@assays$RNA), ]

# Assemble cell data set object
cds.integration <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cells,
                                     gene_metadata = genes)

# transfer dataset information
orig.idents <- sc.obj.integration@meta.data[["orig.ident"]]
names(orig.idents) <- sc.obj.integration@assays[["RNA"]]@data@Dimnames[[2]]

cds.integration@colData$orig.ident <- 
  orig.idents[cds.integration@colData$barcode]
cds.integration@colData$orig.ident <- 
  factor(cds.integration@colData$orig.ident)

# preprocess cds with variable features from seurat
cds.integration <- 
  preprocess_cds(cds.integration, method = "PCA", 
                 use_genes = VariableFeatures(sc.obj.integration))

# align cds based on dataset
cds.integration <- align_cds(cds.integration, alignment_group = "orig.ident")

# perform umap dimensionality reduction
cds.integration <- reduce_dimension(cds.integration, reduction_method = "UMAP")

# cluster cells
cds.integration <- cluster_cells(cds.integration, reduction = "UMAP", 25)

# store clustering information in cds
cds.integration@colData$clusters_low <- 
  cds.integration@clusters$UMAP$clusters[cds.integration@colData$barcode]

# learn trajectory graph
cds.integration <- learn_graph(cds.integration, use_partition = F, close_loop = F, learn_graph_control = 
                                 list(minimal_branch_len = 20))

# plot clustering
plot_cells(cds.integration, 
           color_cells_by = 'clusters_low',
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F
           , show_trajectory_graph = T
           , label_cell_groups = T)

# example for plotting genes
plot_cells(cds.integration, 
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F
           , show_trajectory_graph = F
           , label_cell_groups =F
           , genes = "EGFR"
           , min_expr = 1
           , cell_size = 1) + scale_color_gradient(low = "grey", high = "red")

# generate annotation for media and genotype
annot_table <- cds.integration@colData[,c(1, 3)]
annot_table$media_genotype[annot_table$orig.ident %in% 
                             c("81950.TSC2", "82546.TSC2")] = "LN.TSC2"
annot_table$media_genotype[annot_table$orig.ident %in% 
                             c("81949.WT", "82546.WT")] = "LN.Ctrl"
annot_table$media_genotype[annot_table$orig.ident %in% 
                             "TSC_June_Mutant"] = "HN.TSC2"
annot_table$media_genotype[annot_table$orig.ident %in% 
                             "Ctrl_June_Mutant"] = "HN.Ctrl"
annot_table$media_genotype[annot_table$orig.ident %in% 
                             "TSC_Dissection"] = "Old.Tumor"
cds.integration@colData$media_genotype <- annot_table$media_genotype

cds.integration@colData$media_genotype <- 
  factor(cds.integration@colData$media_genotype
         , levels = c("LN.Ctrl", "LN.TSC2", "HN.Ctrl", "HN.TSC2", "Old.Tumor"))

# plot media and genotype
plot_cells(cds.integration, 
           color_cells_by = 'media_genotype',
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F
           , show_trajectory_graph = F
           , label_cell_groups = F
           , alpha = .5) + scale_color_manual(values = brewer.dark2(5))

saveRDS(cds.integration, file = paste0(OutDir, 'AllOrgs_plusDiss_monocle.Rds'))

# Plot statistics on cell type abundances --------------------------------------

# sums of groups
sums_id <- data.frame(cds.integration@colData) %>%
  count(media_genotype)

# downsampled relative of clusters per genotype
combined_cl_ident_monocle <- data.frame(cds.integration@colData) %>%
  group_by(clusters_low, .drop = F) %>%
  count(media_genotype) %>%
  mutate(downsample = (n / sums_id$n) * min(sums_id$n)) %>%
  mutate(relative = round(downsample/ sum(downsample) *100)) %>%
  group_by(clusters_low, .drop = F) %>%
  mutate(annot = c(sum(relative) - relative[1] /2
                   , sum(relative) - relative[1] - relative[2]/2
                   , sum(relative) - sum(relative[1:2]) - relative[3]/2
                   , sum(relative) - sum(relative[1:3]) - relative[4]/2
                   , sum(relative) - sum(relative[1:4]) - relative[5]/2)/100) %>%
  ungroup()


combined_cl_ident_monocle_HN <- combined_cl_ident_monocle %>%
  filter(media_genotype %in% c("HN.TSC2", "HN.Ctrl")) %>%
  mutate(media = "HN") %>%
  group_by(media_genotype, .drop = FALSE) %>%
  mutate(relative = round(downsample/ sum(downsample) *100)) %>%
  group_by(clusters_low, .drop = F) %>%
  mutate(annot = relative/2) %>%
  ungroup()

combined_cl_ident_monocle_LN <- combined_cl_ident_monocle %>%
  filter(media_genotype %in% c("LN.TSC2", "LN.Ctrl")) %>%
  mutate(media = "LN") %>%
  group_by(media_genotype, .drop = FALSE) %>%
  mutate(relative = round(downsample/ sum(downsample) *100)) %>%
  group_by(clusters_low, .drop = F) %>%
  mutate(annot = relative/2) %>%
  ungroup()

combined_cl_ident_monocle_Diss <- combined_cl_ident_monocle %>%
  filter(media_genotype %in% c("Old.Tumor")) %>%
  mutate(media = "Old.Tumor") %>%
  group_by(media_genotype, .drop = FALSE) %>%
  mutate(relative = round(downsample/ sum(downsample) *100)) %>%
  group_by(clusters_low, .drop = F) %>%
  mutate(annot = relative/2) %>%
  ungroup()

# combine different conditions
combined_plotting <- combined_cl_ident_monocle_HN %>%
  rbind(combined_cl_ident_monocle_LN, combined_cl_ident_monocle_Diss)

# create media annotation
media_label <- c("Old Tumor", "High Nutrient d110", "Low Nutrient d110")
names(media_label) <- c("Old.Tumor", "HN", "LN")
combined_plotting$clusters_low <- as.numeric(combined_plotting$clusters_low)

# plot main clusters, highlighting increased cell types
ggplot(combined_plotting[combined_plotting$clusters_low %in% 1:16,]
       , aes(x = clusters_low, y = relative, fill = media_genotype)) +
  annotate(geom = "rect", xmin = 4.5, xmax = 9.5, ymin = 0, ymax = Inf
           , fill = "red", alpha = .1, position = position_nudge(x = -5)) + 
  annotate(geom = "rect", xmin = 10.5, xmax = 11.5, ymin = 0, ymax = Inf
           , fill = "red", alpha = .1, position = position_nudge(x = -5)) +
  annotate(geom = "rect", xmin = 12.5, xmax = 13.5, ymin = 0, ymax = Inf
           , fill = "red", alpha = .1, position = position_nudge(x = -5)) +
  new_scale_fill()  +
  geom_col( aes(x = clusters_low, y = relative, fill = media_genotype), 
            position = "dodge", width = .8) +
  geom_text(data = combined_plotting[combined_plotting$media_genotype %in% 
                                       c("HN.TSC2", "LN.TSC2")&
                                       combined_plotting$clusters_low %in% 1:16,]
            , aes(x = clusters_low, y = relative +4, label = round(relative)), 
            color = "black", nudge_x = .2) +
  geom_text(data = combined_plotting[combined_plotting$media_genotype %in% 
                                       c("HN.Ctrl", "LN.Ctrl")&
                                       combined_plotting$clusters_low %in% 1:16,]
            , aes(x = clusters_low, y = relative +4, label = round(relative)), 
            color = "black", nudge_x = -.2) +
  geom_text(data = combined_plotting[combined_plotting$media_genotype %in% 
                                       c("Old.Tumor")&
                                       combined_plotting$clusters_low %in% 1:16,]
            , aes(x = clusters_low, y = relative +4, label = round(relative)), 
            color = "black", nudge_x = 0) +
  scale_fill_manual(values = c(rgb(102, 153, 102, maxColorValue = 255)
                               , rgb(51, 102, 153, maxColorValue = 255)
                               , rgb(102, 153, 102, maxColorValue = 255)
                               , rgb(51, 102, 153, maxColorValue = 255)
                               , rgb(51, 102, 153, maxColorValue = 255))) +
  theme_light() +
  theme(strip.text.y = element_text(face = "bold", size = 10, color = "black")
        , strip.background.y = element_rect(fill = 
                                              c(alpha("lightgrey", alpha = .5), 
                                                alpha("blue", alpha = .5), 
                                                alpha("green", alpha = .5)))) +
  xlab("Clusters") + ylab("Percentage of Dataset")+
  scale_x_continuous(breaks = 1:16, sec.axis = dup_axis(labels = waiver())) + 
  scale_y_continuous(limits = c(0,28))+
  facet_grid(rows = vars(media), labeller = labeller(media = media_label))

saveRDS(combined_plotting, file = pate0(OutDir, "statistics_AllOrgs.Rds"))


# calculate top markers plot dotplot -------------------------------------------

top_markers_diss_d110 <- top_markers(cds.integration, group_cells_by = "clusters_low")

top_markers_diss_d110_order <- top_markers_diss_d110 %>%
  arrange(cell_group)

write.simple.xlsx(top_markers_diss_d110_order)  

# select individual genes of interest for plotting
goi <- c("BCAN", "CLU", "CST3", "DLX6.AS1", "GAD2", "DLX2", "HMGB2", "TOP2A"
         , "PCNA", "MKI67", "SCGN", "DLX5", "EOMES", "NEUROD2", "NEUROD6", "EMX2"
         , "SLC1A3", "TTYH1", "LHX2", "SATB2", "FEZF2", "PTGDS", "EDNRB")

plot_genes_by_group(cds.integration, group_cells_by = "clusters_low"
                    , markers = goi, scale_max = 3
                    , max.size = 8, axis_order = "marker_group")

# Correlation plotting of fetal vs organoid cell types -------------------------
# Uses cds.integration_fetal (transferring Organoid + fetal to monocle)

# calculate gene modules in organoids
gene_modules_organoids <- find_gene_modules(cds.integration, resolution = 1e-1)

# generate cluster annotation df
clusters_low_df <- cds_all@cds.integration %>% data.frame() %>% 
  select(barcode, clusters_low)

# aggregate expression of modules per cluster annotation
aggregated_modules_organoids <- 
  aggregate_gene_expression(cds.integration, cell_group_df = clusters_low_df, 
                            gene_group_df = gene_modules_organoids)

# generate subtype annotation df for fetal
subtype_fetal_df <- cds.integration_fetal@colData %>% data.frame() %>% 
  filter(orig.ident %in% "ORC.reanalysis") %>% select(barcode, Subtype)

# subset subtype annotation df of fetal
subtype_fetal_df_subset <- subtype_fetal_df %>% 
  filter(Subtype %in% c("Deep Layer", "Newborn", "Layer IV", "Upper Layer", 
                        "IPC/newborn", "Cajal Retzius", "MGE2", "SST-MGE1", 
                        "IPC-div1", "IPC-div2", "vRG", "OPC", "late", "tRG", 
                        "IPC-new", "Mural", "oRG", "oRG/Astrocyte"))

# aggregate per module and subsetted subtype annotation
aggregated_modules_fetal_subset <- 
  aggregate_gene_expression(cds.integration_fetal, 
                            cell_group_df = subtype_fetal_df_subset, 
                            gene_group_df = gene_modules_organoids)

# bind fetal and organoid aggregated expression
aggregated_modules_bound_subset <- cbind(aggregated_modules_fetal_subset, 
                                         aggregated_modules_organoids)

# correlation aggregated expression of organoids and fetal
cor_aggregated_modules_subset <- 
  cor(x = as.matrix(aggregated_modules_bound_subset))

pl1 <- pheatmap(cor_aggregated_modules_subset[1:18, 19:38]
                , clustering_distance_rows = "manhattan"
                , clustering_distance_cols = "manhattan"
                , clustering_method = "mcquitty"
                , cutree_rows = 6, cutree_cols = 5)

png(paste0(OutDir, "Correlation_Subtypes_fetalOrg.png"), 
    width =750, height = 500)
pl1
dev.off()

pdf(paste0(OutDir, "Correlation_Subtypes_fetalOrg.pdf")
    , width = 10, height = 7.5)
pl1
dev.off()

write.csv(cor_aggregated_modules_subset, 
          file = paste0(OutDir, "Correlation_Subtypes_fetalOrg.csv"))
write.csv(gene_modules_organoids %>% arrange(module), 
          file = paste0(OutDir, "modules_Correlation_Subtypes_fetalOrg.csv"))
