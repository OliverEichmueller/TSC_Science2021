## Transfer to monocle and visualisation of organoid + fetal integration
##------ Tue Sep 28 16:01:01 2021 ------##
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

# set Output Directory ---------------------------------------------------------
OutDir <- 'path_to/OutDir'

# read in seurat object --------------------------------------------------------
# Integration of organoid data (d110 H- and L-medium, d220 tumors)

sc.obj.integration <- readRDS("path_to/LN_HN_d220_fetal_integration.Rds")

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

ages <- sc.obj.integration@meta.data[["Age"]]
names(ages) <- sc.obj.integration@assays[["RNA"]]@data@Dimnames[[2]]

cds.integration@colData$Age <- 
  ages[cds.integration@colData$barcode]
cds.integration@colData$Age <- 
  factor(cds.integration@colData$Age)


annot_age <- data.frame(cds.integration@colData$Age)
colnames(annot_age) = "Age"
annot_age$Age <- factor(annot_age$Age, levels = c(6, 10, 14, 18, 22, "Organoid"))
annot_age$Age[is.na(annot_age$Age) == T] = "Organoid"
cds.integration@colData$Age = annot_age[cds.integration@colData$barcode,1]



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

# cluster cells lower res
cds.integration <- cluster_cells(cds.integration, reduction = "UMAP", 40)

# store clustering information in cds
cds.integration@colData$clusters_superlow <- 
  cds.integration@clusters$UMAP$clusters[cds.integration@colData$barcode]

# learn trajectory graph
cds.integration <- learn_graph(cds.integration, use_partition = F, 
                               close_loop = F, learn_graph_control = 
                                list(minimal_branch_len = 20))

# calculate top markers and save
top_markers_fetal_intall <- 
  top_markers(cds.integration, group_cells_by = "clusters_low")

top_markers_fetal_intall_arranged <- top_markers_fetal_intall %>%
  arrange(cell_group, marker_score)

write.simple.xlsx(top_markers_fetal_intall_arranged)

# plot UMAPs

plot_cells(cds.integration, 
           color_cells_by = 'Age',
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F
           , show_trajectory_graph = F
           , label_cell_groups = F, cell_size = 1
           , alpha = .5) + 
  scale_color_manual(values = c(kovesi.rainbow(5), "lightgrey"))

plot_cells(cds.integration, 
           color_cells_by = 'clusters_superlow',
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F
           , show_trajectory_graph = F
           , label_cell_groups = F
           , cell_size = 1
           , alpha = .5)

# subset only progenitor cells -------------------------------------------------
# select progenitors
cds.integration_subset <- choose_cells(cds.integration, 
                                       reduction_method = "UMAP")

# process subsetted dataset
cds.integration_subset <- reduce_dimension(cds.integration_subset, 
                                           reduction_method = "UMAP")
cds.integration_subset <- cluster_cells(cds.integration_subset, 
                                        reduction = "UMAP", 30)
cds.integration_subset@colData$clusters_low <- 
  cds.integration_subset@clusters$UMAP$clusters[cds.integration_subset@colData$barcode]

cds.integration_subset <- learn_graph(cds.integration_subset, use_partition = F, 
                                      close_loop = F, learn_graph_control = 
                                        list(minimal_branch_len = 20))

# re-annotate age
annot_age <- data.frame(cds.integration_subset@colData$Age)
colnames(annot_age) = "Age"
annot_age$Age <- factor(annot_age$Age, levels = c(6, 10, 14, 18, 22, "Organoid"))
annot_age$Age[is.na(annot_age$Age) == T] = "Organoid"
cds.integration_subset@colData$Age = annot_age[cds.integration_subset@colData$barcode,1]

# plot UMAPs
plot_cells(cds.integration_subset, 
           color_cells_by = 'Age',
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F
           , show_trajectory_graph = F
           , label_cell_groups = F
           , cell_size = 1
           , alpha = .5) + s
cale_color_manual(values = c(kovesi.rainbow(5), "lightgrey"))

plot_cells(cds.integration_subset, 
           color_cells_by = 'cluster',
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F
           , show_trajectory_graph = F
           , label_cell_groups = F
           , cell_size = 1)

# select quiescent cells and calculate modules of quiescent cells

quiescent_cells <- row.names(colData(cds.integration_subset)[colData(
  cds.integration_subset)$clusters_low %in% c(1, 3, 5, 8, 9),])

quiescent_modules <- find_gene_modules(cds.integration_subset[, quiescent_cells],
                                       reduction_method = "UMAP", 
                                       resolution = 1e-2)

age_df <- tibble::tibble(cell = row.names(colData(cds.integration_subset))
                         , age = cds.integration_subset@colData$Age)

age_df_fetal <- age_df %>%
  filter(is.na(age) == F)

org_cells_df <- 
  tibble::tibble(cell = row.names(colData(cds.integration_subset)[
    colData(cds.integration_subset)$orig.ident != "ORC.reanalysis",]), 
    cluster = colData(cds.integration_subset)[colData(
      cds.integration_subset)$orig.ident != "ORC.reanalysis",]$clusters_low)

org_cells_df2 <- org_cells_df %>%
  mutate(cluster = paste0(cluster, "_organoid"))

age_df_fetal2 <- age_df_fetal %>%
  mutate(age = paste0("GW_", age)) %>%
  dplyr::rename(cluster = age)

combined_cell_df <- org_cells_df2 %>%
  bind_rows(age_df_fetal2)

combined_cell_df_quiescent <- combined_cell_df %>%
  filter(cell %in% quiescent_cells)

# calculate aggregated expression per module and age/cell type
agg_mat_combined_quiescent <- 
  aggregate_gene_expression(cds.integration_subset, 
                            quiescent_modules, combined_cell_df_quiescent)

cor_agg_mat_combined_quiescent <- 
  cor(x = data.frame(agg_mat_combined_quiescent), method = "pearson")

pheatmap(cor_agg_mat_combined_quiescent[1:5, 6:10], 
         cluster_rows = T, cluster_cols = T, 
         display_numbers = T, fontsize_number = 15, number_color = "black", 
         color = viridis(100))

pheatmap(cor_agg_mat_combined_quiescent, cluster_rows = T, cluster_cols = T
         , display_numbers = T, fontsize_number = 15, number_color = "black"
         , color = viridis(100))


plot_cells(cds.integration_subset,
           genes = quiescent_modules
           , scale_to_range = T
           , show_trajectory_graph = F
           , label_groups_by_cluster = F
           , label_cell_groups = F)+
  scale_color_gradient(low = "grey", high = "red")

# to save:
# - cor_agg_mat_combined_quiescent
# - agg_mat_combined_quiescent
# - quiecsent_modules
# - cds.integration
# - cds.integration_subset
OutDir <- 'path_to/OutDir'

write.simple.xlsx(cor_agg_mat_combined_quiescent)

cor_agg_mat_combined_quiescent_cutout <- 
  cor_agg_mat_combined_quiescent[1:5, 6:10]
write.simple.xlsx(cor_agg_mat_combined_quiescent_cutout)

write.simple.xlsx(quiescent_modules)

saveRDS(cds.integration, file = paste0(OutDir, 'AllOrg_Diss_Fetal_monocle.Rds'))
saveRDS(cds.integration_subset, 
        file = paste0(OutDir, 'AllOrg_Diss_Fetal_progs_monocle.Rds'))
saveRDS(combined_cl_ident_monocle_age, 
        file = paste0(OutDir, 'AllOrg_Diss_Fetal_progs_stat.Rds'))

# pct fetal ages of progenitors ------------------------------------------------

cds_from_seurat_subset@colData

sums_id <- data.frame(cds.integration_subset@colData) %>%
  count(Age) %>%
  filter(Age != "Organoid")

combined_cl_ident_monocle_age <- data.frame(cds.integration_subset@colData) %>%
  group_by(clusters_low, .drop = F) %>%
  count(Age) %>%
  filter(Age != "Organoid") %>%
  mutate(downsample = (n / sums_id$n) * min(sums_id$n)) %>%
  mutate(relative = round(downsample/ sum(downsample) *100)) %>%
  group_by(clusters_low, .drop = F) %>%
  mutate(annot = c(sum(relative) - relative[1] /2
                   , sum(relative) - relative[1] - relative[2]/2
                   , sum(relative) - sum(relative[1:2]) - relative[3]/2
                   , sum(relative) - sum(relative[1:3]) - relative[4]/2
                   , sum(relative) - sum(relative[1:4]) - relative[5]/2)/100) %>%
  ungroup()

ggplot(combined_cl_ident_monocle_age
       , aes(x = clusters_low, y = relative, fill = Age)) +
  geom_col(position = "fill") + scale_fill_manual(values = c(kovesi.rainbow(5)))+
  geom_text(aes(x = clusters_low, y = annot, label = relative))


