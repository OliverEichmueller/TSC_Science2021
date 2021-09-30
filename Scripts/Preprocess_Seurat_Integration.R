## Pre-processing Reference Seurat objects and integration
##------ Tue Sep 28 12:58:58 2021 ------##
## Oliver Eichmueller

# load required libraries ------------------------------------------------------
library(Seurat)
library(dplyr)


# define Output Directory
OutDir <- "path_to_OutDir"

# read in 10X data and create Seurat Object ====================================
# this is was performed on all 10X experiments
sc.data <- Read10X(data.dir = 'path_to_CellRanger_output')

sc.obj <- CreateSeuratObject(sc.data, min.cells = 10, min.features = 200)


sc.obj[["percent.mt"]] <- PercentageFeatureSet(sc.obj, pattern = "^MT\\-")

saveRDS(sc.obj, file = paste0(OutDir, "preprocessed_obj_min.cells.10.min.features.200.Rds"))

# subset d220 tumor dataset based on MULTI-seq barcodes ------------------------

multi_d220           <- read.csv("path_to/MULTI_barcodes_d220_tum.csv")

multi_d220_filtered  <- multi_d220 %>% filter(MULTI_barcode != "Doublet")

sc.obj.H.d220.tumors <- sc.obj.H.d220.tumors[,multi_d220_filtered$cell_barcode]

# Integrate 10X data ===========================================================
# this was performed with the respective combinations of datasets:
# d110 ctrl + mutant of L- and H-medium + d220 tumors
# d110 ctrl + mutant of L- and H-medium + d220 tumors + fetal dataset
# datasets were transferred to monocle
# alignment was performed in monocle again from raw values
# here shown exemplarily for combination with fetal data

# read in pre-processed objects ------------------------------------------------
sc.obj.fetal <- 
  readRDS('path_to_fetal/preprocessed_obj_min.cells.10.min.features.200.Rds')
sc.obj.L.ctrl1 <- 
  readRDS('path_to_L.ctrl1/preprocessed_obj_min.cells.10.min.features.200.Rds')
sc.obj.L.ctrl2 <- 
  readRDS('path_to_L.ctrl2/preprocessed_obj_min.cells.10.min.features.200.Rds')
sc.obj.L.mut1 <- 
  readRDS('path_to_L.mut1/preprocessed_obj_min.cells.10.min.features.200.Rds')
sc.obj.L.mut2 <- 
  readRDS('path_to_L.mut1/preprocessed_obj_min.cells.10.min.features.200.Rds')
sc.obj.H.ctrl1 <- 
  readRDS('path_to_H.ctrl1/preprocessed_obj_min.cells.10.min.features.200.Rds')
sc.obj.H.mut1 <- 
  readRDS('path_to/preprocessed_obj_min.cells.10.min.features.200.Rds')
sc.obj.H.d220.tumors <- 
  readRDS('path_to/preprocessed_obj_min.cells.10.min.features.200.Rds')

# update gene names ------------------------------------------------------------
map_check <- HGNChelper::getCurrentHumanMap()

Check_fetal   <- HGNChelper::checkGeneSymbols(row.names(sc.obj.fetal), 
                                              unmapped.as.na = F, 
                                              map = map_check)
Check_L.ctrl1 <- HGNChelper::checkGeneSymbols(row.names(sc.obj.L.ctrl1), 
                                              unmapped.as.na = F,
                                              map = map_check)
Check_L.ctrl2 <- HGNChelper::checkGeneSymbols(row.names(sc.obj.L.ctrl2), 
                                              unmapped.as.na = F, 
                                              map = map_check)
Check_L.mut1  <- HGNChelper::checkGeneSymbols(row.names(sc.obj.L.mut1), 
                                              unmapped.as.na = F, 
                                              map = map_check)
Check_L.mut2  <- HGNChelper::checkGeneSymbols(row.names(sc.obj.L.mut2), 
                                              unmapped.as.na = F, 
                                              map = map_check)
Check_H.ctrl1 <- HGNChelper::checkGeneSymbols(row.names(sc.obj.H.ctrl1), 
                                              unmapped.as.na = F, 
                                              map = map_check)
Check_H.mut1  <- HGNChelper::checkGeneSymbols(row.names(sc.obj.H.mut1), 
                                              unmapped.as.na = F, 
                                              map = map_check)
Check_H.d220  <- HGNChelper::checkGeneSymbols(row.names(sc.obj.H.d220.tumors), 
                                              unmapped.as.na = F, 
                                              map = map_check)

Check_fetal$Suggested.Symbol   <- make.names(Check_fetal$Suggested.Symbol, 
                                             unique = T)   
Check_L.ctrl1$Suggested.Symbol <- make.names(Check_L.ctrl1$Suggested.Symbol, 
                                             unique = T)   
Check_L.ctrl2$Suggested.Symbol <- make.names(Check_L.ctrl2$Suggested.Symbol, 
                                             unique = T)   
Check_L.mut1$Suggested.Symbol  <- make.names(Check_L.mut1$Suggested.Symbol, 
                                             unique = T)   
Check_L.mut2$Suggested.Symbol  <- make.names(Check_L.mut2$Suggested.Symbol, 
                                             unique = T)   
Check_H.ctrl1$Suggested.Symbol <- make.names(Check_H.ctrl1$Suggested.Symbol, 
                                             unique = T)   
Check_H.mut1$Suggested.Symbol  <- make.names(Check_H.mut1$Suggested.Symbol, 
                                             unique = T)   
Check_H.d220$Suggested.Symbol  <- make.names(Check_H.d220$Suggested.Symbol, 
                                             unique = T)   

sc.obj.fetal         <- RenameGenesSeurat(obj = sc.obj.fetal, 
                                          newnames = Check_fetal)
sc.obj.L.ctrl1       <- RenameGenesSeurat(obj = sc.obj.L.ctrl1, 
                                          newnames = Check_L.ctrl1)
sc.obj.L.ctrl2       <- RenameGenesSeurat(obj = sc.obj.L.ctrl2, 
                                          newnames = Check_L.ctrl2)
sc.obj.L.mut1        <- RenameGenesSeurat(obj = sc.obj.L.mut1, 
                                          newnames = Check_L.mut1)
sc.obj.L.mut2        <- RenameGenesSeurat(obj = sc.obj.L.mut2, 
                                          newnames = Check_L.mut2)
sc.obj.H.ctrl1       <- RenameGenesSeurat(obj = sc.obj.H.ctrl1, 
                                          newnames = Check_H.ctrl1)
sc.obj.H.mut1        <- RenameGenesSeurat(obj = sc.obj.H.mut1, 
                                          newnames = Check_H.mut1)
sc.obj.H.d220.tumors <- RenameGenesSeurat(obj = sc.obj.H.d220.tumors, 
                                          newnames = Check_H.d220)



# preprocess for integration ---------------------------------------------------

object.list <- list(sc.obj.L.mut1, sc.obj.L.ctrl1, sc.obj.L.mut2, sc.obj.L.ctrl1, 
                    sc.obj.H.mut1, sc.obj.H.ctrl1,
                    sc.obj.H.d220.tumors, sc.obj.fetal)

sample.tree <- matrix(data = c(-1, -3, -5, 1, 3, 5, 6, -2, -4, -6, 2, 4, -7, -8), nrow = 7)


for (i in 1: length(object.list)) {
  object.list[[i]][["percent.mt"]] <- PercentageFeatureSet(object.list[[i]], pattern = "^MT\\.")
  
}


for (i in 1: length(object.list)) {
  object.list[[i]] <- subset(object.list[[i]], subset = nFeature_RNA >1000)
  object.list[[i]] <- subset(object.list[[i]], subset = percent.mt <15)
  
}

for (i in 1: length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], 
                                    normalization.method = "LogNormalize", 
                                    scale.factor = 10000)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], 
                                           selection.method = "vst", 
                                           nfeatures = 2000)
}

for (i in 1: length(object.list)) {
  row.names(object.list[[i]]@assays$RNA@meta.features) <- 
    row.names(object.list[[i]]@assays$RNA)
  
  object.list[[i]]@active.assay <- "RNA"
}

# perform integration ----------------------------------------------------------
object.list_int_anchors <- FindIntegrationAnchors(object.list = object.list, 
                                                  dims = 1:30)


sc.obj.integration <- IntegrateData(anchorset = object.list_int_anchors, 
                                    dims = 1:30, sample.tree = sample.tree)

saveRDS(sc.obj.integration, file = paste0(OutDir, "LN_HN_d220_fetal_integration.Rds"))

# scale data and calc. visualis. -----------------------------------------------
# Data will be transferred to monocle where visualisation will be performed


# sc.obj.integration <- ScaleData(sc.obj.integration, verbose = FALSE)
# sc.obj.integration <- RunPCA(sc.obj.integration, npcs = 30)
# 
# sc.obj.integration <- RunUMAP(sc.obj.integration, dims = 1:30)
# sc.obj.integration <- CellCycleScoring(sc.obj.integration, 
#                                        s.features = cc.genes$s.genes, 
#                                        g2m.features = cc.genes$g2m.genes)
# sc.obj.integration <- FindNeighbors(sc.obj.integration,dims = 1:10)
# sc.obj.integration <- FindClusters(sc.obj.integration, resolution = .6)




