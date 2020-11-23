#' vIntegrate Matrix
#' vIntegrate scRNA-seq from host and pathogen
#' @title vIntegrate
#' @param x Tab-delimited sample information file, such as sample, path
#' @param path Path that holds scRNA matrix
#' @param column Column of genes.tsv or features.tsv to use for gene names, default column=2
#' @param project Project name for the Seurat object
#' @param group Group name for the Seurat object
#' @param min.nFeature Minimum number of features in a cell
#' @param max.nFeature Maximum number of features in a cell
#' @param percent.mito Maximum percent of mito in a cell
#' @param min.nCount_RNA Minimum number of UMIs in a cell
#' @param max.nCount_RNA Maximum number of UMIs in a cell
#' @param nVariable Number of features to select as top variable features
#' @param nPC Total Number of PCs to compute
#' @param res Resolution parameter to set. Default res=0.5
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return Seurat object.
#'
#' @author rstatistics
#' @export
vIntegrate <- function(x="/home/scrna/Lung/Pathogen-Track/scNSCLC/SampleInfo.txt", min.nFeature=500, max.nFeature=4000,
                       percent.mito=10, min.nCount_RNA=1000, max.nCount_RNA=Inf, nVariable=2000, nPC=30, res=0.5){
  requireNamespace("Seurat")
  requireNamespace("Matrix")
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  requireNamespace("cowplot")
  requireNamespace("harmony")
  # read sample information
  SampleInfo <- read.table(x, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "", comment.char = "")
  SampleInfoDict <- as.list(SampleInfo$group)
  SamplePathDict <- as.list(SampleInfo$path)
  GeneColumnDict <- as.list(SampleInfo$gene_column)
  # store sample infor into SampleInfoDict
  names(SampleInfoDict) <- SampleInfo$sample
  names(SamplePathDict) <- SampleInfo$sample
  names(GeneColumnDict) <- SampleInfo$sample
  objectList <- list()
  pb <- txtProgressBar(min = 1, max = length(SampleInfo$sample), initial = 1, style=3, file = stderr())

  for (sample in SampleInfo$sample){
    if (!file.exists(paste0(SamplePathDict[[sample]], "/", sample))){
      stop(paste0("Host single cell matrix ", sample, " was not found!"))
    }
    if (!file.exists(paste0(SamplePathDict[[sample]], "/", sample, "_matrix.txt"))){
      stop(paste0("Pathogen single cell matrix ", sample, "_matrix.txt", " was not found!"))
    }
    warning(paste0("Load ", sample, " matrix ... "))
    setTxtProgressBar(pb, sample)
    suppressMessages(objectList[[sample]] <- vReadTenX(sample=sample, path=SamplePathDict[[sample]], project=sample, column=GeneColumnDict[[sample]],
                                                       group=SampleInfoDict[[sample]], min.nFeature=min.nFeature, max.nFeature=max.nFeature,
                                                       percent.mito=percent.mito, min.nCount_RNA=min.nCount_RNA, max.nCount_RNA=max.nCount_RNA))
    warning("done.\n")
  }

  SeuratNames <- SampleInfo$sample

  # Merge each seurat object
  object <- merge(objectList[[SeuratNames[1]]],  y=objectList[SeuratNames[2:length(SeuratNames)]], add.cell.ids = as.character(SeuratNames), project = "object")

  object <- NormalizeData(object, verbose = TRUE)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nVariable)
  object <- ScaleData(object, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = TRUE)
  object <- RunPCA(object, pc.genes = VariableFeatures(object), npcs = nPC, verbose = TRUE)

  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = object, reduction = "pca", pt.size = .1, group.by = "sample")
  p2 <- VlnPlot(object = object, features = "PC_1", pt.size = .1, group.by = "sample")
  cowplot::plot_grid(p1,p2)

  object <- object %>% harmony::RunHarmony("sample", plot_convergence = TRUE)

  harmony_embeddings <- Embeddings(object, 'harmony')
  harmony_embeddings[1:5, 1:5]
  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = object, reduction = "harmony", pt.size = .1, group.by = "sample")
  p2 <- VlnPlot(object = object, features = "harmony_1", pt.size = .1, group.by = "sample")
  cowplot::plot_grid(p1,p2)

  object <- RunUMAP(object, reduction = "harmony", dims = 1:nPC)
  object <- RunTSNE(object, reduction = "harmony", dims = 1:nPC)
  object <- FindNeighbors(object, reduction = "harmony", dims = 1:nPC)
  object <- FindClusters(object, resolution = res)
  object <- identity(object)

  cols <- NA
  if (length(levels(Idents(object))) <= 36){
    cols = c("#1660A7","#FF6A00","#219418","#CD0C18","#814BB2","#794339","#DC59B6","#CC79A7","#FF0000","#11B3C6",
             "#AFB400","#00FFFF", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00",
             "#CC79A7", "#00AFBB", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#4477AA",
             "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
    object@misc$cols = cols
  }
  DimPlot(object, reduction = "umap", pt.size = 0.1, cols = cols)
  return(object)
}
