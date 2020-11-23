#' vReadTenX
#' vReadTenX reads host and pathogen matrix simultaneously
#' @title xReadTenX
#' @param sample Tab-delimited sample information file, such as sample, path
#' @param path Path that holds scRNA matrix
#' @param column Column of genes.tsv or features.tsv to use for gene names, default column=2
#' @param project Project name for the Seurat object
#' @param group Group name for the Seurat object
#' @param min.nFeature Minimum number of features in a cell
#' @param max.nFeature Maximum number of features in a cell
#' @param percent.mito Maximum percent of mito in a cell
#' @param min.nCount_RNA Minimum number of UMIs in a cell
#' @param max.nCount_RNA Maximum number of UMIs in a cell
#' @return Seurat object.

#' @author rstatistics
#' @export
vReadTenX <- function(sample=sample, path=path, column=column, project=project, group=group, min.nFeature=min.nFeature,
                      max.nFeature=max.nFeature, percent.mito=percent.mito, min.nCount_RNA=min.nCount_RNA, max.nCount_RNA=max.nCount_RNA){
  Dict <- list()
  Dict[["sample"]] = sample
  Dict[["path"]] = path
  Dict[["column"]] = column
  Dict[["project"]] = project
  Dict[["group"]] = group
  Dict["min.nFeature"] = min.nFeature
  Dict[["max.nFeature"]] = max.nFeature
  Dict[["percent.mito"]] = percent.mito
  Dict[["min.nCount_RNA"]] = min.nCount_RNA
  Dict[["max.nCount_RNA"]] = max.nCount_RNA
  object <- Read10X(paste(Dict[["path"]], "/", Dict[["sample"]], sep=""), gene.column = Dict[["column"]])
  object <- as.matrix(object)
  colnames(object) <- gsub("-[0-9]+", "", colnames(object))
  Microbe <- read.table(paste(Dict[["path"]], "/", Dict[["sample"]], "_matrix.txt", sep=""), header = TRUE, row.names = 1,
                        stringsAsFactors = FALSE, sep = "\t", quote = "", comment.char = "")
  object <- rbind(object, Microbe)
  object <- CreateSeuratObject(counts = object, min.cells = 5, min.features = 200, project = Dict[["project"]])
  object$sample <- Dict[["sample"]]
  object$group <- Dict[["group"]]
  object[["percent.mito"]]<-PercentageFeatureSet(object, pattern = "^MT-")
  VlnPlot(object=object, features=c("nFeature_RNA","nCount_RNA","percent.mito"), ncol=3, pt.size=0.01)
  FeatureScatter(object=object, feature1="nCount_RNA", feature2="nFeature_RNA")
  # QC and selecting cells for further analysis
  object <- subset(x=object, subset=nFeature_RNA > Dict["min.nFeature"] & nFeature_RNA < Dict["max.nFeature"] &
                     percent.mito < Dict[["percent.mito"]] & nCount_RNA > Dict[["min.nCount_RNA"]] & nCount_RNA < Dict[["max.nCount_RNA"]])
  return(object)
}
