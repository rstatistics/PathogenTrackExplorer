#' xReadTenX
#' xReadTenX reads host and pathogen matrix simultaneously
#' @title xReadTenX
#' @param sample Tab-delimited sample information file, such as sample, path
#' @param path Minimum number of features in a cell.
#' @param column Maximum number of features in a cell.
#' @param project Maximum percent of mito in a cell.
#' @param group Number of variable features.
#' @param min.nFeature Minimum number of features in a cell.
#' @param max.nFeature Maximum number of features in a cell.
#' @param percent.mito Maximum percent of mito in a cell.
#' @param nCount_RNA Minimum number of UMIs in a cell.
#' @return a seurat object.

#' @author rstatistics
#' @export
xReadTenX <- function(sample=sample, path=path, column=column, project=project, group=group,
                      min.nFeature=min.nFeature, max.nFeature=max.nFeature, percent.mito=percent.mito, nCount=nCount){
  Dict <- list()
  Dict[["sample"]] = sample
  Dict[["path"]] = path
  Dict[["column"]] = column
  Dict[["project"]] = project
  Dict[["group"]] = group
  Dict["min.nFeature"] = min.nFeature
  Dict[["max.nFeature"]] = max.nFeature
  Dict[["percent.mito"]] = percent.mito
  Dict[["nCount_RNA"]] = nCount_RNA
  object <- Read10X(paste(Dict[["path"]], "/", Dict[["sample"]], sep=""), gene.column = Dict[["column"]])
  object <- CreateSeuratObject(counts = object, min.cells = 5, min.features = 200, project = Dict[["project"]])
  object$sample <- sample
  object$group <- group
  object[["percent.mito"]]<-PercentageFeatureSet(object, pattern = "^MT-")
  VlnPlot(object=object, features=c("nFeature_RNA","nCount_RNA","percent.mito"), ncol=3, pt.size=0.01)
  FeatureScatter(object=object, feature1="nCount_RNA", feature2="nFeature_RNA")
  # QC and selecting cells for further analysis
  object <- subset(x=object, subset=nFeature_RNA > Dict["min.nFeature"] & nFeature_RNA < Dict["max.nFeature"] & percent.mito < Dict[["percent.mito"]] & nCount_RNA > Dict[["nCount_RNA"]])
  return(object)
}
