#' Average HeatMap for specific feature
#' Average HeatMap for specific feature in all clusters
#' @title vHeatMap visualization
#' @param object Seurat object
#' @param key The DE matrix are stored in object@misc$key
#' @param cluster Cluster specified by the user
#' @param feature.by Feature which is used to divid cells into two groups
#' @param markers A vector of features to plot, default is top_n features
#' @param colors User defined colours
#' @param top_n Top n markers to plot
#' @param font.size Font size
#' @param assay Assay to pull from
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @return
#' @example
#' objec@misc$vDE <- vDE(object = object, features.by = features.by, min.cells = 20, clusters = NULL)
#'
#' vHeatMap(object=object, key="vDE", cluster=cluster, feature.by = feature.by, markers = NULL)
#'
#'
#' @author rstatistics
#' @export
vHeatMap <- function(object = object, key = "vDE", cluster = cluster, feature.by = feature.by, markers = NULL,
                     colors = NULL, top_n = 5, font.size = 8, assay = NULL, slot = NULL){
  if (is.null(colors)){
    colors <- c("white", "red")
  }
  if (is.null(assay)){
    assay = "RNA"
  }
  if (is.null(slot)){
    slot = "data"
  }
  DefaultAssay(object = object) <- assay
  cluster = as.character(cluster)
  object <- subset(object, seurat_clusters %in% cluster)
  DE <- object@misc[[key]]
  Expr <- as.matrix(GetAssayData(object=object, assay="RNA", slot="data"))
  object$group <- ifelse(Expr[feature.by, ] > 0, "Pos", "Neg")
  if (is.null(markers)){
    markers.up <- rownames(subset(object@misc[[key]][[cluster]][[feature.by]], avg_logFC > 0 & p_val_adj < 1e-2))[1:top_n]
    markers.dn <- rownames(subset(object@misc[[key]][[cluster]][[feature.by]], avg_logFC < 0 & p_val_adj < 1e-2))[1:top_n]
    markers <- unique(c(markers.up, markers.dn))
    if (length(markers) < 1){
      warning("Some markers does not exist!")
    }else if (length(markers)==0){
      stop("No markers left for heatmap visualization!")
    }
  }else{
    markers_len <- length(unique(markers))
    markers = unique(intersect(markers, rownames(object)))
    if (length(markers) < markers_len){
      warning("Some markers does not exist!")
    }else if (length(markers)==0){
      stop("No markers left for heatmap visualization!")
    }
  }
  #markers <- intersect(markers, VariableFeatures(object = object))
  if(length(markers)==0){
    stop("All markers do not exist in the variable features!\n")
  }
  return(suppressMessages(DoHeatmap(object, features = markers, group.by = "group", assay = "RNA", slot = "data") +
                            scale_fill_gradientn(colors = colors)))
}
