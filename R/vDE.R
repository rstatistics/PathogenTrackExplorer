#' Differential gene expression analysis for all identity classes
#' Find differentially expressed genes in features for each clusters in a dataset
#' @param object Seurat object
#' @param features.by A vector of features which is used to divid cells into two groups
#' @param min.cells Minimum number of cells to do DGE
#' @param clusters Clusters selected to do DGE
#'
#' @param return
#' @example
#' object@misc$vDE <- vDE(object = object, features = c("TP53", "FLT3"), min.cells = 20, clusters = clusters)
#'
#'
#' @author rstatistics
#' @export
vDE <- function(object = object, features.by = features.by, min.cells = 20, clusters = NULL){
  results <- list()
  if (is.null(clusters))
    clusters = levels(object)

  for (cluster in clusters){
    cat(paste0("### ", "Analysis Cluster-", cluster, " ...\n"))
    if (is.null(results[[cluster]])){ results[[cluster]] = list() }
    Object <- subset(object, seurat_clusters==cluster)

    Expr <- as.matrix(GetAssayData(Object))

    for (feature in features.by){
      cat(paste0("====== ", "Analysis feature ", feature, " ... "))
      if (!feature %in% rownames(Expr)){
        cat("all cells with 0 reads.\n")
        next
      }
      Object$group <- ifelse(Expr[feature, ] > 0, "Pos", "Neg")
      if (length(Object$group[Object$group=="Pos"]) < min.cells){
        cat("too few cells to process.\n")
        next
      }
      results[[cluster]][[feature]] <- suppressMessages(suppressWarnings(FindMarkers(object = Object, assay = 'RNA', ident.1 = "Pos", ident.2 = "Neg", group.by = "group", only.pos = FALSE, min.pct = 0.25, test.use = 'MAST', verbose = FALSE)))
      cat("done.\n")
    }
  }
  # clean up results
  for (cluster in clusters){
    if (length(results[[cluster]]) == 0){ results[[cluster]] <- NULL; next }
    for (feature in features.by){
      if (length(results[[cluster]][[feature]]) == 0){ next }
      results[[cluster]][[feature]] <- subset(results[[cluster]][[feature]], p_val_adj <= 1)
      if (nrow(results[[cluster]][[feature]]) <= 1){ results[[cluster]][[feature]] <- NULL; next }
    }
    if (length(results[[cluster]]) == 0){ results[[cluster]] <- NULL; next }
  }
  return(results)
}
