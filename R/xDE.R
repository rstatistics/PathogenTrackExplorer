#' Differential gene expression analysis for each cluster
#' Perform differential analysis (DE) between two groups for each clusters in a dataset
#' @param object Seurat object.
#' @param comparison Vector of nominal variables from group.by. Eg., comparison=c("Normal", "Tumor")
#' @param min.cells Minimum number of cells to perform DE
#' @param group.by Regroup cells before performing DE
#' @param clusters Vector of Clusters to perform DE
#'
#' @param return A list.
#' @example
#' object@misc$xDE <- xDE(object = object, comparison = c("Normal", "Tumor"), group.by="group", min.cells = 20, clusters = NULL)
#'
#' @author rstatistics
#' @export
xDE <- function(object = object, comparison = c("condA", "condB"), group.by=NULL, min.cells = 20, clusters = NULL){
  results <- list()
  if (is.null(clusters))
    clusters = levels(object)
  if (length(comparison) != 2){
    stop("Comparison must have 2 elements!")
  }
  condA <- comparison[1]
  condB <- comparison[2]
  if (is.null(group.by)){
    group.by <- "group"
  }
  results[[paste0(condA, "_vs_", condB)]] <- list()
  clusters <- as.character(clusters)
  for (cluster in clusters){
    cat(paste0("### ", "Comparing Cluster-", cluster, " in ", condA, " and ", condB, " ...\n"))
    Object <- subset(object, seurat_clusters==cluster)
    results[[paste0(condA, "_vs_", condB)]][[cluster]] = list()
    results[[paste0(condA, "_vs_", condB)]][[cluster]] <- suppressWarnings(
      FindMarkers(object = Object, assay = 'RNA', ident.1 = condA, ident.2 = condB, group.by = group.by,
                  min.cells.group = min.cells, test.use = 'MAST', ...))
    cat("done.\n")
  }
  return(results)
}
