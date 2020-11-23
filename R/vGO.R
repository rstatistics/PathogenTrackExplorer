#' vGO
#' vGO annotation for specific cluster
#' @title vGO annotation
#' @param object Seurat object
#' @param key The DEGs are stored in object@misc$vDE$key
#' @param cluster Cluster specified by user
#' @param feature.by Feature which is used to divid cells into two groups
#' @return a list object that stores vGO results
#' @example
#' object@misc$vDE <- vDE(object = object, features.by = features.by, min.cells = 20, clusters = clusters)
#'
#' object@misc$vGO <- vGO(object = object, key="vDE", cluster=cluster, feature.by=feature.by)
#'
#' @author rstatistics
#' @export

vGO <- function(object = object, key="vDE", cluster=cluster, feature.by=feature.by){
  requireNamespace("topGO")
  requireNamespace("org.Hs.eg.db")

  if (is.null(key)){
    stop("Parameter \"key\" must be specified!\n")
  }
  DEdata = object@misc[[key]]
  x <- as.list(org.Hs.egALIAS2EG)
  geneList <- rep(0, length(rownames(object)))
  names(geneList) <- rownames(object)
  geneList <- geneList[intersect(names(geneList), names(x))]
  TotalGenes <- names(geneList)
  for (ii in 1:length(geneList)){
    names(geneList)[ii] = x[[names(geneList)[ii]]][1]
  }
  go_enrichment_results = list()

  cluster=as.character(cluster)
  if (! cluster %in% names(DEdata)){ stop(paste0("Cluster-", cluster, " does not exist!\n")) }
  if (! feature.by %in% names(DEdata[[cluster]])) { stop(paste0(feature.by, " does not exist in Cluster- ", cluster, "!\n")) }
  go_enrichment_results[[cluster]] = list()
  go_enrichment_results[[cluster]][[feature.by]] = list()

  print(paste0("Running ", feature.by, " in Cluster-", cluster, " ... "))
  rawMatrix = DEdata[[cluster]][[feature.by]]

  for (direction in c("up", "dn")){
    if (direction == "up"){
      queryMatrix = subset(rawMatrix[order(rawMatrix$pct.1 / rawMatrix$pct.2, decreasing = TRUE),], avg_logFC > 0.585 & p_val_adj < 1e-2)
    }else{
      queryMatrix = subset(rawMatrix[order(rawMatrix$pct.2 / rawMatrix$pct.1, decreasing = TRUE),], avg_logFC < -0.585 & p_val_adj < 1e-2)
    }
    # only use top 50 genes to do GO
    queryGene = rownames(queryMatrix)
    # next if queryGene with only one Gene
    if(length(queryGene) < 2){ next }
    # Run against topGO ####
    queryGeneList = geneList
    queryGeneList[which(TotalGenes %in% queryGene)] = 1
    if(length(table(queryGeneList)) < 2){ next }
    queryMatrix$gene = rownames(queryMatrix)
    queryGene = queryMatrix$gene
    go_enrichment_results[[cluster]][[feature.by]][[direction]] = list()

    tab = list()
    for (ont in c("BP", "CC", "MF")){
      GOdata <- suppressMessages(new("topGOdata", ontology = ont, allGenes = as.factor(queryGeneList),
                                     nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "entrez"))
      resultTopGO.elim <- suppressMessages(runTest(GOdata, algorithm = "elim", statistic = "Fisher"))
      tab[[ont]] <- rbind(tab[[ont]], GenTable(GOdata, Fisher.elim = resultTopGO.elim,
                                               orderBy = "Fisher.elim", topNodes = 200))
    }
    go_enrichment_results[[cluster]][[feature.by]][[direction]] = tab
  }
  print("done.")
  return(go_enrichment_results)
}
