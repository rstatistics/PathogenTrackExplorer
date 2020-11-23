#' Cluster GO annotation
#' GO annotation for each cluster
#' @title GO annotation
#' @param object Seurat object
#' @param key The DEGs are stored in object@misc$key
#' @param clusters A vector of clusters, default clusters=NULL
#' @param only.pos Only use Upregulated genes. Default only.pos=TRUE
#' @return a list object that stored GO results.
#' @example
#' object@misc$DE <- FindAllMarkers(object = object, assay = 'RNA', only.pos = TRUE, test.use = 'MAST')
#'
#' object@misc$GO <- xGO(object = object, key="DE", clusters=NULL)
#'
#'
#' @author rstatistics
#' @export

xGO <- function(object = object, key="DE", clusters=NULL, only.pos=TRUE){
  requireNamespace("topGO")
  requireNamespace("org.Hs.eg.db")
  if (is.null(key)){
    stop("Parameter \"key\" must be specified!\n")
  }
  markers <- object@misc[[key]]
  if (!is.null(clusters)){
    markers <- subset(markers, cluster %in% clusters)
  }
	x <- as.list(org.Hs.egALIAS2EG)
	geneList <- rep(0, length(rownames(object)))
	names(geneList) <- rownames(object)
	geneList <- geneList[intersect(names(geneList), names(x))]
	TotalGenes <- names(geneList)
	for (ii in 1:length(geneList)){
	  names(geneList)[ii] = x[[names(geneList)[ii]]][1]
	}
	go_enrichment_results = list()

	for (c in as.character(unique(markers$cluster))){
	  print(paste0("Running cluster ", c))
	  if (only.pos==TRUE){
	    queryMatrix = subset(markers, cluster == c & pct.2<=pct.1 & avg_logFC > 0.585 & p_val_adj < 1e-2)
	  }else{
	    queryMatrix = subset(markers, cluster == c & abs(avg_logFC) > 0.585 & p_val_adj < 1e-2)
	  }
	  queryGene = rownames(queryMatrix)
	  # next if queryGene with less than 10 genes
	  if(length(queryGene) < 10){ next }
	  # order by pct.2/pct.1
	  queryMatrix = queryMatrix[order(queryMatrix$pct.2 / queryMatrix$pct.1, decreasing = FALSE),]
	  # only use top 50 genes to do GO
	  # queryGene = queryGene[1:ifelse(length(queryMatrix$gene)>=50, 50, length(queryMatrix$gene))]
	  go_enrichment_results[[c]] = list()
	  # Run against topGO ####
	  queryGeneList = geneList
	  queryGeneList[which(TotalGenes %in% queryGene)] = 1
	  # skip if no query gene left
	  if (all(queryGeneList==0)){ next }
	  tab = list()
	  for (ont in c("BP", "CC", "MF")){
	    GOdata <- suppressMessages(new("topGOdata", ontology = ont, allGenes = as.factor(queryGeneList),
	                  nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "entrez"))
	    resultTopGO.elim <- suppressMessages(runTest(GOdata, algorithm = "elim", statistic = "Fisher"))
	    tab[[ont]] <- rbind(tab[[ont]], GenTable(GOdata, Fisher.elim = resultTopGO.elim,
	                               orderBy = "Fisher.elim", topNodes = 200))
	  }
	  go_enrichment_results[[c]] = tab
	}
	return(go_enrichment_results)
}
