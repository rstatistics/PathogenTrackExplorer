#' GO BarPlot visualization
#' GO BarPlot visualization for each cluster
#' @title GOBarPlot visualization
#' @param object Seurat object
#' @param key The vGO results are stored in object@misc$key
#' @param cluster Cluster specified by the user
#' @param ont Ontology to use, only "BP", "CC", or "MF".
#' @param top_n Only top_n entries to plot
#' @return
#' @example
#' object@misc$DE <- FindAllMarkers(object = object, assay = 'RNA', only.pos = TRUE, test.use = 'MAST')
#'
#' object@misc$GO <- xGO(object=object, key="DE")
#'
#' xGOBarPlot(object=object, key="GO", cluster=cluster, ont="BP", top_n=20)
#'
#'
#' @author rstatistics
#' @export


xGOBarPlot <- function(object=object, key="GO", cluster=cluster, ont="BP", top_n=20){
  requireNamespace("ggplot2")
  if (is.null(key)){
    stop("Parameter \"key\" must be specified!\n")
  }
  GOdata = object@misc[[key]]
  cluster = as.character(cluster)
  if (! cluster %in% names(GOdata)){
    stop("cluster ", cluster, " does not exist!")
  }
  if (is.null(ont)){
    ont = "BP"
  }
  if (! ont %in% c("BP", "CC", "MF")){
    stop(paste0("ont must be \"BP\", \"CC\" or \"MF\"."))
  }
  z <- GOdata[[cluster]][[ont]]
  if (is.null(z)){ stop("No data available!\n") }
  z <- head(z, top_n)
  x <- data.frame(Term=z$Term, FDR=as.numeric(gsub("< ", "", z$Fisher.elim)))
  p <- ggplot2::ggplot(x,aes(x=reorder(Term,-log10(FDR)),y=-log10(FDR))) +
    geom_bar(stat="identity", fill="red", color="NA",width = 0.8) +
    geom_hline(yintercept = 2, linetype = "dashed", size = 0.5, color = "grey50") +
    coord_flip() + theme(panel.grid.major =element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.title = element_text(family = "ArialMT", size = 12, colour = "black"),
                         axis.text = element_text(family = "ArialMT", colour = "black"),
                         axis.line.x = element_line(colour = "black"),
                         axis.line.y = element_line(colour = "black"),
                         axis.ticks.y = element_blank()) + scale_y_continuous(expand = c(0, 0.015)) +
    labs(x = "GO Terms", y = "-Log10(Pvalue)")
  return(p)
}
