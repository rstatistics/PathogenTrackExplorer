#' vGOBarPlot
#' vGO BarPlot visualization for specific cluster
#' @title GOBarPlot visualization
#' @param object A seurat object
#' @param key Key name. The vGO results are stored in object@misc$key
#' @param cluster Cluster specified by the user
#' @param feature.by Feature which is used to divid cells into two groups
#' @param ont Ontology to use, only "BP", "CC", or "MF".
#' @param top_n Only top_n entries to plot
#' @return
#' @example
#' object@misc$vGO <- vGO(object=object, key="vDE", cluster=cluster, feature.by=feature.by)
#'
#' vGOBarPlot(object=object, key="vGO", cluster=cluster, feature.by=feature.by, ont="BP", top_n=20)
#'
#' @author rstatistics
#' @export
vGOBarPlot <- function(object=object, key="vGO", cluster=cluster, feature.by=feature.by, ont="BP", top_n=20){
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
  if (! cluster %in% names(GOdata)) { stop(paste0("Cluster-", cluster, " does not exist!\n")) }
  if (! feature.by %in% names(GOdata[[cluster]])) { stop(paste0(feature.by, " does not exist in Cluster- ", cluster, "!\n")) }
  z_up <- GOdata[[cluster]][[feature.by]][["up"]][[ont]]
  z_up$color <- "#FF0000"
  z_up$Log10FDR <- -log10(as.numeric(gsub("< ", "", z_up$Fisher.elim)))
  z_dn <- GOdata[[cluster]][[feature.by]][["dn"]][[ont]]
  z_dn$color <- "#0000FF"
  z_dn$Log10FDR <- log10(as.numeric(gsub("< ", "", z_dn$Fisher.elim)))
  z_duplicate <- intersect(z_up$Term, z_dn$Term)
  z <- rbind(head(z_up, top_n), head(z_dn, top_n))
  z <- z[!z$Term %in% z_duplicate, ]
  x <- data.frame(Term=z$Term, Log10FDR=z$Log10FDR, color=z$color)
  ymax <- ceiling(max(abs(z$Log10FDR)))
  p <- ggplot(x, aes(x=reorder(Term, Log10FDR), y=Log10FDR)) + labs(x = "GO Terms", y = "Log10(P-adjust)") +
    ggtitle(label = paste0(feature.by, " in Cluster-", cluster)) +
    theme(plot.title = element_text(color = "red", size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(color = "blue", face = "bold", hjust = 0.5)) +
    geom_bar(stat="identity", fill=z$color, color="NA", width = 0.8) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, color = "black") +
    geom_hline(yintercept = 2, linetype = "dashed", size = 0.5, color = "gray50") +
    geom_hline(yintercept = -2, linetype = "dashed", size = 0.5, color = "gray50") +
    coord_flip(ylim = c(-ymax,ymax)) + theme_bw()
  p + theme(panel.grid.major =element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title = element_text(family = "Times", size = 12, colour = "black"),
            axis.text = element_text(family = "Times", colour = "black"),
            axis.line.x = element_line(colour = "black"),
            axis.ticks.y = element_blank()) + scale_y_continuous(expand = c(0, 0.015))

  return(p)
}
