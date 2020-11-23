#' vGODotPlot
#' vGO DotPlot visualization for specific cluster
#' @title GODotPlot visualization
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
#' vGODotPlot(object=object, key="vGO", cluster=cluster, feature.by=feature.by, ont="BP", top_n=20)
#'
#' @author rstatistics
#' @export
vGODotPlot <- function(object=object, key="vGO", cluster=cluster, feature.by=feature.by, ont="BP", top_n=20){
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
  z_up$Group <- "Up"
  z_up$Pvalue <- as.numeric(gsub("< ", "", z_up$Fisher.elim))
  z_up <- z_up[order(z_up$Pvalue, decreasing=FALSE), ]
  z_dn <- GOdata[[cluster]][[feature.by]][["dn"]][[ont]]
  z_dn$Group <- "Dn"
  z_dn$Pvalue <- as.numeric(gsub("< ", "", z_dn$Fisher.elim))
  z_dn <- z_dn[order(z_dn$Pvalue, decreasing=FALSE), ]
  z_duplicate <- intersect(z_up$Term, z_dn$Term)
  z <- rbind(head(z_up, top_n), head(z_dn, top_n))
  z <- z[!z$Term %in% z_duplicate, ]
  z$Term <- factor(z$Term, levels=z$Term)
  z$Group <- factor(z$Group, levels=c("Up", "Dn"))
  x <- data.frame(Term=z$Term, Group=z$Group, GeneNumber=z$Annotated, Pvalue=z$Pvalue)

  p <- ggplot(x, aes(x=Group, y=Term, size=GeneNumber)) + geom_point(aes(colour=-log10(Pvalue))) +
    scale_color_gradient(low = "orange", high = "red", na.value = NA) + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.title = element_text(family = "Times", size = 12, colour = "black"),
                 axis.text = element_text(family = "Times", colour = "black"),
                 axis.line.x = element_line(colour = "black"),
                 axis.ticks.y = element_blank()) + labs(x = "Direction", y = "GO Terms") +
    ggtitle(label = paste0(feature.by, " in Cluster-", cluster)) +
    theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(color = "blue", face = "bold", hjust = 0.5)) +
    scale_y_discrete(limits=rev(levels(x$Term)))
  return(p)
}
