#' Population Plot
#' A stacked bar plot to show the fraction composition of each sample or cluster
#' @param object Seurat object
#' @param by String used to separate cluster/sample, only "cluster" or "sample" is accepted
#' @param cols Vector of colors, each color corresponds to an identity class
#' @return
#' @export
xPopulationPlot <- function(object, by, cols){
  if(is.null(by)){
    by = "sample"
  }
  data <- table(Idents(object), object$orig.ident)
  ClusterNames <- rownames(data)
  data <- as.data.frame(data)
  names(data) <- c("Cluster", "Sample", "Cells")
  data$Cluster <- as.character(data$Cluster)
  data$Cluster <- factor(data$Cluster, levels = ClusterNames)
  data$Sample <- factor(data$Sample, levels = unique(object$sample))
  if (by=="sample"){
    p <- ggplot(data=data, aes(x = Sample, y = Cells, fill = Cluster)) +
      geom_bar(stat = "identity", width = 0.6, position = position_fill(reverse = FALSE), size = 0.3, colour = NA) +
      labs(x = 'Sample', y = 'Fraction of cells (%)') + cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(), axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 12), panel.grid.minor = element_blank()) +
      scale_fill_manual(values = cols, guide = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=1)) +
      scale_x_discrete(expand = c(0.05, 0.05)) + scale_y_continuous(expand = c(0.01, 0.01), labels = c(0,25,50,75,100))
  }else if (by=="cluster"){
    p <- ggplot(data=data, aes(x = Cluster, y = Cells, fill = Sample)) +
      geom_bar(stat = "identity", width = 0.6, position = position_fill(reverse = FALSE), size = 0.3, colour = NA) +
      labs(x = 'Cluster', y = 'Fraction of cells (%)') + cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(), axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 12), panel.grid.minor = element_blank()) +
      scale_fill_manual(values = cols, guide = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=1)) +
      scale_x_discrete(expand = c(0.05, 0.05)) + scale_y_continuous(expand = c(0.01, 0.01), labels = c(0,25,50,75,100))
  }else{
    stop("Parameter \"by\" must be \"cluster\" or \"sample\"!\n")
  }
  return(p)
}


