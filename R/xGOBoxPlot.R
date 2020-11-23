#' GO BoxPlot visualization
#' GO BoxPlot visualization for all clusters
#' @title GOBoxPlot visualization
#' @param object Seurat object
#' @param key The DEGs are stored in object@misc$key
#' @param ont Ontology to use, only "BP", "CC", and "MF".
#' @param top_n Only top_n entries to plot.
#' @return
#' @example
#' object@misc$DE <- FindAllMarkers(object = object, assay = 'RNA', only.pos = TRUE, test.use = 'MAST')
#'
#' object@misc$GO <- xGO(object=object, key="DE")
#'
#' xGOBoxPlot(object, key="GO", ont="BP", top_n=20)
#'
#'
#' @author rstatistics
#' @export

xGOBoxPlot <- function(object = object, key = "GO", ont = "BP", top_n = 5){
  requireNamespace("ggplot2")
  if (is.null(key)){
    stop("Parameter \"key\" must be specified!\n")
  }
  GOdata = object@misc[[key]]
  if (is.null(top_n)){
    top_n = 5
  }
  GO.terms = vector()
  for (i in 1:length(GOdata)){
    term <- GOdata[[i]]$BP$Term[1:top_n]
    GO.terms <- unique(c(GO.terms, term))
  }

  GO.data <- data.frame(Cluster=character(), Term=character(), Value=numeric())
  for (i in names(GOdata)){
    z = GOdata[[i]]$BP
    if (is.null(z)){ next }
    for (GO.term in GO.terms){
      if(nrow(subset(z, Term==GO.term))==0){
        GO.value = NA
      }else{
        GO.value = -log10(as.numeric(gsub("< ", "", subset(z, Term == GO.term)[1,]$Fisher.elim)))
      }
      GO.value <- ifelse(is.na(GO.value), 0, GO.value)
      GO.row <- data.frame(Cluster=i, Term=GO.term, Value=GO.value)
      GO.data <- rbind(GO.data, GO.row)
    }
  }

  # fix the order of go terms
  GO.data$Term <- factor(GO.data$Term, levels = GO.terms)
  # fix the order of clusters
  GO.data$Cluster <- factor(GO.data$Cluster, levels = names(GOdata))
  # normalize the significant value
  GO.data$Value <- ifelse(GO.data$Value > 10, 10, GO.data$Value)
  GO.data$Value <- ifelse(GO.data$Value < 2, 0, GO.data$Value)
  p <- ggplot(GO.data, aes(x=Cluster, y=Term)) + geom_tile(aes(fill=Value)) + coord_equal()
  p <- p + xlab("Cluster") + ylab("GO Terms") + labs(fill = "-log10(Pvalue)") +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                  axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, colour="black"),
                  axis.text = element_text(size = 6), axis.line = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.grid = element_blank(), axis.ticks = element_blank(),
                  panel.border = element_rect(colour="black",fill=NA,size=1),
                  legend.direction = "vertical") +
            scale_fill_gradientn(colours=colorRampPalette(c("white", "red"))(n = 100)) +
            scale_y_discrete(limits=rev(levels(GO.data$Term)), )
  return(p)
}
