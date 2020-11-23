#' GO DotPlot visualization
#' GO DotPlot visualization for all clusters
#' @title GODotPlot visualization
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
#' xGODotPlot(object=object, key="GO", ont="BP", top_n=5)
#'
#'
#' @author rstatistics
#' @export

xGODotPlot <- function(object = object, key = "GO", ont = "BP", top_n = 5){
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
    term <- GOdata[[i]][[ont]]$Term[1:top_n]
    GO.terms <- unique(c(GO.terms, term))
  }

  GO.data <- data.frame(Cluster=character(), Term=character(), GeneNum=character(), Value=numeric())
  for (i in names(GOdata)){
    z = GOdata[[i]][[ont]]
    if (is.null(z)){ next }
    for (GO.term in GO.terms){
      if(nrow(subset(z, Term==GO.term))==0){
        GO.value = NA
        GeneNum = 0
      }else{
        GO.value = -log10(as.numeric(gsub("< ", "", subset(z, Term == GO.term)[1,]$Fisher.elim)))
        GeneNum = as.integer(subset(z, Term == GO.term)[1,]$Annotated)
      }
      GO.value <- ifelse(is.na(GO.value), 0, GO.value)
      GO.row <- data.frame(Cluster=i, Term=GO.term, GeneNum=GeneNum, Value=GO.value)
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

  breaks = c(seq(0,floor(max(GO.data$Value)),length.out = 5))
  p <- ggplot(GO.data, aes(x=reorder(Cluster, Value, median), y=Term)) + geom_point(aes(colour=Value, size=GeneNum), stat="identity") +
    theme(plot.title=element_text(face="bold",vjust=1.0), axis.title.x=element_text(face="bold",vjust=-0.2),
          axis.title.y=element_text(face="bold"), axis.text.y=element_text(hjust=1.0,colour="black"),
          axis.text.x=element_text(angle=45, vjust = 1, hjust = 1, colour="black"), panel.background=element_blank(),
          axis.ticks = element_blank(), axis.text = element_text(size = 7), panel.border=element_rect(colour="black",fill=NA,size=1),
          panel.grid.minor=element_line(colour="grey", linetype="dotted", size=0.1),
          panel.grid.major=element_line(colour="grey", linetype="dotted", size=0.1)) +
    scale_colour_gradient(low="white",high="red",breaks=breaks,labels=as.character(breaks),limits=c(0,max(breaks))) +
    labs(x = "Cluster", y = "GO Terms", color = "-Log10(Pvalue)") + scale_y_discrete(limits=rev(levels(GO.data$Term)))
  return(p)
}
