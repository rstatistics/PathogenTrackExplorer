#' xFeaturePlot
#' wrapper for Seurat::FeaturePlot
#' @title xFeaturePlot
#' @param object Seurat object
#' @param reduction Dimensional reduction to use
#' @param features Vector of features to plot
#' @param pt.size Point size to plot
#' @param font.size Font size
#' @param cols Vector of two colors to form the gradient over
#' @param ncols Number of columns to combine features
#' @param combine Combine plots into a single patchworked ggplot object
#' @param order Determine whether to plot cells in order of expression
#'
#' @return
#' @example
#' xFeaturePlot(object=object, reduction="umap", features=features)
#'
#' @author rstatistics
#' @export

xFeaturePlot <- function(object=object, reduction=NULL, features=features, pt.size=NULL, font.size=NULL,
                         cols=NULL, ncol=NULL, combine=FALSE, order=FALSE, ...){
  if(is.null(reduction)){
    reduction = "umap"
  }
  if(is.null(pt.size)){
    pt.size = 0.1
  }
  if(is.null(cols)){
    cols = c("lightgrey","#FF0000")
  }
  if(is.null(font.size)){
    font.size = 9
  }
  features <- intersect(unique(features), rownames(object))
  if(is.null(ncol)){
    ncol = ceiling(sqrt(length(features)))
  }
  pp = Seurat::FeaturePlot(object = object, reduction = reduction, features = features, pt.size = pt.size,
                           cols = cols, combine = combine, order = order, ...)
  plots <- lapply(X = pp, FUN = function(p){
    p + theme(axis.title = element_text(size = font.size), axis.text = element_text(size = font.size),
              plot.title = element_text(family = 'sans',face='italic',size=10), legend.text = element_text(size = 10),
              legend.key.height = unit(0.9,"line"), legend.key.width = unit(0.6,"line"))})
  return(patchwork::wrap_plots(plots, ncol = ncol))
}
