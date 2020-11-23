#' xDotPlot
#' wrapper for Seurat::DotPlot
#' @title xDotPlot
#' @param object Seurat object
#' @param features Vector of features to plot
#' @param cols Vector of two colors to form the gradient over
#' @param dot.scale Scale the size of the dot
#' @param font.size Font size
#'
#' @return
#' @example
#' object@misc$DE %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC) -> top3
#' xDotPlot(object, features)
#'
#' @author rstatistics
#' @export
xDotPlot <- function(object, features, cols=NULL, dot.scale=NULL, font.size=NULL, ...){
  features <- rev(features)
  features <- intersect(unique(features), rownames(object))
  if (is.null(cols)){
    cols <- c('#FFFFFF','#16388E')
  }
  if (is.null(dot.scale)){
    dot.scale <- 3
  }
  if (is.null(font.size)){
    font.size <- 8
  }
  p <- Seurat::DotPlot(object=object, features = features, cols = cols, dot.scale = dot.scale) +
    RotatedAxis() + theme(axis.text.x = element_text(size = font.size), axis.text.y = element_text(size = font.size)) +
    guides(color = guide_colorbar(title = 'Scaled Expression'), size = guide_legend(title = 'Percent Expressed')) +
    theme(axis.line = element_line(size = 0.6)) + labs(x='', y='')
  return(p)
}
