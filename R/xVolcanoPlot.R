#' VolcanoPlot for specific group
#' @title VolcanoPlot visualization
#' @param object A seurat object.
#' @param key Key name. The DE matrix are stored in object@misc$key
#' @param cluster Cluster specified by the user
#' @param top_n Number of top up/down genes to show.
#' @example
#' object@misc$DE <- FindAllMarkers(object = object, assay = 'RNA',
#'                                  only.pos = TRUE, test.use = 'MAST')
#'
#' xVolcanoPlot(object=object, key=key, cluster=cluster, top_n=NULL)
#'
#' @author rstatistics
#' @export
xVolcanoPlot <- function(object = object, key = key, cluster = cluster, top_n = NULL){
  if (is.null(top_n)){ top_n = 5 }
  x <- object@misc[[key]][[cluster]]
  if(is.null(x)){ stop("No differentially expressed genes in cluster ", cluster, "!\n") }
  x$gene <- rownames(x)
  logFCcut <- 0.585
  pvalCut <- 1e-3
  adjPcut <- 1e-5
  logFCcut2 <- 1.0
  logFCcut3 <- 2.0
  pvalCut2 <- 1e-20
  pvalCut3 <- 1e-10
  xrange <- range(x$avg_logFC)
  xrange <- xrange[is.finite(xrange)]
  xmin <- floor(min(xrange))
  xmax <- ceiling(max(xrange))
  xmin <- -max(abs(xmin), abs(xmax))
  xmax <- max(abs(xmin), abs(xmax))
  yrange <- -log10(x$p_val_adj)
  yrange <- yrange[is.finite(yrange)]
  ymin <- 0
  ymax <- ceiling(max(yrange)) * 1.1
  mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  cols <- rep("grey40", n1)
  names(cols)<- rownames(x)
  # set colors for each threshold
  cols[x$p_val_adj < pvalCut & x$avg_logFC>logFCcut]<- "#FB9A99"
  cols[x$p_val_adj < pvalCut2 & x$avg_logFC> logFCcut2]<- "#ED4F4F"
  cols[x$p_val_adj < pvalCut & x$avg_logFC< -logFCcut]<- "#B2DF8A"
  cols[x$p_val_adj < pvalCut2 & x$avg_logFC< -logFCcut2]<- "#329E3F"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent

  n1 <- length(x[, 1])
  size <- rep(1, n1)

  # set dot size for each threshold
  size[x$p_val_adj < pvalCut & x$avg_logFC> logFCcut]<- 1
  size[x$p_val_adj < pvalCut2 & x$avg_logFC> logFCcut2]<- 2
  size[x$p_val_adj < pvalCut3 & x$avg_logFC> logFCcut3]<- 3
  size[x$p_val_adj < pvalCut & x$avg_logFC< -logFCcut]<- 1
  size[x$p_val_adj < pvalCut2 & x$avg_logFC< -logFCcut2]<- 2
  size[x$p_val_adj < pvalCut3 & x$avg_logFC< -logFCcut3]<- 3

  # build the plot object
  p <- ggplot(data=x, aes(avg_logFC, -log10(p_val_adj), label = gene)) +
    geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
    labs(x = "Log2(Fold Change)", y = "-Log10(P-value)", title = paste0(key, " in ", "Cluster-", cluster)) +
    ylim(c(ymin, ymax)) +
    scale_x_continuous(
      breaks = c(-3, -2, -1, -logFCcut, 0, logFCcut, 1, 2, 3),
      labels = c(-3, -2, -1, -logFCcut, 0, logFCcut, 1, 2, 3),
      limits = c(xmin, xmax)
    ) +
    geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey60",
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCut), color="grey60",
               linetype="longdash", lwd = 0.5) +

    theme_bw(base_size = 12, base_family = "Times") +
    theme(panel.grid=element_blank()) +
    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey20",
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCut2), color="grey20",
               linetype="longdash", lwd = 0.5)

  x_up <- subset(x, avg_logFC>0) %>% head(top_n)
  x_dn <- subset(x, avg_logFC<0) %>% head(top_n)
  x <- rbind(x_up, x_dn)

  p <- p +
    geom_point(data = x, alpha = 1, size = 3, shape = 1,
               stroke = 0.5, color = "black") +
    scale_color_manual(values = mycol) +
    geom_text_repel(data = x,
                    show.legend = FALSE,
                    size = 4, box.padding = unit(0.25, "lines"),
                    point.padding = unit(0.3, "lines")) +
    guides(color=guide_legend(title = NULL))

  return(p)
}
