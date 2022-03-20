#' A template for plotting functions documented with rOxygen
#'
#' @param x Data matrix for the Heatmap
#' @return The heatmap as a \code{plotly} interactive object \code{p}.
#' @examples
#' # Not run
#' # niceHeatMap(mat)
#' @export
niceInteractiveHeatMap <- function(x) {
  require(heatmaply)
  p <- heatmaply(data.matrix(x),
                 #dendrogram = "row",
                 xlab = "", ylab = "",
                 main = "",
                 colorRampPalette(c("blue", "cyan", "white", "orange", "red"))(100),
                 scale = "none",
                 margins = c(60,100,40,20),
                 #grid_color = "white",
                 #grid_width = 0,
                 titleX = FALSE,
                 hide_colorbar = FALSE,
                 branches_lwd = 0.1,
                 label_names = c("Feature", "Sample", "Value"),
                 fontsize_row = 5, fontsize_col = 10,
                 labCol = colnames(data.matrix(mat)),
                 labRow = rownames(data.matrix(mat)),
                 heatmap_layers = theme(axis.line=element_blank())
  )
  p
  # save the widget
  # library(htmlwidgets)
  # saveWidget(p, file= "~/.../heatmap.html")  
}

