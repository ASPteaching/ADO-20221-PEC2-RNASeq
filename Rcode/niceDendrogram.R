#' A template for plotting Dendrograms
#' Preparat per Jordi Cortés
#' @param hcObj The results of a hierarchical clustering
#' @return The (gg)plot of the Dendrogram
#' @examples
#' niceDendrogram (hc)
#' @export
niceDendrogram <- function (hcObj){
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols = gg_color_hue(2)
  gg2 <- fviz_dend(hc, k = 2,rect=FALSE,lwd=2,k_colors = c("grey","grey"),
                   label_cols =  ifelse(substr(colnames(log2count_norm)[hc$order],1,1)=='C',
                                        cols[1],cols[2]))
  gg2 <- gg2 +
    theme(axis.title = element_text(face='bold'),
          title = element_blank())
  gg2
  # ggsave('figures/dendograma.png',gg2) 
}
