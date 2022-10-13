#' A template for plotting Dendrograms
#' Preparat per Jordi Cortés
#' @param hcObj The results of a hierarchical clustering
#' @return The (gg)plot of the Dendrogram
#' @examples
#' niceDendrogram (hc)
#' @export

# niceDendrogram <- function (hcObj){
#   gg_color_hue <- function(n) {
#     hues = seq(15, 375, length = n + 1)
#     hcl(h = hues, l = 65, c = 100)[1:n]
#   }
#   cols = gg_color_hue(2)
#   gg2 <- fviz_dend(hc, k = 2,rect=FALSE,lwd=2,k_colors = c("grey","grey"),
#                    label_cols =  ifelse(substr(colnames(log2count_norm)[hc$order],1,1)=='C',
#                                         cols[1],cols[2]))
#   gg2 <- gg2 +
#     theme(axis.title = element_text(face='bold'),
#           title = element_blank())
#   gg2
#   # ggsave('figures/dendograma.png',gg2) 
# }

niceDendrogram <- function (hcObj){
  hc <- hclust(sampleDists,method = 'ward.D2')
  require("ggdendro", "dendextend")
  dend <- as.dendrogram(hc)
  dend_data <- dendro_data(dend, type = "rectangle")
  dend_data$labels$color <- ifelse(grepl('COV',dend_data$labels$label),'COVID','HEALTHY')
  head(dend_data$segments)
  gg2 <- ggplot(dend_data$segments) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend),size=1.2)+
    geom_text(data = dend_data$labels, aes(x, y-1, label = label,color=color,fontface='bold'),
              hjust = 1, angle = 90, size = 8)+
    ylim(-20, 150) + ylab('Dissimilarity') + xlab('') +
    theme(axis.ticks.length.x = unit(0,units = "cm"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10),
          legend.position = "none",
          text=element_text(face='bold'),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"))
  gg2
}