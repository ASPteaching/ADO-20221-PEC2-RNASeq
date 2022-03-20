#' Visualització de MDS amb ggplot2
#' Preparat per Jordi Cortés
#' 
#' @param aDGEObj An object of class DGEList created by edgeR
#' @return The (gg)plot of the Dendrogram
#' @examples
#' niceDendrogram (hc)
#' @export
# Function to plot a MDS from 
niceMDS <- function (aDGEObj){
  require(magrittr)
  require(tidyverse)
  require(ggplot2)
  require(metafolio)
  
  mds <- plotMDS.DGEList(aDGEObj, plot = FALSE)
  
  gg3 <- data.frame(x=mds$x,y=mds$y,
                    co=ifelse(substr(colnames(log2count_norm),1,3)=='COV','COVID','HEALTHY'),
                    ind=substr(colnames(log2count_norm),4,6)) %>%
    set_colnames(c("Dim1", "Dim2","co","ind")) %>%
    rownames_to_column("SampleID") %>%
    ggplot(aes(x = Dim1, y = Dim2,col=co,label=ind)) +
    geom_point(size = 3) +
    geom_text(hjust=-0.3, vjust=0,check_overlap = TRUE, show.legend = FALSE) +
    scale_colour_manual(values = gg_color_hue(2)) +
    xlab(paste0("Dimension 1 (",
                formatC(100*mds$var.explained[1],di=1,fo='f'),'%)')) + 
    ylab(paste0("Dimension 2 (",
                formatC(100*mds$var.explained[2],di=1,fo='f'),'%)')) +
    theme(legend.title=element_blank(),
          legend.position = "bottom",
          axis.title = element_text(face='bold'))
  gg3  
  # ggplot2::ggsave('figures/MDS.png',gg3) 
}
