#' Visualització de MDS amb ggplot2
#' Preparat per Jordi Cortés
#' 
#' @param aDGEObj An object of class DGEList created by edgeR
#' @return The (gg)plot of the Dendrogram
#' @examples
#' niceDendrogram (hc)
#' @export
# Function to plot a MDS from 
# niceMDS <- function (aDGEObj){
#   require(magrittr)
#   require(tidyverse)
#   require(ggplot2)
#   require(metafolio)
#   
#   mds <- plotMDS.DGEList(aDGEObj, plot = FALSE)
#   
#   gg3 <- data.frame(x=mds$x,y=mds$y,
#                     co=ifelse(substr(colnames(log2count_norm),1,3)=='COV','COVID','HEALTHY'),
#                     ind=substr(colnames(log2count_norm),4,6)) %>%
#     set_colnames(c("Dim1", "Dim2","co","ind")) %>%
#     rownames_to_column("SampleID") %>%
#     ggplot(aes(x = Dim1, y = Dim2,col=co,label=ind)) +
#     geom_point(size = 3) +
#     geom_text(hjust=-0.3, vjust=0,check_overlap = TRUE, show.legend = FALSE) +
#     scale_colour_manual(values = gg_color_hue(2)) +
#     xlab(paste0("Dimension 1 (",
#                 formatC(100*mds$var.explained[1],di=1,fo='f'),'%)')) + 
#     ylab(paste0("Dimension 2 (",
#                 formatC(100*mds$var.explained[2],di=1,fo='f'),'%)')) +
#     theme(legend.title=element_blank(),
#           legend.position = "bottom",
#           axis.title = element_text(face='bold'))
#   gg3  
#   # ggplot2::ggsave('figures/MDS.png',gg3) 
# }

niceMDS <- function (aDGEObj){
  require(magrittr)
  require(emojifont)
  require(ggpubr)
  
  mds <- plotMDS.DGEList(aDGEObj,plot = FALSE)
  ddd <- data.frame(x=mds$x,y=mds$y,
                    co=ifelse(substr(colnames(log2count_norm),1,3)=='COV','COVID','HEALTHY'),
                    ind=substr(colnames(log2count_norm),4,6)) %>%
    set_colnames(c("Dim1", "Dim2","co","ind")) %>%
    rownames_to_column("SampleID")
  ddd <- ddd %>%
    mutate(label=if_else(co=='COVID',
                         emoji('sneezing_face'),
                         emoji('smile')))
  
  gg4_1 <- ggplot(ddd,aes(x = Dim1, y = Dim2,col=label,label=label)) +
    geom_text(family="EmojiOne", size=6) +
    scale_colour_manual(values = rev(gg_color_hue(2))) +
    xlab(paste0("Dimension 1 (",
                formatC(100*mds$var.explained[1],di=1,fo='f'),'%)')) + 
    ylab(paste0("Dimension 2 (",
                formatC(100*mds$var.explained[2],di=1,fo='f'),'%)')) +
    theme(legend.title=element_blank(),
          legend.position = "bottom",
          axis.title = element_text(face='bold'),
          legend.text=element_text(family='EmojiOne'))
  
  ##-- awesome fav icons -------------------------------------------
  # https://stackoverflow.com/questions/56605100/unable-to-load-fontawesome-icons-in-shiny-app
  # https://fontawesome.com/v4/icons/
  
  ddd <- ddd %>%
    mutate(label=if_else(co=='COVID',
                         fontawesome('fa-ambulance'),
                         fontawesome('fa-heart-o')))
  gg4_2 <- ggplot(ddd,aes(x = Dim1, y = Dim2, col=label,label=label)) +
    geom_text(family="fontawesome-webfont", size=20) +
    scale_colour_manual(values = rev(gg_color_hue(2))) +
    xlab(paste0("Dimension 1 (",
                formatC(100*mds$var.explained[1],di=1,fo='f'),'%)')) + 
    ylab(paste0("Dimension 2 (",
                formatC(100*mds$var.explained[2],di=1,fo='f'),'%)')) +
    theme(legend.title=element_blank(),
          legend.position = "bottom",
          axis.title = element_text(face='bold',size=30),
          legend.text=element_text(family='fontawesome-webfont'),
          panel.grid = element_line(color = "grey"),# ,size = 0.75,linetype = 2)
          panel.background = element_rect(fill = "transparent"),            # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.key.size = unit(x = 0.01,units = "cm"),
          legend.background = element_rect(fill = "transparent"),           # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"),
          axis.text = element_text(size=20)
    ) + guides(color="none")
  ggarrange(gg4_2, nrow=1) 
}
