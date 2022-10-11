#' Creacio de Boxplots amb ggplot2 per matrius de comptatges
#' Preparat per Jordi Cortes
#' @param aCPMObj Data matrix with Counts Per Million 
#' @return A ggplot2-based Boxplot.
#' @examples
#' # Not run
#' # niceHeatMap(mat)
#' @export
# niceBoxPlot <- function (aCPMObj){
#   require(tidyverse)
#   d <- tibble(y=as.numeric(aCPMObj),
#               x=rep(colnames(aCPMObj),each=nrow(aCPMObj)),
#               # x=rep(sub(pattern = 'HEA','IND_',
#               #           sub(pattern = 'COV',replacement = 'IND_',colnames(aCPMObj))),
#               #       each=nrow(aCPMObj)),
#               co=factor(rep(ifelse(
#                 substr(colnames(aCPMObj),1,3)=='COV',
#                 'COVID','HEALTHY'),
#                 each=nrow(aCPMObj))))
#   gg1 <- ggplot(d,aes(x,y,fill=co)) + geom_boxplot() +
#     xlab('') + ylab(expression(bold('Normalized'~'counts'~(log[2])))) +
#     theme(axis.text.x = element_text(angle=90,face='bold',vjust = 0.5),
#           axis.title = element_text(face='bold'),
#           legend.title = element_blank(),
#           legend.position = 'bottom')
#   gg1
#   #ggsave(filename = 'figures/boxplot.png',gg1)
# }

niceBoxPlot <- function (aCPMObj){
  require(tidyverse)
  d <- tibble(y=as.numeric(aCPMObj),
              x=rep(colnames(aCPMObj),each=nrow(aCPMObj)),
              co=factor(rep(ifelse(
                substr(colnames(aCPMObj),1,3)=='COV',
                'COVID','HEALTHY'),
                each=nrow(aCPMObj))))
  gg1 <- ggplot(d,aes(x,y,fill=co)) + geom_boxplot() +
    xlab('') + ylab(expression(bold('Normalized'~'counts'~(log[2])))) +
    theme(axis.text.x = element_text(angle=90,face='bold',vjust = 0.5, size=13),
          axis.text.y = element_text(size=13),
          axis.title.y = element_text(size=15),
          axis.title = element_text(face='bold'),
          legend.title = element_blank(),
          legend.text = element_text(size=13),
          legend.position = 'bottom',
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"))
  gg1
}

niceBoxPlotly <- function (aCPMObj){
  require(plotly)
  cols = gg_color_hue(2)
  d <- tibble(y=as.numeric(aCPMObj),
              x=rep(colnames(aCPMObj),each=nrow(aCPMObj)),
              co=factor(rep(ifelse(
                substr(colnames(aCPMObj),1,3)=='COV',
                'COVID','HEALTHY'),
                each=nrow(aCPMObj))))
  plot_ly(d, x = ~x, y = ~y, type = "box", #color=~cols,
          colors = gg_color_hue(2)) %>%
    layout(xaxis = list(title = ''),
           yaxis = list(title = 'Normalized counts (log2)',
                        hoverformat = '.2f'),
           legend = list(x = 0.5, y = 17, orientation = 'h',xanchor = "center"))
}