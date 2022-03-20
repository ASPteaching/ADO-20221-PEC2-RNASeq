#' Creacio de Boxplots amb ggplot2 per matrius de comptatges
#' Preparat per Jordi Cortes
#' @param aCPMObj Data matrix with Counts Per Million 
#' @return A ggplot2-based Boxplot.
#' @examples
#' # Not run
#' # niceHeatMap(mat)
#' @export
niceBoxPlot <- function (aCPMObj){
  require(tidyverse)
  d <- tibble(y=as.numeric(aCPMObj),
              x=rep(colnames(aCPMObj),each=nrow(aCPMObj)),
              # x=rep(sub(pattern = 'HEA','IND_',
              #           sub(pattern = 'COV',replacement = 'IND_',colnames(log2count_norm))),
              #       each=nrow(log2count_norm)),
              co=factor(rep(ifelse(
                substr(colnames(aCPMObj),1,3)=='COV',
                'COVID','HEALTHY'),
                each=nrow(aCPMObj))))
  gg1 <- ggplot(d,aes(x,y,fill=co)) + geom_boxplot() +
    xlab('') + ylab(expression(bold('Normalized'~'counts'~(log[2])))) +
    theme(axis.text.x = element_text(angle=90,face='bold',vjust = 0.5),
          axis.title = element_text(face='bold'),
          legend.title = element_blank(),
          legend.position = 'bottom')
  gg1
  #ggsave(filename = 'figures/boxplot.png',gg1)
}


