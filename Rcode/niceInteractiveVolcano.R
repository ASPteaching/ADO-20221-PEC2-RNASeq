#' An interactive Volcano Plot
#' Preparat per Ferran Reverter
#' @param fittedObj A limma object produced after fitting and shrinking a linear model
#' @return An interactive Volcano Plot
#' @examples
#' # niceInteractiveVolcano (fit.cont)
#' @export
niceInteractiveVolcano <- function(fittedObj){
  require(gridExtra)
  require(plotly)
  diseased_vs_healthy<-data.frame(fittedObj$genes, fittedObj$coefficients,fittedObj$p.value)
  colnames(diseased_vs_healthy)<-c("geneid","foldchange","adjpval")
  vol_plot <- diseased_vs_healthy %>%
    ggplot(aes(x = foldchange,
               y = -log10(adjpval))) + 
    geom_point() 
  
  vol_plot + 
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed") + xlim(-6.5, 6.5) 
  
  diseased_vs_healthy <- diseased_vs_healthy %>%
    mutate(gene_type = case_when(foldchange >= 1.5 & adjpval <= 0.05 ~ "up",
                                 foldchange <= -1.5 & adjpval <= 0.05 ~ "down",
                                 TRUE ~ "ns"))   
  
  # Obtain gene_type counts ------------------------------------------------------           
  diseased_vs_healthy %>%
    count(gene_type)
  
  diseased_vs_healthy %>%
    distinct(gene_type) %>%
    pull()  
  #> [1] "down" "up"   "ns"    
  
  # Add colour, size and alpha (transparency) to volcano plot --------------------
  cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("up" = 3, "down" = 3, "ns" = 1) 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
  
  diseased_vs_healthy %>%
    ggplot(aes(x = foldchange,
               y = -log10(adjpval),
               fill = gene_type,    
               size = gene_type,
               alpha = gene_type)) + 
    geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
               colour = "black") + 
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(-2, 2),
               linetype = "dashed") +
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas) + # Modify point transparency
    scale_x_continuous(breaks = c(seq(-6.5, 6.5, 1)),       
                       limits = c(-7, 7))  
  
  # add a grouping column; default value is "not significant"
  
  diff_df<-diseased_vs_healthy
  colnames(diff_df)<-c("geneid","foldchange","adjpval","genetype")
  diff_df["group"] <- "NotSignificant"
  
  # for our plot, we want to highlight 
  # FDR < 0.05 (significance level)
  # Fold Change > 1.5
  
  # change the grouping for the entries with significance but not a large enough Fold change
  diff_df[which(diff_df['adjpval'] < 0.05 & abs(diff_df['foldchange']) < 1.5 ),"group"] <- "Significant"
  
  # change the grouping for the entries a large enough Fold change but not a low enough p value
  diff_df[which(diff_df['adjpval'] > 0.05 & abs(diff_df['foldchange']) > 1.5 ),"group"] <- "FoldChange"
  
  # change the grouping for the entries with both significance and large enough fold change
  diff_df[which(diff_df['adjpval'] < 0.05 & abs(diff_df['foldchange']) > 1.5 ),"group"] <- "Significant&FoldChange"
  
  # Find and label the top peaks..
  top_peaks <- diff_df[with(diff_df, order(foldchange, adjpval)),][1:5,]
  top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-foldchange, adjpval)),][1:5,])
  
  # Add gene labels for all of the top genes we found
  # here we are creating an empty list, and filling it with entries for each row in the dataframe
  # each list entry is another list with named items that will be used by Plot.ly
  a <- list()
  for (i in seq_len(nrow(top_peaks))) {
    m <- top_peaks[i, ]
    a[[i]] <- list(
      x = m[["foldchange"]],
      y = -log10(m[["adjpval"]]),
      text = m[["geneid"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40
    )
  }
  
  # make the Plot.ly plot
  p <- plot_ly(data = diff_df, x=~foldchange, y=~-log10(adjpval), text = ~geneid, mode = "markers", color = ~group) %>% 
    layout(title ="") %>% layout(annotations = a)
  
  p
  
}
