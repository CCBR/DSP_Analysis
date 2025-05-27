library(dplyr)
library(PCAtools)
library(RColorBrewer)

PCA_plot <- function(counts, 
                     annotation.df, 
                     annotation.field){
  
  counts.df <- as.data.frame(counts) %>% 
    mutate_all(~ log2(.)) %>% 
    rename_all(~ gsub("\\.dcc", "", .))
  
  
  # Remove the negative controls from the log counts
  control.probes <- c("NegProbe-WTX")
  counts.df <- counts.df[!(rownames(counts.df) %in% control.probes), ]
  
  # Order of rownames of annotation need to match columns of count data
  ordered.annotation.df <- annotation.df[order(rownames(annotation.df)), ]
  
  ordered.counts.df <- counts.df[order(colnames(counts.df))]
  
  # Remove .dcc from Sample ID row names
  ordered.annotation.df <- ordered.annotation.df %>% 
    `rownames<-`(sub("\\.dcc", "", rownames(.)))
  
  pca.table <- pca(ordered.counts.df, 
                   metadata = ordered.annotation.df, 
                   removeVar = 0.1)
  
  # Get the annotation values
  annotation.values <- unique(annotation.df[[annotation.field]])
  
  # Get as many colors as there are annotation values
  n.values <- length(annotation.values)
  paired.colors <- brewer.pal(n = max(n.values, 3), 
                              name = "Paired")[1:n.values]
  
  annotation.colors <- setNames(paired.colors, 
                                annotation.values)
  
  
  pca.plot <- biplot(pca.table, 
                     colby = annotation.field, 
                     colkey = annotation.colors, 
                     legendPosition = "right", 
                     legendLabSize = 6, 
                     legendIconSize = 3, 
                     lab = NULL,
                     title = "Quantile Normalization", 
                     subtitle = "NTCs removed")
    
  
  
  return(pca.plot)
  
}