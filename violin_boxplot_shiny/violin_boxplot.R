library(dplyr)
library(ggplot2)
library(tibble)
library(ggpubr)
library(tidyr)

violin_boxplot <- function(counts, 
                           annotation.df, 
                           gene.list, 
                           annotation.field, 
                           display.summary.stat = FALSE, 
                           compare = FALSE, 
                           compare.groups = NULL){
  
  # Set up the annotation df
  annotation.fields <- c("SampleID", annotation.field)
  annotation.df$SampleID <- rownames(annotation.df)
  
  subset.annotation <- as.data.frame(annotation.df)[, annotation.fields, 
                                                    drop = FALSE]
  
  # Convert counts to data frame 
  counts <- as.data.frame(counts)
  
  # Check if goi are found in counts
  for(gene in gene.list){
    
    if(!(gene %in% rownames(counts))){
      
      print(paste0(gene, " not found in counts file"))
      
      gene.list <- gene.list[-which(gene.list == gene)]
      
    }
  }
  
  # Convert gene counts to log2 for only genes of interest
  gene.counts <- counts %>% 
    filter(rownames(counts) %in% gene.list) %>% 
    mutate(across(where(is.numeric), ~.+1)) %>% 
    mutate(across(where(is.numeric), log2))
  
  # Set up counts for merge with annotation
  gene.counts.transpose <- as.data.frame(t(gene.counts)) %>% 
    rownames_to_column(var = "SampleID")
  
  # Create master annotation/counts df
  counts.anno.df <- merge(gene.counts.transpose, 
                          subset.annotation, 
                          by = "SampleID")
  
  # Set up the annotation/counts df for ggplot2
  counts.anno.df.melt <- counts.anno.df %>% 
    pivot_longer(cols = all_of(gene.list), 
                 names_to = "gene", 
                 values_to = "log_counts")
  
  counts.anno.df.melt$log_counts <- as.numeric(counts.anno.df.melt$log_counts)
  
  max.value <- max(counts.anno.df.melt$log_counts)
  
  # Create a combined boxplot and violin plot
  if(display.summary.stat == TRUE){
    
    field.violin.plot <- ggplot(counts.anno.df.melt, aes(x = !!sym(annotation.field), 
                                                         y = log_counts, 
                                                         fill = !!sym(annotation.field))) +
      geom_violin() + 
      geom_boxplot(width = 0.2, fill = "white") + 
      facet_wrap(~ gene) + 
      labs(x = paste(gsub("_", " ", annotation.field)), 
           y = "Log2 Counts", 
           title = paste0("Log2 Counts for ", gsub("_", " ", annotation.field))) + 
      stat_summary(
        fun = mean, 
        geom = "text", 
        aes(label = paste("mean:", round(after_stat(y), 2))), 
        vjust = -0.5, 
        color = "darkblue"
      ) + 
      stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "darkblue", 
                   fatten = 0.5) +
      stat_summary(
        fun.min = function(x) quantile(x, 0.25),
        fun = function(x) quantile(x, 0.25), 
        geom = "text", 
        aes(label = paste("Q1:", round(after_stat(y), 2))),
        vjust = 1.5, 
        color = "black"
      ) +
      stat_summary(
        fun.max = function(x) quantile(x, 0.75),
        fun = function(x) quantile(x, 0.75), 
        geom = "text", 
        aes(label = paste("Q3:", round(after_stat(y), 2))),
        vjust = -1.5, 
        color = "black"
      ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else if(compare == TRUE) {
    
    field.violin.plot <- ggplot(counts.anno.df.melt, aes(x = !!sym(annotation.field), 
                                                         y = log_counts, 
                                                         fill = !!sym(annotation.field))) +
      geom_violin() + 
      geom_boxplot(width = 0.2, fill = "white") + 
      facet_wrap(~ gene) + 
      labs(x = paste(gsub("_", " ", annotation.field)), 
           y = "Log2 Counts", 
           title = paste0("Log2 Counts for ", gsub("_", " ", annotation.field))) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.title.x = element_blank())  + 
      stat_compare_means(comparisons = list(compare.groups), 
                         label = "p.signif", 
                         label.y = max.value*1.01)
    
    #field.violin.plot <- ggplot(counts.anno.df.melt, aes(x = !!sym(annotation.field), 
    #                                                     y = log_counts, 
    #                                                     fill = !!sym(annotation.field))) +
    #  geom_violin() + 
    #  geom_boxplot(width = 0.2, fill = "white") + 
    #  stat_compare_means(method = "wilcox.test", 
    #                     label.y = max.value*1.01) + 
    #  facet_wrap(~ gene) + 
    #  labs(x = paste(gsub("_", " ", annotation.field)), 
    #       y = "Log2 Counts", 
    #       title = paste0("Log2 Counts for ", gsub("_", " ", annotation.field))) + 
    #  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
    #  stat_compare_means(label = "p.signif", 
    #                     label.y = max.value + (max.value*0.01))
    
  } else(
    
    field.violin.plot <- ggplot(counts.anno.df.melt, aes(x = !!sym(annotation.field), 
                                                         y = log_counts, 
                                                         fill = !!sym(annotation.field))) +
      geom_violin() + 
      geom_boxplot(width = 0.2, fill = "white") + 
      facet_wrap(~ gene) + 
      labs(x = paste(gsub("_", " ", annotation.field)), 
           y = "Log2 Counts", 
           title = paste0("Log2 Counts for ", gsub("_", " ", annotation.field))) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  )
  
  
  
  
  return(field.violin.plot)
  
}