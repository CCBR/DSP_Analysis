
initialize_object <- function(dcc.files,
                              pkc.files,
                              annotation.file,
                              annotation.sheet.name = "Template",
                              sample.id.field.name = "Sample_ID",
                              roi.field.name = "roi",
                              panel.field.name = "panel",
                              slide.field.name = "slide name", 
                              class.field.name = "class", 
                              region.field.name = "region", 
                              segment.field.name = "segment",
                              area.field.name = "area",
                              nuclei.field.name = "nuclei", 
                              segment.id.length = 4){
    
    # load all input data into a GeoMX object
    object <-
      readNanoStringGeoMxSet(
        dccFiles = dcc.files,
        pkcFiles = pkc.files,
        phenoDataFile = annotation.file,
        phenoDataSheet = annotation.sheet.name,
        phenoDataDccColName = sample.id.field.name, 
        experimentDataColNames = panel.field.name
      )
    
    # Check the column names for required fields exist in the annotation
    
    required.field.names = c(slide.field.name, 
                             class.field.name, 
                             region.field.name, 
                             segment.field.name, 
                             roi.field.name)
    given.field.names = colnames(sData(object))
    
    # Check each of the required fields for correct naming
    for (field in required.field.names) {
      if (!(field %in% given.field.names)) {
        stop(
          paste0(
            field,
            " is not found in the annotation sheet field names.\n"
          )
        )
      }
    }
    
    # Check for the optional fields
    optional.field.names = c("area", "nuclei")
    for (field in optional.field.names) {
      if (!(field %in% given.field.names)) {
        warning(
          paste0(
            field,
            " is not found in the annotation and will not be considered \n"
          )
        )
      }
    }
    
    # Rename all of the required columns based on user parameters in data
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == slide.field.name] = "slide_name"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == class.field.name] = "class"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == region.field.name] = "region"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == segment.field.name] = "segment"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == roi.field.name] = "roi"
    
    # Rename all of the required columns based on user parameters in metadata
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == slide.field.name] = "slide_name"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == class.field.name] = "class"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == region.field.name] = "region"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == segment.field.name] = "segment"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == roi.field.name] = "roi"
    
    # Rename optional columns if they are present
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == area.field.name] = "area"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == nuclei.field.name] = "nuclei"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == area.field.name] = "area"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == nuclei.field.name] = "nuclei" 
    
    # Reformat to remove spaces and dashes in the main annotation columns
    annotation.columns <- c("class", "region", "segment", "slide_name")
    
    for(column in annotation.columns){
      pData(object)[[column]] <- gsub("\\s+", "", pData(object)[[column]])
      pData(object)[[column]] <- gsub("-", "", pData(object)[[column]])
    }
    
    # Establish the segment specific IDs
    pData(object)$segmentID <- paste0(substr(pData(object)$class, 1, segment.id.length),
                                      "|",
                                      substr(pData(object)$region, 1, segment.id.length),
                                      "|",
                                      substr(pData(object)$segment, 1, segment.id.length),
                                      "|", 
                                      substr(pData(object)$slide_name, 1, segment.id.length), 
                                      "|", 
                                      sData(object)$roi)

    
    return(object)
  
}

# Set up the MA plot table
make_MA <- function(contrast.field, 
                    condition.label, 
                    reference.label, 
                    log.counts, 
                    raw.log.counts, 
                    annotation){
  
  # Gather the sample IDs for condition and reference groups
  condition.samples <- rownames(annotation[annotation[[contrast.field]] == condition.label, ])
  reference.samples <- rownames(annotation[annotation[[contrast.field]] == reference.label, ])
  
  # Gather normalized and raw counts for both groups
  condition.counts <- as.data.frame(log.counts[, condition.samples])
  reference.counts <- as.data.frame(log.counts[, reference.samples])
  
  condition.raw.counts <- as.data.frame(raw.log.counts[, condition.samples])
  reference.raw.counts <- as.data.frame(raw.log.counts[, reference.samples])  
  
  # Get the mean log score for each gene for both 
  # normalized counts
  condition.row.order <- rownames(condition.counts)
  condition.counts <- as.data.frame(sapply(condition.counts, as.numeric))
  condition.counts$cond_mean <- rowMeans(condition.counts)
  condition.counts$gene <- condition.row.order
  
  reference.row.order <- rownames(reference.counts)
  reference.counts <- as.data.frame(sapply(reference.counts, as.numeric))
  reference.counts$ref_mean <- rowMeans(reference.counts)
  reference.counts$gene <- reference.row.order
  
  # raw counts
  condition.row.order <- rownames(condition.raw.counts)
  condition.raw.counts <- as.data.frame(sapply(condition.raw.counts, as.numeric))
  condition.raw.counts$cond_raw_mean <- rowMeans(condition.raw.counts)
  condition.raw.counts$gene <- condition.row.order
  
  reference.row.order <- rownames(reference.raw.counts)
  reference.raw.counts <- as.data.frame(sapply(reference.raw.counts, as.numeric))
  reference.raw.counts$ref_raw_mean <- rowMeans(reference.raw.counts)
  reference.raw.counts$gene <- reference.row.order
  
  
  # Create a new data frame of the gene and group means with M and A values
  normalized.counts <- merge(condition.counts, reference.counts, by = "gene") %>% 
    select(gene, cond_mean, ref_mean) %>% 
    mutate(M.value = cond_mean - ref_mean) %>% 
    mutate(A.value = (cond_mean + ref_mean)/2)
  
  raw.counts <- merge(condition.raw.counts, reference.raw.counts, by = "gene") %>% 
    select(gene, cond_raw_mean, ref_raw_mean) %>% 
    mutate(M.raw.value = cond_raw_mean - ref_raw_mean) %>% 
    mutate(A.raw.value = (cond_raw_mean + ref_raw_mean)/2)
  
  # Add the DE results and log counts together
  ma.plot.counts <- merge(normalized.counts, raw.counts, by = "gene")
  
  # Set the bounds for the y axix so that they are aligned
  min.y <- min(c(min(ma.plot.counts$M.value),min(ma.plot.counts$M.raw.value)))
  max.y <- max(c(max(ma.plot.counts$M.value),max(ma.plot.counts$M.raw.value)))
  
  ma.plot.norm <- ggplot(ma.plot.counts, aes(x = A.value, y = M.value)) +
    geom_point(alpha = 0.5, col = "black") + 
    geom_smooth(method=loess, col="steelblue1") + 
    geom_hline(yintercept = 0, lty = "dashed") + 
    labs(x = "Average log expression",
         y = paste0("log(", condition.label, ") - log(", reference.label, ")"), 
         title = "Post-normalization") + 
    ylim(min.y, max.y) + 
    theme_classic()
  
  ma.plot.raw <- ggplot(ma.plot.counts, aes(x = A.raw.value, y = M.raw.value)) + 
    geom_point(alpha = 0.5, col = "black") + 
    geom_smooth(method=loess, col="steelblue1") + 
    geom_hline(yintercept = 0, lty = "dashed") + 
    labs(x = "Average log expression",
         y = paste0("log(", condition.label, ") - log(", reference.label, ")"), 
         title = "Pre-normalization") + 
    ylim(min.y, max.y) + 
    theme_classic()
  
  combined.MA.plots <- arrangeGrob(ggplotGrob(ma.plot.raw), 
                               ggplotGrob(ma.plot.norm), 
                               nrow = 1, ncol = 2)
  
  return(combined.MA.plots)
  
}


plot_distribution <- function(object, annotation.fields){
  
  # run reductions
  color.variable <- Value <- Statistic <- NegProbe <- Q3 <- Annotation <- NULL
  
  # Start Function
  neg.probes<- "NegProbe-WTX"
  
  # Set up a list of annotation fields and values
  annotation.list <- list()
  for(field in annotation.fields){
    annotation.list[[field]] <- unique(pData(object)[[field]])
  }
  
  count.data <- t(exprs(object))
  
  annotation.data <- pData(object)
  
  stat.data <- base::data.frame(row.names = colnames(exprs(object)),
                                AOI = colnames(exprs(object)),
                                Annotation = Biobase::pData(object)[, annotation.fields],
                                Q3 = unlist(apply(exprs(object), 2,
                                                  quantile, 0.75, na.rm = TRUE)),
                                NegProbe = exprs(object)[neg.probes, ])
  
  
  
  stat.data <- stat.data %>% 
    mutate(sig2noise = Q3 / NegProbe)
  
  
  stat.data.melt <- melt(stat.data, measures.vars = c("Q3", "NegProbe"),
                      variable.name = "Statistic", value.name = "Value")
  
  stat.data.melt <- melt(stat.data, 
                         measure.vars = annotation.fields, 
                         variable.name = "field", 
                         value.name = "annotation")
  
  distribution.plot <- ggplot(stat.data.melt, aes(x=Value, 
                                               color=Annotation, 
                                               fill=Annotation)) + 
    geom_density(alpha=0.6) + 
    scale_x_continuous(limits = c(0, max(stat.data.melt$Value) + 10), 
                       expand = expansion(mult = c(0, 0))) + 
    labs(title=" Distribution per AOI of All Probes vs Negative", 
         x="Probe Counts per AOI", 
         y = "Density from AOI Count", 
         color = "Annotation", 
         fill = "Annotation") +
    theme_bw()
  
  #stat.data.mean <- stat.data.m %>% 
  #  mutate(group = paste0(Annotation, Statistic)) %>% 
  #  group_by(group) %>% 
  #  mutate(group_mean = mean(Value)) %>% 
  #  ungroup() %>% 
  #  select(Annotation, Statistic, group_mean) %>% 
  #  distinct()
  
  distribution.plot <- ggplot(stat.data.melt, aes(x=Value, 
                                               color=Statistic, 
                                               fill=Statistic)) + 
    geom_density(alpha=0.6) +
    geom_vline(data=stat.data.mean, aes(xintercept=group_mean, color=Statistic),
               linetype="dashed") +
    scale_color_manual(values = c("#56B4E9", "#E69F00")) +
    scale_fill_manual(values=c("#56B4E9", "#E69F00")) + 
    scale_x_continuous(limits = c(0, max(stat.data.melt$Value) + 10), 
                       expand = expansion(mult = c(0, 0))) +  
    facet_wrap(~Annotation, nrow = 1) + 
    labs(title=" Distribution per AOI of All Probes vs Negative", 
         x="Probe Counts per AOI", 
         y = "Density from AOI Count", 
         color = "Statistic", 
         fill = "Statistic") +
    theme_bw()
  
}  


normalize_counts <- function() {}

top_variable_heatmap <- function(log2.counts, 
                                 top.x.genes = 500, 
                                 annotation.column, 
                                 annotation.row = NULL, 
                                 anno.colors, 
                                 cluster.rows = FALSE, 
                                 cluster.columns = FALSE, 
                                 main.title, 
                                 row.gaps = NULL, 
                                 column.gaps = NULL, 
                                 show.rownames = FALSE, 
                                 show.colnames = FALSE){
  
  # create Coefficient of Variation (CV) function and apply to the log counts
  calc_CV <- function(x) {sd(x) / mean(x)}
  cv.df <- data.frame(CV = apply(log2.counts, 1, calc_CV))
  
  # Take the top X most variable genes by CV score
  cv.df.top <- cv.df %>% arrange(desc(CV)) %>% slice(1:top.x.genes)
  
  # Get the list of top CV genes
  top.cv.gene.list <- rownames(cv.df.top)
  
  # Subset the counts for the top CV genes
  top.cv.heatmap.counts <- log2.counts[rownames(log2.counts) %in% top.cv.gene.list, ]
  
  # Order the counts by top CV
  top.cv.heatmap.counts <- top.cv.heatmap.counts[match(top.cv.gene.list, rownames(top.cv.heatmap.counts)), ]
  
  # Subset the annotation and arrange the order
  annotation.column.fields <- names(anno.colors)
  
  annotation.row.order <- gsub("\\.dcc", "", rownames(annotation.column))
  
  # Order the samples in counts the same as the annotation
  top.cv.heatmap.counts <- top.cv.heatmap.counts[, annotation.row.order]
  
  heatmap.plot <- pheatmap(top.cv.heatmap.counts, 
                           main = main.title, 
                           show_rownames = show.rownames, 
                           scale = "row",   
                           show_colnames = show.colnames,
                           border_color = NA, 
                           cluster_rows = cluster.rows, 
                           cluster_cols = cluster.columns, 
                           clustering_method = "average", 
                           clustering_distance_rows = "correlation", 
                           clustering_distance_cols = "correlation", 
                           color = colorRampPalette(c("blue", "white", "red"))(120), 
                           annotation_row = annotation.row, 
                           annotation_col = annotation.column,  
                           annotation_colors = anno.colors, 
                           gaps_row = row.gaps, 
                           gaps_col = column.gaps, 
                           fontsize_row = 4)
  
  
  return(heatmap.plot)
  
}

plot_umap <- function(log.counts, 
                      annotation, 
                      group.field, 
                      roi.field, 
                      slide.field){
  
  # Set up the counts and order by sample ID
  log.counts.transpose <- as.data.frame(t(log.counts))
  log.counts.transpose <- log.counts.transpose[order(rownames(log.counts.transpose)), ]
  
  # Order the annotation by sample ID
  annotation <- annotation[order(rownames(annotation)), ]
  
  # Run 2D UMAP and select PCs
  umap <- umap(log.counts.transpose, 
               n_components = 2, 
               random_state = 15) 
  layout <- umap[["layout"]] 
  layout <- data.frame(layout) 
  
  # Merge the annotation and UMAP
  layout$sampleID <- rownames(layout)
  annotation$sampleID <- rownames(annotation)
  umap.df <- merge(layout, annotation, by = "sampleID") 
  
  # Use the correct column names in mutate and select
  umap.df <- umap.df %>% 
    mutate(segmentID = paste({{ roi.field }}, {{ slide.field }}, sep = "|")) %>% 
    select(segmentID, X1, X2, {{ group.field }})
  
  # Create the UMAP plot
  umap.plot <- ggplot(umap.df, 
                         aes(x = X1, 
                             y = X2, 
                             color = !!sym(group.field), 
                             fill = !!sym(group.field))) +
    geom_point() + 
    geom_encircle(inherit.aes = TRUE, 
                  alpha = 0.2)
  
  return(umap.plot)
}

gene_detect_plot <- function(object, 
                             facet.column = NULL, 
                             loq.mat = NULL){
  
  # Create the plot for the all genes
  gene.stacked.bar.plot.total <- ggplot(fData(object),
                                        aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = Module)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate (Detected AOIs/Total AOIs)",
         y = "Genes, #",
         fill = "Probe Set")
  
  
  # If a facet has been selected also make a faceted bar plot
  if(!is.null(facet.column)) {
    
    # Gather the facet annotation information
    annotation.data <- pData(object)
    facet.values <- unique(annotation.data[[facet.column]])
    
    # A master df to hold all feature (gene) detection for facet values
    feature.detect.facet.df <- data.frame(feature = rownames(fData(object)))
    
    
    # Gather the IDs for each facet value
    for(value in facet.values){
      
      # Gather the sample IDs for only the current facet value
      value.df <- annotation.data %>% 
        filter(!!sym(facet.column) == value)
      
      value.IDs <- rownames(value.df)
      
      total.AOIs <- length(value.IDs)
      
      # Gather the detection per gene for value Sample IDs
      loq.mat.value <- loq.mat[, value.IDs]
      
      # Compute the detection for each feature
      value.feature.df <- data.frame(feature = rownames(fData(object)))
      
      value.feature.df[[value]] <- 100*(rowSums(loq.mat.value, na.rm = TRUE)/total.AOIs)
      
      # Add the detection per feature for this value to the master df
      feature.detect.facet.df <- merge(feature.detect.facet.df, 
                                       value.feature.df, 
                                       by = "feature")
    }
    
    # Melt the feature detect facet df for easier ggplot faceting
    
    facet.df.melt <- feature.detect.facet.df %>% 
      pivot_longer(cols = -feature, 
                   names_to = "class", 
                   values_to = "detection")
    
    # Create bins for the boxplot
    detection.bins <- c("0", 
                        "<1", 
                        "1-5", 
                        "5-10", 
                        "10-20", 
                        "20-30", 
                        "30-40", 
                        "40-50", 
                        ">50")
    
    # Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
    facet.df.melt$detection_bin <- 
      cut(facet.df.melt$detection,
          breaks = c(-1, 0, 1, 5, 10, 20, 30, 40, 50, 100),
          labels = detection.bins)
    
    facet.table <- table(facet.df.melt$detection_bin,
                         facet.df.melt$class)
    
    max.count.facet <- max(facet.table)
    
    gene.stacked.bar.plot.facet <- ggplot(facet.df.melt,
                                          aes(x = detection_bin, 
                                              fill = class)) +
      geom_bar(position = "dodge") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                         breaks = seq(0, max(max.count.facet), by = 500)) +
      labs(x = "Gene Detection Rate (Detected AOIs/Total AOIs)",
           y = "Number of Genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  }
  
  
  return(list("total.plot" = gene.stacked.bar.plot.total, 
                 "facet.plot" = gene.stacked.bar.plot.facet, 
                 "facet.table" = facet.table))
}

plot_sankey <- function(object, 
                        lane.1, 
                        lane.2, 
                        lane.3, 
                        lane.4, 
                        fill.lane){
  
  #Rename the slide name column for formatting
  pData(object) <- pData(object) %>% 
    mutate(slide = gsub("slide_", "", slide_name))
  
  lanes <- c(lane.1, lane.2, lane.3, lane.4)
  
  
  #Establish variables for the Sankey plot
  x <- id <- y <- n <- NULL
  
  # select the annotations we want to show, use `` to surround column
  # names with spaces or special symbols
  
  # Create a count matrix
  count.mat <- count(pData(object), 
                     !!as.name(lane.1), 
                     !!as.name(lane.2), 
                     !!as.name(lane.3), 
                     !!as.name(lane.4))
  
  # Remove any rows with NA values
  na.per.column <- colSums(is.na(count.mat))
  na.total.count <- sum(na.per.column)
  
  if(na.total.count > 0){
    count.mat <- count.mat[!rowSums(is.na(count.mat)),]
    rownames(count.mat) <- 1:nrow(count.mat)
  }
  
  
  # Gather the data and plot in order: lane 1, lane 2, ..., lane n
  # gather_set_data creates x, id, y, and n fields within sankey.count.data
  # Establish the levels of the Sankey
  sankey.count.data <- gather_set_data(count.mat, 1:4)
  
  # Define the annotations to use for the Sankey x axis labels
  sankey.count.data$x[sankey.count.data$x == 1] <- lane.1
  sankey.count.data$x[sankey.count.data$x == 2] <- lane.2
  sankey.count.data$x[sankey.count.data$x == 3] <- lane.3
  sankey.count.data$x[sankey.count.data$x == 4] <- lane.4
  
  sankey.count.data$x <-
    factor(
      sankey.count.data$x,
      levels = c(as.name(lane.1), 
                 as.name(lane.2), 
                 as.name(lane.3), 
                 as.name(lane.4)))
  
  # For position of Sankey 100 segment scale
  adjust.scale.pos = -1.1
  
  # plot Sankey diagram
  sankey.plot <-
    ggplot(sankey.count.data,
           aes(
             x,
             id = id,
             split = y,
             value = n
           )) +
    geom_parallel_sets(aes(fill = !!as.name(fill.lane)), 
                       alpha = 0.5, 
                       axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.2, 
                            fill = "seashell", 
                            color = "seashell4") +
    geom_parallel_sets_labels(color = "black",
                              size = 3,
                              angle = 0) + 
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      axis.text.y = element_blank()
    ) + 
    scale_y_continuous(expand = expansion(0)) +
    scale_x_discrete(expand = expansion(0)) +
    labs(x = "", y = "") +
    annotate(
      geom = "segment",
      x = (3.25 - adjust.scale.pos),
      xend = (3.25 - adjust.scale.pos),
      y = 20,
      yend = 120,
      lwd = 2
    ) +
    annotate(
      geom = "text",
      x = (3.19 - adjust.scale.pos),
      y = 70,
      angle = 90,
      size = 5,
      hjust = 0.5,
      label = "100 AOIs"
    )
  
  
  # Make the annotation bar plot
  
  AOI.counts <- sankey.count.data
  
  # Gather the counts for each annotation
  AOI.counts$AOI_count <- as.numeric(AOI.counts$n)
  AOI.counts$type <- as.character(AOI.counts$x)
  AOI.counts$annotation <- AOI.counts$y
  
  # Create a sum for each annotation
  AOI.annotation.sum <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(AOI.annotation.sum) <- c("annotation", "AOI_sum")
  
  # Create a data frame of AOI sums per annotation 
  for(anno in unique(AOI.counts$annotation)){
    
    # Filter for a specific annotation
    anno.subset <- AOI.counts %>% 
      filter(annotation == anno)
    
    # Add together the AOI counts
    anno.sum.row <- data.frame(AOI_sum = sum(anno.subset$AOI_count), annotation = anno)
    
    # Append to the master AOI sum df
    AOI.annotation.sum <- rbind(AOI.annotation.sum, anno.sum.row)
    
  }
  
  # Creare a final df for plotting
  AOI.counts.all <- merge(AOI.annotation.sum, AOI.counts, by = "annotation")
  
  AOI.counts.all  <-  AOI.counts.all %>% 
    select(all_of(c("AOI_sum", "type", "annotation"))) %>% 
    distinct()
  
  # Create the bar plots
  AOI.bar.plot <- ggplot(AOI.counts.all, aes(x = annotation, y = AOI_sum)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~ type, ncol = 2, scales = "free_x") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
    geom_text(aes(label = AOI_sum), vjust = -0.3, size = 3.5) +
    labs(x = NULL, y = "AOI Count") + 
    ylim(0, max(AOI.counts.all$AOI_sum) + 30)
  
  return(list("sankey.plot" = sankey.plot, 
              "AOI.bar.plot" = AOI.bar.plot, 
              "sankey.count.data" =  sankey.count.data))
  
}

upsetr_plot <- function(object, 
                        annotation.groups){
  
  # To hold all annotation values for each annotation of interest
  all.group.values <- c()
  
  # Gather all of the values for the upsetr plot
  for(group in annotation.groups){ 
    
    group.values <- unique(pData(object)[[group]])
    
    all.group.values <- c(all.group.values, group.values)
    
  }
  
  # Create the upset df with all FALSE values
  upset.df <- as.data.frame(matrix(FALSE, nrow = nrow(pData(object)), 
                                   ncol = length(all.group.values)))
  
  # Rename the columns to be all possible values for the upsetr plot
  colnames(upset.df) <- all.group.values
  
  # Subset the annotation for only the relevant columns for upsetr
  anno.subset <- pData(object) %>% select(all_of(annotation.groups))
  
  # For each row in the annotation data, if it contains the value of a column in the upsetr plot mark as TRUE
  for (i in 1:nrow(anno.subset)) {
    row.values <- as.character(unlist(anno.subset[i, ]))
    upset.df[i, row.values] <- TRUE
  }
  
  # Create the UpSetR Plot
  AOI.inter.count.plot <- upset(upset.df,  
                                intersect = all.group.values, 
                                width_ratio = 0.4, 
                                min_size = 4, 
                                set_sizes=(upset_set_size() + 
                                             geom_text(aes(label=..count..),
                                                       hjust=1.1, stat='count') +
                                             expand_limits(y=nrow(upset.df)) +
                                             theme(axis.text.x=element_text(angle=90))))
  
  
  return(AOI.inter.count.plot)
  
}

loq_detection <- function(object, 
                          pkc.file.names){
  
  # Set up lists of segment IDs
  segment.list.total <- pData(object)$segmentID
  
  # Define Modules
  modules <- gsub(".pkc", "", pkc.file.names)
  
  # Calculate limit of quantification (LOQ) in each segment
  # LOQ = geomean(NegProbes) * geoSD(NegProbes)^(LOQ cutoff)
  # LOQ is calculated for each module (pkc file)
  loq <- data.frame(row.names = colnames(object))
  
  loq.min <- 2
  loq.cutoff <- 2
  
  for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                   module)
    if(all(vars[1:2] %in% colnames(pData(object)))) {
      
      neg.geo.mean <- vars[1]
      neg.geo.sd <- vars[2]
      
      loq[, module] <-
        pmax(loq.min,
             pData(object)[, neg.geo.mean] * 
               pData(object)[, neg.geo.sd] ^ loq.cutoff)
    }
  }
  
  # Store the loq df in the annotation df
  pData(object)$loq <- loq
  
  # Setup a master loq matrix
  loq.mat <- c()
  
  
  for(module in modules) {
    # Gather rows with the given module
    ind <- fData(object)$Module == module
    
    # Check if each feature has counts above the LOQ
    mat.i <- t(esApply(object[ind, ], MARGIN = 1,
                       FUN = function(x) {
                         x > loq[, module]
                       }))
    
    # Store results in the master loq matrix
    loq.mat <- rbind(loq.mat, mat.i)
  }
  
  # ensure ordering since this is stored outside of the geomxSet
  loq.mat <- loq.mat[fData(object)$TargetName, ]
  
  # Evaluate and Filter Segment Gene Detection Rate
  # Save detection rate information to pheno data
  pData(object)$GenesDetected <- colSums(loq.mat, na.rm = TRUE)
  pData(object)$GeneDetectionRate <- 100*(pData(object)$GenesDetected / nrow(object))
  
  # Establish detection bins
  detection.bins <- c("<1", "1-5", "5-10", "10-15", ">15")
  
  # Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
  pData(object)$DetectionThreshold <- 
    cut(pData(object)$GeneDetectionRate,
        breaks = c(0, 1, 5, 10, 15, 100),
        labels = detection.bins)
  
  return(object)
  
}

aoi_flag_table <- function(aoi.flags){
  
  flag.column.detect <- sapply(aoi.flags, is.logical)
  flag.column.names <- names(aoi.flags[flag.column.detect])
  
  # A function for coloring TRUE flags as red
  red.flag <- function(x) {
    x <- as.logical(x)
    ifelse(x, "red", "white") 
  }
  
  # Create the table using the flag coloring function
  aoi.flag.table <- qc.output$segment.flags %>% 
    gt() %>% 
    data_color(columns = flag.column.names, 
               fn = red.flag, 
               alpha = 0.7)
  
  return(aoi.flag.table)
  
}

probe_flag_table <- function(probe.flags, 
                             object){
  
  # Create the table for probe flags
  probe.flags.df <- probe.flags %>% separate_rows(LocalFlag, sep = ",")
  
  # Rename the dcc file name column
  probe.flags.df$Sample_ID <- probe.flags.df$LocalFlag
  
  # Grab the annotation for only the columns to map
  annotation <- pData(object)
  annotation$Sample_ID <- rownames(annotation)
  
  annotation.subset <- annotation %>% 
    select(Sample_ID, segmentID)
  
  # Map the AOI names in the flags to the segmentID
  probe.flags.df <- merge(probe.flags.df, annotation.subset, by = "Sample_ID")
  
  # Remove the dcc file name column 
  probe.flags.table <- probe.flags.df %>% 
    select(TargetName, RTS_ID, segmentID, FlagType) %>% 
    gt()
  
  # For a summary of only probe names
  probe.flag.summary <- qc.output$probe.flags %>% 
    select(TargetName, RTS_ID, FlagType) %>% 
    gt()
  
  return(list("probe.flag.table" = probe.flags.table, 
              "probe.flag.summary" = probe.flag.summary))
  
}
