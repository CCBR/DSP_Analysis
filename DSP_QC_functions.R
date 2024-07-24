
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
        protocolDataColNames = roi.field.name,
        experimentDataColNames = panel.field.name
      )

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