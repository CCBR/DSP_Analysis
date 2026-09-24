# Example YAML:
#image_folder: /data/CCBR/spitr/spitr10_tosato/slide_images
#annotation_file: /data/CCBR/spitr/spitr10_tosato/annotations/annotation_tosato_CPTR10.xlsx
#tiff_file_list:
#  12:
#    tiff.file: /data/CCBR/spitr/spitr10_tosato/slide_images/12/1_2_WTA_050124.ome.tiff
#    slide.name: 1-2 WTA 050124
#  1A1:
#    tiff.file: /data/CCBR/spitr/spitr10_tosato/slide_images/1A1/1A_1_WTA_042324.ome.tiff
#    slide.name: 1A-1 WTA 042324
#  1B1:
#    tiff.file: /data/CCBR/spitr/spitr10_tosato/slide_images/1B1/1B_1_WTA_042624.ome.tiff
#    slide.name: 1B-1 WTA 042624


## To run:
## Rscript write_xml_for_halo.R project_params.yaml

## Setup ##
library(SpatialOmicsOverlay)
library(ggplot2)
library(xml2)
library(EBImage)
library(dplyr)
library(readxl)
library(yaml)

# Read the project specific parameters
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) stop("Usage: Rscript write_xml_for_halo.R project_params.yaml")
config.path <- args[1]
config <- yaml::read_yaml(config.path)

project.folder <- config[["image_folder"]]
annotation.file <- config[["annotation_file"]]
tiff.file.list <- config[["tiff_file_list"]]

if (is.null(project.folder) || is.null(annotation.file) || is.null(tiff.file.list)) {
    stop("Project config is missing one or more required fields: image_folder, annotation_file, or tiff_file_list")
}

## Functions ##
annotMatchingFixed <- function(annots, ROInum, maskText_segment, segCol = NULL) {
  if (!"ROILabel" %in% colnames(annots)) {
    stop("The column ROILabel is not in annots.")
  }
  if (suppressWarnings(!is.na(as.numeric(ROInum)))) {
    ROInum <- as.numeric(ROInum)
    annots$ROILabel <- as.numeric(annots$ROILabel)
  }
  w2kp <- which(annots$ROILabel == ROInum)
  if (length(w2kp) == 0) return(NULL)
  annots <- annots[w2kp, ]
  
  if (!"SegmentID" %in% colnames(annots)) {
    stop("SegmentID column required for matching")
  }
  
  # Match directly on the mask's own Text attribute, not position
  match_row <- which(annots$SegmentID == maskText_segment)
  if (length(match_row) == 0) return(NULL)
  annots <- annots[match_row, , drop = FALSE]
  
  if ("SegmentDisplayName" %in% colnames(annots)) sampCol <- "SegmentDisplayName"
  if ("Sample_ID" %in% colnames(annots)) sampCol <- "Sample_ID"
  
  annots <- annots[, c(sampCol, "SegmentID")]
  colnames(annots) <- c("Sample_ID", "SegmentID")
  return(annots)
}

parseOverlayAttrsFixed <- function(omexml, annots, ...) {
  ROIs <- omexml[which(names(omexml) == "ROI")]
  names(ROIs) <- paste0(names(ROIs), 1:(length(ROIs)))
  AOIattrs <- NULL
  
  for (ROI_name in names(ROIs)) {
    ROInum <- ROIs[[ROI_name]]$Union$Label[["Text"]]
    ROInum <- gsub("\\W", "", ROInum)
    
    ROI <- ROIs[[ROI_name]]$Union
    masks <- which(names(ROI) == "Mask")
    
    for (mask in masks) {
      mask.attrs <- ROI[[mask]]$.attrs
      segmentation <- ifelse(length(masks) == 1, "Geometric", "Segmented")
      
      # Use the mask's own Text field directly, if present; else NA (geometric ROI, no segment text)
      maskText_segment <- if ("Text" %in% names(mask.attrs)) mask.attrs[["Text"]] else NA
      
      if (is.na(maskText_segment)) {
        # Geometric case — only one mask, no segment disambiguation needed
        ROIannot <- annots[annots$ROILabel == ROInum, ]
        if (nrow(ROIannot) == 0) next
        sampCol <- if ("Sample_ID" %in% colnames(ROIannot)) "Sample_ID" else "SegmentDisplayName"
        ROIannot <- ROIannot[, c(sampCol, "SegmentID")]
        colnames(ROIannot) <- c("Sample_ID", "SegmentID")
      } else {
        ROIannot <- annotMatchingFixed(annots, ROInum, maskText_segment)
      }
      
      if (is.null(ROIannot) || nrow(ROIannot) == 0) next
      
      AOIattr <- as.data.frame(c(ROILabel = ROInum,
                                 ROIannot,
                                 mask.attrs[c("Height", "Width", "X", "Y")],
                                 Segmentation = segmentation))
      AOIattr$Height <- as.numeric(AOIattr$Height)
      AOIattr$Width  <- as.numeric(AOIattr$Width)
      AOIattr$X <- as.numeric(AOIattr$X)
      AOIattr$Y <- as.numeric(AOIattr$Y)
      
      AOIattrs <- rbind(AOIattrs, cbind(AOIattr, Position = ROI[[mask]]$BinData$text))
    }
  }
  
  return(SpatialPosition(position = AOIattrs))
}

# Identify holes in the AOI for drawing inner boundaries
find_holes <- function(mask) {
  inverted <- 1 - mask
  inv_labeled <- bwlabel(inverted)
  n_inv <- max(inv_labeled)
  
  hole_labels <- c()
  for (i in seq_len(n_inv)) {
    coords_i <- which(inv_labeled == i, arr.ind = TRUE)
    # If this "background" component touches any edge of the mask's bounding box,
    # it's the true exterior background, not an enclosed hole — skip it
    touches_border <- any(coords_i[,1] == 1 | coords_i[,1] == nrow(mask) |
                           coords_i[,2] == 1 | coords_i[,2] == ncol(mask))
    if (!touches_border) {
      hole_labels <- c(hole_labels, i)
    }
  }
  
  if (length(hole_labels) == 0) return(list())
  
  hole_mask <- inv_labeled
  hole_mask[!(hole_mask %in% hole_labels)] <- 0
  hole_labeled <- bwlabel(hole_mask > 0)
  
  ocontour(hole_labeled)
}

trace_all_aois <- function(coords_split) {
  results <- list()
  aoi_ids <- names(coords_split)
  n <- length(aoi_ids)
  
  for (i in seq_along(aoi_ids)) {
    id <- aoi_ids[i]
    aoi_coords <- coords_split[[id]]
    
    x_min <- min(aoi_coords$xcoor); y_min <- min(aoi_coords$ycoor)
    x_max <- max(aoi_coords$xcoor); y_max <- max(aoi_coords$ycoor)
    
    width  <- x_max - x_min + 1
    height <- y_max - y_min + 1
    
    mask <- matrix(0L, nrow = width, ncol = height)
    mask[cbind(aoi_coords$xcoor - x_min + 1, aoi_coords$ycoor - y_min + 1)] <- 1L
    
    labeled <- bwlabel(mask)
    outer_contours_local <- ocontour(labeled)
    hole_contours_local  <- find_holes(mask)
    
    # Shift both outer and hole contours back to real slide-wide pixel coordinates
    shift_contours <- function(contour_list) {
      lapply(contour_list, function(cont) {
        cont[,1] <- cont[,1] + x_min - 1
        cont[,2] <- cont[,2] + y_min - 1
        cont
      })
    }
    
    results[[id]] <- list(
      outer = shift_contours(outer_contours_local),
      holes = shift_contours(hole_contours_local)
    )
    
    if (i %% 10 == 0 || i == n) message(sprintf("Processed %d / %d AOIs", i, n))
    
    rm(mask, labeled, outer_contours_local, hole_contours_local)
  }
  
  results
}

### Testing ###
plot_slide_qc_by_segment <- function(contour_results, segment_lookup, image_dest, max_dim = 2000) {
  
  all_x <- unlist(lapply(contour_results, function(cs) {
    c(unlist(lapply(cs$outer, function(c) c[,1])),
      unlist(lapply(cs$holes, function(c) c[,1])))
  }))
  all_y <- unlist(lapply(contour_results, function(cs) {
    c(unlist(lapply(cs$outer, function(c) c[,2])),
      unlist(lapply(cs$holes, function(c) c[,2])))
  }))
  
  x_range <- range(all_x); y_range <- range(all_y)
  scale <- max_dim / max(diff(x_range), diff(y_range))
  
  seg_colors <- c("Segment 1" = "red", "Segment 2" = "blue")
  
  png(image_dest, width = round(diff(x_range) * scale), height = round(diff(y_range) * scale))
  par(mar = c(0, 0, 0, 0))
  plot(NA, xlim = x_range, ylim = rev(y_range), asp = 1,
       xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  
  total_holes <- 0
  
  for (id in names(contour_results)) {
    seg <- segment_lookup[[id]]
    col <- seg_colors[[seg]]
    
    for (cont in contour_results[[id]]$outer) {
      polygon(cont[,1], cont[,2], border = col)
    }
    
    # Holes now match their AOI's segment color, but stay dashed so they're still distinguishable
    for (cont in contour_results[[id]]$holes) {
      polygon(cont[,1], cont[,2], border = col, lty = 1)
      total_holes <- total_holes + 1
    }
  }
  
  legend("topright", legend = names(seg_colors), col = seg_colors, lty = 1, bg = "white")
  
  dev.off()
  
  segments <- names(seg_colors)
  message(sprintf("Plotted %d AOIs, %d segment types, %d total holes detected", 
                 length(contour_results), length(segments), total_holes))
}

# Write XML output file
write_halo_xml_by_segment <- function(contour_results, segment_lookup, aoi_name_lookup, outfile) {
  doc <- xml_new_root("Annotations")
  
  segments <- sort(unique(segment_lookup[names(contour_results)]))
  aoi_ids <- names(contour_results)
  
  seg_linecolors <- c("Segment 1" = "16711680", "Segment 2" = "65280")  # adjust to match your explicit mapping
  
  # Helper: write outer + hole regions for a list of contour sets into a given Regions node
  write_regions <- function(regions_node, aoi_contours) {
    for (cont in aoi_contours$outer) {
      region_node <- xml_add_child(regions_node, "Region", Type = "Polygon", 
                                     HasEndcaps = "0", NegativeROA = "0")
      verts_node <- xml_add_child(region_node, "Vertices")
      for (j in seq_len(nrow(cont))) {
        xml_add_child(verts_node, "V", X = as.character(round(cont[j,1])), 
                                        Y = as.character(round(cont[j,2])))
      }
    }
    for (cont in aoi_contours$holes) {
      region_node <- xml_add_child(regions_node, "Region", Type = "Polygon", 
                                     HasEndcaps = "0", NegativeROA = "1")
      verts_node <- xml_add_child(region_node, "Vertices")
      for (j in seq_len(nrow(cont))) {
        xml_add_child(verts_node, "V", X = as.character(round(cont[j,1])), 
                                        Y = as.character(round(cont[j,2])))
      }
    }
  }
  
  # --- Segment layers ---
  for (seg in segments) {
    seg_aoi_ids <- names(contour_results)[segment_lookup[names(contour_results)] == seg]
    
    annot_node <- xml_add_child(doc, "Annotation", Name = seg, 
                                  LineColor = as.character(seg_linecolors[[seg]]), 
                                  Visible = "1")
    regions_node <- xml_add_child(annot_node, "Regions")
    
    for (id in seg_aoi_ids) {
      write_regions(regions_node, contour_results[[id]])
    }
  }
  
  # --- Individual AOI layers ---
  for (aoi in aoi_ids) {
    aoi_alt_name <- aoi_name_lookup[[aoi]]
    
    annot_node <- xml_add_child(doc, "Annotation", Name = aoi_alt_name, 
                                  LineColor = "16777215",  # White
                                  Visible = "1")
    regions_node <- xml_add_child(annot_node, "Regions")
    
    write_regions(regions_node, contour_results[[aoi]])
  }
  
  write_xml(doc, outfile)
  message(sprintf("Wrote HALO XML with %d segment layers + %d individual AOI layers", 
                   length(segments), length(aoi_ids)))
}

## Main ##

# --- Step 1: Read and filter annotations (same as before) ---
annotation.df <- read_excel(annotation.file, sheet = "SegmentProperties")

for (tiff_id in names(tiff.file.list)) {
  tiff_info <- tiff.file.list[[tiff_id]]
  tiff.file <- tiff_info$tiff.file
  slide.name <- tiff_info$slide.name
  annotation.df.slide <- annotation.df |> filter(`slide name` == slide.name)

  # --- Step 2: Extract XML (same as readSpatialOverlay does internally) ---
  xml <- xmlExtraction(ometiff = tiff.file, saveFile = FALSE)

  # --- Step 3: Parse scan metadata (same as readSpatialOverlay does internally) ---
  scan_metadata <- SpatialOmicsOverlay:::parseScanMetadata(omexml = xml)

  # --- Step 4: Parse overlay data using OUR FIXED function instead of the buggy one ---
  AOIattrs <- parseOverlayAttrsFixed(omexml = xml, annots = annotation.df.slide)

  # --- Step 5: Determine overall segmentation status (same logic as readSpatialOverlay) ---
  if (any(meta(AOIattrs)$Segmentation == "Segmented")) {
    scan_metadata[["Segmentation"]] <- "Segmented"
  } else {
    scan_metadata[["Segmentation"]] <- "Geometric"
  }

  # --- Step 6: Build the SpatialOverlay object (same constructor readSpatialOverlay uses) ---
  slide.image.overlay <- SpatialOmicsOverlay:::SpatialOverlay(
    slideName = slide.name,
    scanMetadata = scan_metadata,
    overlayData = AOIattrs,
    workflow = list(labWorksheet = FALSE, outline = TRUE, scaled = FALSE),
    image = list(filePath = NULL, imagePointer = NULL, resolution = NULL)
  )

  # --- Step 7: Generate coordinates, same as readSpatialOverlay does when image = FALSE ---
  slide.image.overlay <- createCoordFile(overlay = slide.image.overlay, outline = TRUE)

  roi_meta <- meta(overlay(slide.image.overlay))
  roi_coords <- coords(slide.image.overlay)

  output_dir <- file.path(project.folder, tiff_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  contour_rds <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(tiff.file)), "_contour_results.rds")
  )

  if (file.exists(contour_rds)) {
    contour_results <- readRDS(contour_rds)
    message(sprintf("Loaded existing contour results for %s (%s) -> %s", basename(tiff.file), tiff_id, contour_rds))
  } else {
    coords_split <- split(roi_coords, roi_coords$sampleID)
    #rm(roi_coords); gc()

    contour_results <- trace_all_aois(coords_split)

    saveRDS(contour_results, contour_rds)
    message(sprintf("Saved contour results for %s (%s) -> %s", basename(tiff.file), tiff_id, contour_rds))
  }

  roi_meta <- roi_meta %>%
    mutate(ROI_Segment = paste0("ROI", ROILabel, "_", SegmentID))

  segment_lookup <- setNames(roi_meta$SegmentID, roi_meta$Sample_ID)
  aoi_name_lookup <- setNames(roi_meta$ROI_Segment, roi_meta$Sample_ID)

  if (!all(names(contour_results) %in% names(segment_lookup))) {
    warning(sprintf("Some traced AOIs do not have a matching Sample_ID in %s; check ROI metadata.", tiff_id))
  }
  print(table(segment_lookup[names(contour_results)]))

  qc_png <- file.path(output_dir, paste0("full_slide_", tiff_id, "_qc.png"))
  plot_slide_qc_by_segment(contour_results = contour_results,
                           segment_lookup = segment_lookup,
                           image_dest = qc_png,
                           max_dim = 32000)

  xml_outfile <- file.path(output_dir, paste0("roi_overlay_segment_aoi_", tiff_id, ".xml"))
  write_halo_xml_by_segment(contour_results = contour_results,
                           segment_lookup = segment_lookup,
                           aoi_name_lookup = aoi_name_lookup,
                           outfile = xml_outfile)
}



