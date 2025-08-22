library(circlize)
library(grid)
library(gridBase)
library(ComplexHeatmap)

#' Generate Circular Visualization for Enrichment Analysis Results
#'
#' @param EnrichDat Dataframe of enrichment results (must contain: ID, Category, Gene, Total, Pvalue, adjustPvalue)
#' @param gene_diff Dataframe of differential genes (must contain: avg_log2FC with row names as gene symbols)
#' @param adjust.pvalue Logical, whether to use adjusted p-value (default: FALSE)
#' @param show_RichFactor Logical, whether to show RichFactor values (default: TRUE)
#' @param title_override Optional custom title for the plot
#' @param color_palette Optional custom color palette for categories
#' @return Generates circos plot (no return value)
PlotEnrichCircos <- function(EnrichDat, gene_diff, adjust.pvalue = FALSE,rankby="pvalue",topn=20,
                            show_RichFactor = TRUE, title_override = NULL,
                            color_palette = NULL) {
  
  # Validate input data
  if (!all(c("ID", "GeneRatio","BgRatio", "geneID",  "pvalue") %in% colnames(EnrichDat))) {
    stop("EnrichDat must contain columns: ID, geneID, Total, pvalue")
  }
  if ("Subontologies" %in% colnames(EnrichDat)) {
    EnrichDat$Category <- EnrichDat$Subontologies
	}else if ("level_a" %in% colnames(EnrichDat)) {
	  EnrichDat$Category <- EnrichDat$level_a
		EnrichDat$geneID = EnrichDat$geneSymbol
	}else if (!"Category" %in% colnames(EnrichDat)) {
	  stop("EnrichDat must contain 'Category' column")
	}

  if (!"avg_log2FC" %in% colnames(gene_diff) || is.null(rownames(gene_diff))) {
    stop("gene_diff must have 'avg_log2FC' column and gene symbols as row names")
  }
  
  # Set colors
  default_colors <-c(
        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
        "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
        "#CCEBC5", "#FFED6F", "#71A99F", "#CCCC8F", "#9895AE",
        "#C9665B", "#668EA9", "#CA904E", "#8FB254", "#CAA4B7"
    )

  color <- if (is.null(color_palette)) default_colors else color_palette

  second_col <- c('#C51B7D', '#F8F0F4', '#4D9221')
  upcol <- "#F6756F"
  downcol <- "#1BBFC3"
  
  # Data preprocessing
  diff_gene <- setNames(gene_diff$avg_log2FC, gene_diff$gene)
  EnrichDat <- EnrichDat[order(EnrichDat[[rankby]])[1:topn],]
	EnrichDat$Total = as.numeric(gsub("(\\d+)\\/\\d+","\\1",EnrichDat$BgRatio))
	EnrichDat$RichFactor <- EnrichDat$Count/EnrichDat$Total
  circlize_df <- EnrichDat[order(EnrichDat$Category), ]
  circlize_df$Category_color <- color[as.numeric(factor(circlize_df$Category))]
  rownames(circlize_df) <- circlize_df$ID
  
  # Prepare track data
  circlize_df$gene_num.min <- 0
  circlize_df$gene_num.rich <- circlize_df$Total
  circlize_df$gene_num.max <- max(circlize_df$gene_num.rich)
  
  # Handle p-values
  if (!adjust.pvalue) {
    circlize_df$stat_val <- -log10(circlize_df$pvalue)
    title_second_data <- "-Log10(Pvalue)"
  } else {
    circlize_df$stat_val <- -log10(circlize_df$p.adjust)
    title_second_data <- "-Log10(adjust.Pvalue)"
  }
  circlize_df$stat_val[is.infinite(circlize_df$stat_val)] <- 325
  
  # Initialize plot
  circle_size <- unit(1, "snpc")
  circos.par(gap.degree = 0.5, start.degree = 90)
  plot.new()
  pushViewport(viewport(
    x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")
  ))
  par(omi = gridOMI(), new = TRUE)
  
  # Track 1: ID labels
  plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.max")]
  circos.genomicInitialize(plot_data, plotType = NULL)
  circos.track(
    ylim = c(0, 1), track.margin = c(0, 0), track.height = 0.08, bg.border = NA, 
    bg.col = circlize_df$Category_color,
    panel.fun = function(x, y) {
      ylim <- get.cell.meta.data("ycenter")
      xlim <- get.cell.meta.data("xcenter")
      sector.name <- get.cell.meta.data("sector.index")
      circos.axis(h = "top", labels.cex = 0.6, labels.pos.adjust = TRUE, 
                 labels.col = "#000000", labels.niceFacing = FALSE)
      circos.text(xlim, ylim, sector.name, cex = 0.6, col = "#000000", 
                 niceFacing = FALSE)
    }
  )
  
  # Track 2: Enrichment significance
  plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.rich", "stat_val")]
  label_data <- circlize_df[, "gene_num.rich", drop = FALSE]
  p_max <- round(max(circlize_df[["stat_val"]])) + 1
  colorsChoice <- colorRampPalette(second_col)
  color_assign <- colorRamp2(0:p_max, colorsChoice(p_max + 1))
  
  circos.genomicTrackPlotRegion(
    plot_data,
    track.margin = c(0, 0), track.height = 0.1, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), 
                        border = NA, ...)
      ylim <- get.cell.meta.data("ycenter")
      xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
      sector.name <- label_data[get.cell.meta.data("sector.index"), 1]
      circos.text(xlim, ylim, sector.name, cex = 0.8, col = "#000000", 
                 niceFacing = FALSE)
    }
  )
  
  # Track 3: Up/down regulated genes
  circlize_df$Up <- unlist(lapply(strsplit(circlize_df$geneID,'/'), function(x) sum(diff_gene[x]>0)))
  circlize_df$Down <- unlist(lapply(strsplit(circlize_df$geneID,'/'), function(x) sum(diff_gene[x]<=0)))
  
  if(any(circlize_df$Up > 0 & circlize_df$Down > 0)) {
    circlize_df$up.proportion <- circlize_df$Up / (circlize_df$Up+circlize_df$Down)
    circlize_df$down.proportion <- circlize_df$Down / (circlize_df$Up+circlize_df$Down)
    circlize_df$up_v <- circlize_df$up.proportion * circlize_df$gene_num.max
    plot_data_up <- circlize_df[, c("ID", "gene_num.min", "up_v")]
    names(plot_data_up) <- c("id", "start", "end")
    plot_data_up$type <- 1
    
    circlize_df$down_v <- circlize_df$down.proportion * circlize_df$gene_num.max + circlize_df$up_v
    plot_data_down <- circlize_df[, c("ID", "up_v", "down_v")]
    names(plot_data_down) <- c("id", "start", "end")
    plot_data_down$type <- 2
    
    plot_data <- rbind(plot_data_up, plot_data_down)
    label_data <- circlize_df[, c("up_v", "down_v", "Up", "Down")]
    label_data <- as.data.frame(label_data)
    rownames(label_data) <- circlize_df$ID
    color_assign <- colorRamp2(breaks = c(1, 2), col = c(upcol, downcol))
    
    suppressWarnings(suppressMessages(circos.genomicTrackPlotRegion(
      plot_data,
      track.margin = c(0, 0), track.height = 0.1, bg.border = NA, stack = TRUE,
      panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
        ylim <- get.cell.meta.data("cell.bottom.radius") - 0.5
        xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
        sector.name <- ifelse(label_data[get.cell.meta.data("sector.index"), 3] > 0, 
                             label_data[get.cell.meta.data("sector.index"), 3], " ")
        circos.text(xlim, ylim, sector.name, cex = 0.6, col = "#000000", niceFacing = FALSE)
        xlim <- (label_data[get.cell.meta.data("sector.index"), 2] + 
                label_data[get.cell.meta.data("sector.index"), 1]) / 2
        sector.name <- ifelse(label_data[get.cell.meta.data("sector.index"), 4] > 0, 
                             label_data[get.cell.meta.data("sector.index"), 4], " ")
        circos.text(xlim, ylim, sector.name, cex = 0.6, col = "#000000", niceFacing = FALSE)
      }
    )))
  } else {
    plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.rich")]
    label_v <- ifelse(sum(circlize_df$Up)>0, "Up", "Down")
    label_col <- ifelse(sum(circlize_df$Up)>0, upcol, downcol)
    label_data <- as.data.frame(cbind(circlize_df[, "gene_num.rich"], circlize_df[[label_v]]))
    rownames(label_data) <- circlize_df$ID
    
    suppressWarnings(suppressMessages(circos.genomicTrackPlotRegion(
      plot_data,
      track.margin = c(0, 0), track.height = 0.1, bg.border = NA, stack = TRUE,
      panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = label_col, border = NA, ...)
        ylim <- get.cell.meta.data("cell.bottom.radius") - 0.6
        xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
        sector.name <- label_data[get.cell.meta.data("sector.index"), 2]
        circos.text(xlim, ylim, sector.name, cex = 0.6, col = "#000000", niceFacing = FALSE)
      }
    )))
  }
  
  # Track 4: Rich factor
  plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.max", "RichFactor")]
  label_data <- circlize_df[, c("Category", "RichFactor")]
  label_data$RichFactor <- round(label_data$RichFactor, 2)
  label_data <- as.data.frame(label_data)
  rownames(label_data) <- circlize_df$ID
  
  circos.genomicTrack(
    plot_data,
    track.margin = c(0.01, 0.04), track.height = 0.3, ylim = c(0.05, 0.95), 
    bg.col = "gray95", bg.border = NA,
    panel.fun = function(region, value, ...) {
      sector.name <- get.cell.meta.data("sector.index")
      for (i in seq(0.1, 0.9, by = 0.1)) {
        circos.lines(c(0, max(region)), c(i, i), col = "gray", lwd = 0.3)
      }
      circos.genomicRect(region, value, col = circlize_df[sector.name,'Category_color'], 
                        border = NA, ytop.column = 1, ybottom = 0, ...)
      if (show_RichFactor) {
        ylim <- label_data[get.cell.meta.data("sector.index"), 2] - 0.05
        xlim <- get.cell.meta.data("xcenter")
        sector.name <- label_data[sector.name, 2]
        circos.text(xlim, ylim, sector.name, cex = 0.8, col = "#000000", niceFacing = FALSE)
      }
    }
  )
  
  circos.clear()
  
  # Draw legends
  updown_legend <- Legend(
    labels = c("Number of Genes", "Up-regulated", "Down-regulated", "Rich Factor(0-1)"),
    graphics = list(
      function(x = 0, y, w, h) {
        grid.draw(gTree(
          children = gList(
            rectGrob(x = 0, y, w * 2.2, h * 1, gp = gpar(fill = "gray80", col = "gray80")),
            textGrob("100", x = 0, y)
          ),
          gp = gpar(col = "black", fontsize = 10, fontface = "bold")
        ))
      },
      function(x = 0, y, w, h) {
        grid.rect(x = 0, y, w * 2.2, h * 1, gp = gpar(fill = upcol, col = upcol))
      },
      function(x = 0, y, w, h) {
        grid.rect(x = 0, y, w * 2.2, h * 1, gp = gpar(fill = downcol, col = downcol))
      },
      function(x = 0, y, w, h) {
        grid.polygon(y = c(0.09, -0.05, -0.12, 0.16), x = c(-0.13, -0.13, 0.13, 0.13), 
                    gp = gpar(fill = "gray85", col = "gray85"))
      }
    ),
    labels_gp = gpar(fontsize = 10),
    row_gap = unit(2, "mm")
  )
  
  category_legend <- Legend(
    labels = unique(circlize_df$Category),
    type = "points", pch = NA, background = unique(circlize_df$Category_color),
    labels_gp = gpar(fontsize = 10), row_gap = unit(1, "mm")
  )
  
  pvalue_legend <- Legend(
    col_fun = colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1)),
    legend_height = unit(3, "cm"), labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10), title_position = "topleft", 
    title = paste0(title_second_data, "\n"), row_gap = unit(3, "mm")
  )
  
  lgd_list_vertical <- packLegend(updown_legend)
  lgd_list_vertical2 <- packLegend(category_legend, pvalue_legend)
  
  draw(lgd_list_vertical, just = "center")
  upViewport()
  draw(lgd_list_vertical2, x = circle_size, just = "left")
}