#' Enhanced ggplot2 Saver with Multi-Format Support
#'
#' An extension of \code{ggplot2::ggsave} that supports batch export to multiple formats
#' with intelligent defaults for scientific publishing.
#'
#' @param filename Character. Output file path with or without extension. If no extension
#'                 is provided, both PDF and PNG will be generated.
#' @param plot ggplot object. Plot to save, defaults to last plot displayed.
#' @param ... Additional arguments passed to \code{ggplot2::ggsave}:
#'            \describe{
#'              \item{width, height}{Plot dimensions in inches (default: 7)}
#'              \item{dpi}{Resolution for raster formats (default: 300)}
#'              \item{bg}{Background color (default: transparent)}
#'              \item{limitsize}{Logical. Prevent large outputs? (default: FALSE)}
#'            }
#'
#' @return Invisibly returns the input plot. Primary side effect is files written to disk.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(mpg, hp)) + geom_point()
#'
#' # Save as PDF+PNG with default settings
#' ggsave("output/plot")
#'
#' # Custom dimensions and single format
#' ggsave("results/figure1.png", width = 5, height = 4, dpi = 600)
#' }
#'
#' @export
ggsave <- function(filename, plot = last_plot(), ...) {
  dots <- list(...)
#   if (!"filename" %in% names(dots)) stop("必须提供 filename 参数")
#   filename <- dots$filename
  ext <- tools::file_ext(filename)
	if (!ext %in% c("pdf", "png", "svg", "jpeg", "jpg", "tiff", "bmp")) {ext=""}
  
  if(ext == ""){
    ext = c("pdf","png")
		base_name = filename
  }else{
		base_name <- tools::file_path_sans_ext(filename)
	}
  # 自定义设备函数
  # 设置默认设备
  for (e in ext){
      device <- switch(
        tolower(e),
        pdf = Cairo::CairoPDF,
        png = Cairo::CairoPNG,
        svg = Cairo::CairoSVG,
        jpeg = Cairo::CairoJPEG,
        jpg = Cairo::CairoJPEG,
        tiff = Cairo::CairoTIFF,
        bmp = Cairo::CairoBMP,
        stop(paste0("不支持的文件格式: ", e))
      )
      if(tolower(e) == "png"){
        dpi = ifelse("dpi" %in% names(dots), dots$dpi, 300)
        # 用户是否传入了 width/height？
        width <- ifelse("width" %in% names(dots), dots$width, 7)
        height <- ifelse("height" %in% names(dots), dots$height, 7)

        dots$width <- round(width * dpi)
        dots$height <- round(height * dpi)
        dots$dpi <- dpi
      }
      dots$limitsize <- ifelse("limitsize" %in% names(dots), dots$limitsize, FALSE)
      # 使用 do.call 调用 ggsave，避免参数冲突
      args <- c(list(filename = paste0(base_name,".",e), plot = plot), dots, list(device = device))

      do.call(ggplot2::ggsave, args)
      #message("Saved: ", paste0(base_name,".",e), " ; height:",dots$height, " | width:",dots$width)
  }
}

#' Plot Cell Type Abundance as Stacked Barplot
#'
#' Visualize cell type proportions across samples/groups using stacked barplots.
#'
#' @param x Seurat object containing cell metadata
#' @param prop.by Character. Metadata column name for cell types/clusters (default: "res.1")
#' @param group.by Character. Metadata column name for sample/group grouping (default: "sampleid")
#' @param split.by Character. Optional metadata column for faceting (default: NULL)
#' @param cols Character vector. Color palette for cell types (default: auto-generated)
#' 
#' @return ggplot2 object with stacked barplot visualization
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' PlotBarplot(seurat_obj)
#' 
#' # Custom grouping and colors
#' PlotBarplot(seurat_obj, 
#'                prop.by = "celltype", 
#'                group.by = "patient",
#'                cols = c("B" = "blue", "T" = "red"))
#' }
#' @export
PlotBarplot <- function(
  x,
  prop.by = "seurat_clusters",
  group.by = "sampleid",
  split.by = NULL,
  cols = NULL
) {
  # Input validation
  if (!inherits(x, "Seurat")) {
    stop("Input must be a Seurat object", call. = FALSE)
  }
  
  md <- x@meta.data
  required_cols <- c(prop.by, group.by)
  missing_cols <- setdiff(required_cols, colnames(md))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in metadata:",
               paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }
  
  # Calculate proportions
  counts <- table(md[, prop.by], md[, group.by])
  df <- as.data.frame(
    prop.table(counts, margin = 2) * 100,
    responseName = "freq"
  )
  colnames(df)[1:2] <- c(prop.by, group.by)
  
  # Add metadata for faceting
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(md)) {
      stop(paste("split.by column", split.by, "not found in metadata"),
           call. = FALSE)
    }
    
    facet_data <- unique(md[, c(group.by, split.by)])
    df <- merge(df, facet_data, by = group.by)
  }
  
  # Create plot
  p <- ggplot(df, aes_string(x = group.by, y = "freq", fill = prop.by)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
    labs(x = NULL, y = "Proportion [%]") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_rect(fill = NA, color = NA),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  
  # Add faceting if specified
  if (!is.null(split.by)) {
    p <- p + facet_grid(reformulate(split.by), 
                       scales = "free_x", 
                       space = "free_x")
  }
  
  return(p)
}


#' Calculates The Optimum Dimensionality From An Elbow Plot.
#'
#' A function that calculates the optimum dimensionality from an elbow plot. Involves calculating
#' slopes between two consecutive points. The optimum dimensionality will be determined as the point at
#' which slopes no longer change or no visible change occurs to the slopes thereafter. In this function, a dimension
#' is determined optimal when a slope connecting that point to another is flatter than 10 times the flattest slope
#' in the plot and the stdev at that point is proximal to the lower limit.
#' The plot data must be plotted by dims (x-axis) and stdev (y-axis).
#'
#' @param object An elbow plot data used for detemining dimensionality, generated by calling ElbowPlot() on a Seurat object.
#' @param reduction the reduction method to use, default "pca".
#' @param method the method used find elbow point on the screen plot, choice can be "slope"(default), "Kerstin".
#' @param ndims the maximium number of components used to determine the elbow point, 30 as default.
#' @return Returns an integer that represents the optimal dimensionality that was determined by calculation.
#' @examples
#' library(Seurat)
#' FindElbow(p$data)
#' @export
FindElbow <- function(
  object,
  reduction = "pca",
  method = "slope",
  ndims = 30
) {
  method <- match.arg(method, c("slope", "Kerstin"))
  data.use <- SeuratObject::Stdev(object, reduction = reduction)
  if (length(data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }
  if (ndims > length(data.use)) {
    warning("The object only has information for ", length(x = data.use), " reductions")
    ndims <- length(x = data.use)
  }

  optimal_dims <- switch(method,
          slope = {
            data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])

            # Prepare all required parameters from plot_data.
            dimensions <- data$dims
            total_dims <- length(dimensions)
            stdev <- data$stdev
            slopes_so_far <- c() # Will keep track of all the slopes (between two consecutive points) calculated.
            last_dim <- length(dimensions)

            # Calculates slopes between every pair of consecutive points.
            for (dim in 2:last_dim) {
              slopes_so_far <- c(slopes_so_far, (stdev[dim] - stdev[dim - 1]))
            }

            # Default dimensionality to return is 1.
            dimensionality <- 1
            nth_point <- 2
            for (slope in slopes_so_far) {
              # After trial and error, these conditions seem to work well to determine the elbow.
              # Note: slopes are negative in this plot, so using max(slopes_so_far).
              # The first condition looks to see if we are approaching a lower limit, but this wouldn't be enough to
              # confidently say that it is, so we need the second condition.
              # The second condition makes sure that the point used to determine the slope is very close to the
              # smallest value, which means that it increases the confidence that this slope is indeed approaching
              # a lower limit.
              if (slope > 10 * max(slopes_so_far) && stdev[nth_point] < min(stdev) + ((max(stdev) - min(stdev)) * 0.05)) {
                dimensionality <- nth_point
                return(dimensionality - 1) # Subtract 1 because we want the point that leads into the lower limit or the
                # flat part of the graph.
              }
              nth_point <- nth_point + 1
            }
            dimensionality
          },
          Kerstin = {
            data <- Embeddings(object, reduction = reduction)[, 1:ndims]
            dims <- round(as.numeric(intrinsicDimension::maxLikGlobalDimEst(data, k = 20)))
          }
  )

  # The function could not determine dimensionality or elbow from the plot.
  return(optimal_dims)
}

#' 严格模式列名验证
#' @param object Seurat对象
#' @param columns 用户输入的列名（可含大小写/前后空格）
#' @return 标准化后的有效列名
validate_metadata_columns_strict <- function(object, columns) {
  # 分割多参数（处理 splitby="sample, group" 情况）
  input_cols <- trimws(unlist(strsplit(columns, ",\\s*")))
  
  # 获取实际元数据列名（保留原始大小写）
  meta_cols <- colnames(object@meta.data)
  found_cols <- character(0)
  
  # 大小写不敏感匹配
  for (col in input_cols) {
    # 标准化：仅处理前后空格，保留中间空格（如"seurat clusters"视为整体）
    clean_col <- trimws(col)
    
    # 精确匹配（区分大小写）
    if (clean_col %in% meta_cols) {
      found_cols <- c(found_cols, clean_col)
      next
    }
    
    # 大小写不敏感匹配
    matched <- meta_cols[tolower(meta_cols) == tolower(clean_col)]
    if (length(matched) > 0) {
      found_cols <- c(found_cols, matched[1])
      warning("Case-insensitive match: ", col, " -> ", matched[1])
    } else {
      stop("Column not found: ", col, 
           "\nAvailable columns: ", paste(meta_cols, collapse = ", "))
    }
  }
  
  return(unique(found_cols))
}




VisualizeCluster <- function(single_ob,reductions="umap,tsne",groupby="seurat_clusters",splitby="sample,group",outdir="./",colors,pointsize=NULL){
    if ( is.null(pointsize) ){
      if (dim(single_ob)[2] < 500){
          pointsize = 1.5
      } else pointsize = 0.5
		} else {
				pointsize = pointsize
		}
    #容错判断
    splitby = validate_metadata_columns_strict(single_ob,splitby)
		groupby = validate_metadata_columns_strict(single_ob,groupby)
		reductions = unlist(strsplit(reductions,","))
		if(!is.factor(single_ob@meta.data[,groupby])){single_ob@meta.data[,groupby] = as.factor(single_ob@meta.data[,groupby])}
    color_use = colors[1:length(levels(single_ob@meta.data[,groupby]))]
		names(color_use) = levels(single_ob@meta.data[,groupby])
		for (reduction in reductions){
			output_dir <- file.path(outdir,reduction)
		  if (!dir.exists(output_dir)) {dir.create(output_dir,recursive = TRUE,mode = "777")}

			# 总体绘图的图片（基础款）
			for (tag in c(groupby,splitby)){
				tag_vis = Seurat::DimPlot(object = single_ob, dims = c(1,2),
																		reduction = reduction,
																		pt.size = as.numeric(pointsize),
																		group.by = tag)+
										ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) + theme(aspect.ratio = 1/1)
				if (tag =="seurat_clusters"){
					cell_count = table(single_ob$seurat_clusters)
    			#cell_count_labels = paste(paste(names(cell_count),cell_count,sep="-")," cells")
					cell_count_labels = paste0("C",paste(names(cell_count),cell_count,sep=":")," cells")
					tag_vis = tag_vis+ ggplot2::scale_colour_manual( values = colors,
                                  breaks=levels(single_ob$seurat_clusters),
                                  labels = cell_count_labels)

				}else{
					tag_vis = tag_vis+ggplot2::scale_colour_manual( values = colors)
				}

				width = ifelse(length(single_ob@meta.data[,tag]) > 20, 7.4, 9)
				ggsave(file.path(output_dir,paste0(reduction,"_reduction_groupby_",tag)),tag_vis,width =width,height = 7)
			}
      # 总体绘图的图片（升级款，只支持groupby参数）
			p=dimplot_integrated(seurat_obj=single_ob,reduct =reduction,groupby =groupby,mode='level',color_use = colors,label=T,label_color=T,pt.size=pointsize)
      ggsave(file.path(output_dir,paste0(reduction,"_reduction_density_isobars_groupby_",groupby)),p,width = 7,height = 7)
		}
		# splitby参数的图片
		for (reduction in reductions){
			output_dir <- file.path(outdir,reduction)
			for (tag in splitby){
				#----------------
				if ( !is.null(tag) ){
             nrow = ceiling(length(unique(single_ob@meta.data[,tag]))/2)
         }

				tagvis_split = Seurat::DimPlot(object = single_ob,
                          dims = c(1,2),reduction = reduction,
                          pt.size = as.numeric(pointsize), ncol = 2,
                          order = rev(names(color_use)),
                          group.by = groupby, split.by = tag)+
             ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5))
         tagvis_split = tagvis_split + theme(aspect.ratio = 1/1, legend.position = "right")+
                          ggplot2::scale_colour_manual( values = color_use,
                                               breaks = names(color_use),
                                               labels = names(color_use) )

         nlevel = length(unique(single_ob@meta.data[,groupby]))
         max.nchar = max(nchar(as.vector(unique(single_ob@meta.data[,groupby]))))
         width_char = ifelse(nlevel <= 20,1 + max.nchar/10, 2 + max.nchar/5)
         ggsave(file.path(output_dir,paste0(reduction,"_reduction_groupby_",groupby,"_splitby_",tag)),
                                    tagvis_split, limitsize = FALSE, 
                                    width = 8 + width_char, height = 2 + 4*nrow, bg="white" )

        #-----------------------------
			}
		  
		}
    
}

Plot_dodge_barplot=function(data,x.by = "seurat_clusters", group.by = "group",cols= colors){
  df <- as.data.frame(data)
  df[[x.by]] = factor(as.character(df[[x.by]]), levels= levels(data[[x.by]]))
  df[[group.by]]=factor(as.character(df[[group.by]]),levels=levels(data[[group.by]]))
  if (x.by == "seurat_clusters") {
  angle_size = 0
  hjust_size = 0.5
  vjust_size = 0
  } else {
  angle_size = 45
  hjust_size = 0.85
  vjust_size = 0.85
  }

  p <- ggplot(df, aes(x = !!sym(x.by), y = cell_number, fill = !!sym(group.by), color = !!sym(group.by))) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.85),alpha = 0.5,width = 0.7) +
  # 误差线（均值标准误）
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, position = position_dodge(0.85), color = "black",  # 统一误差线颜色
    size = 0.8 ) + geom_point( position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),size = 2, alpha = 0.8) +
  scale_fill_manual(values = cols) + scale_color_manual(values = cols) + theme_bw() +
  labs( x = gsub("_", " ", x.by),  # 替换下划线为空格
    y = "Cell Number", fill = group.by, color = group.by) +
  theme(axis.text = element_text(size = 10, color = "black", face = "bold"), axis.text.x = element_text(
      angle = angle_size, hjust = hjust_size,vjust = vjust_size ),strip.background = element_rect(fill = NA,color = NA), 
          axis.title.x=element_text(size=18),axis.title.y=element_text(size=14),panel.grid =element_blank())
	#将x,y轴的原点改为c(0,0)	
  p = p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  return(p)

}




Plot_dodge_barplot_withP=function(data,x.by = "seurat_clusters", group.by = "group",cols= summary_colors,output_dir){
  df <- as.data.frame(data)
  df[[x.by]] = factor(as.character(df[[x.by]]), levels= levels(data[[x.by]]))
  df[[group.by]]=factor(as.character(df[[group.by]]),levels=levels(data[[group.by]]))
  stat.test <- df %>%  group_by(.dots=x.by) %>% wilcox_test(as.formula(paste("cell_number ~", group.by))) %>% add_significance("p")  ### wilcox_test
  if ("p.adj" %in% names(stat.test) & "p.adj.signif" %in% names(stat.test)){
  stat.test=stat.test %>% select(-c("p.adj","p.adj.signif"))
  }
  #输出统计结果到表格
  output_data=stat.test %>% select(-c(".y."))
  write.table(output_data,file.path(output_dir,paste0("groupby_",group.by,"_cell_proportion_statistic_pvalue.xls")),sep="\t",row.names=FALSE)
  #都是ns时全部ns都展示。不都是ns时，隐藏掉ns值
  if (all(stat.test[["p.signif"]] == "ns")){
  stat.test <- stat.test %>% add_xy_position(x = x.by, dodge = 0.85,step.increase=0.03)
  }else{
  stat.test <- subset(stat.test,p.signif!="ns")
  max_freq <- df %>%
  group_by(.dots=x.by) %>%
  summarise(max_freq = max(cell_number, na.rm = TRUE)+0.02*max(df$cell_number))
  stat.test <- stat.test %>%
    left_join(max_freq, by = x.by) %>%
    rename(y.position = max_freq) %>% add_x_position(x = x.by, dodge = 0.85)
  }

  if (x.by == "seurat_clusters") {
  angle_size = 0
  hjust_size = 0.5
  vjust_size = 0
  } else {
  angle_size = 45
  hjust_size = 0.85
  vjust_size = 0.85
  }
  p<- ggplot(df, aes_string(x = x.by, y = "cell_number",fill = group.by)) +  #下面fun后边跟的median或者mean一定要加引号！！否则会报错
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(0.85),alpha=0.5,width=0.7) +
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width = 0.4,position = position_dodge(0.85),color="black")+
  geom_point(aes_string(fill = group.by),position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),size = 2,alpha = 0.8) +
   scale_fill_manual(values = cols )+ theme_bw() + scale_color_manual(values = cols )+
   labs(x = gsub("_"," ",x.by) ,y = "Cell Number") +
   theme(axis.text=element_text(size=10, angle = 0,color="black",face = "bold"),
         axis.text.x=element_text(size=10, hjust = hjust_size, angle = angle_size, vjust = vjust_size),
         strip.background = element_rect(fill = NA,color = NA), 
         axis.title.x=element_text(size=18),axis.title.y=element_text(size=14),panel.grid =element_blank()) +
    stat_pvalue_manual(
    stat.test, label = "p.signif", tip.length = 0.01,step.increase=0.05,step.group.by=x.by,label.size = 6) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
	return(p)
}





SummaryCluster <- function(single_ob,groupby="seurat_clusters",splitby="sample.name,group",outdir="./",colors){
	#容错判断**-
    splitby = validate_metadata_columns_strict(single_ob,splitby)
		groupby = validate_metadata_columns_strict(single_ob,groupby)
    library(ggalluvial) 
#普通比例柱状图和ggalluvial图
		for (split in splitby){
			p<- PlotBarplot(single_ob,prop.by = groupby, 
                group.by = split, cols =colors)
			ggsave(file.path(outdir,paste0("SummaryCluster_",groupby,"_splitby_",split)),p,limitsize = FALSE, width = 7.5, height = 7, bg="white" )
			DATA <- suppressMessages(single_ob@meta.data[,c(groupby,split)] %>%
                dplyr::group_by( .dots= c(groupby,split)) %>%
                dplyr::summarize(cell_number = dplyr::n()))
			df = merge(p$data,DATA,by=c(groupby,split),all.x=TRUE)
			write.table(df,file.path(outdir,paste0("SummaryCluster_",groupby,"_splitby_",split,".xls")),sep="\t",row.names=FALSE)

      rel_allu <- ggplot(data = df,
							aes(x = !!sym(split), y = freq, fill = !!sym(groupby), 
									stratum =!!sym(groupby),
									alluvium = !!sym(groupby))) +
					geom_alluvium()+
						geom_stratum(width=0.45, size=0.1) +
					geom_bar(aes(fill = !!sym(groupby)),stat='identity', width=0.45) +
					scale_y_continuous(expand=c(0, 0))+
					theme_bw() +
					#facet_grid( ~ grazing,scales = "fixed")+
					scale_fill_manual(values = colors,name=groupby) +
					scale_color_manual(values = colors) +
					theme(legend.position = "top",
								panel.grid=element_blank(),
								panel.spacing.x = unit(0,units = "cm"),
								strip.background = element_rect(
						color="white", fill="white", 
						linetype="solid",size=0.8),
						strip.placement = "outside",
						#axis.line.y.left = element_line(color="black",size=0.8),
						#axis.line.x.bottom = element_line(color="black",size=0.8),
						strip.text.x = element_text(size = 14,face = "bold"),
						axis.text = element_text(face = "bold", 
																				size = 12,color="black"),
						axis.title = element_text(face = "bold", 
																					size = 14,colour = "black"),
							legend.title = element_text(face = "bold", 
																					size =12,color="black"),
						legend.text = element_text(face = "bold", size =12,color="black"),
						axis.ticks.x = element_blank(),
						axis.ticks.y = element_line(size = 0.3),
						
								)+
					labs(x = split,y= paste0("Relative Abundance of ", groupby,"(%)"))
        ggsave(file.path(outdir,paste0("SummaryCluster_",groupby,"_splitby_",split,"_alluvial")),rel_allu,limitsize = FALSE, width = 7.5, height = 7, bg="white" )
		}

#绘dogebarplot
#基于sample.name绘制加p值的图
		splitby <- splitby[!splitby %in% c("sample.name","seurat_clusters")]
		for (i in splitby) {
        group.by=i
        result <- single_ob@meta.data %>%
          dplyr::group_by(sample.name) %>%
          dplyr::summarise(unique_num = dplyr::n_distinct(!!sym(group.by)))
        if (!all(result$unique_num)==1){
            splitby <- setdiff(splitby, i)
        }
    }
    DATA <- suppressMessages(single_ob@meta.data[,c("sample.name", groupby,splitby)] %>%
                dplyr::group_by( .dots= c("sample.name", groupby,splitby)) %>%
                dplyr::summarize(cell_number = dplyr::n()))

		for (split.by in splitby) {
      x_by=groupby
     #判断是否每个横坐标有2个及2个以上的分组没有样本
      total_sub_num=c()
      ##首先，计算一共有多少个分组
      num_group=length(unique(DATA[[split.by]]))
      ##计算每个横坐标中有多少个分组存在
      for (i in unique(DATA[[groupby]])){
      sub_data=subset(DATA,DATA[[groupby]]==i)
      ##计算差值，算出各横坐标中缺少几个分组
      sub_num=num_group-length(unique(sub_data[[split.by]]))
      ##统计各个横坐标中缺少分组的数值，后续会判断是否都小于2
      total_sub_num=c(total_sub_num,sub_num)
      }
      #填充每个cluster缺失的group的cellnumber为0
      # 步骤1: 为每个sampleid生成完整的clusters序列
      full_clusters <- expand.grid(sample.name = unique(DATA$sample.name), x_by = unique(DATA[[x_by]]))
      colnames(full_clusters)[which(colnames(full_clusters) == "x_by" )] <- x_by
      # 步骤2: 合并原始数据和完整的clusters序列  
      df_full <- merge(full_clusters, DATA, by = c("sample.name", x_by), all.x = TRUE)
      # 步骤3: 填充缺失的group和cellnumber  
      df_full[[split.by]][is.na(df_full[[split.by]])] <-   
        sapply(df_full$sampleid[is.na(df_full[[split.by]])], function(sid) {  
          DATA[[split.by]][DATA$sampleid == sid][1]
        })  
      df_full$cell_number[is.na(df_full$cell_number)] <- 0
      #df_full$freq[is.na(df_full$freq)] <- 0 
      DATA=df_full
      if (is.factor(single_ob@meta.data[[groupby]]) == FALSE) {
          if (groupby=="seurat_clusters"){
            single_ob@meta.data[[groupby]] = factor(single_ob@meta.data[[groupby]], levels = unique(sort(as.numeric(single_ob@meta.data[[groupby]]))))
          } else {
            single_ob@meta.data[[groupby]] = as.factor(single_ob@meta.data[[groupby]])
          }
      }
      if (is.factor(single_ob@meta.data[[split.by]]) == FALSE) {
          single_ob@meta.data[[split.by]] = as.factor(single_ob@meta.data[[split.by]])
      }

			DATA[[groupby]] = factor(DATA[[groupby]], levels = levels(single_ob@meta.data[[groupby]]))
      DATA[[split.by]] = factor(DATA[[split.by]], levels = levels(single_ob@meta.data[[split.by]]))
      dodge_barplot=Plot_dodge_barplot(DATA,x_by,split.by,cols= colors)

      ggsave(file.path(outdir,paste0(split.by,"_dodge_barplot_groupby",x_by)),
          dodge_barplot,width = 0.25*length(unique(DATA[[x_by]]))*length(unique(DATA[[split.by]]))+3, height = 6, bg="white" ,limitsize = FALSE)
      #判断是否组内有至少3个重复样本。若有，才出带p值的图
      n_sampleid <- as.data.frame(DATA) %>%
      group_by(.dots=split.by) %>%
      summarize(count_distinct_sampleid = n_distinct(sample.name))
      
      if (length(unique(DATA[[split.by]])) >= 2) {
				if (all(n_sampleid$count_distinct_sampleid >= 3)){
					#每个横坐标有2个及2个以上的分组中没有细胞，不绘制带p值的dodge barplot
					if (all(total_sub_num < 2)) {
					dodge_barplot_withP=Plot_dodge_barplot_withP(DATA,groupby,split.by,cols= colors,outdir)
					ggsave(file.path(outdir,paste0(split.by,"_dodge_barplot_by",groupby,"_withpvalue",collapse="")),
						dpi=300, dodge_barplot_withP,width = 0.25*length(unique(DATA[[groupby]]))*length(unique(DATA[[split.by]]))+3, height = 6, bg="white",limitsize = FALSE )
					}else{
								print(paste0(split.by,"某个横坐标中，所有样本cell_number都为0的分组数大于二，不出带p值的图"))
								}
				}else{
							print(paste0(split.by,"不满足每个组内有至少3个重复样本，不出带p值的图"))
				}
			}else {
				print(paste0(split.by,"只有1个分组，不出带p值的图"))
			}
    }
}


Topn_ave_marker_heatmap<-function(immune.combined,collapseby = collapseby,Topn = 30,seurat_exp_cluster_dir=seurat_exp_cluster_dir,markers=markers,colors =colors ){
        suppressPackageStartupMessages(library(ComplexHeatmap))
        cluster_topn_markers=markers %>% group_by(cluster)%>% top_n(n = Topn, wt = avg_log2FC)
        out_dir = seurat_exp_cluster_dir
        #write.table(cluster_topn_markers, file.path(out_dir,glue::glue("cluster_top{Topn}_markers_avg_heatmap.xls")),quote = F, row.names = F, sep = "\t")

        markers2vis = cluster_topn_markers$gene
        count = as.matrix(GetAssayData(immune.combined, slot = "data")[markers2vis,])
        metadata = immune.combined@meta.data
        metadata$id = rownames(metadata)
        collapsed_count = vector()
        if ( !collapseby %in% colnames(metadata) ){
            stop("NO specified column found!")
        }
        collapsed_group = metadata %>% group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))
        for ( cells in collapsed_group$cellid ){
            samplex = unlist(strsplit(cells, ",", perl =T))
            collapsed_count= cbind(collapsed_count,rowMeans( count[,samplex,drop=F] ))
        }
        collapsed_count = as.matrix( collapsed_count )
        collapsed_group = as.data.frame(collapsed_group)
        colnames(collapsed_count) = as.matrix(collapsed_group[,1])
        data = tibble::rownames_to_column(as.data.frame(collapsed_count),"GeneID")
        #out_dir = file.path(seurat_exp_cluster_dir,glue::glue('Top{Topn}_Marker_heatmap'))
        out_dir = seurat_exp_cluster_dir
        if (!file.exists(out_dir)){
                dir.create(out_dir,recursive = TRUE)
        }
        #write.table(data, file.path(out_dir,glue::glue("Top{Topn}_avgExpression_marker_heatmap_data.xls")),quote = F, row.names = F, sep = "\t")
        if (dim(collapsed_count)[2]>2) {
            df <- t(scale(t(collapsed_count)))
                scale_data = tibble::rownames_to_column(as.data.frame(df),"GeneID")
                write.table(scale_data, file.path(out_dir,glue::glue("Top{Topn}_avgExpression_marker_heatmap_scaledata.xls")),quote = F, row.names = F, sep = "\t")
        }else{
                df <- collapsed_count
        }
        #增加Annobar
        if(is.factor(metadata[,collapseby]) ){
            Anno_df = levels(metadata[,collapseby])
        }else{
            Anno_df = as.vector(colnames(collapsed_count))
        }
				#collapseby="celltype"
        Anno_df = factor(Anno_df,levels=Anno_df)
        col <- setNames(colors[1:length(Anno_df)], Anno_df)
                ncol = ceiling(length(unique(metadata[[collapseby]]))/10)
        Ha = HeatmapAnnotation(CellType=Anno_df,col=list(CellType =col ),simple_anno_size = unit(0.3, "cm"),show_annotation_name=F,annotation_legend_param =list(ncol=ncol))
        palette <-colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))(299)

        plot2 =Ha %v% Heatmap(df,
                row_names_gp = gpar(fontsize = 4),
                column_names_gp = gpar(fontsize = 8),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names= FALSE,
                col = palette,
                heatmap_legend_param = list(title = "Expression-level",title_position ="leftcenter-rot"),
                column_names_rot = 45,
                name=" ")
        ggheatmap <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot2, padding = grid::unit(c(2, 3, 2, 3), "mm"))))
        ggsave(plot=ggheatmap,filename=file.path(out_dir,glue::glue("cluster_top{Topn}_markers_avg_heatmap.pdf")),width=3.8+ncol*0.2,height=4)
        ggsave(plot=ggheatmap,filename=file.path(out_dir,glue::glue("cluster_top{Topn}_markers_avg_heatmap.png")),width=3.8+ncol*0.2,height=4,dpi=1000)
}


RunDiffexp <-function(object, test = "wilcox", contrast = NULL, fdr = NULL, 
    fc.thres = 2, pval.thres = 0.05, outputdir = NULL, logfc_threshold = logfc_threshold , ...)
{
    quiet_library("openxlsx")
    contrasts = unlist(strsplit(contrast, ":", perl = T))
    numerator_subset = unlist(strsplit(contrasts[2], ",", perl = T))
    denominator_subset = unlist(strsplit(contrasts[3], ",", perl = T))
		if( contrasts[1] !="sample.name"){
			cellmeta = Seurat::FetchData(object, vars = c("sample.name",contrasts[1])) %>% tibble::rownames_to_column(var = "barcode")
		}else{
			cellmeta = Seurat::FetchData(object, vars = contrasts[1]) %>% tibble::rownames_to_column(var = "barcode")
			
		}
    
    numerator = cellmeta %>% dplyr::filter(!!rlang::sym(contrasts[1]) %in%
        numerator_subset) %>% dplyr::pull(barcode)
    denominator = cellmeta %>% dplyr::filter(!!rlang::sym(contrasts[1]) %in%
        denominator_subset) %>% dplyr::pull(barcode)
    #如果numerator或者denominator 小于3，则不进行差异表达分析
    if (length(numerator) < 3 || length(denominator) < 3) {
        warning(paste0("细胞数小于3 ", contrast,
            sep = ""))
        return(0)
    }
    res = Seurat::FindMarkers(object, ident.1 = numerator_subset,
        ident.2 = denominator_subset, group.by = contrasts[1],
        min.pct = 0, test.use = test, only.pos = F, logfc.threshold = logfc_threshold ,...)
    if ("auc" %in% names(res)){res = res %>% dplyr::select(-auc)}
    res = res %>% tibble::rownames_to_column(var = "gene")
    if (nrow(res) == 0) {
        warning(paste0("NO Differential Genes Identified for ",
            contrast, sep = ""))
        return(0)
    }
		expr <- SeuratObject::GetAssayData(object, slot = "data")[res$gene, , drop=FALSE]
    numerator_means = Matrix::rowMeans(expr[res$gene, numerator])
    denominator_means = Matrix::rowMeans(expr[res$gene, denominator])
    res = res %>%  dplyr::mutate(
													FoldChange = 2^avg_log2FC,
													!!paste0(contrasts[1], "_mean") := numerator_means,
													!!paste0(contrasts[2], "_mean") := denominator_means
													)
		if (contrasts[1] != "sample.name") {
			
			sample_list <- split(cellmeta$barcode, cellmeta$sample.name)

			sample_means_list <- lapply(sample_list, function(cells) {
				Matrix::rowMeans(expr[, cells, drop=FALSE])
			})

			sample_means_df <- as.data.frame(sample_means_list)
			colnames(sample_means_df) = paste0(colnames(sample_means_df), "_mean")
			sample_means_df <- tibble::rownames_to_column(sample_means_df, "gene")

			res <- res %>%
				left_join(sample_means_df, by="gene")
		}
    
    if (!is.null(pval.thres)) {
        res_Significant = dplyr::filter(res, p_val < pval.thres,
            abs(avg_log2FC) > log2(fc.thres))
    }
    if (!is.null(fdr)) {
        res_Significant = dplyr::filter(res, p_val_adj < fdr, abs(avg_log2FC) >
            log2(fc.thres))
    }
    if (nrow(res_Significant) == 0) {
        warning(paste0("NO Significant Differential Genes Identified for ",
            contrast, sep = ""))
        return(0)
    }
    res_Significant[which(res_Significant$avg_log2FC > 0),
        "Regulation"] <- "Up"
    res_Significant[which(res_Significant$avg_log2FC < 0),
        "Regulation"] <- "Down"
    res[match(res_Significant$gene, res$gene), "Regulation"] <- res_Significant$Regulation
		res$Regulation[is.na(res$Regulation)] <- "NO_DEG"
    res =res %>% dplyr::select(gene,Regulation,everything()) 
    if (!is.null(outputdir)) {
        base_prefix = paste0(contrasts[1], "_", contrasts[2],
            "-vs-", contrasts[3])
				#xlsx的标签设定
				readme_title = createStyle(
					fontSize = 12, textDecoration = "bold", 
					halign = "center", fgFill = "#D9D9D9", 
					border = "bottom"
				)
        wb <- createWorkbook()
				sheet <- addWorksheet(wb,"all_diffexp_genes")
        writeData(wb, sheet, res, startCol = 1, startRow = 1, colNames = TRUE)
				addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:8, gridExpand = TRUE)
        if (!is.null(pval.thres)) {
		#colnames(res_Significant) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange","Regulation")
            # sheet <- addWorksheet(wb,paste0("diff-", "p_val-", pval.thres,
            #       "-FC-", fc.thres))
						# writeData(wb, sheet, res_Significant, startCol = 1, startRow = 1, colNames = TRUE)
						# addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:9, gridExpand = TRUE)
						# setColWidths(wb, sheet, cols = 1:9, widths = 20)
            stat <- matrix(c("Case", "Control", "Up_diff", "Down_diff",
                paste("Total_diff(", "p_val<", pval.thres, "&FoldChange>",
                  fc.thres, ")", sep = "")), ncol = 5)
        }
        if (!is.null(pval.thres)) {
		#colnames(res_Significant) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange","Regulation")
            # sheet <- addWorksheet(wb,paste0("diff-", "p_val-", pval.thres,
            #       "-FC-", fc.thres))
						# writeData(wb, sheet, res_Significant, startCol = 1, startRow = 1, colNames = TRUE)
						# addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:9, gridExpand = TRUE)
						# setColWidths(wb, sheet, cols = 1:9, widths = 20)
            stat <- matrix(c("Case", "Control", "Up_diff", "Down_diff",
                paste("Total_diff(", "p_val<", pval.thres, "&FoldChange>",
                  fc.thres, ")", sep = "")), ncol = 5)
        }
        if (!is.null(fdr)) {
            # sheet <- addWorksheet(wb,paste0("diff-", "p_val_adj-", fdr,
            #       "-FC-", fc.thres))
						# writeData(wb, sheet, res_Significant, startCol = 1, startRow = 1, colNames = TRUE)
						# addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:9, gridExpand = TRUE)
						# setColWidths(wb, sheet, cols = 1:9, widths = 20)
            # write.table(res_Significant, file.path(outputdir,
            #     paste0(base_prefix, "-diff-", "p_val_adj-", fdr, "-FC-",
            #       fc.thres, ".xls")), sep = "\t", col.names = TRUE,
            #     row.names = FALSE, quote = FALSE, na = "")
            stat <- matrix(c("Case", "Control", "Up_diff", "Down_diff",
                paste("Total_diff(", "p_val_adj<", fdr, "&FoldChange>",
                  fc.thres, ")", sep = "")), ncol = 5)
        }
				#readme页面设置
				sheet = addWorksheet(wb,"README")
				writeData(wb, sheet, t(as.data.frame(c("标题", "说明"))), startCol = 1, startRow = 1, colNames = FALSE)
				addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:2, gridExpand = TRUE)
				readme_list <- list(
					"gene"        = "基因名称",
					"p_val"    = "显著性p值",
					"avg_log2FC"       = "差异倍数的log2值",
					"pct.1"       = "基因在当前cluster中有表达的细胞比例",
					"pct.2"           = "基因在除当前cluster以外所有的cluster有表达的细胞比例",
					"p_val_adj" = "校正后的p值",
					"FoldChange" = "差异倍数",
					"Regulation" = "上调或下调",
					"*_mean" = "基因在该群细胞中的表达均值"
				)
				setColWidths(wb, sheet, cols = 2, widths = 80)
				readme_list_filter <- rownames_to_column(as.data.frame(t(as.data.frame(readme_list))),var = "tittle")
        writeData(wb, sheet, readme_list_filter, startCol = 1, startRow = 2, colNames = FALSE)
				saveWorkbook(wb, file = file.path(outputdir, paste0(base_prefix,"_diff_result.xlsx")), overwrite = TRUE)
        
				#计算差异基因数量
				up_num <- length(which(res_Significant$Regulation == "Up"))
        down_num <- length(which(res_Significant$Regulation == "Down"))
        total <- sum(up_num, down_num)
        stat1 <- c(paste(contrasts[2], "(", contrasts[1], ")",
            sep = ""), paste(contrasts[3], "(", contrasts[1],
            ")", sep = ""), up_num, down_num, total)
        
        #stat <- rbind(stat, stat1)
       	if (file.exists(file.path(outputdir, "diffexp_results_stat.xls"))) {
            stat <- rbind(stat1)
        }
        else {
            stat <- rbind(stat, stat1)
        }
        write.table(stat, file.path(outputdir, "diffexp_results_stat.xls"),
            append = file.exists(file.path(outputdir, "diffexp_results_stat.xls")),     
            quote = F, row.names = F, col.names = F, sep = "\t",
            na = "")
    }else{
        return(res_Significant)
    }
}

volcano_plot <- function(DEG, outputdir, base_prefix,foldchange,pvalue) {
	DEG_label = DEG$gene
  up = DEG %>% filter(Regulation =="Up") %>% arrange(desc(avg_log2FC ))  %>% top_n(5,avg_log2FC )
  down =  DEG %>% filter(Regulation =="Down") %>% arrange(desc(avg_log2FC ))  %>% top_n(-5,avg_log2FC )
  DEG$significant <- DEG$gene %in% c(as.character(up$gene),as.character(down$gene))
  DEG$label <-''
  DEG[DEG$significant,"label"] = as.character(DEG[DEG$significant,"gene"])
  ## 极值处理
	DEG$p_val[which(DEG$p_val < 5E-300 )] = 5E-300
	color_map = c("#54BCD1FF", "grey", "#FC6882FF")
	#c("up" = "#FC6882FF", "down" = "#54BCD1FF", "NO_DEG" = "grey")
	names=c("Down","NO_DEG","Up")
  color_lable = table(DEG$Regulation)
	color_lable=paste0(names(color_lable),' ',color_lable," genes")
	names(color_map) <- names
	## visualization
	p= ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val))) +
  # 绘制基础散点图，并根据 gene_type 对点的颜色进行分类，设置点的透明度 (alpha=0.6)，形状 (shape = 16)，大小 (size = 1)
  geom_point(aes(color = Regulation), alpha = 0.6, shape = 16, size = 2) +
  # 从 up_genes 数据框中绘制特定形状的散点图，填充颜色为红色，边框颜色为黑色，大小为 2
  geom_point(data = up, shape = 21, size = 3, fill = "#C70E7BFF", colour = "#C70E7BFF") +
  # 从 down_genes 数据框中绘制特定形状的散点图，填充颜色为钢蓝色，边框颜色为黑色，大小为 2
  geom_point(data = down, shape = 21, size = 3, fill = "#007BC3FF", colour = "#007BC3FF") +
  # 添加水平虚线，y 轴截距为 -log10(0.05),表示显著性阈值为 0.05
  geom_hline(yintercept = -log10(pvalue), linetype = "dashed") +
  # 添加垂直虚线，x 轴截距为 log2(0.5) 和 log2(2)，表示折叠变化范围为 0.5 到 2
  geom_vline(xintercept = c(-log2(foldchange),log2(foldchange)), linetype = "dashed") +
  # 在图中显示 sig_genes 数据框中基因符号的标签
  geom_label_repel(data = DEG[DEG$significant=="TRUE",], aes(label = gene), force = 2, nudge_y = 1) +
  # 设置 gene_type 对应的颜色映射
  scale_color_manual(values =color_map,
                     labels = color_lable) +
  # 设置 x 轴的刻度和范围
  scale_y_continuous(
  expand = expansion(mult = c(0, 0.2)),  # 下边界不留空隙，顶端多留20%
  limits = c(0, NA)                      # 强制从0开始，最大值自动
) +
  # 设置 x 轴和 y 轴的标签
  labs(x = "log2(fold change)", y = "-log10( P-value)", colour = "Expression change") +
  # 调整图例外观，将图例大小设为 5，位置设置为右上角
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_bw() + #  # 设置图的主题为白色背景
  # 设置图的主题样式，包括边框、网格线、背景等
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text( color = "black", size = 10),
        axis.text = element_text(color = "black", size = 9),
        legend.background = element_blank(),
        legend.title = element_text( color = "black", size = 10),
        legend.text = element_text( color = "black", size = 9),
        legend.spacing.x = unit(0, "cm"),
        legend.position = c(0.88, 0.89)  # 设置图例位置为右上角
  )
  ggsave(file.path(outputdir, paste0(base_prefix,"_volcano")), p, width = 8, dpi = 300,bg="white")
	# ggsave(file.path(outputdir, paste0(base_prefix,"_volcano.pdf")), p, width = 8)

}
downsample <- function(data_ob, limits_num, groupby) {
  #如果细胞数太多，采用随机抽样到limits_num
	#downsample = opt$downsample
	if(ncol(data_ob) > limits_num){
	# library(sampling) need install,instead with manual function
		ratio <- as.numeric(downsample) / ncol(data_ob)
		metadata_temp <- as.data.frame(data_ob@meta.data)
		cells_sample <- c()
		for (i in unique(data_ob@meta.data[[groupby]])) {
			cells_temp <- rownames(metadata_temp)[which(metadata_temp[[groupby]] == i)]
			#设置随机种子，确保结果可复原
			set.seed(2024)
			cells_temp_sample <- sample(cells_temp, ceiling(length(cells_temp) * ratio), replace = FALSE, prob = NULL)
			cells_sample <- append(cells_sample, cells_temp_sample)
		}
		data_ob <- subset(data_ob, cells = cells_sample)
		print("单组细胞数过多，已进行降采样。")
		print(paste0(c("基因数是：","细胞数是："),dim(data_ob)))
    gc()
	}
	return(data_ob)	
}




diff_heatmap = function(DEG,single_ob, colors,groupby,outputdir, base_prefix) {
	#groupby=diff_group[1]
	up = DEG %>% filter(Regulation =="Up") %>% arrange(desc(avg_log2FC ))  %>% top_n(20,avg_log2FC )
  down =  DEG %>% filter(Regulation =="Down") %>% arrange(desc(avg_log2FC ))  %>% top_n(-20,avg_log2FC )
	topn_markers = rbind(up, down)
	markers2vis4heatmap = unique(as.vector(topn_markers$gene))
	single_ob<-Seurat::ScaleData(single_ob,
                           features = markers2vis4heatmap, ## 更改为只针对marker基因进行scale,以减小内存压力
                           verbose = T )
   sc = ggplot2::scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                                mid = "white",
                                high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                                midpoint = 0,
                                na.value="white",
                                guide = "colourbar",
                                aesthetics = "fill")
  ggheat = Seurat::DoHeatmap( object = single_ob,
                        features = markers2vis4heatmap,
                        group.colors = colors,
                        group.by = groupby, group.bar = T, label = F, draw.lines = F) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6-log2(length(markers2vis4heatmap)/80),face = "bold"))+
      ggplot2::guides(fill = ggplot2::guide_colorbar( title.position = "top", order = 1)  ) +#, color = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 1)))+
      sc+ scale_color_discrete(name = "Identity", labels = levels(single_ob@meta.data[[groupby]]))
  ggsave(file.path(outputdir, paste0(base_prefix,"_Top20_diff_gene_heatmap")), ggheat, width =7, height = 7, bg="white")
}

box_fecth = function(DEG,single_ob,groupby,outputdir, base_prefix) {
	quiet_library("tidyr")
	up = DEG %>% filter(Regulation =="Up") %>% arrange(desc(avg_log2FC ))  %>% top_n(5,avg_log2FC )
  down =  DEG %>% filter(Regulation =="Down") %>% arrange(desc(avg_log2FC ))  %>% top_n(-5,avg_log2FC )
	topn_markers = rbind(up, down)
  genes <- topn_markers$gene
	# 提取表达矩阵（slot = "data" 表示归一化后的表达量）
	expr_mat <- GetAssayData(single_ob, assay = "RNA", slot = "data")[genes, , drop = FALSE]
	expr_df <- as.data.frame(t(as.matrix(expr_mat)))
	# 合并进 meta.data
	single_ob@meta.data <- cbind(single_ob@meta.data, expr_df)
  metadata <- single_ob@meta.data
  metadata <- metadata[, c(groupby, genes)]
	
	df_long <- metadata %>%
  pivot_longer(
    cols = -!!rlang::sym(groupby),         # 除 group 外所有列
    names_to = "gene", 
    values_to = "expression"
  )
	df_long$gene = factor(df_long$gene, levels = genes)
	colnames(df_long)[1] = "group"
	p <- ggplot(df_long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol =5) +  # 每个基因一个子图
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  ) +
  labs(y = "Expression", x = "Group")
  ggplot2::ggsave(file.path(outputdir,paste0(base_prefix,"_Top5_diff_gene_boxplot.pdf")),p,width=14)
	ggplot2::ggsave(file.path(outputdir,paste0(base_prefix,"_Top5_diff_gene_boxplot.png")),p,width=14)
}
