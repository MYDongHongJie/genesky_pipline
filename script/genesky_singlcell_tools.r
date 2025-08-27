#!/usr/bin/env Rscript
# Autor：donghongjie
# PSN singlecell Analysis tool
# Version 1.0
# Date 2025-07-21
# =======================GLOBAL Installation department loading=================


if (!requireNamespace('argparse', quietly = TRUE)) {
  stop("Please conda activate")
}
quiet_library <- function(pkg) {
  suppressWarnings(suppressPackageStartupMessages(library(pkg, character.only = TRUE)))
}
print("正在加载分析环境")
quiet_library("Seurat") 
quiet_library("argparse") 
quiet_library("magrittr") 
quiet_library("rlang") 
quiet_library("ggplot2") 
quiet_library("harmony")
quiet_library("cowplot")
quiet_library("stringr")
#quiet_library("SingleCellExperiment"))
quiet_library("dplyr")
quiet_library("parallel")
quiet_library("tibble")
quiet_library("pbapply")
quiet_library("ggrepel")
quiet_library("tidyr")
#颜色代码
source("/home/donghj/scRNA_snakemake/extra_dir/color.r") #TODO:这个路径需要修改
source("/home/donghj/scRNA_snakemake/extra_dir/public_func.r")
# 额外函数
print("加载完成")
# ======================= COMMAND LINE PARAMETERS SETTING =======================
# ======================= GLOBAL  parameters =================
parser = ArgumentParser(description = "GLOBAL  parameters",
                        usage = "%(prog)s [global options]" )
parser$add_argument("-i", "--input", type = "character",
             help = "The input exprssion matrix in several possible format.")
parser$add_argument("-f", "--informat", type = "character", default = 'rds',
             help = "The format of input file, maybe:h5seurat,(seurat)rds,(sce)rds, loom.[default: %(default)s]")
parser$add_argument("-o", "--output", type = "character", default = "./",
             help = "the output directory of results."  )
parser$add_argument("-j", "--ncores", type="integer", default = 10,
             help="the number of CPUs used to improve the performace.[default: %(default)s]")
parser$add_argument("--prefix", type = "character", default = "data_ob",
             help = "the prefix of output file without file extension.[default %(default)s]")
parser$add_argument("--assay", type = "character", default = "RNA",
             help = "Assay used for analysis,maybe RNA or SCT")
parser$add_argument("--subset", type = "character", default = NULL,
             help = "The conditional expression to subset cells used for subtyping.[default: %(default)s]")
parser$add_argument("--update", default= TRUE, type="logical",
             help="whether update the data in the object on disk using the newly produced results. Set this to FALSE when you use this script for subclustering! Only availiable for h5seurat input.[default TRUE]")
parser$add_argument("--report", default= TRUE, type="logical",
             help="Is it the result of the report outfile.[default TRUE]")
subparsers = parser$add_subparsers(help = "subcommands:")

# ======================= create  parameters =================
# create the parsers for subcommand create
sub_create = subparsers$add_parser("create",
     help = "make the specified format data object out of the supported input data sources.")
sub_create$add_argument("-s", "--source", type = "character", default = "02_cellranger",
     help = "THE load of cellrang martrix")
sub_create$add_argument("-m", "--metadata", type = "character", default = "sample_info.txt",
     help = "sample metadata ")
sub_create$add_argument("-t", "--type", type = "character", default ="SCRNA",
     help = "Sequencing data type,maybe VDJ SCRNA other")
sub_create$add_argument("--genecolumn", type = "double", default =2,
     help = "If features.tsv.gz only has one column of data, change this parameter to 1")
sub_create$add_argument("--minCells", type = "double", default = 3,
     help = "Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff")
sub_create$add_argument("--minfeatures", type = "double", default = 200,
     help = "Include cells where at least this many features are detected")

# ======================= QC  parameters =================
## create the parsers for subcommand QC
sub_filterQC=subparsers$add_parser("QC",help = "Filter low-quality cells.")
sub_filterQC$add_argument("-s","--species",type="character",default='human',
																	help = "The species of the reference to use. ")
sub_filterQC$add_argument("-g","--gmt",type="character",default=NULL,help = "List of genes that need to be filtered,The column name will be used as the name for filtering, corresponding to the -- filter parameter,eg:percent.mt")
sub_filterQC$add_argument("-f","--filter",type="character",default="percent.mt,nFeature_RNA,nCount_RNA",help = "Filtered name,split by',',eg:percent.mt")
sub_filterQC$add_argument("-m","--high",type="character",default="20,7000,50000",help = "Filter cells that exceed this value,split by ',' and Keep the length consistent with -- filter")
sub_filterQC$add_argument("-l","--low",type="character",default="0,400,0",help = "Filter cells that below this value,split by ',' and Keep the length consistent with -- filter")
sub_filterQC$add_argument("-d","--doublet",type="character",default="FALSE",help = "Whether to filter doublet cells,default FALSE")
sub_filterQC$add_argument("--chem_version",type="character",default="other",help = "The version of chemistry used,choices can be: v4,other")
sub_filterQC$add_argument("--pattle",type="character",default="customecol2_light",help = "Color palette")

# ======================= Clustering  parameters =================
## create the parsers for subcommand Clustering
sub_Clustering=subparsers$add_parser("Clustering",help = "Clustering cells.")
sub_Clustering$add_argument("-g","--gather",type="character",default="harmony",help="The method of dimensionality reduction clustering")
sub_Clustering$add_argument("-c","--cycle",type="logical",default=TRUE,help = "Whether to filter Cycle gene,default TRUE")
sub_Clustering$add_argument("-n","--normmeth",type="character",default = "LogNormalize",help="the method to normalize the raw counts. Choices can be:  LogNormalize, CR, CLR. For feature-barcoding/cite-seq, CLR is recommended. For scRNA-seq, sctransform is recommended.")
sub_Clustering$add_argument("-b","--batchid",type="character",default=NULL,help="If a sample batch needs to be specified, a two column file can be passed in, with the sample name in the first column. The second column contains batch information")
sub_Clustering$add_argument("-r","--rely",type="character",default="sample.name",help="If batch information(-b optional) is not provided, which information will be used to batch,eg:sample.name,group")
sub_Clustering$add_argument("-s","--res",type="double",default=0.6,help="resolution ratio")
sub_Clustering$add_argument("-k","--knn",type="integer",default=30,help="The number of nearest neighbors to be used for graph construction")
sub_Clustering$add_argument("-d","--dim",type="integer",default=NULL,help="Dimensionality reduction,Automatically calculated or artificially specified")
#sub_Clustering$add_argument("-t","--thread",type="integer",default=1,help="The number of nearest neighbors to be used for graph construction")
sub_Clustering$add_argument("--regression",type="character",default=NULL,help="The parameters to be included in the regression analysis must all be in the metadata column,Multiple options available, separated by commas[optional]")
sub_Clustering$add_argument("-f","--scFeature",type="character",default="nCount_RNA,nFeature_RNA,percent.mt",help="Select data features for plotting")
sub_Clustering$add_argument("--pattle",type="character",default="customecol2_light",help = "Color palette")
#======================= Visualize  parameters =======================
## create the parsers for subcommand Clustering
sub_Visualize=subparsers$add_parser("Visualize",help = "Clustering Visualize.")
sub_Visualize$add_argument("-g","--groupby",type="character",default="seurat_clusters",help="Select the RDS tag for visualization")
sub_Visualize$add_argument("-s","--splitby",type="character",default="sample.name,group",help="Select the RDS tag for split visualization")
sub_Visualize$add_argument("-z","--pointsize",type="double",default=NULL,help="dimplot pointsize")
sub_Visualize$add_argument("--pattle",type="character",default="customecol2_light",help = "Color palette")
sub_Visualize$add_argument("-t","--type",type="character",default="all",help = "Dimplot,summary or all,splitby ','")
sub_Visualize$add_argument("-r","--reductions",type="character",default="tsne,umap",help = "tsne or umap,splitby ','")

#======================= Marker  parameters =======================
sub_Marker=subparsers$add_parser("marker",help = "Find Marker")
sub_Marker$add_argument("-m","--min_pct",type="double",default="0.1",help="Only retain genes with a proportion greater than this parameter in two groups of cells [Findallmarker parameter]")
sub_Marker$add_argument("-l","--logfcthreshold",type="double",default="0.25",help="Perform restriction testing on genes with an average difference of at least x times (logarithmic scale) between two groups of cells. The default value is 0.25 for log threshold acceleration function, but weaker signals may be missed [Findallmarker parameter]")
sub_Marker$add_argument("-t","--testuse",type="character",default="wilcox",help="Test method, including DESeq2, negbinom, Poisson,bimod,t,roc, LR,MAST,and wilcox_limma [Findallmarker parameter]")
sub_Marker$add_argument("-k","--slot",type="character",default="data",help="Note that if the test uses negbinom, poisson, or DESeq2, the slot will be set to counts [Findallmarker parameter]")
sub_Marker$add_argument("--plot",type="logical",default=TRUE,help="whether to plot the data")
sub_Marker$add_argument("-g","--groupby",type="character",default="seurat_clusters",help="Still the cell label of findamarker ")
sub_Marker$add_argument("-n","--topn",type="integer",default="10",help="Number of top markers selected")
sub_Marker$add_argument("-a","--avetopn",type="integer",default="30",help="The top number of marker mean heatmaps")
#sub_Marker$add_argument("-s","--species",type="character",default='human',help = "The species of the reference to use. ")
sub_Marker$add_argument("--pattle",type="character",default="customecol2_light",help = "Color palette")
#----------------------------------------------findmarker的检验方法参数介绍-------------------------------------------------------
# "wilcox" : Identifies differentially expressed genes
              # between two groups of cells using a Wilcoxon Rank Sum
              # test (default); will use a fast implementation by Presto
              # if installed
#"wilcox_limma" : Identifies differentially expressed
        
#"bimod" : Likelihood-ratio test for single cell gene
            #expression, (McDavid et al., Bioinformatics, 2013)							
#"roc" : Identifies 'markers' of gene expression using ROC
              # analysis.  For each gene, evaluates (using AUC) a
              # classifier built on that gene alone, to classify between
              # two groups of cells. An AUC value of 1 means that
              # expression values for this gene alone can perfectly
              # classify the two groupings (i.e. Each of the cells in
              # cells.1 exhibit a higher level than each of the cells in
              # cells.2). An AUC value of 0 also means there is perfect
              # classification, but in the other direction. A value of
              # 0.5 implies that the gene has no predictive power to
              # classify the two groups. Returns a 'predictive power'
              # (abs(AUC-0.5) * 2) ranked matrix of putative
              # differentially expressed genes.
#"t" : Identify differentially expressed genes between two
              #groups of cells using the Student's t-test.
# "negbinom" : Identifies differentially expressed genes
#               between two groups of cells using a negative binomial
#               generalized linear model.  Use only for UMI-based
#               datasets
# "poisson" : Identifies differentially expressed genes
#               between two groups of cells using a poisson generalized
#               linear model.  Use only for UMI-based datasets
# "LR" : Uses a logistic regression framework to determine
#               differentially expressed genes. Constructs a logistic
#               regression model predicting group membership based on
#               each feature individually and compares this to a null
#               model with a likelihood ratio test.
# "MAST" : Identifies differentially expressed genes
#               between two groups of cells using a hurdle model tailored
#               to scRNA-seq data. Utilizes the MAST package to run the
#               DE testing.
#  "DESeq2" : Identifies differentially expressed genes
#               between two groups of cells based on a model using DESeq2
#               which uses a negative binomial distribution (Love et al,
#               Genome Biology, 2014).This test does not support
#               pre-filtering of genes based on average difference (or
#               percent detection rate) between cell groups. However,
#               genes may be pre-filtered based on their minimum
#               detection rate (min.pct) across both cell groups. To use
#               this method, please install DESeq2, using the
#               instructions at
#               https://bioconductor.org/packages/release/bioc/html/DESeq2.html														

#======================= Enrichment  parameters =======================
sub_Enrichment=subparsers$add_parser("Enrichments",help = "gene enrichment analysis")
sub_Enrichment$add_argument("-f","--file",type="character",default=NULL,help="marker gene file or diffexp file")
sub_Enrichment$add_argument("-s","--species",type="character",default="human",help="The species of the reference to use. ")
sub_Enrichment$add_argument("-n","--topn",type="integer",default="10",help="Number of KEGG pathway to map")
sub_Enrichment$add_argument("-r","--rankby",type="character",default="pvalue",help="Sort top based on this column")
sub_Enrichment$add_argument("-c","--CirchordColor",type="character",default=NULL,help="Color palette of the circhord")
#======================= GSEA  parameters =======================
sub_GSEA=subparsers$add_parser("GSEA",help = "gene enrichment analysis")
sub_GSEA$add_argument("-f","--file",type="character",default=NULL,help="marker gene file or diffexp file")
sub_GSEA$add_argument("-s","--species",type="character",default="human",help="The species of the reference to use. ")
#======================= singleR  parameters =======================
sub_singleR=subparsers$add_parser("singleR",help = "singleR analysis")
sub_singleR$add_argument("-r","--ref",type="character",default=NULL,help="marker gene file or diffexp file")
sub_singleR$add_argument("-g","--groupby",type="character",default='celltype',help="ref celltype")
sub_singleR$add_argument("--downsample",type="integer",default=30000,help="downsample")
sub_singleR$add_argument("--pattle",type="character",default="customecol2_light",help = "Color palette")
#======================= diff  parameters =======================
sub_diffexp=subparsers$add_parser("diffexp",help = "diffexp analysis")
sub_diffexp$add_argument("-e","--pvalue",type="double",default=0.05,help="选择显著性差异基因的pvalue阈值")
sub_diffexp$add_argument("-f","--fdr",type="double",default=NULL,help="选择显著性差异基因的矫正后的pvalue阈值，与pvalue互斥")
#sub_diffexp$add_argument("-g","--groupby",type="character",default='celltype',help="ref celltype")
sub_diffexp$add_argument("-g","--logfc",type="double",default=0.25,help="筛选差异基因的阈值，如果想要全部差异基因就改为0")
sub_diffexp$add_argument("--FC",type="integer",default=2,help="选择显著性差异基因的FC值的阈值")
sub_diffexp$add_argument("--pattle",type="character",default="customecol2_light",help = "Color palette")
sub_diffexp$add_argument("--test",type="character",default="wilcox",help = "计算差异基因的方法")
sub_diffexp$add_argument("-s","--splitby",type="character",default="celltype",help = "选择某类细胞进行拆分去计算差异基因")
sub_diffexp$add_argument("-c","--contrasts",type="character",default=NULL,help = "差异分组比较信息，eg：group:groupA:groupB")

#========================loupeR parameters =========================
sub_loupeR=subparsers$add_parser("loupeR",help = "rds to loupeR analysis")
sub_loupeR$add_argument("-n","--cloupe_name",type="character",default="dimensional_reduction",help="loupe name")

#========================Celltyping parameters =========================
sub_Celltyping=subparsers$add_parser("Celltyping",help = "add celltype to rds")
sub_Celltyping$add_argument("-s","--re_sample_diff",type="character",default="re_sample_diff.txt",help="rename sample and define group file")
sub_Celltyping$add_argument("-c","--cluster",type="character",default="anno.txt",help="anno file")
sub_Celltyping$add_argument("--pattle",type="character",default="customecol2_light",help = "Color palette")
opt = parser$parse_args()
#------------------------------------- GLOBAL Process ------------------------------------------
output = opt$output
if(!file.exists(output)){
	dir.create(output,recursive = TRUE)
}
invisible(gc(full = T, verbose = F))
options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
future::plan("multicore", workers = min(future::availableCores(), opt$ncores)) # parallization using specified CPUs start from here
# ------------------------------------ SUB Create ----------------------------------------------
args<-commandArgs(TRUE)
if ("create" %in% args ){
	quiet_library("future.apply")
	if (!is.null(opt$metadata)){
    sample = read.delim(opt$metadata,header=F)
		sample_list=c(sample[[2]])
    sample_group=c(sample[[1]])
    ifnb_list=list()
	}else{
		print('please input the metadata of sample')
		quit()
	}
	if (is.null(opt$source)){
		print('please input the load of cellranger matrix')
		quit()
	}

	if (opt$type =="SCRNA"){
		data_structure = 'outs/filtered_feature_bc_matrix'
	}else if (opt$type == 'VDJ') {
		data_structure = "count/sample_filtered_feature_bc_matrix"
	}else{
		data_structure='filtered_feature_bc_matrix'
	}

	
	ifnb_list <- lapply(sample_list, function(x) {
		group_X = sample[sample[,2] == x, 1]
		if (opt$type == "other"){
			datadir=file.path(opt$source,x)
		}else{
			datadir=file.path(opt$source,x,data_structure)
		}
		single_data<-Read10X(datadir,gene.column =opt$genecolumn)
		single_ob<-CreateSeuratObject(counts = single_data, project = x, min.cells = opt$minCells, min.features = opt$minfeatures)
		single_ob[['sample.name']] = x
		single_ob[['group']] = group_X
		return(single_ob)
	})

	if(length(sample_list)>1){
        single_ob=merge(ifnb_list[[1]],ifnb_list[2:length(ifnb_list)],add.cell.ids=sample_list)}else{
        single_ob=ifnb_list[[1]]
  }
	saveRDS(single_ob,file=file.path(output,paste0(opt$prefix,".rds")))
	quit()
}

# ------------------------------------ SUB   filterQC ----------------------------------------------
if ( "QC" %in% args ){
	options(warn = -1)
	#设置一下配色
	colors = discrete_palette[[opt$pattle]]

  if(is.null(opt$input)){
		print('No input file or rds')
		quit()
	}else {
		single_ob = readRDS(opt$input)
		if (grepl('5',single_ob@version)){
			 single_ob <- JoinLayers(single_ob)
		}
	}
	if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
  }
	###解析输入的过滤参数
	filter = unlist(strsplit(opt$filter,","))
	nfilter = length(filter)
	fiter_max =  unlist(strsplit(opt$high,","))
	nfiter_max = length(fiter_max)
	filter_low = unlist(strsplit(opt$low,","))
	nfilter_low = length(filter_low)

  # if(nfilter!=nfiter_max & nfilter!=nfilter_low){
  #   print('The number of filter parameters is not consistent')
  #   quit()
	# }

 if(nfilter!=nfiter_max | nfilter!=nfilter_low){
    print('The number of filter parameters is not consistent')
    quit()
  }
	filter_list = setNames(Map(c, filter_low, fiter_max), filter)
 	#var_parameter = filter
	gmt_filter=c()
	if(!is.null(opt$gmt)){
		file_path <- opt$gmt
		#
		lines <- readLines(file_path)
		header <- strsplit(lines[1], "\t")[[1]]
		data <- vector("list", length(header))
		names(data) <- header

		# 遍历剩余的行，拆分并添加到相应的列
		for (line in lines[-1]) {
			values <- strsplit(line, "\t")[[1]]
			for (i in seq_along(values)) {
				if (values[i] != "") {
					data[[i]] <- c(data[[i]], values[i])
				}
			}
		}
		gmt_filter = intersect(filter,names(data))
		for (filter_name in gmt_filter){
			features = CaseMatch(data[[filter_name]],rownames(single_ob))
			single_ob[[filter_name]] = PercentageFeatureSet(single_ob, features = unname(features))
    }
		
  }
	if (!"percent.mt" %in% gmt_filter & "percent.mt" %in% filter){
		quiet_library("rjson")
		quiet_library("jsonlite")
		
		json_load <- "/home/donghj/scRNA_snakemake/extra_dir/ref.json" #TODO:替换成正式环境的路径
		json_dic<-read_json(json_load)
		mtgene_file <- json_dic[[opt$species]]$mito
		if(file.exists(mtgene_file)){
			mtdf <- read.table(mtgene_file)
			gene = CaseMatch( mtdf$V1,rownames(single_ob))
			single_ob[["percent.mt"]] <- PercentageFeatureSet(single_ob, features = gene)
		}else{
			print("提供的表格无线粒体信息，基因组数据也无线粒体信息，请核查")
      quit()
    }
  }
  
  #绘制原始的数据结果图片和图表
	cell_num_raw=table(single_ob$sample.name)
	cellqc_raw_list=list()
	cellqc_filter_list = list()
	for (i in filter){
			cellqc_raw_list[[i]] = VlnPlot(single_ob, features = i, ncol = 1,pt.size = 0,group.by="sample.name",cols=colors)+
			geom_boxplot(width=.2,col="black",fill="white") +theme(plot.margin = margin(t = 0, r = 0, b = 0,l = 0,unit = "cm"))
	}
	cellqc.raw <- do.call(plot_grid, c(plot = lapply(cellqc_raw_list, function(p) p + theme(legend.position = 'none')), ncol = length(filter)))
	#过滤
	print("开始按照规定的标准过滤阈值,过滤阈值如下:")
	if (nfilter != 0) {
    for (i in names(filter_list)) {
    # 获取当前过滤条件的上下限值
      lower_bound <- as.numeric(filter_list[[i]][1])
      upper_bound <- as.numeric(filter_list[[i]][2])
    # 构建过滤表达式
      filter_expr <- sprintf("%s > %f & %s < %f", i, lower_bound, i, upper_bound)
    # 打印过滤表达式用于调试
		  
      print(paste("Applying filter:", filter_expr))
    # 直接使用过滤表达式进行子集化
      subset_condition1 <- eval(parse(text = filter_expr), envir = single_ob@meta.data)
      single_ob <- single_ob[,subset_condition1]
      cell_num = table(single_ob$sample.name)
      cell_num_raw = rbind(cell_num_raw,cell_num)
		
  	}
	}
	filter_data = as.data.frame(t(cell_num_raw))
	colnames(filter_data) = c("raw",paste0(filter,"_after_filter"))
	if(opt$doublet == TRUE){
		print("开始去除双细胞")  
		quiet_library("DoubletFinder")
		source("/home/donghj/scRNA_snakemake/extra_dir/doubletfinder.r") #TODO:替换成正式环境的路径
#------------------------
    get_doublet_rates_table <- function(chem_version) {
			if (chem_version == "v4") {
				data.frame(
					upper_bound = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 
													10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 
													18000, 19000, 20000, Inf),
					rate = c(0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.036,
									0.040, 0.044, 0.048, 0.052, 0.056, 0.060, 0.064, 0.068, 0.072,
									0.076, 0.080, 0.1)
				)
			} else {
				data.frame(
					upper_bound = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 
													10000, Inf),
					rate = c(0.008, 0.015, 0.023, 0.030, 0.038, 0.046, 0.053, 0.061, 0.068,
									0.080, 0.1)
				)
			}
		}
		get_doublet_rate <- function(cell_count, rate_table) {
  		rate_table$rate[findInterval(cell_count, c(0, rate_table$upper_bound))]
		}
			# 主处理逻辑
		if (opt$chem_version %in% c("v4", "other")) {
			# 设置并行
			future::plan("multicore", workers = min(future::availableCores(), 2))
			rate_table <- get_doublet_rates_table(opt$chem_version)
			# 获取当前版本对应的双峰率配置
			ifnb.list =  SplitObject(single_ob, split.by = "sample.name")
			# 批量处理所有数据
			ifnb.list <- lapply(ifnb.list, function(obj) {
				cell_count <- length(Seurat::Cells(obj))
				dbl_rate <- get_doublet_rate(cell_count, rate_table)
				# 运行双峰检测
				obj <- RunDoubletFinder(obj, doublet.rate = dbl_rate)
			})
		}

#------------
		predicted_doublets= unlist(lapply(ifnb.list,FUN=function(x) x = x$predicted_doublets))
		rm(ifnb.list);gc()
		predicted_doublets=as.data.frame(predicted_doublets)
		rownames(predicted_doublets) = gsub("^[^.]*\\.", "",rownames(predicted_doublets))
		single_ob = Seurat::AddMetaData(single_ob, predicted_doublets)
		single_ob = subset(single_ob, subset = predicted_doublets == FALSE)
		doublet =as.data.frame((table(single_ob$sample.name)))
		rownames(doublet) <- doublet[,1]
		filter_data=cbind(filter_data,doublet[,2])
		colnames(filter_data)[ncol(filter_data)] = "doublet_after_filter"
		print("双细胞计算完成")
	}
	
	total <- colSums(filter_data)
	filter_data = rbind(filter_data,total)
	rownames(filter_data)[nrow(filter_data)] <- "total"
	filter_data$Percentage <- sprintf("%.2f%%",(filter_data[, ncol(filter_data)] / filter_data[, 1]) * 100)
	filter_data <- rownames_to_column(filter_data, var = "Sample")
	quiet_library("openxlsx")
	wb <- createWorkbook()
	sheet <- addWorksheet(wb,"Cell Quality Summary")
	writeData(wb, sheet, filter_data, startCol = 1, startRow = 1, colNames = TRUE)
	sheet = addWorksheet(wb,"README")
	readme_title = createStyle(
      fontSize = 12, textDecoration = "bold", 
      halign = "center", fgFill = "#D9D9D9", 
      border = "bottom"
    )
	writeData(wb, sheet, t(as.data.frame(c("标题", "说明"))), startCol = 1, startRow = 1, colNames = FALSE)
	addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:2, gridExpand = TRUE)
	readme_list <- list(
    "Sample"        = "样本名称",
    "raw"      = "原始细胞数量",
    "nFeature_RNA_after_filter"    = "nFeatrure_RNA过滤后细胞数量",
    "nCount_RNA_after_filter"       = "nCount_RNA过滤后细胞数量",
    "doublet_after_filter"       = "双细胞过滤后细胞数量",
    "Percentage"           = "过滤后的细胞占原始细胞数量占比",
		"percent.mt_after_filter" = "线粒体基因过滤后细胞数量"
  )
	setColWidths(wb, sheet, cols =1:2, widths = 40)
	readme_list_filter = readme_list[names(readme_list) %in% colnames(filter_data)]
	readme_list_filter <- rownames_to_column(as.data.frame(t(as.data.frame(readme_list_filter))),var = "tittle")
	#write.table(filter_data,file.path(output,"QC_cell_counts_per_step.xls"),quote = F,row.names = T,col.names = NA,sep='\t')
	writeData(wb, sheet, readme_list_filter, startCol = 1, startRow = 2, colNames = FALSE)
	saveWorkbook(wb, file = file.path(output, "QC_cell_counts_per_step.xlsx"), overwrite = TRUE)
	#绘制过滤后的结果
	for (i in filter){
			cellqc_filter_list[[i]] = VlnPlot(single_ob, features = i, ncol = 1,pt.size = 0,group.by="sample.name",cols=colors)+
			geom_boxplot(width=.2,col="black",fill="white") +theme(plot.margin = margin(t = 0, r = 0, b = 0,l = 0,unit = "cm"))
	}
	cellqc.filter <- do.call(plot_grid, c(plot = lapply(cellqc_filter_list, function(p) p + theme(legend.position = 'none')), ncol = length(filter)))
	cellqc.raw/cellqc.filter
  width <- 8+ceiling(length(unique(single_ob$sample.name))/3)+ceiling(length(filter)/1.5)
  ggsave(file.path(output,"QC_metrics_violin_before_vs_after_filtering"),width = width,height = 7,bg='white')
  #ggsave(file.path(output,"QC_metrics_violin_before_vs_after_filtering.pdf"),width = width,height = 7,bg='white')

	saveRDS(single_ob,file.path(output,paste0(opt$prefix,".rds")))
	quit()
}

# ------------------------------------ SUB Clustering ----------------------------------------------
if ( "Clustering" %in% args ){
	print("开始进行降维聚类")
	start_time = Sys.time()
	colors = discrete_palette[[opt$pattle]]
  if(is.null(opt$input)){
		print('No input file or rds')
		quit()
	}else {
		single_ob = readRDS(opt$input)
		if (grepl('5',single_ob@version)){
			 single_ob <- JoinLayers(single_ob)
		}
	}
	if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
  }
	DefaultAssay(single_ob) = opt$assay
	gather = opt$gather
	dim = opt$dim
	#细胞周期分析
	#在某些情况下，我们发现这会对下游的分析产生一定的负面影响，尤其是在细胞分化过程中(如鼠类造血
	# 过程)。在此过程中干细胞处于静止状态，而分化的细胞正在增殖。在这种情况下，清除所有细胞周期效
	# 应也会使干细胞和祖细胞之间的区别模糊。作为替代方案，建议可以逐步消除G2M和S期评分
	# 之间的差异。这意味着将保持非周期细胞和周期细胞的组分差异，但是增殖细胞之间的细胞周
	# 期阶段的差异将从数据中去除。
	if ( opt$cycle ){
		  #library(dittoSeq)
			single_ob <- NormalizeData(single_ob,normalization.method=opt$normmeth,scale.factor = 10000,margin = 1,verbose = TRUE)
      s.genes = Seurat::CaseMatch(search = Seurat::cc.genes$s.genes, match = rownames(single_ob) )
      g2m.genes = Seurat::CaseMatch(search = Seurat::cc.genes$g2m.genes, match = rownames(single_ob) )
      single_ob = Seurat::CellCycleScoring(single_ob, s.features = s.genes, g2m.features = g2m.genes,set.ident = F)
      cycleS = Seurat::FetchData(single_ob, vars = c("S.Score", "G2M.Score") )
      single_ob = AddMetaData(single_ob, col.name = "CC.Difference",
                           metadata = cycleS[["S.Score"]] - cycleS[["G2M.Score"]])
			
			
  }

	if (!is.null(opt$batchid) ){
		batchid_data = read.delim(opt$batchid)
		if (all((unique(single_ob$sample.name %in% batchid_data$sample.name)))){
			if("batchid" %in% colnames(single_ob@meta.data)){
				single_ob$batchid =NULL
			}
			single_ob@meta.data$Barcode = rownames(single_ob@meta.data)
			single_ob@meta.data <- single_ob@meta.data %>%
  		left_join(batchid_data, by = "sample.name") 
			rownames(single_ob@meta.data) =single_ob@meta.data$Barcode
			print(head(single_ob@meta.data))
		}else{
			print("样本名与rds的样本名不匹配")
      quit()
		}
		rely='batchid'
	}else if(is.null(opt$rely)){
		 print("没有提供去批次的指标，默认使用样本作为批次")
		 rely ="sample.name"
	}else{
		rely =opt$rely
	}


##判断加入回归分析的参数
	if(!is.null(opt$regression)){
		regression =unlist(strsplit(opt$regression,","))
	}else{
		regression=NULL
	}

	#去批次
	if (gather == 'CCA'){
		message("Started CCA integration analysis,next DefaultAssay use integrated")
		ifnb.list <-SplitObject(single_ob,split.by = rely)
		ifnb.list <- lapply(ifnb.list, function(x) SCTransform(x, vars.to.regress = regression))
		features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)
		ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
		single_ob <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",anchor.features = features)
		single_ob <- IntegrateData(anchorset = single_ob, normalization.method = "SCT")
		DefaultAssay(single_ob) <- "integrated"
		if(is.null(dim)){dim=30}
	}else{
		if (opt$cycle ==FALSE){single_ob <- NormalizeData(single_ob,normalization.method=opt$normmeth,scale.factor = 10000,margin = 1,verbose = TRUE)}
		single_ob <- Seurat::FindVariableFeatures(single_ob, loess.span = 0.3,
                              clip.max = "auto", mean.function = "FastExpMean",
                              dispersion.function = "FastLogVMR", num.bin = 20,
                              nfeature = 2000, binning.method = "equal_width" )
		single_ob <- ScaleData(single_ob, vars.to.regress = regression,features = VariableFeatures(single_ob))#只针对高变基因进行scale，降低存储和运行时间，只有绘制热图的时候需要再对绘图基因进行scale 
		#https://www.jianshu.com/p/c4fcc51be2f4 回归方法的说明
	#PCA的时候还是要选择高变基因，否则会引入噪声。低丰度，变化低的基因
	#updata:20250722 新增自动判断最佳的npcs值,且pca降维只会对高变基因进行主成分分析
    #进行第一步的降维，目前只使用pca，ica看情况添加
		message("Start the first step of dimensionality reduction, the dimensionality reduction method is pca, and the dim value is set to 30")
		single_ob=single_ob %>% RunPCA(verbose = FALSE,npcs =30,features =VariableFeatures(single_ob))
    
    #通过肘形图的落差最大的拐点确定最佳的dim值，也可以自己指定dim值
		if(is.null(dim)){
			dim = FindElbow(single_ob,reduction="pca")
		}
		message("Start the sceond step of dimensionality reduction, the dimensionality reduction method is ", gather, " and the dim value is set to ",dim)
		if (gather == 'harmony') {
		 single_ob <- single_ob %>% RunHarmony( group.by.vars = rely,dims.use=1:dim)
		 }else if(gather == "pca"){
			single_ob <- single_ob %>% RunPCA(verbose = FALSE,npcs =dim)
		} 
	}
		
	#降维聚类
	reduction = ifelse(gather =="pca" ,'pca',gather)
	single_ob <- RunUMAP(single_ob, reduction = reduction, dims = 1:dim)
	single_ob <- RunTSNE(single_ob, reduction = reduction, dims = 1:dim,check_duplicates = FALSE)
	single_ob <- FindNeighbors(single_ob, reduction = reduction, dims = 1:dim)
	single_ob <- FindClusters(single_ob,resolution = opt$res)
	
	#保存数据
	saveRDS(single_ob,file.path(output,paste0(opt$prefix,".rds")))
	if (opt$report == TRUE){
		preprocess=file.path(output,'pre_clustering')
		if(!file.exists(preprocess)){dir.create(preprocess,recursive = T)}
		if (opt$cycle){
			p=single_ob@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+theme_minimal()
			ggsave(paste(preprocess,"cells_cell_cycle",sep="/"),p,width = 8,height = 7)
			#$ggsave(paste(preprocess,"cells_cell_cycle.png",sep="/"),width = 8,height = 7)
			#p<-dittoBarPlot(single_ob,"Phase", group.by = "sample",color.panel=colors)
			p<- PlotBarplot(single_ob,prop.by = "Phase", 
                group.by = "sample.name", cols =colors)
			ggsave(p,filename=paste(preprocess,"cell_cycle_per_sample",sep="/"))
			#ggsave(p,filename=paste(preprocess,"cell_cycle_per_sample.png",sep="/"))

			cycle_p<-PlotBarplot(single_ob,prop.by="Phase", group.by = "seurat_clusters",cols = colors)
      ggsave(paste(preprocess,"cell_cycle_per_cluster_barplot",sep="/"),cycle_p)     
			cycle_p2 <- DimPlot(single_ob,group.by="Phase",cols=colors,reduction='umap',label=F)+theme(axis.text.x = element_text(size=14,face="bold"),
		axis.text.y = element_text(size=14,face="bold"),
		axis.title.x = element_text(size=16,face="bold"),
		axis.title.y = element_text(size=16,face="bold"),
		aspect.ratio = 1/1)
			ggsave(paste(preprocess,"cell_cycle_per_cluster_umap",sep="/"),cycle_p2,bg = 'white')
		}
		pc_use=dim  #可修改，默认20
    Elbowplot = ElbowPlot(single_ob,ndims = 30)$data %>% ggplot() +geom_point(aes(x=dims,y=stdev)) + 
	geom_vline(xintercept=pc_use,color='darkred')+
	  ggtitle('Elbow plot-quantitative approach')	+xlab("PCs")+ ylab("Standard Deviation")+
		theme(axis.text.x = element_text(size=14,face="bold"),
			axis.text.y = element_text(size=14,face="bold"),
			axis.line = element_line(colour = "black"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			panel.background = element_blank(),
			plot.title=element_text(hjust=0.5))  
			ggsave(file.path(preprocess,'PCA_ElbowPlot'),Elbowplot,width=10,bg='white')

			##髙变基因绘图
			assay <- ifelse(gather == 'CCA', 'SCT', 'RNA')
			vargene_plot <- VariableFeaturePlot(single_ob,assay=assay,cols=c('#999999','#007ACC')) 
			vargene_top10<- head(VariableFeatures(single_ob,assay=assay), 10)
			if ( all(nchar(vargene_top10) < 10)){
				vargene_plot <- LabelPoints(plot = vargene_plot, points = vargene_top10, repel = TRUE)
			}

			vargene_plot =vargene_plot+theme_bw()+theme(
				panel.grid=element_blank(),
				axis.text.x = element_text(size=14,face="bold"),
				axis.text.y = element_text(size=14,face="bold"),
				axis.line = element_line(colour = "black"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),aspect.ratio = 1/1) 

			ggsave(file.path(preprocess,'VariableFeatures_distribution'),vargene_plot,bg='white',width=8)

			###数据特征绘图
			quiet_library("RColorBrewer")
			scFeature=unlist(str_split(opt$scFeature,","))
			parallel::mclapply(scFeature, function(feature) {
  			if (feature %in% colnames(single_ob@meta.data)) {
    			nfeatureplot = FeaturePlot(single_ob, reduction = "umap", label = FALSE, features = feature) +
      			scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
    			ggsave(file.path(preprocess,paste0(feature,"_umap")),nfeatureplot,bg='white',width=8)
				} else {
					print(glue::glue("The {feature} is not in the meta.data"))
				}}, mc.cores = detectCores() - 1)
			}
	Anal_time = Sys.time() - start_time
	print(paste("降维聚类 分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
	quit()
}
#------------------------------------SUB cluster visualize--------------------------------------------
if ( "Visualize" %in% args ){
	quiet_library("rstatix")
	quiet_library("ggpubr")
	source("/home/donghj/scRNA_snakemake/special_umap.r") #TODO:路径修改
	groupby = opt$groupby
	message("开始针对单细胞数据的",groupby,"标签进行可视化")
	#公用参数
	start_time = Sys.time()
	colors = discrete_palette[[opt$pattle]]
	color_use = colors
  if(is.null(opt$input)){
		print('No input file or rds')
		quit()
	}else {
		single_ob = readRDS(opt$input)
		if (grepl('5',single_ob@version)){
			 single_ob <- JoinLayers(single_ob)
		}
	}
	if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
  }
	DefaultAssay(single_ob) = opt$assay
	#umap图和Tsne图
	#文件夹创建
	splitby = opt$splitby #默认sample,group
	reductions = opt$reductions #默认umap,tsne
	pointsize = opt$pointsize #默认NULL
	if (opt$report) {
		output_dir =file.path(output,"02_Cluster")
		summary_result = file.path(output,paste0("03_",groupby,"_proportion_visualization"))
	}else{
		output_dir = file.path(output,"01_Cluster")
		summary_result = file.path(output,paste0("02_",groupby,"_proportion_visualization"))
	}
	if(!file.exists(output_dir)){dir.create(output_dir,recursive = T)}
	#将metadata表格导出来
	splitbys = validate_metadata_columns_strict(single_ob,splitby)
	groupbys = validate_metadata_columns_strict(single_ob,groupby)
  clustering_df = single_ob@meta.data %>% mutate(Barcode = rownames(single_ob@meta.data)) %>%
    dplyr::select( Barcode, sample.name, seurat_clusters, group,setdiff(c(splitbys,groupbys), c("Barcode", "sample.name", "seurat_clusters", "group")))
  write.table(clustering_df, quote = F,sep =",",row.names = F,
              file.path(output_dir,paste0("clustering_results.csv",collapse = "")))
  
	if (toupper(opt$type) =="ALL"){
		type = c("Dimplot","Summary")
	}else {
		type = unlist(strsplit(type,","))
	}
  #将type改为首字母大写
  type = tolower(type)
	#可视化
  if ("dimplot" %in% type){
		VisualizeCluster(single_ob,reductions=reductions,groupby=groupby,splitby=splitby,outdir=output_dir,colors=color_use,pointsize=pointsize)
  }
	if ("summary" %in% type) {
		if(!file.exists(summary_result)){dir.create(summary_result)}
		SummaryCluster(single_ob,groupby=groupby,splitby=splitby,outdir=summary_result,colors)
	}
	
	Anal_time = Sys.time() - start_time
	print(paste("cluster可视化 分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
}
# ------------------------------------SUB findallmarker--------------------------------------------
if ( "marker" %in% args ){
	print("开始进行findmarker")
	start_time = Sys.time()
	colors = discrete_palette[[opt$pattle]]
	heat_colors = opt$heat_colors #TODO:增加热图的色板等 还没写到代码里面，现在还是默认的颜色
	#source("enrichment.r")
  if(is.null(opt$input)){
		print('No input file or rds')
		quit()
	}else {
		single_ob = readRDS(opt$input)
		if (grepl('5',single_ob@version)){
			 single_ob <- JoinLayers(single_ob)
		}
	}
	if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
  }
	DefaultAssay(single_ob) = opt$assay
	if (opt$groupby %in% colnames(single_ob@meta.data)){
		groupby = opt$groupby
		Idents(single_ob)<-groupby
	}else{
		stop("The groupby is not in the meta.data,please check your input")
	}

	#findmarker 
	if (opt$testuse %in% c("negbinom", "poisson","DESeq2")){
		slot ="counts"
		print(glue::glue("您选择的方法是{opt$testuse},slot将强制采用counts数据"))
	}else{
		print(glue::glue("使用{opt$slot}数据进行下游的findmarker运算"))
		slot =opt$slot
	}


	markers <- FindAllMarkers(single_ob, only.pos = T, min.pct = opt$min_pct, logfc.threshold = opt$logfcthreshold,test.use = opt$testuse,slot = slot)
	markers = markers %>%filter(p_val < 0.05 )
	markers = markers %>% select(gene,cluster,p_val,avg_log2FC,pct.1,pct.2,p_val_adj)

	quiet_library("openxlsx")
	readme_title = createStyle(
      fontSize = 12, textDecoration = "bold", 
      halign = "center", fgFill = "#D9D9D9", 
      border = "bottom"
    )
	wb <- createWorkbook()
	sheet <- addWorksheet(wb,"Marker Gene Summary")
	writeData(wb, sheet, markers, startCol = 1, startRow = 1, colNames = TRUE)
	addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:7, gridExpand = TRUE)
	sheet = addWorksheet(wb,"README")
	writeData(wb, sheet, t(as.data.frame(c("标题", "说明"))), startCol = 1, startRow = 1, colNames = FALSE)
	addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:2, gridExpand = TRUE)
	readme_list <- list(
    "gene"        = "基因名称",
    "cluster"      = "细胞标签信息",
    "p_val"    = "显著性p值",
    "avg_log2FC"       = "差异倍数的log2值",
    "pct.1"       = "基因在当前cluster中有表达的细胞比例",
    "pct.2"           = "基因在除当前cluster以外所有的cluster有表达的细胞比例",
		"p_val_adj" = "校正后的p值"
  )
	setColWidths(wb, sheet, cols = 2, widths = 100)
	readme_list_filter <- rownames_to_column(as.data.frame(t(as.data.frame(readme_list))),var = "tittle")
	#write.table(filter_data,file.path(output,"QC_cell_counts_per_step.xls"),quote = F,row.names = T,col.names = NA,sep='\t')
	writeData(wb, sheet, readme_list_filter, startCol = 1, startRow = 2, colNames = FALSE)
	saveWorkbook(wb, file = file.path(output, "allmarkers.xlsx"), overwrite = TRUE)

	#write.table(markers,paste(output,"allmarkers.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)


  if (opt$plot){
		topn=opt$topn
		lapply(unique(as.character(markers$cluster)), function(clust_num) {
			cluster_dir = file.path(output, paste0("Each_",groupby,"_marker"), paste(groupby, clust_num, sep = "_"))
			
			if (!file.exists(cluster_dir)) {
					dir.create(cluster_dir, recursive = TRUE)
			}
			# 子集提取
			cluster_markers <- subset(markers, cluster == clust_num)
			rownames(cluster_markers) <- cluster_markers$gene
			index <- paste0(groupby,"_", clust_num)
				
			write.table(
					cluster_markers, 
					file.path(cluster_dir, paste("cluster", clust_num, "markers.xls", sep = "_")), 
					sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
				)
			

			if(nrow(cluster_markers) > 1 ){
					top10_markers = cluster_markers %>% top_n(n =topn , wt = avg_log2FC)
					height= ceiling(length(top10_markers$gene)/5)*3
					vln_p = VlnPlot(single_ob, features = top10_markers$gene, pt.size = 0, ncol = 5, cols = colors)
					ggsave(paste(cluster_dir, paste0(index, glue::glue("_top{topn}_vlnplot")), sep = "/"),vln_p, width = 20, height = height, bg = '#ffffff')
		
					Fea_p=FeaturePlot(single_ob, features = top10_markers$gene, min.cutoff = "q9", ncol = 5, order = TRUE, cols = c("lightgrey", "red"))
					ggsave(paste(cluster_dir, paste0(index,  glue::glue("_top{topn}_umap")), sep = "/"),Fea_p,width = 20, height = height, bg = '#ffffff')
			}
	})	
		#heatmap
		avetopn = opt$avetopn
		Topn_ave_marker_heatmap(single_ob,collapseby = groupby,Topn = avetopn,
		seurat_exp_cluster_dir=output,markers=markers,colors =colors )
		cluster_top10_markers=markers %>% group_by(cluster)%>% top_n(n = topn, wt = avg_log2FC)
		if (!all(cluster_top10_markers$gene %in% VariableFeatures(single_ob))){
			single_ob <- ScaleData(single_ob, feature=cluster_top10_markers$gene,verbose = FALSE)
		}
		if (length(cluster_top10_markers$gene) > 135){
			sz = 4-log(length(cluster_top10_markers$gene)/100)
			height = 5+log2(length(cluster_top10_markers$gene)/10)
			width = 7.5
		}else if (length(cluster_top10_markers$gene) > 75){
			sz = 4-log2(length(cluster_top10_markers$gene)/120);height = 7;width = 7
		}else{
			sz = 6-log2(length(cluster_top10_markers$gene)/80);height = 7;width = 7
		}

		ggheat = Seurat::DoHeatmap( object = single_ob,
                        features = cluster_top10_markers$gene,
                        group.colors = colors,
                        group.by = groupby, group.bar = T, label = F, draw.lines = F) +ggplot2::theme(axis.text.y = ggplot2::element_text(size = sz,face = "bold"))+ggplot2::scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white",
                        high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0,
                        na.value="white",
                        guide = "colourbar",
                        aesthetics = "fill")
    ggsave(paste(output,glue::glue("cluster_top{topn}_markers_heatmap"),sep="/"),ggheat,width = width,height = height)
	}
	if (opt$report){
		#marker_number
		marker_number<-table(markers$cluster)%>% reshape2::melt()
    marker_number$Cluster <- as.factor(marker_number$Var1)
		number_of_bar <- nrow(marker_number)
		marker_number$id <- seq(1, number_of_bar)
		angle <- 90 - 360 * (marker_number$id-0.5) /number_of_bar

		marker_number$hjust <- ifelse( angle < -90, 1, 0)
		marker_number$angle <- ifelse(angle < -90, angle+180, angle)

		p <- ggplot(marker_number, aes(x=as.factor(id), y=value, fill=Cluster)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
		geom_bar(stat="identity", alpha=0.5)+scale_fill_manual(values=colors) +
		ylim(-(min(marker_number$value)-10),NA) +
		theme_minimal() +
		theme(
			#legend.position = "none",
			axis.text = element_blank(),
			axis.title = element_blank(),
			panel.grid = element_blank(),
			plot.margin = unit(rep(0.1,4), "cm") 
		)  +
		coord_polar() + 
		geom_text(data=marker_number, aes(x=id, y=value+10, label=as.vector(value), hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3, angle= marker_number$angle, inherit.aes = FALSE ) 

    ggsave(file.path(output,"marker_number"),p,width = 8,height = 8,bg='white')

		#dotplot,如何clusters的数量超过30个，则只取前20个cluster进行绘图
		if (length(unique(markers$cluster))<31){
			eachclusters_top5_markers=markers %>% filter(avg_log2FC != Inf) %>%group_by(cluster)  %>%  top_n(n = 5, wt = avg_log2FC) %>% as.data.frame() %>% dplyr::distinct(.,gene,.keep_all = T)
			data_ob =single_ob
		}else {
				#markers$cluster = as.numeric(markers$cluster)
				eachclusters_top5_markers=markers %>% filter( cluster %in% c(0:20), avg_log2FC != Inf) %>%group_by(cluster)  %>%  top_n(n = 5, wt = avg_log2FC) %>% as.data.frame() %>% dplyr::distinct(.,gene,.keep_all = T)
				data_ob = subset(single_ob,seurat_clusters %in% c(0:20))
		}
		 
		
		if(length(eachclusters_top5_markers$gene) > 9){
				direction ="vertical"
				box = "vertical"
		}else{
				direction ="vertical"
				box = "horizontal"
		}

		ggdots = Seurat::DotPlot(object = data_ob,dot.scale=4,
										features = eachclusters_top5_markers$gene ) +							
										Seurat::RotatedAxis() +
										guides(color = guide_colorbar(title = "Exp avg"), 
														size = guide_legend(title = "Pct Exp %")) +		
										#ggplot2::coord_flip() + 
										ggplot2::scale_colour_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100) ) + 
										theme( legend.position = "right", # 设置图例位置在底部
														legend.direction = direction, # 设置图例排列方向为水平
														legend.box = box, # 设置图例框的方向为水平
														legend.box.just = "center", # 设置图例框居中显示
														legend.spacing.x = unit(0.3, "cm"), # 设置图例条目之间的水平间距
														axis.text.x = element_text(size=9)
												)

	  ggsave(file.path(output,"cluster_top5_dotplot"),ggdots,width = length(eachclusters_top5_markers$gene)*0.16+1.2,height = 0.3*length(unique(data_ob@meta.data$seurat_clusters)),limitsize = FALSE,bg="white")
   
	}
	Anal_time = Sys.time() - start_time
	print(paste("findmarker 分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
	quit()
}
# ------------------------------------SUB Enrichement ----------------------------------------------
if ("Enrichments" %in% args) {
	source("/home/donghj/scRNA_snakemake/extra_dir/Enrichment_circor.r")
	quiet_library("future.apply")
	quiet_library("openxlsx")
	species_mapping <- list(
  "human" = "Homo_sapiens", # human
  "mouse" = "Mus_musculus", # mouse
  "rat" = "Rattus_norvegicus", # rat
  "dme" = "fruit_fly", # fruit_fly
  "dre" = "zebrafish_danio_rerio", # zebrafish
  "ath" = "Oryza_sativa_Japonica_Group", # Arabidopsis
  "sce" = "Saccharomyces_cerevisiae_S288C", # yeast
  "cel" = "Caenorhabditis_elegans", # C.elegans
  "bta" = "Bos_taurus", # Bovine
  "mcc" = "macaca_mulatta", # monkey
  "cfa" = "canis_lupus_familiaris", # dog
  "ssc" = "sus_scrofa", # pig
  "gga" = "Gallus_gallus" # chicken
	)

	species_package <- species_mapping[[opt$species]]
	if(is.null(species_package)){message("species not found, please check your input");quit()}
	start_time = Sys.time()
	print("富集分析开始")

	if (endsWith(opt$file, ".xlsx")){
		quiet_library("openxlsx")
		gene_file = read.xlsx(opt$file, sheet = 1)
	}else{
		gene_file = read.delim(opt$file, header = T, stringsAsFactors = F)
	}
  
	if (!is.null(opt$CirchordColor)){
		CirchordColor = discrete_palette[[opt$CirchordColor]]
	}else{
		CirchordColor = c(
        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
        "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
        "#CCEBC5", "#FFED6F", "#71A99F", "#CCCC8F", "#9895AE",
        "#C9665B", "#668EA9", "#CA904E", "#8FB254", "#CAA4B7"
    )
	}
  readme_title = createStyle(
					fontSize = 12, textDecoration = "bold", 
					halign = "center", fgFill = "#D9D9D9", 
					border = "bottom"
	)
	readme_Kegg <- list(
  "ID" = " KEGG编号",
  "Description" = "KEGG描述信息",
  "GeneRatio" = "和该功能相关的差异基因与差异基因总数的比值",
  "BgRatio" = "该功能相关的基因总数与基因组所有基因总数的比值",
  "pvalue" = "P值",
  "p.adjust" = "校正P值",
  "qvalue" = "Q值（Q方法校正后的P值）",
  "geneID" = "与该功能相关的差异基因ID",
  "Count" = "与该功能相关的差异基因的数目",
  "geneSymbol" = "与该功能相关的差异基因",
  "keggLink" = "KEGG网页绘制通路链接",
  "level_a" = "KEGG第一层级",
  "level_b" = "KEGG第二层级"
  )
  readme_Go <- list(
    "ID" = "GO编号",
    "Description" = "GO描述信息",
    "GeneRatio" = "和该功能相关的差异基因与差异基因总数的比值",
    "BgRatio" = "该功能相关的基因总数与基因组所有基因总数的比值",
    "pvalue" = "P值",
    "p.adjust" = "校正P值",
    "qvalue" = "Q值（Q方法校正后的P值）",
    "geneID" = "与该功能相关的差异基因",
    "Count" = "与该功能相关的差异基因的数目",
    "Subontologies" = "GO分类"
  )
	style <- createStyle(
					fontSize = 14,     # 字体大小
					border = "TopBottomLeftRight" # 给单元格加上框线
	)
	
	if ("cluster" %in% colnames(gene_file)) {
		results <- future_lapply(unique(as.character(gene_file$cluster)), function(clust_num) {
			tryCatch({
				# 创建临时和输出目录
				temp <- file.path(output, "temp", clust_num)
				dir.create(temp, recursive = TRUE, showWarnings = FALSE)
				cluster_dir <- file.path(output, paste0(clust_num, "_marker_enrichment"))
				dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
				
				# 筛选该 cluster 的基因
				temp_maker <- gene_file[gene_file$cluster == clust_num, ]
				genes <- temp_maker$gene
				temp_file <- file.path(temp, paste0(clust_num, "_marker_gene.txt"))
				write.table(genes, temp_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
				
				# 构造 Perl 命令，并重定向所有输出到 /dev/null
				cmd <- glue::glue(
					"perl /home/genesky/pipeline/tools_script/gene_enrichment/latest/gene_enrichment.pl ",
					"-g {temp_file} -a enrich -t {temp} -o {cluster_dir} ",
					"-s {species_package} -mode local > /dev/null 2>&1"
				)
				
				wb <- createWorkbook()
				sheet <- addWorksheet(wb,"Gene")
				writeData(wb, sheet, genes, startCol = 1, startRow = 1, colNames = TRUE)
				# 给数据区域加样式（包括表头 + 数据）
				addStyle(
					wb, sheet,
					style = style,
					rows = 1:(length(genes) + 1),   # 行数：1 行表头 + n 行数据
					cols = 1,                     # 只有第 1 列
					gridExpand = TRUE
				)
				
				#saveWorkbook(wb, file = file.path(output, "allmarkers.xlsx"), overwrite = TRUE)
				# 执行命令，不捕获输出，也不生成日志
				system(cmd, intern = FALSE)
				
				# 读取 GO 富集结果并绘图
				EnrichDat_file <- file.path(cluster_dir, "enrichGO.txt")
				if (file.exists(EnrichDat_file)) {
					EnrichDat <- read.delim(EnrichDat_file, header = TRUE, stringsAsFactors = FALSE)
					sheet <- addWorksheet(wb,"GO Summary")
					writeData(wb, sheet, EnrichDat, startCol = 1, startRow = 1, colNames = TRUE)
					addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:10, gridExpand = TRUE)
					setColWidths(wb, sheet, cols = 1:10, widths = 15)
					pdf(file.path(cluster_dir, paste0("GO_Top", opt$topn, "_enrichment_pathway_circos.pdf")),
							width = 8, height = 7)
					PlotEnrichCircos(EnrichDat, gene_diff = temp_maker, adjust.pvalue = FALSE,
													show_RichFactor = TRUE, title_override = NULL,
													color_palette = CirchordColor, topn = opt$topn, rankby = opt$rankby)
					dev.off()
					
					png(file.path(cluster_dir, paste0("GO_Top", opt$topn, "_enrichment_pathway_circos.png")),
							width = 800, height = 700)
					PlotEnrichCircos(EnrichDat, gene_diff = temp_maker, adjust.pvalue = FALSE,
													show_RichFactor = TRUE, title_override = NULL,
													color_palette = CirchordColor, topn = opt$topn, rankby = opt$rankby)
					dev.off()
					file.remove(EnrichDat_file)
				}
				
				# 读取 KEGG 富集结果并绘图
				EnrichDat_file_kegg <- file.path(cluster_dir, "enrichKEGG.txt")
				if (file.exists(EnrichDat_file_kegg)) {
					EnrichDat_kegg <- read.delim(EnrichDat_file_kegg, header = TRUE, stringsAsFactors = FALSE)
					sheet <- addWorksheet(wb,"KEGG Summary")  
					writeData(wb, sheet, EnrichDat_kegg, startCol = 1, startRow = 1, colNames = TRUE)
					addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:13, gridExpand = TRUE)
					setColWidths(wb, sheet, cols = 1:13, widths = 15)
					pdf(file.path(cluster_dir, paste0("KEGG_Top", opt$topn, "_enrichment_pathway_circos.pdf")),
							width = 9, height = 7)
					PlotEnrichCircos(EnrichDat_kegg, gene_diff = temp_maker, adjust.pvalue = FALSE,
													show_RichFactor = TRUE, title_override = NULL,
													color_palette = CirchordColor, topn = opt$topn, rankby = opt$rankby)
					dev.off()
					
					png(file.path(cluster_dir, paste0("KEGG_Top", opt$topn, "_enrichment_pathway_circos.png")),
							width = 900, height = 700)
					PlotEnrichCircos(EnrichDat_kegg, gene_diff = temp_maker, adjust.pvalue = FALSE,
													show_RichFactor = TRUE, title_override = NULL,
													color_palette = CirchordColor, topn = opt$topn, rankby = opt$rankby)
					dev.off()
					file.remove(EnrichDat_file_kegg)
				}
				#生成一个xlsx
				for (sheet_name in c("README GO", "README KEGG")) {
					sheet <- addWorksheet(wb, sheet_name)
					writeData(wb, sheet, t(as.data.frame(c("标题", "说明"))), 
										startCol = 1, startRow = 1, colNames = FALSE)
					addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:2, gridExpand = TRUE)
					# 选择对应的数据
					if (sheet_name == "README GO") {
						readme_data <- readme_Go
					} else {
						readme_data <- readme_Kegg
					}
					# 设置列宽和换行
					setColWidths(wb, sheet, cols = 1, widths = 15)
					setColWidths(wb, sheet, cols = 2, widths = 40)
					# addStyle(wb, sheet, createStyle(wrapText = TRUE), 
					# 				rows = 2:(nrow(readme_data)+1), cols = 2, gridExpand = TRUE)
					readme_list_filter <- rownames_to_column(as.data.frame(t(as.data.frame(readme_data))), var = "tittle")
					writeData(wb, sheet, readme_list_filter, startCol = 1, startRow = 2, colNames = FALSE)
				}
				
				saveWorkbook(wb, file = file.path(cluster_dir, "Gene_Enrich_Summary.xlsx"), overwrite = TRUE)

				return(TRUE)  # 标记成功
			}, error = function(e) {
				message("ERROR in cluster ", clust_num, ": ", e$message)
				return(FALSE) # 标记失败
			})
		}, future.seed = TRUE)
	}
  #删除临时文件
	if (dir.exists(file.path(output, "temp"))) {
  # 删除里面的文件
  #file.remove(list.files(file.path(output, "temp"), full.names = TRUE))
  # 删除目录本身
  unlink(file.path(output, "temp"), recursive = TRUE, force = TRUE)
	}
	Anal_time = Sys.time() - start_time
	print(paste("富集分析 分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
	quit()
}
#------------------------------------sub singleR ----------------------------------------------------
if ( "singleR" %in% args ){
	print("开始使用singleR进行自动细胞类型注释")
	start_time = Sys.time()
	colors = discrete_palette[[opt$pattle]]
	#heat_colors = opt$heat_colors #TODO:增加热图的色板等 还没写到代码里面，现在还是默认的颜色
	#source("enrichment.r")
  if(is.null(opt$input)){
		print('No input file or rds')
		quit()
	}else {
		single_ob = readRDS(opt$input)
		if (grepl('5',single_ob@version)){
			 single_ob <- JoinLayers(single_ob)
		}
	}
	if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
  }
	DefaultAssay(single_ob) = opt$assay
	
  #加载额外的包
	quiet_library("SingleR")
	quiet_library("RColorBrewer")
	quiet_library("SingleCellExperiment")
	#参考的rds路径
	ref_file = opt$ref
	if (!file.exists(ref_file)){
		print("singleR ref file not exist")
		stop()
	}else {
		 ref_sc = readRDS(ref_file)
	}
  
	ref_celltype = opt$groupby

	
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	sc <- scale_colour_gradientn(colours = myPalette(100))

	#如果细胞数太多，采用随机抽样到30000
	downsample = opt$downsample
	if(ncol(single_ob) > downsample){
	# library(sampling) need install,instead with manual function
		ratio <- as.numeric(downsample) / ncol(single_ob)
		metadata_temp <- as.data.frame(single_ob@meta.data)
		# strata(metadata_temp,stratanames="clusters",ratio,description=FALSE)
		cells_sample <- c()
		for (i in unique(single_ob$seurat_clusters)) {
			cells_temp <- rownames(metadata_temp)[which(metadata_temp$seurat_clusters == i)]
			#设置随机种子，确保结果可复原
			set.seed(2024)
			cells_temp_sample <- sample(cells_temp, ceiling(length(cells_temp) * ratio), replace = FALSE, prob = NULL)
			cells_sample <- append(cells_sample, cells_temp_sample)
		}
		single_ob <- subset(single_ob, cells = cells_sample)
		print("单组细胞数过多，已进行降采样。")
		print(paste0(c("基因数是：","细胞数是："),dim(single_ob)))
    gc()
	}
	data_SingleR <- GetAssayData(single_ob, slot="data")
  DefaultAssay(ref_sc)="RNA"
	if (ref_celltype %in% colnames(ref_sc@meta.data)){
		ref_sc = as.SingleCellExperiment(ref_sc)
		ref_lab = ref_sc[[ref_celltype]]
	}else{
		stop("The celltype is not in the ref meta.data,please check your input")
	}
  anno_result <- SingleR(test = data_SingleR, ref = ref_sc, labels = ref_lab,  num.threads =  opt$ncores)

##初步展示那些细胞类型显著标记出来
	p1 = plotScoreHeatmap(anno_result)
	if( max(nchar(p1$tree_row$labels)) > 30){
		resize_w = 20
	}else if(max(nchar(p1$tree_row$labels)) > 15 & max(nchar(p1$tree_row$labels)) <=30){
		resize_w = 10 + 0.3* (max(nchar(p1$tree_row$labels))/10)
	} else {
		resize_w = 10 
	}

	ggsave(paste0(output, "/SingleR_Preview.heatmap"), plot = p1, height = 10, width = resize_w)

	raw_anno = as.data.frame(anno_result['labels'])
	single_ob = AddMetaData(single_ob,raw_anno,col.name="SingleR_raw")

	#2.按照clusters 注释
	single_ob@meta.data$seurat_clusters = factor(single_ob@meta.data$seurat_clusters,levels = sort(unique(single_ob@meta.data$seurat_clusters)))

	anno_result2 = SingleR(test = data_SingleR, 
								ref = ref_sc, 
					method = "cluster",
								clusters = single_ob$seurat_clusters, 
								labels = ref_lab)

	single_ob@meta.data['SingleR_anno'] = anno_result2$labels[match(single_ob$seurat_clusters,
																						rownames(anno_result2))]

	anno_df = as.data.frame(single_ob@meta.data) %>% tibble::rownames_to_column(.,var="barcodes")
	anno_df = anno_df %>% dplyr::select(barcodes,sample.name,group,seurat_clusters,SingleR_raw,SingleR_anno)
#------------
  quiet_library("openxlsx")
	readme_title = createStyle(
      fontSize = 12, textDecoration = "bold", 
      halign = "center", fgFill = "#D9D9D9", 
      border = "bottom"
    )
	wb <- createWorkbook()
	sheet <- addWorksheet(wb,"Metadata")
	setColWidths(wb, sheet, cols = 1:6, widths = 20)
	writeData(wb, sheet, anno_df, startCol = 1, startRow = 1, colNames = TRUE)
	addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:6, gridExpand = TRUE)
  sum = table(single_ob@meta.data$sample.name,single_ob@meta.data$SingleR_anno)
	percent =   sum / rowSums(sum)
  sum = rownames_to_column(as.data.frame.matrix(sum),var="Sample\\Cell_Type")
  sheet = addWorksheet(wb,"Sample Count")
	setColWidths(wb, sheet, cols = 1:ncol(sum), widths = 20)
  writeData(wb, sheet, sum, startCol = 1, startRow = 1, colNames = TRUE)
  addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:ncol(sum), gridExpand = TRUE)
#计算每个样本里面的细胞类型比例
	sheet = addWorksheet(wb,"Sample Percent")
	setColWidths(wb, sheet, cols = 1:ncol(percent), widths = 20)
	percent = rownames_to_column(as.data.frame.matrix(percent),var="Sample\\Cell_Type")
  writeData(wb, sheet, percent, startCol = 1, startRow = 1, colNames = TRUE)
  addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:ncol(percent), gridExpand = TRUE)

	sheet = addWorksheet(wb,"README")
	writeData(wb, sheet, t(as.data.frame(c("标题", "说明"))), startCol = 1, startRow = 1, colNames = FALSE)
	addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:2, gridExpand = TRUE)
	readme_list <- list(
    "barcodes"        = "细胞标签",
    "sample.name"      = "样本名",
    "group"    = "组名",
    "seurat_clusters"       = "降维聚类分组编号",
    "SingleR_raw"       = "singleR初步注释结果",
    "SingleR_anno"           = "最终注释结果"
  )
	setColWidths(wb, sheet, cols = 1:2, widths = 20)
	readme_list_filter <- rownames_to_column(as.data.frame(t(as.data.frame(readme_list))),var = "tittle")
	#write.table(filter_data,file.path(output,"QC_cell_counts_per_step.xls"),quote = F,row.names = T,col.names = NA,sep='\t')
	writeData(wb, sheet, readme_list_filter, startCol = 1, startRow = 2, colNames = FALSE)
	saveWorkbook(wb, file = file.path(output, "SingleR_anno.meta.xlsx"), overwrite = TRUE)
#------------------------
	# write.table(anno_df,file = file.path(output,"SingleR_anno.meta.xls"),
	# 						col.names = T,row.names = F,sep = "\t",quote=F)

	anno_result2$labels <- paste0("C",rownames(anno_result2)," ",anno_result2$labels)
	p2 = plotScoreHeatmap(anno_result2)

	if( max(nchar(p2$tree_row$labels)) > 30){
		resize_w = 20
	}else if(max(nchar(p2$tree_row$labels)) > 15 & max(nchar(p2$tree_row$labels)) <=30){
		resize_w = 10 + 0.3* (max(nchar(p2$tree_row$labels))/10)
	} else {
		resize_w = 10 
	}
	ggsave(paste0(output, "/SingleR_clusters_anno.heatmap"), plot = p2, height = 10, width = resize_w)

	p1 = DimPlot(single_ob, group.by = "seurat_clusters",reduction = "umap") + 
				scale_colour_manual( values = colors[1:length(unique(single_ob$seurat_clusters))])
	p2 = DimPlot(single_ob, group.by = "SingleR_anno",reduction = "umap") + 
				scale_colour_manual( values = colors[1:length(unique(single_ob$SingleR_anno))])
	p = p1 + p2
	ggsave(paste0(output, "/compare_clusters.anno"), plot = p, height = 7, width = 15)

	anno_score = as.data.frame(anno_result2$scores)
	anno_score$celltype <- anno_result2$labels
	anno_score$seurat_clusters <- rownames(anno_result2)
	write.table(anno_score, file.path(output,paste0("SingleR_clusters_anno.score.xls",collapse = "")),sep="\t",quote=F,row.names = F,col.names=T)
  saveRDS(single_ob,file.path(output,paste0(opt$prefix,".rds")))
	Anal_time = Sys.time() - start_time
	print(paste("singleR自动化分析鉴定 分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
	quit()
}
#------------------------------------sub loupe ----------------------------------------------------
if ( "loupeR" %in% args ){
	print("开始将RDS转成loupeR文件")
	start_time = Sys.time()
  if(is.null(opt$input)){
		print('No input file or rds')
		quit()
	}else {
		single_ob = readRDS(opt$input)
		if (grepl('5',single_ob@version)){
			 single_ob <- JoinLayers(single_ob)
		}
	}
	if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
  }
	DefaultAssay(single_ob) = opt$assay
	
  #加载额外的包
	quiet_library("loupeR")
  #loupeR::setup(executable_path = opt$loupe)
	output_name=opt$cloupe_name
  single_ob@meta.data <- single_ob@meta.data %>%
  select(sample.name, seurat_clusters, starts_with("RNA_snn_res"))
	fit <- try(create_loupe_from_seurat(single_ob, output_dir = output, output_name = output_name),silent=TRUE)
	if('try-error' %in% class(fit)){
		print('Some barocdes in your dataset seems not from 10X, so we use some whitelist barcode replace them')
		sample <- single_ob@meta.data$sample.name
		white <- read.table(gzfile("/home/genesky/software/cellranger/7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"), header = F)
		single_ob <- RenameCells(single_ob,new.names=paste(sample,white$V1[1:length(Cells(single_ob))],sep="_"))
		Idents(single_ob) <- single_ob@meta.data$seurat_clusters
		clusters <- select_clusters(single_ob)
		Wrong_col <- which(lapply(clusters, nlevels) >= 32768)
		if(length(Wrong_col)>0){
			print(paste(colnames(single_ob@meta.data)[Wrong_col],'have too many factors, auto remove them',collapse = "\n"))
			single_ob@meta.data <- single_ob@meta.data[,-Wrong_col]
		}
		create_loupe_from_seurat(single_ob, output_dir = output, output_name = output_name)
	}
	Anal_time = Sys.time() - start_time
	print(paste("RDS文件转loupe 分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
	quit()
}
# ------------------------------------SUB celltype -------------------------------------------------
if ( "Celltyping" %in% args ){
	print("开始根据注释文件对细胞进行注释")
	start_time = Sys.time()
	colors = discrete_palette[[opt$pattle]]
	#heat_colors = opt$heat_colors #TODO:增加热图的色板等 还没写到代码里面，现在还是默认的颜色
	#source("enrichment.r")
  if(is.null(opt$input)){
		print('No input file or rds')
		quit()
	}else {
		single_ob = readRDS(opt$input)
		if (grepl('5',single_ob@version)){
			 single_ob <- JoinLayers(single_ob)
		}
	}
	if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
  }
	DefaultAssay(single_ob) = opt$assay
#读取重命名文件
  if (!is.null(opt$re_sample_diff)){
		re_sample_diff = read.delim(opt$re_sample_diff)
		if (!all(unique(re_sample_diff$raw) %in% unique(single_ob$sample.name))){
			print("重命名文件中的样本名不在rds文件中")
			quit()}
		
		single_ob$group <- re_sample_diff$group[match(single_ob$sample.name,re_sample_diff$raw)]
    other_groups = setdiff(colnames(re_sample_diff),c("raw","sample","group"))
    if(length(other_groups)>0){
			for(other_group in other_groups){
				single_ob@meta.data[[other_group]] <- re_sample_diff[[other_group]][match(single_ob$sample.name,re_sample_diff$raw)]
			}
      #single_ob@meta.data[other_group] <- re_sample_diff[other_group][match(single_ob$sample.name,re_sample_diff$raw)]
		}
		single_ob$sample.name <- re_sample_diff$sample[match(single_ob$sample.name,re_sample_diff$raw)]
	}
  
#读取注释文件
	seurat_data <- read.delim(opt$cluster)
	colnames(seurat_data) <- tolower(colnames(seurat_data))

	if (!all(c("anno","cluster") %in% colnames(seurat_data))){
		print("anno 和cluster 必须要存在于注释表格中")
		quit()
	}else {
		selected_cols <- c("cluster","anno")
		seurat_anno <- seurat_data[,selected_cols]
		#colnames(seurat_anno) <- c("cluster","anno")
		seurat_anno <-seurat_anno %>% as_tibble() %>% 
  separate_rows(cluster, sep = ",")
		if(length(seurat_anno$cluster)!=length(unique(seurat_anno$cluster))){return(message("错误，cluster重复"))}
		new.cluster.ids <-as.vector(seurat_anno$anno)
		names(new.cluster.ids) <-seurat_anno$cluster
		Idents(single_ob) = "seurat_clusters"
		single_ob <- subset(single_ob,idents=as.vector(seurat_anno$cluster))
		single_ob <- RenameIdents(single_ob, new.cluster.ids)
		single_ob$celltype <- Idents(single_ob)
		}
		saveRDS(single_ob,file.path(output,paste0(opt$prefix,".rds")))
		#-----------------------------plot----------------------------------------
		rds = file.path(output,paste0(opt$prefix,".rds"))
		plot_out = file.path(output,"Ident_CellType_Markers")
		dir.create(plot_out,recursive = T)
    if("marker" %in% colnames(seurat_data)){
			cmd = glue::glue("source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && Rscript /home/donghj/script/plot_gene_unified.R -i {rds} -g {opt$cluster} -t dot  -b celltype -c customecol2_light -f F -o {plot_out}  &&  Rscript /home/donghj/script/plot_gene_unified.R -i {rds} -g {opt$cluster} -t violin -b celltype -c customecol2_light -f F -o {plot_out} --stack T && Rscript /home/donghj/script/plot_gene_unified.R -i {rds} -g {opt$cluster} -t heatmap -b celltype -c customecol2_light -f F -o {plot_out} --averageexp T")
			print("marker基因在加标签之后的rds的表达情况，该类型图片作为细胞类型鉴定的证据！")
			#print(cmd)
			system(cmd,ignore.stdout=T)
		}

    Anal_time = Sys.time() - start_time
		print(paste("细胞加标签已经分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
		quit()
}
# ------------------------------------SUB different gene--------------------------------------------
if ( "diffexp" %in% args ){
	print("差异分析开始")
	colors = discrete_palette[[opt$pattle]]

  if(is.null(opt$input)){
		print('No input file or rds')
		quit()
	}else {
		single_ob = readRDS(opt$input)
		if (grepl('5',single_ob@version)){
			 single_ob <- JoinLayers(single_ob)
		}
	}
	if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
  }
	DefaultAssay(single_ob) = opt$assay


	start_time = Sys.time()
	
	# -------------------------------findmarker--------------------------------------------
  if ( is.null(opt$pvalue) && is.null(opt$fdr) ){
      stop("None of P-value or FDR is AVAILABLE! Filtering can be performed using any one of (--pvalue), (--fdr) at one time.", call.=FALSE)
  	}else if ( !is.null(opt$fdr)){
      fdr = opt$fdr
      pvalue = NULL
  	}else{
      pvalue = opt$pvalue
      fdr = NULL
    }
  	logfc_threshold = opt$logfc
  	if (is.null(opt$FC)){
      fc.thres = 0
  	}else{
      fc.thres = opt$FC
  	}
  	if (is.null(opt$test)){
      test = "wilcox"
  	}else{
      test = opt$test
		}
		contrasts_list = unlist(strsplit(opt$contrasts,","))
  	groupby_list = unlist(lapply(contrasts_list, function(x) strsplit(x,":")[[1]][1]))
		group_diff = setdiff(groupby_list, colnames(single_ob@meta.data))
		if (length(group_diff)!=0){
			stop(paste("选择的差异分组中的", group_diff, "不在 rds 的 meta.data 中"), call. = TRUE)
		}
		splitby=opt$splitby
		if(!is.null(splitby)) {
			if (!splitby %in% colnames(single_ob@meta.data)){
				stop(paste("选择的差异分组中的", group_diff, "不在 rds 的 meta.data 中"), call. = TRUE)
			}else{
				print(glue::glue("将按照{splitby}拆分做两组细胞之间的差异分析"))
			}
			for (split in unique(single_ob@meta.data[[splitby]])){
				cells_keep <- rownames(single_ob@meta.data[single_ob@meta.data[[splitby]] == split, ])
				sub_ob <- subset(single_ob, cells = cells_keep)
				output_dir = file.path(output,split)
				dir.create(output_dir,recursive = T)
				future.apply::future_lapply( contrasts_list, function( contrast){
				RunDiffexp( object = sub_ob,
								test = opt$test,fdr = fdr,
								fc.thres = opt$FC,
								pval.thres = pvalue,
								contrast = contrast,
								outputdir = output_dir,
								logfc_threshold = logfc_threshold)
				}, future.seed = 2020)
				#绘图，热图和火山图，boxplot
				#火山图
				suppressWarnings(suppressMessages(future.apply::future_lapply( dir(output_dir,pattern = ".*-vs-.*\\.xlsx"), function( diff_xlsx){
					base_prefix = gsub("_diff_result.xlsx","",diff_xlsx)
					diff_gene_file = file.path(output_dir,diff_xlsx)
					diff_gene = openxlsx::read.xlsx(diff_gene_file)
					volcano_out = file.path(output_dir,"volcano")
					dir.create(volcano_out,recursive = T)
					volcano_plot(diff_gene, volcano_out, base_prefix,opt$FC,pvalue)
				}, future.seed = 2020)))
        #热图，这里涉及到scaledata计算，不适用并行
				for(diff_group in contrasts_list){
					heatmap_out = file.path(output_dir,"heatmap")
					dir.create(heatmap_out,recursive = T)
					diff_group = unlist(strsplit(diff_group,":"))
          diff_xlsx = paste0(diff_group[1],"_",diff_group[2],"-vs-",diff_group[3],"_diff_result.xlsx")
					diff_gene_file = file.path(output_dir,diff_xlsx)
					if (file.exists(diff_gene_file)){
						base_prefix = gsub("_diff_result.xlsx","",diff_xlsx)
						diff_gene = openxlsx::read.xlsx(diff_gene_file)
						cells_keep <- rownames(sub_ob@meta.data[sub_ob@meta.data[[diff_group[1]]] %in% c(diff_group[2],diff_group[3]), ])
						sub_ob2 <- subset(sub_ob, cells = cells_keep)
						if(is.factor(sub_ob2@meta.data[,diff_group[1]] )){
							sub_ob2@meta.data[,diff_group[1]] = droplevels(sub_ob2@meta.data[,diff_group[1]])
						}else{
							sub_ob2@meta.data[,diff_group[1]] = factor(sub_ob2@meta.data[,diff_group[1]])
						}
						box_dir = file.path(output_dir,"boxplot")
						dir.create(box_dir,recursive = T)
						suppressWarnings(suppressMessages(box_fecth(diff_gene,sub_ob2,diff_group[1],box_dir, base_prefix)))
						sub_ob2=downsample(sub_ob2,8000,diff_group[1])
						suppressWarnings(suppressMessages(diff_heatmap(diff_gene,sub_ob2,colors,diff_group[1],heatmap_out,base_prefix)))

					}
				}
			}
			}else{
					output_dir = file.path(output,"all_cell")
					dir.create(output_dir,recursive = T)
					future.apply::future_lapply( contrasts_list, function( contrast){
					RunDiffexp( object = single_ob,
									test = opt$test,fdr = fdr,
									fc.thres = opt$FC,
									pval.thres = pvalue,
									contrast = contrast,
									outputdir = output_dir,
									logfc_threshold = logfc_threshold)
					}, future.seed = 2020)
					#火山图
					suppressWarnings(suppressMessages(future.apply::future_lapply( dir(output_dir,pattern = ".*-vs-.*\\.xlsx"), function( diff_xlsx){
						base_prefix = gsub("_diff_result.xlsx","",diff_xlsx)
						diff_gene_file = file.path(output_dir,diff_xlsx)
						diff_gene = openxlsx::read.xlsx(diff_gene_file)
						volcano_out = file.path(output_dir,"volcano")
						dir.create(volcano_out,recursive = T)
						volcano_plot(diff_gene, volcano_out, base_prefix,opt$FC,pvalue)
					}, future.seed = 2020)))

					for(diff_group in contrasts_list){
					heatmap_out = file.path(output_dir,"heatmap")
					dir.create(heatmap_out,recursive = T)
					diff_group = unlist(strsplit(diff_group,":"))
          diff_xlsx = paste0(diff_group[1],"_",diff_group[2],"-vs-",diff_group[3],"_diff_result.xlsx")
					diff_gene_file = file.path(output_dir,diff_xlsx)
					if (file.exists(diff_gene_file)){
						base_prefix = gsub("_diff_result.xlsx","",diff_xlsx)
						diff_gene = openxlsx::read.xlsx(diff_gene_file)
						cells_keep <- rownames(single_ob@meta.data[single_ob@meta.data[[diff_group[1]]] %in% c(diff_group[2],diff_group[3]), ])
						sub_ob2 <- subset(single_ob, cells = cells_keep)
						if(is.factor(sub_ob2@meta.data[,diff_group[1]] )){
							sub_ob2@meta.data[,diff_group[1]] = droplevels(sub_ob2@meta.data[,diff_group[1]])
						}else{
							sub_ob2@meta.data[,diff_group[1]] = factor(sub_ob2@meta.data[,diff_group[1]])
						}
						suppressWarnings(suppressMessages(box_dir = file.path(output_dir,"boxplot")))
						dir.create(box_dir,recursive = T)
						box_fecth(diff_gene,sub_ob2,diff_group[1],box_dir, base_prefix)
						sub_ob2=downsample(sub_ob2,8000,diff_group[1])
						suppressWarnings(suppressMessages(diff_heatmap(diff_gene,sub_ob2,colors,diff_group[1],heatmap_out,base_prefix)))
					}
				}
			}
		Anal_time = Sys.time() - start_time
		print(paste("差异分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
		quit()
}
#-------------------------------------sub GSEA------------------------------------------------
if ("GSEA" %in% args) {
	quiet_library("future.apply")
	quiet_library("openxlsx")
  species_mapping <- list(
  "human" = "Homo_sapiens", # human
  "mouse" = "Mus_musculus", # mouse
  "rat" = "Rattus_norvegicus", # rat
  "dme" = "fruit_fly", # fruit_fly
  "dre" = "zebrafish_danio_rerio", # zebrafish
  "ath" = "Oryza_sativa_Japonica_Group", # Arabidopsis
  "sce" = "Saccharomyces_cerevisiae_S288C", # yeast
  "cel" = "Caenorhabditis_elegans", # C.elegans
  "bta" = "Bos_taurus", # Bovine
  "mcc" = "macaca_mulatta", # monkey
  "cfa" = "canis_lupus_familiaris", # dog
  "ssc" = "sus_scrofa", # pig
  "gga" = "Gallus_gallus" # chicken
	)
	species_package <- species_mapping[[opt$species]]
	if(is.null(species_package)){message("species not found, please check your input");quit()}
	start_time = Sys.time()
	print("GSEA分析开始")
	if (endsWith(opt$file, ".xlsx")){
		quiet_library("openxlsx")
		gene_file = read.xlsx(opt$file, sheet = 1)
	}else{
		gene_file = read.delim(opt$file, header = T, stringsAsFactors = F)
	}
	readme_title = createStyle(
					fontSize = 12, textDecoration = "bold", 
					halign = "center", fgFill = "#D9D9D9", 
					border = "bottom"
	)
	readme_Kegg <- list(
    "ID" = "KEGG编号",
    "Description" = "KEGG描述信息",
    "setSize" = "每次置换检验的数据集大小",
    "enrichmentScore" = "富集分数",
    "NES" = "归一化之后的富集分数",
		"pvalue" = "P值",
    "p.adjust" = "校正P值（BH方法校正后的P值）",
    "qvalue" = "Q值（Q方法校正后的P值）",
    "rank" = "核心基因的最大排序值",
    "leading_edge" = "tags表示核心基因占该基因集基因总数的比例，list表示核心基因占所有基因总数的比例，singal通过公式tag*(1-list)*(N/(N-Nh))计算得到,N代表该基因集下的基因总数，Nh代表核心基因数",
    "core_enrichment" = "富集的核心基因的基因ID，用“/”分割",
		"geneSymbol" = "富集的核心基因的基因名，用“/”分割"
  )
  readme_Go <- list(
    "ID" = "GO编号",
    "Description" = "GO描述信息",
    "setSize" = "每次置换检验的数据集大小",
    "enrichmentScore" = "富集分数",
    "NES" = "归一化之后的富集分数",
		"pvalue" = "P值",
    "p.adjust" = "校正P值（BH方法校正后的P值）",
    "qvalue" = "Q值（Q方法校正后的P值）",
    "rank" = "核心基因的最大排序值",
    "leading_edge" = "tags表示核心基因占该基因集基因总数的比例，list表示核心基因占所有基因总数的比例，singal通过公式tag*(1-list)*(N/(N-Nh))计算得到,N代表该基因集下的基因总数，Nh代表核心基因数",
    "core_enrichment" = "富集的核心基因的基因ID，用“/”分割",
		" Subontologies" = "GO分类",
		"geneSymbol" = "富集的核心基因的基因名，用“/”分割"
  )
  if ("cluster" %in% colnames(gene_file)) {
		results <- future_lapply(unique(as.character(gene_file$cluster)), function(clust_num) {
			tryCatch({
				# 创建临时和输出目录
				temp <- file.path(output, "temp", clust_num)
				dir.create(temp, recursive = TRUE, showWarnings = FALSE)
				cluster_dir <- file.path(output, paste0(clust_num, "_GSEA"))
				dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
				# 筛选该 cluster 的基因
				temp_maker <- gene_file[gene_file$cluster == clust_num, ]
				temp_maker <- temp_maker %>% select(gene,avg_log2FC) %>% arrange(desc(avg_log2FC))
				temp_file <- file.path(temp, paste0(clust_num, "_marker_gene.txt"))
				write.table(temp_maker, temp_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
				
				# 构造 Perl 命令，并重定向所有输出到 /dev/null
				cmd <- glue::glue(
					"perl /home/genesky/pipeline/tools_script/gene_enrichment/latest/gene_enrichment.pl ",
					"-g {temp_file} -a gsea -t {temp} -o {cluster_dir} ",
					"-s {species_package} -mode local > /dev/null 2>&1"
				)
        system(cmd, intern = FALSE)
				EnrichDat_file <- file.path(cluster_dir, "enrichGO.txt")
				if (file.exists(EnrichDat_file)) {
					EnrichDat <- read.delim(EnrichDat_file, header = TRUE, stringsAsFactors = FALSE)
					wb <- createWorkbook()
					sheet <- addWorksheet(wb,"GO Summary")
					writeData(wb, sheet, EnrichDat, startCol = 1, startRow = 1, colNames = TRUE)
					addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:13, gridExpand = TRUE)
					setColWidths(wb, sheet, cols = 1:13, widths = 15)
					setColWidths(wb, sheet, cols = 2, widths = 30)
					sheet <- addWorksheet(wb, "README")
					writeData(wb, sheet, t(as.data.frame(c("标题", "说明"))), 
										startCol = 1, startRow = 1, colNames = FALSE)
					addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:2, gridExpand = TRUE)
					setColWidths(wb, sheet, cols = 1, widths = 15)
					setColWidths(wb, sheet, cols = 2, widths = 40)
					# addStyle(wb, sheet, createStyle(wrapText = TRUE), .
					# 				rows = 2:(nrow(readme_data)+1), cols = 2, gridExpand = TRUE)
					readme_list_filter <- rownames_to_column(as.data.frame(t(as.data.frame(readme_Go))), var = "tittle")
					writeData(wb, sheet, readme_list_filter, startCol = 1, startRow = 2, colNames = FALSE)
					saveWorkbook(wb, file = file.path(cluster_dir, "GSEA_GO_Summary.xlsx"), overwrite = TRUE)
					file.remove(EnrichDat_file)
				}
				EnrichDat_file <- file.path(cluster_dir, "enrichKEGG.txt")
				if (file.exists(EnrichDat_file)) {
					EnrichDat <- read.delim(EnrichDat_file, header = TRUE, stringsAsFactors = FALSE)
					wb <- createWorkbook()
					sheet <- addWorksheet(wb,"KEGG Summary")
					writeData(wb, sheet, EnrichDat, startCol = 1, startRow = 1, colNames = TRUE)
					addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:15, gridExpand = TRUE)
					setColWidths(wb, sheet, cols = 1:15, widths = 15)
					setColWidths(wb, sheet, cols = 2, widths = 30)
					sheet <- addWorksheet(wb, "README")
					writeData(wb, sheet, t(as.data.frame(c("标题", "说明"))), 
										startCol = 1, startRow = 1, colNames = FALSE)
					addStyle(wb, sheet, style = readme_title, rows = 1, cols = 1:2, gridExpand = TRUE)
					setColWidths(wb, sheet, cols = 1, widths = 15)
					setColWidths(wb, sheet, cols = 2, widths = 40)
					# addStyle(wb, sheet, createStyle(wrapText = TRUE), .
					# 				rows = 2:(nrow(readme_data)+1), cols = 2, gridExpand = TRUE)
					readme_list_filter <- rownames_to_column(as.data.frame(t(as.data.frame(readme_Kegg))), var = "tittle")
					writeData(wb, sheet, readme_list_filter, startCol = 1, startRow = 2, colNames = FALSE)
					saveWorkbook(wb, file = file.path(cluster_dir, "GSEA_KEGG_Summary.xlsx"), overwrite = TRUE)
					file.remove(EnrichDat_file)
				}
				return(TRUE)  # 标记成功
			}, error = function(e) {
				message("ERROR in cluster ", clust_num, ": ", e$message)
				return(FALSE) # 标记失败
			})
		}, future.seed = TRUE)
	}
  #删除临时文件
	# if (dir.exists(file.path(output, "temp"))) {
  # # 删除里面的文件
  # #file.remove(list.files(file.path(output, "temp"), full.names = TRUE))
  # # 删除目录本身
  # unlink(file.path(output, "temp"), recursive = TRUE, force = TRUE)
	# }
	Anal_time = Sys.time() - start_time
	print(paste("GSEA 分析已结束，耗时：", round(as.numeric(Anal_time, units = "mins"),2),"分钟"))
	quit()
}