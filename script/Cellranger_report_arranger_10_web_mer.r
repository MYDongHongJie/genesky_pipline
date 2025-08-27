suppressWarnings(suppressMessages({
  library(readr)
  library(openxlsx)
  library(fs)
	library(optparse)
}))
option_list <- list(
  make_option(c("--webs"), type = "character", help = "多个 web_summary.html 路径（逗号分隔）", metavar = "FILES"),
  make_option(c("--metrics"), type = "character", help = "多个 metrics_summary.csv 路径（逗号分隔）", metavar = "FILES"),
  make_option(c("--samples"), type = "character", help = "多个样本名（逗号分隔，需与文件一一对应）", metavar = "NAMES"),
  make_option(c("-o", "--output_dir"), type = "character", help = "输出目录", metavar = "DIR"),
	make_option(c("-p", "--platforms"), type = "character",default = "10X", help = "建库平台，10X or huadaC4")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ========================
#     参数解析与检查
# ========================
web_files <- unlist(strsplit(opt$webs, " "))
metrics_files <- unlist(strsplit(opt$metrics, " "))
sample_names <- unlist(strsplit(opt$samples, " "))
output_dir <- opt$output_dir

header_style <- createStyle(
  fontSize = 11,
  fontColour = "#000000",
  halign = "center",
  fgFill = "#DCE6F1",  # 浅蓝背景
  textDecoration = "bold"
)

if (!(length(web_files) == length(metrics_files) && length(metrics_files) == length(sample_names))) {
  stop("❌ 错误：webs、metrics 和 samples 参数的数量不一致！")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
html_output <- file.path(output_dir, "CellRanger_Html")
dir_create(html_output)

# 初始化 Excel
xlsx_file <- file.path(output_dir, "CellRanger_Summary.xlsx")
wb <- createWorkbook()

# ======= 写入 Summary sheet =======
addWorksheet(wb, "CellRanger Summary")
writeData(wb, "CellRanger Summary", "Metrics\\Sample", startCol = 1, startRow = 1)

row_names <- NULL

for (i in seq_len(length(sample_names))) {
  html_file <- web_files[i]
  metrics_file <- metrics_files[i]

  # 提取样本名
  sample <- sample_names[i]
  writeData(wb, "CellRanger Summary", sample, startCol = i + 1, startRow = 1)

  # 拷贝 HTML 文件
  dst <- file.path(html_output, paste0(sample, "_web.html"))
  file_copy(html_file, dst, overwrite = TRUE)

  # 写入 metrics
	if (opt$platforms == "10X") {
  	df <- read_csv(metrics_file, col_names = TRUE, show_col_types = FALSE)
	} else if (opt$platforms == "huadaC4") {
		df <- read.delim(metrics_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	}
  if (is.null(row_names)) {
    row_names <- names(df)
    writeData(wb, "CellRanger Summary", row_names, startCol = 1, startRow = 2)
  }
  writeData(wb, "CellRanger Summary", as.character(df[1, ]), startCol = i + 1, startRow = 2)

}
setColWidths(wb, "CellRanger Summary", cols = 1, widths = 40)
setColWidths(wb, "CellRanger Summary", cols = 2, widths = 20)
# ======= 写入 README sheet =======
addWorksheet(wb, "README")
readme_10X <- list(
  "Estimated Number of Cells" = "样本细胞数",
  "Mean Reads per Cell" = "细胞平均 Reads 数",
  "Median Genes per Cell" = "细胞中基因数量的中位数",
  "Number of Reads" = "样本 Reads 数",
  "Valid Barcodes" = "有效 Barcodes 比例",
  "Sequencing Saturation" = "测序饱和度",
  "Q30 Bases in Barcode" = "Barcode 序列 Q30 比例",
  "Q30 Bases in RNA Read" = "RNA 序列 Q30 比例",
  "Q30 Bases in UMI" = "UMI 序列 Q30 比例",
  "Reads Mapped to Genome" = "Reads 比对到参考基因组比例",
  "Reads Mapped Confidently to Genome" = "Reads 高质量比对到参考基因组的比例",
  "Reads Mapped Confidently to Intergenic Regions" = "Reads 高质量比对到基因之间区域的比例",
  "Reads Mapped Confidently to Intronic Regions" = "Reads 高质量比对到内含子的比例",
  "Reads Mapped Confidently to Exonic Regions" = "Reads 高质量比对到外显子的比例",
  "Reads Mapped Confidently to Transcriptome" = "Reads 高质量比对到转录本的比例",
  "Reads Mapped Antisense to Gene" = "Reads 比对到基因反义链的比例",
  "Fraction Reads in Cells" = "包含在细胞内的 Reads 比例",
  "Total Genes Detected" = "检测到的基因数量",
  "Median UMI Counts per Cell" = "细胞中 UMI 数量的中位数"
)
readme_huadaC4 <- list(
  "SampleName" = "样本名",
  "species" = "物种信息",
  "Estimated.number.of.cell" = "样本细胞数",
  "Mean.reads.per.cell" = "细胞平均 Reads 数",
  "Mean.UMI.count.per.cell" = "细胞中 UMI 数量的平均数",
  "Median.UMI.counts.per.cell" = "细胞中 UMI 数量的中位数",
  "Total.genes.detected" = "检测到的基因数量",
  "Mean.genes.per.cell" = "细胞中基因数量的平均数",
  "Median.genes.per.cell" = "细胞中基因数量的中位数",
  "Sequencing.saturation" = "测序饱和度",
  "Fraction.Reads.in.cell" = "包含在细胞内的 Reads 比例",
  "cDNA.Number.of.reads" = "cDNA的测序reads数",
  "cDNA.Reads.pass.QC" = "通过QC的Reads比例",
  "cDNA.Adapter.Reads" = "Reads 中接头的比例",
  "cDNA.Q30.bases.in.reads" = "Q评分≥30的read数比例",
  "index.Number.of.reads" = "index的测序reads数",
  "index.Reads.pass.QC" = "通过QC的index的Reads比例",
  "index.Q30.bases.in.reads" = "Q评分≥30的index的read数比例",
  "Mitochondria.ratio" = "线粒体比率",
	"Reads.mapped.to.genome" = "Reads 比对到参考基因组比例",
	"Reads.mapped.to.exonic.regions" = "Reads 比对到外显子的比例",
	"Reads.mapped.to.intronic.regions" = "Reads 比对到内含子的比例",
	"Reads.mapped.antisense.to.gene" = "Reads 比对到基因反义链的比例",
	"Reads.mapped.to.intergenic.regions" = "Reads 比对到基因之间区域的比例"
)
if (opt$platforms == "10X") {
  readme = readme_10X
}else{
	readme = readme_huadaC4
}


readme_df <- data.frame(
  标题 = names(readme),
  说明 = unname(unlist(readme)),
  stringsAsFactors = FALSE
)

writeData(wb, "README", readme_df, startRow = 1, startCol = 1)

addStyle(wb, sheet = "README", style = header_style, rows = 1, cols = 1:2, gridExpand = TRUE)

addStyle(wb, sheet = "CellRanger Summary", style = header_style, rows = 1, cols = 1:length(sample_names)+1, gridExpand = TRUE)
setColWidths(wb, sheet="README", cols = 1:2, widths = 32)
# 保存 Excel 文件
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
