suppressPackageStartupMessages({
  library(fs)
  library(optparse)
	library(openxlsx)
	library(dplyr)
})

# ========================
# 命令行参数定义
# ========================
option_list <- list(
  make_option(c("--r1_dir_raw"), type="character", help="R1 质量文件路径"),
  make_option(c("--r2_dir_raw"), type="character", help="R2 质量文件路径"),
  make_option(c("--r1_dir_final"), type="character", help="R1最终质量文件路径"),
	make_option(c("--r2_dir_final"), type="character", help="R2最终质量文件路径"),
  make_option(c("-s", "--sample"), type="character", help="输出文件前缀"),
	make_option(c("-t", "--stat"), type="character", help="stat的路径"),
  make_option(c("-o", "--output"), type="character", help="输出文件路径")
  #make_option(c("-t", "--title"), type="character", default="Sample", help="图表标题前缀")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
r1_dir_raw = unlist(strsplit(opt$r1_dir_raw," "))
r2_dir_raw = unlist(strsplit(opt$r2_dir_raw," "))
r1_dir_final = unlist(strsplit(opt$r1_dir_final," "))
r2_dir_final = unlist(strsplit(opt$r2_dir_final," "))
stat = unlist(strsplit(opt$stat," "))
sample = unlist(strsplit(opt$sample," "))

output = opt$output
Quality_Images = paste0(output,"/","Quality_Images")
if(!file.exists(Quality_Images)){dir.create(Quality_Images,recursive = T)}



if (!(length(r1_dir_raw) == length(r2_dir_raw) && length(sample) == length(r1_dir_raw))) {
  stop("❌ 错误：r1_dir_raw、r2_dir_raw 和 sample 参数的数量不一致！")
}
if (!(length(r1_dir_final) == length(r2_dir_final) && length(sample) == length(r2_dir_final))) {
  stop("❌ 错误：r1_dir_final、r2_dir_final 和 sample 参数的数量不一致！")
}

for (i in seq_len(length(sample))) {
	raw_file1 = paste0(r1_dir_raw[i],"/","Images/per_base_quality.png")
	cp_file1 = paste0(Quality_Images,"/",sample[i],"_r1_base_quality.png")
	file_copy(raw_file1, cp_file1, overwrite = TRUE)

	raw_file2 = paste0(r2_dir_raw[i],"/","Images/per_base_quality.png")
	cp_file2 = paste0(Quality_Images,"/",sample[i],"_r2_base_quality.png")
	file_copy(raw_file2, cp_file2, overwrite = TRUE)

	raw_file3 = paste0(r1_dir_final[i],"/","Images/per_base_quality.png")
	cp_file3 = paste0(Quality_Images,"/",sample[i],"_final_r1_base_quality.png")
	file_copy(raw_file3, cp_file3, overwrite = TRUE)

	raw_file4 = paste0(r2_dir_final[i],"/","Images/per_base_quality.png")
	cp_file4 = paste0(Quality_Images,"/",sample[i],"_final_r2_base_quality.png")
	file_copy(raw_file4, cp_file4, overwrite = TRUE)
}


#------------------------------------------------------
# 定义样式
create_excel_styles <- function(wb) {
  list(
    title = createStyle(
      textDecoration = "bold", fontColour = "white", 
      halign = "center", fgFill = "#4F81BD", 
      border = "bottom", fontSize = 11
    ),
    header = createStyle(textDecoration = "bold", halign = "center"),
    normal = createStyle(halign = "center"),
    percent = createStyle(numFmt = "0.0%", halign = "center"),
    readme_title = createStyle(
      fontSize = 12, textDecoration = "bold", 
      halign = "center", fgFill = "#D9D9D9", 
      border = "bottom"
    ),
    readme_bold = createStyle(textDecoration = "bold"),
    readme_desc = createStyle(wrapText = TRUE)
  )
}

# 创建数据工作表
create_data_sheet <- function(wb, styles, sample, stat) {
  sheet <- addWorksheet(wb, "Quality Summary")
  
 # 定义样式
	title_style <- createStyle(
		textDecoration = "bold", 
		fontColour = "white", 
		halign = "center",
		fgFill = "#4F81BD", 
		border = "bottom",
		fontSize = 11
	)
  header_style <- createStyle(
  textDecoration = "bold",
  halign = "center"
	)
  # 写入表头
  # 写入第一行标题
	first_row <- c("Sample", "Raw Data", "", "", "", "", "Clean Data", "", "", "", "", "Clean Reads%")
	writeData(wb, "Quality Summary", 
          matrix(first_row, nrow = 1), 
          startCol = 1, startRow = 1, 
          colNames = FALSE)
  # 写入第二行标题
	second_row <- c("Sample", "# reads", "# bases", "Q20", "Q30", "GC%",
									"# reads", "# bases", "Q20", "Q30", "GC%", "Clean Reads%")
	writeData(wb, "Quality Summary", 
						matrix(second_row, nrow = 1), 
						startCol = 1, startRow = 2, 
						colNames = FALSE)
	# 合并单元格
	mergeCells(wb, "Quality Summary", cols = 1, rows = 1:2)  # 合并Sample列
	mergeCells(wb, "Quality Summary", cols = 12, rows = 1:2) # 合并Clean Reads%列
	mergeCells(wb, "Quality Summary", cols = 2:6, rows = 1)   # 合并Raw Data
	mergeCells(wb, "Quality Summary", cols = 7:11, rows = 1)  # 合并Clean Data

 addStyle(wb, "Quality Summary", style = title_style, rows = 1, cols = 1:12, gridExpand = TRUE)
	addStyle(wb, "Quality Summary", style = header_style, rows = 2, cols = 1:12, gridExpand = TRUE)

	# 设置列宽
	setColWidths(wb, "Quality Summary", cols = 1, widths = 15)    # Sample列
	setColWidths(wb, "Quality Summary", cols = 12, widths = 12)   # Clean Reads%列
	setColWidths(wb, "Quality Summary", cols = c(2,7), widths = 10)  # reads列
	setColWidths(wb, "Quality Summary", cols = c(3,8), widths = 15)  # bases列
	setColWidths(wb, "Quality Summary", cols = c(4:6,9:11), widths = 8)
  
  # 写入数据
  for (i in seq_along(sample)) {
    data <- read_stat_file(stat[i], sample[i])
    writeData(wb, sheet, data, startCol = 1, startRow = i + 2, colNames = FALSE)
    
    # 应用数据样式
    addStyle(wb, sheet, style = styles$normal, rows = i + 2, cols = c(1:3,7:8), gridExpand = TRUE)
    addStyle(wb, sheet, style = styles$percent, rows = i + 2, cols = c(4:6,9:12), gridExpand = TRUE)
  }
}

# 读取统计文件并格式化数据
read_stat_file <- function(stat_file, sampleid) {
  stat_data <- read.delim(stat_file, header = FALSE, stringsAsFactors = FALSE)
  values <- setNames(stat_data$V2, stat_data$V1)
  
  data.frame(
    Sample = sampleid,
    Raw_reads = format(as.numeric(values["raw_reads"]), big.mark = ","),
    Raw_bases = format(as.numeric(values["raw_bases"]), scientific = TRUE, digits = 3),
    Raw_Q20 = as.numeric(values["raw_q20"]),
    Raw_Q30 = as.numeric(values["raw_q30"]),
    Raw_GC = as.numeric(values["raw_gc"]),
    Clean_reads = format(as.numeric(values["clean_reads"]), big.mark = ","),
    Clean_bases = format(as.numeric(values["clean_bases"]), big.mark = ","),
    Clean_Q20 = as.numeric(values["clean_q20"]),
    Clean_Q30 = as.numeric(values["clean_q30"]),
    Clean_GC = as.numeric(values["clean_gc"]),
    Clean_Ratio = as.numeric(values["clean_reads"]) / as.numeric(values["raw_reads"]),
    stringsAsFactors = FALSE
  )
}

# 创建README工作表
create_readme_sheet <- function(wb, styles) {
  sheet <- addWorksheet(wb, "README")
  setColWidths(wb, sheet, cols = 1, widths = 30)
  setColWidths(wb, sheet, cols = 2, widths = 100)
  
  # 写入标题
  writeData(wb, sheet, t(as.data.frame(c("标题", "说明"))), startCol = 1, startRow = 1, colNames = FALSE)
  addStyle(wb, sheet, style = styles$readme_title, rows = 1, cols = 1:2, gridExpand = TRUE)
  
  # README 内容以 list 定义
  readme_list <- list(
    "Sample"        = "样本名称",
    "Raw Data"      = "原始数据",
    "Clean Data"    = "质控后数据",
    "# reads"       = "累计 Reads 数量",
    "# bases"       = "累计 Reads 碱基数量",
    "Q20"           = "Q20 碱基所占比例",
    "Q30"           = "Q30 碱基所占比例",
    "GC%"           = "GC 碱基所占比例",
    "Clean Reads%"  = "质控后数据所占比例"
  )
  
  # 转换为数据框
  readme_df <- data.frame(
    标题 = names(readme_list),
    说明 = unname(unlist(readme_list)),
    stringsAsFactors = FALSE
  )
  
  # 写入内容
  writeData(wb, sheet, readme_df, startCol = 1, startRow = 2, colNames = FALSE)
  
  # 应用样式
  data_rows <- 2:(nrow(readme_df) + 1)
  addStyle(wb, sheet, style = styles$readme_bold, rows = data_rows, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet, style = styles$readme_desc, rows = data_rows, cols = 2, gridExpand = TRUE)
}


# 主函数
generate_quality_report <- function(output, sample, stat) {
  wb <- createWorkbook()
  styles <- create_excel_styles(wb)
  
  create_data_sheet(wb, styles, sample, stat)
  create_readme_sheet(wb, styles)
  
  saveWorkbook(wb, file = file.path(output, "Quality_Summary.xlsx"), overwrite = TRUE)
}
generate_quality_report(output, sample, stat)
