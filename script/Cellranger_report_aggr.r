suppressPackageStartupMessages({
  library(optparse)
})
option_list <- list(
  make_option(c("--barcodes"), type="character", help="aggr barcodes文件路径"),
  make_option(c("--csv"), type="character", help="aggrter csv文件路径"),
  make_option(c("--outfile"), type="character", help="输出文件")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

make_cell_sample <- function(barcodes, csv, outfile) {
  # 读取CSV文件并创建哈希映射
  csv_data <- read.csv(csv, header = FALSE, stringsAsFactors = FALSE)
  hash <- setNames(csv_data$V1, seq_len(nrow(csv_data)) - 1)  # 从0开始索引
  
  # 处理压缩的barcodes文件
  con <- gzfile(barcodes, "r")
  lines <- readLines(con)
  close(con)
  
  # 处理每一行并写入输出文件
  output <- sapply(lines, function(line) {
    line <- gsub("[\r\n]", "", line)
    i <- strsplit(line, "-")[[1]][2]
    paste(line, hash[[i]], sep = "\t")
  })
  
  writeLines(output, outfile)
}

make_cell_sample(opt$barcodes, opt$csv, opt$outfile)