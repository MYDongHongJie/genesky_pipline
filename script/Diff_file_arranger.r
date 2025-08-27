suppressPackageStartupMessages({
  library(optparse)
	#library(readxl)
	library(dplyr)
	library(purrr)
	library(stringr)
	library(openxlsx)
})
option_list <- list(
  make_option(c("--input"), type="character", help="差异文件路径"),
  make_option(c("--output"), type="character", help="输出路径")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input <- opt$input
output_dir <- opt$output

sub_dirs <- list.dirs(input, recursive = FALSE)

# 提取子文件夹中所有xlsx文件
all_xlsx <- map(sub_dirs, function(subdir) {
  files <- list.files(subdir, pattern = "\\.xlsx$", full.names = TRUE)
  if(length(files) == 0) return(NULL)
  data.frame(file = files, cluster = basename(subdir), stringsAsFactors = FALSE)
}) %>% bind_rows()

# 获取所有文件的前缀（假设前缀是文件名去掉扩展名）
all_xlsx <- all_xlsx %>%
  mutate(prefix = str_remove(basename(file), "\\.xlsx$"))

# 对每个前缀合并
prefix_list <- unique(all_xlsx$prefix)

for(pre in prefix_list){
  files_to_merge <- all_xlsx %>% filter(prefix == pre)
  
  merged_data <- map_df(files_to_merge$file, function(f){
    df <- read.xlsx(f)
		df <- df %>% select(gene,Regulation,avg_log2FC) %>% filter(Regulation != "NO_DEG")
    df$cluster <- files_to_merge$cluster[files_to_merge$file == f]
    return(df)
  })
	pre = gsub("_diff_result","",pre)
  output_dir_pre = file.path(output_dir, pre)
  dir.create(output_dir_pre, showWarnings = FALSE, recursive = TRUE)
  # 写出合并文件
  write.table(merged_data, file.path(output_dir_pre,  "all_gene_merged.xls"),sep="\t",quote = F,row.names = F)
  merged_data_up <- merged_data %>% filter(Regulation == "Up")
  # 写出合并文件
	write.table(merged_data_up, file.path(output_dir_pre, "Up_gene_merged.xls"),sep="\t",quote = F,row.names = F)
  merged_data_down <- merged_data %>% filter(Regulation == "Down")
  # 写出合并文件
	write.table(merged_data_down, file.path(output_dir_pre  , "Down_gene_merged.xls"),sep="\t",quote = F,row.names = F)
}