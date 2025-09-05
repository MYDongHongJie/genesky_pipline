#!/usr/bin/env Rscript

# ========================
# 加载必要的包
#install.packages("extrafont")
# ========================
suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
	library(extrafont)
})

# ========================
# 命令行参数定义
# ========================
option_list <- list(
  make_option(c("--read1"), type="character", help="R1 质量文件路径"),
  make_option(c("--read2"), type="character", help="R2 质量文件路径"),
  make_option(c("-s", "--sample"), type="character", help="输出文件前缀"),
	make_option(c("-f", "--final"), type="logical", help="输出文件前缀"),
  make_option(c("-o", "--output"), type="character", help="输出文件路径")
  #make_option(c("-t", "--title"), type="character", default="Sample", help="图表标题前缀")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$read1) || is.null(opt$read2) ) {
  print_help(opt_parser)
  stop("请提供 --read1, --read2 和 --output_prefix 参数", call.=FALSE)
}
# font_import(paths = "/home/donghj/snakemake/fonts", pattern = "Times",prompt = FALSE)
# loadfonts()
output = opt$output
if(!file.exists(output)){
	dir.create(output, recursive = TRUE)
}
read1 = unlist(strsplit(opt$read1," "))
read2 = unlist(strsplit(opt$read2," "))
sample = unlist(strsplit(opt$sample," "))
if(opt$final){sample = paste0(sample,"_final")}

if (!(length(read1) == length(read2) && length(read2) == length(sample))) {
  stop("❌ 错误：read1、read2 和 samples 参数的数量不一致！")
}

# ========================
# 数据读取和合并
# ========================
for (i in seq_len(length(sample))) {
  read_file1 = read1[i]
  read_file2 = read2[i]
	r1 <- read.table(read_file1, header=FALSE, sep="\t", skip=3)
	r2 <- read.table(read_file2, header=FALSE, sep="\t", skip=3)
	data <- rbind(r1, r2)

	# ========================
	# 画 Base Composition 图
	# ========================
	atcgn <- data[,3:7] / data[[2]] * 100
	colnames(atcgn) <- c('A', 'C', 'G', 'T', 'N')
	atcgn2 <- stack(atcgn)
	atcgn2$group <- rep(1:nrow(atcgn), ncol(atcgn))
	atcgn2$ind <- factor(atcgn2$ind, level=unique(atcgn2$ind))
	colnames(atcgn2) <- c('values', 'type', 'pos')
  
	result_file_pdf <- file.path(output, paste0(sample[i], "_base_content.pdf"))
	result_file_png <- file.path(output, paste0(sample[i], "_base_content.png"))

	p <- ggplot(atcgn2, aes(x=pos, y=values, group=type)) +
  geom_line(aes(colour=type)) +
  labs(x = 'Position along reads', y = 'Percent of bases',
       title = paste(sample[i], "Raw Base content along reads")) +
  geom_vline(xintercept = nrow(data) / 2, linetype=2, colour='skyblue') +
  scale_x_continuous(expand=c(0,0), limits=c(0, nrow(data)), breaks=seq(0, nrow(data), nrow(data)/2)) +
  theme(plot.title = element_text(hjust = 0.5, family = "Times New Roman"),
        text = element_text(family = "Times New Roman"))

	# 保存 PDF
	ggsave(filename = result_file_pdf, plot = p, width = 8, height = 7, units = "in")
	# 保存 PNG
	ggsave(filename = result_file_png, plot = p, width = 8, height = 7, units = "in", dpi = 1000)

# ========================
# 画 Error Rate 图
# ========================
	err <- data[, c(1, 8)]
	colnames(err) <- c('position', 'error')
	err$position <- 1:nrow(err)
	err$error <- 10 ^ ((err$error - 0) / -10) * 100
	# 设置文件名
	result_file_pdf <- file.path(output, paste0(sample[i], "_error_rate.pdf"))
	result_file_png <- file.path(output, paste0(sample[i], "_error_rate.png"))

	# 构建 ggplot 图对象，字体设置在 theme(text = element_text(family = "Times"))
	p <- ggplot(err, aes(x = position, y = error)) +
	geom_bar(fill = '#00FF7F', color = '#32CD32', stat = 'identity') +
	ylim(0, 1) +
	labs(
		x = 'Position along reads',
		y = 'Error rate% (limit 1%)',
		title = paste(sample[i], "Raw Error rate distribution along reads")
	) +
	geom_vline(xintercept = nrow(err) / 2, linetype = 2, colour = 'skyblue') +
	scale_x_continuous(
		expand = c(0, 0),
		limits = c(0, nrow(err)),
		breaks = seq(0, nrow(err), nrow(err) / 2)
	) +
	theme(
		plot.title = element_text(hjust = 0.5),
		text = element_text(family = "Times")  # 指定字体
	)

	# 保存 PDF 和 PNG 图像
	ggsave(result_file_pdf, plot = p, width = 8, height = 7, device = cairo_pdf)
	ggsave(result_file_png, plot = p, width = 8, height = 7, dpi = 300, device = "png")
}
