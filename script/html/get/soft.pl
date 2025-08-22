use strict;
use warnings;

my @array = (
	["序列质控", [
		["FastQC",	"0.11.9",	"评估每个位点的测序质量和碱基组成，并采用图表呈现"],
		["seqtk",	"1.3-r106",	"统计碱基质量和测序错误率"],
		["fastp",	"0.23.2",	"去除带接头和低质量序列"],
		["R",		"4.2",		"统计和画图"],
	]],
	["CellRanger分析", [
		["CellRanger",	"7.1.0",	"单细胞数据分析"],
	]],
	["细胞质控", [
		["Seurat",	"5.0.3",	"单细胞数据分析"],
	]],
	["双细胞检测", [
		["scDblFinder",	"1.16.0",	"双细胞检测"],
	]],
	["降维聚类", [
		["Seurat",	"5.0.3",	"单细胞数据分析"],
	]],
	["Cluster标记基因分析", [
		["Seurat",	"5.0.3",	"单细胞数据分析"],
		["R",		"4.2",		"统计和画图"],
	]],
	["细胞类型注释", [
		["Seurat",	"5.0.3",	"单细胞数据分析"],
		["SingleR",	"2.4.1",	"细胞类型预测"],
	]],
	["分组重聚类", [
		["Seurat",	"5.0.3",	"单细胞数据分析"],
	]],
	["细胞类型重注释", [
		["Seurat",	"5.0.3",	"单细胞数据分析"],
	]],
	["Cell_Type标记基因分析", [
		["Seurat",	"5.0.3",	"单细胞数据分析"],
		["R",		"4.2",		"统计和画图"],
	]],
	["差异表达分析", [
		["R",	"4.2",	"统计和画图"],
	]],
	["基因富集分析", [
		["clusterProfiler",	"4.6.1",	"基因富集分析"],
	]],
	["基因富集GSEA分析", [
		["clusterProfiler",	"4.6.1",	"基因富集分析"],
	]],
	["拟时轨迹Monocle分析", [
		["Monocle",	"2.30.1",	"细胞排序分析"],
	]],
	["拷贝数InferCNV分析", [
		["InferCNV",	"1.14.2",	"单细胞CNV分析"],
	]],
	["细胞通讯CellChat分析", [
		["CellChat",	"1.6.1",	"细胞通讯分析"],
	]],
	["差异细胞通讯CellChat分析", [
		["CellChat",	"1.6.1",	"细胞通讯分析"],
	]],
	["RNA速率scVelo分析", [
		["velocyto",	"0.17.17",	"RNA速率分析"],
		["scVelo",		"0.2.5",	"RNA速率分析"],
	]],
	["转录调控SCENIC分析", [
		["pySCENIC",	"0.12.1",	"转录调控分析"],
		["rSCENIC",		"1.1.2",	"转录调控分析"],
	]],
);

print "
	<table border=0 style='margin:auto;text-align:left'>
		<tr>
";
foreach("分析内容", "软件", "版本", "说明"){
	print "<td style='border-top:2px #000000 solid;padding:15px'>$_</td>";
}
print "
		</tr>
";
my $count_1 = 0;
my $count_2 = 0;
foreach my $key(@array){
	my @softs = @{$key->[1]};
	$count_1++;
	$count_2++;
	my $color_1 = ($count_1%2 == 1)?"background-color:#d6e3bc;":"";
	my $color_2 = ($count_2%2 == 1)?"background-color:#d6e3bc;":"";
	print "<tr>
		<td rowspan=".(@softs+0)." style='border-top:1px #000000 solid;padding:15px;$color_1'>$key->[0]</td>
		".(join "", map{"<td style='border-top:1px #000000 solid;padding:15px;$color_2'>$softs[0]->[$_]</td>"} 0..2)."
	</tr>";
	for my $i(1..@softs-1){
		$count_2++;
		my $color_2 = ($count_2%2 == 1)?"background-color:#d6e3bc;":"";
		print "<tr>".(join "", map{"<td style='padding:15px;$color_2'>$softs[$i]->[$_]</td>"} 0..2)."</tr>";
	}
}
print "<tr>".(join "", map{"<td style='border-top:2px #000000 solid;'></td>"} 0..3)."</tr>";
print "</table>";
