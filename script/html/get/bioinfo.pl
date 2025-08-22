# 导入 -> 系统 package   
use warnings;
use strict;

die "Usage perl $0 inputdir group" if(@ARGV != 2);

###################################################################### 主程序

my $dir = shift @ARGV;
my $group = shift @ARGV;

die "[ERROR] 路径不存在, $dir\n" if(not -e $dir);

my @array = (
	["序列质控", [
		["数据质量评估",	"SCGES_*/*.Quality_Statistics/Quality_Summary.xlsx"],
		["质量数据绘图",	"SCGES_*/*.Quality_Statistics/Quality_Images/*_base_content.pdf"],
	]],
	["CellRanger分析", [
		["细胞Barcode分析",	"SCGES_*/*.CellRanger/CellRanger_Html/*.html"],
		["分析结果统计",	"SCGES_*/*.CellRanger/CellRanger_Summary.xlsx"],
		["合并样本",		"SCGES_*/*.CellRanger/filtered_feature_bc_matrix/web_summary.html"],
	]],
	["细胞质控", [
		["细胞表达水平质控",	"SCGES_*/*.Cell_QC/Cell_Quality_Summary.xlsx"],
		["差异基因选取",		"SCGES_*/*.Cell_QC/Cell_Quality_Images/variable_feature.pdf"],
	]],
	["双细胞检测", [
		["双细胞检测",		"SCGES_*/*.Doublet_Cell/Doublet_Cell_Summary.xlsx"],
		["双细胞分布绘图",	"SCGES_*/*.Doublet_Cell/doublet_umap.pdf"],
	]],
	["降维聚类", [
		["PCA特征选取",		"SCGES_*/*.Dimensional_Reduction/jackstraw.pdf"],
		["UMAP降维聚类",	"SCGES_*/*.Dimensional_Reduction/clusters_umap.pdf"],
	]],
	["Cluster标记基因分析", [
		["特征基因分析",	"SCGES_*/*.Cluster_Marker_Gene/Marker_Gene_Summary.xlsx"],
		["表达分布绘图",	"SCGES_*/*.Cluster_Marker_Gene/Gene_Images/*_feature.png"],
	]],
	["细胞类型注释", [
		["细胞类型注释",	"SCGES_*/*.Cell_Type/Cell_Type_Summary.xlsx"],
	]],
	["分组重聚类", [
		["分组重聚类",	"SCGES_*/*.Re_Cluster/*/clusters_umap.pdf"],
	]],
	["细胞类型重注释", [
		["细胞类型重注释",	"SCGES_*/Statistical_Analysis/$group/*.Re_Cell_Type/Re_Cell_Type_Summary.xlsx"],
	]],
	["Cell_Type标记基因分析", [
		["特征基因分析",	"SCGES_*/Statistical_Analysis/$group/*.Cell_Type_Marker_Gene/Marker_Gene_Summary.xlsx"],
		["表达分布绘图",	"SCGES_*/Statistical_Analysis/$group/*.Cell_Type_Marker_Gene/Gene_Images/*_feature.png"],
	]],
	["差异表达分析", [
		["Wilcoxon检验",		"SCGES_*/Statistical_Analysis/$group/*.Differential_Gene_Expression/*/wilcox/*/Differential_Gene_Summary.xlsx"],
		["Kruskal-Wallis检验",	"SCGES_*/Statistical_Analysis/$group/*.Differential_Gene_Expression/*/kruskal/Differential_Gene_Summary.xlsx"],
		["ANOVA方差分析",		"SCGES_*/Statistical_Analysis/$group/*.Differential_Gene_Expression/*/anova/Differential_Gene_Summary.xlsx"],
	]],
	["基因富集分析", [
		["富集分析",		"SCGES_*/Statistical_Analysis/$group/*.Gene_Enrich/*/*/Gene_Enrich_Summary.xlsx"],
		["KEGG通路图绘制",	"SCGES_*/Statistical_Analysis/$group/*.Gene_Enrich/*/*/KEGG_Images"],
	]],
	["基因富集GSEA分析", [
		["富集分析",		"SCGES_*/Statistical_Analysis/$group/*.Gene_GSEA/*/*/Gene_GSEA_Summary.xlsx"],
		["KEGG通路图绘制",	"SCGES_*/Statistical_Analysis/$group/*.Gene_GSEA/*/*/KEGG_Images"],
		["GSEA统计绘图",	"SCGES_*/Statistical_Analysis/$group/*.Gene_GSEA/*/*/GSEA_Images"],
	]],
	["拟时轨迹Monocle分析", [
		["Monocle分析",		"SCGES_*/Statistical_Analysis/$group/*.Trajectory_Monocle/*/Trajectory_Monocle_Summary.xlsx"],
		["基因表达量分布",	"SCGES_*/Statistical_Analysis/$group/*.Trajectory_Monocle/*/Gene_Images/*.pdf"],
	]],
	["拷贝数InferCNV分析", [
		["InferCNV分析",	"SCGES_*/Statistical_Analysis/$group/*.CNV_InferCNV/CNV_InferCNV_Summary.xlsx"],
	]],
	["细胞通讯CellChat分析", [
		["CellChat分析",	"SCGES_*/Statistical_Analysis/$group/*.Cell_Communication_CellChat/CellChat_*_Summary.xlsx"],
		["通讯网络分析",	"SCGES_*/Statistical_Analysis/$group/*.Cell_Communication_CellChat/*/Network_Images/*.pdf"],
		["贡献度分析",		"SCGES_*/Statistical_Analysis/$group/*.Cell_Communication_CellChat/*/Contribution_Images/*.pdf"],
		["模式分析",		"SCGES_*/Statistical_Analysis/$group/*.Cell_Communication_CellChat/*/Pattern_Images/*.pdf"],
		["相似性分析",		"SCGES_*/Statistical_Analysis/$group/*.Cell_Communication_CellChat/*/Similarity_Images/*.pdf"],
		["基因表达分析",	"SCGES_*/Statistical_Analysis/$group/*.Cell_Communication_CellChat/*/Gene_Violin/*.pdf"],
	]],
	["差异细胞通讯CellChat分析", [
		["通讯网络分析",	"SCGES_*/Statistical_Analysis/$group/*.Differential_Cell_Communication_CellChat/*/Network_Images/*.pdf"],
		["贡献度分析",		"SCGES_*/Statistical_Analysis/$group/*.Differential_Cell_Communication_CellChat/*/Contribution_Images/*.pdf"],
		["相似性分析",		"SCGES_*/Statistical_Analysis/$group/*.Differential_Cell_Communication_CellChat/*/Similarity_Images/*.pdf"],
		["基因表达分析",	"SCGES_*/Statistical_Analysis/$group/*.Differential_Cell_Communication_CellChat/*/Gene_Violin/*.pdf"],
	]],
	["RNA速率scVelo分析", [
		["scVelo分析",		"SCGES_*/Statistical_Analysis/$group/*.RNA_Velocity_scVelo/RNA_Velocity_Summary.xlsx"],
		["基因速率绘图",	"SCGES_*/Statistical_Analysis/$group/*.RNA_Velocity_scVelo/Gene_Images/*.pdf"],
	]],
	["转录调控SCENIC分析", [
		["SCENIC分析",	"SCGES_*/Statistical_Analysis/$group/*.Transcription_Factors_SCENIC/*/Transcription_Factors_SCENIC_Summary.xlsx"],
		["Rank绘图",	"SCGES_*/Statistical_Analysis/$group/*.Transcription_Factors_SCENIC/*/RankPlot_Images/*.pdf"],
		["AUC绘图",		"SCGES_*/Statistical_Analysis/$group/*.Transcription_Factors_SCENIC/*/AUC_Images/*.pdf"],
		["RSS绘图",		"SCGES_*/Statistical_Analysis/$group/*.Transcription_Factors_SCENIC/*/RSS_Images/*.pdf"],
	]],
);

print "<table border=0 style='margin:auto;text-align:left'>";
foreach my $key(@array){
	print "<tr><td colspan=4 style='border:1px #000000 solid;padding:8px;background-color:#d6e3bc;text-align:center'>".$key->[0]."</td></tr>";
	my @list = @{$key->[1]}; push @list, "" if(@list%2 != 0);
	for my $i(0..@list-1){
		print "<tr>" if($i%2 == 0);
		my ($v1, $v2) = ("", "");
		if($list[$i] ne ""){
			my @ls = `ls $dir/$list[$i]->[1] 2>/dev/null`;
			$v1 = $list[$i]->[0];
			$v2 = (@ls > 0)?"<strong style='color:#01d800'>√</strong>":"<strong style='color:#d80613'>×</strong>";
		}
		print "<td style='border:1px #000000 solid;padding:8px;padding-right:300px;'>$v1</td><td style='border:1px #000000 solid;padding:8px'>$v2</td>";
		print "</tr>" if($i%2 == 1);
	}
}
print "</table>";



