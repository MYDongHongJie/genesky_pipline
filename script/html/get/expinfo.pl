# 导入 -> 系统 package   
use warnings;
use strict;

die "Usage perl $0 type" if(@ARGV != 1);

###################################################################### 主程序

my $type = shift @ARGV;

my @array;
if($type eq 1){
	@array = (["器材名称", "器材供应商"],
		["Chromium<sup>TM</sup> Single Cell Controller",		"10x Genomics, USA"],
		["Countess<sup>&reg;</sup> II Automated Cell Counter",	"Thermo Fisher Scientific, USA"],
		["Vortex-genie2",										"Scientific Industries, USA"],
		["Qubit Spectrophotometer",								"Thermo Fisher Scientific, USA"],
		["Agilent 2100 bioanalyzer",							"Agilent Technologies, USA"],
		["Illumina PE150 sequencing platform",					"Illumina, USA"],
	);
}
if($type eq 2){
	@array = (["试剂名称", "试剂提供商"],
		["Chromium<sup>TM</sup> Next GEM Single Cell 3' GEM, Library & Gel Bead Kit v3.1",	"10x Genomics, USA"],
		["Chromium<sup>TM</sup> Next GEM Chip G Single Cell Kit", "10x Genomics, USA"],
		["Chromium i7 multiplex kit", "10x Genomics, USA"],
		["Countess<sup>&reg;</sup> II Automated Cell Counting Chamber Slides", "Thermo Fisher Scientific, USA"],
		["Agencourt SPRIselect Reagent Kit", "Beckman Coulter, USA"],
	);
}
die "[ERROR] 输入参数 type 的值未定义, $type\n" if(@array == 0);
print "
	<table border=0 style='margin:auto;text-align:left'>
	<tr>
		<td style='border-top:2px #000000 solid;border-bottom:1px #000000 solid;padding:15px'>".$array[0]->[0]."</td>
		<td style='border-top:2px #000000 solid;border-bottom:1px #000000 solid;padding:15px'>".$array[0]->[1]."</td>
	</tr>
";
for my $i(1..@array-1){
	my $color = ""; $color = "background-color:#d6e3bc;" if($i%2 == 1);
	my $style = ""; $style = "border-bottom:2px #000000 solid;" if($i == @array-1);
	print "
		<tr style='$color'><td style='padding:15px;$style'>".$array[$i]->[0]."</td><td style='padding:15px;$style'>".$array[$i]->[1]."</td></tr>
	";
}
print "</table>";

