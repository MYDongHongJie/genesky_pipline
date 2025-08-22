# 导入 -> 系统 package   
use warnings;
use strict;

die "Usage perl $0 inputfile" if(@ARGV != 1);

###################################################################### 主程序

my $file = shift @ARGV;
die "[ERROR] 文件不存在, $file\n" if(not -e $file);
my %hash;
open FILE, $file;
while(my $line = <FILE>){
	$line =~ s/^\s+//g;
	$line =~ s/\s+$//g;
	if($line =~ /(.+)\s*=\s*(.+)/){
		my $a = $1;
		my $b = $2;
		$a =~ s/\s+$//;
		$b =~ s/^\s+//;
		$hash{$a} = $b;
	}
}
close FILE;
die "[ERROR] 文件 $file 中 Ref 值未定义\n" if(not exists $hash{"Ref"});
die "[ERROR] 文件 $file 中 sequencing 值未定义\n" if(not exists $hash{"sequencing"});
my $border = "border:1px solid #000000;padding:10px;";
my $left = "text-align:left";
my $bg = "background-color:#d6e3bc";
print "
	<table border=0 style='margin:auto;'>
		<tr height='40px'><td colspan=4 style='text-align:center;$border'>项目参考信息</td></tr>
		<tr height='40px'>
			<td width='150px' style='$left;$border;$bg'>参考数据库</td><td colspan=3 style='$left;$border'>".$hash{"Ref"}."</td>
		</tr>
		<tr height='40px'>
			<td style='$left;$border;$bg'>项目类型</td><td style='$left;$border'>10x Genomics 3’单细胞测序</td>
			<td style='$left;$border;$bg'>测序信息</td><td style='$left;$border'>".$hash{"sequencing"}."</td>
		</tr>
		<tr height='40px'><td colspan=4 style='text-align:center;$border'>联系人信息</td></tr>
		
		<tr height='40px'>
			<td style='$left;$border;$bg' rowspan=2>项目售后支持</td><td style='$left;$border' rowspan=2>技术支持</td>
			<td style='$left;$border;'>电话</td><td style='$left;$border'>021-50802060-135</td>
		</tr><tr height='40px'><td style='$left;$border;'>邮箱</td><td style='$left;$border'>gc_service\@cwbio.com.cn</td></tr>
		
		<tr height='40px'><td colspan=4 style='text-align:center;$border'>项目审批</td></tr>
		<tr><td colspan=4 style='text-align:center;$border'>
		<br><br>
		确认报告显示的内容与合同要求完全一致，同意项目结束，批准本项目总结报告发送。
		<br><br><br><br><br><br>
		<div style='float:right'>签名：<img src='html/project/sign.png'><br><br>".time_to_datetime(time)."</div>
		</td></tr>
	</table>
";

###################################################################### 子函数
sub time_to_datetime{
	my $time = shift;
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
	$year += 1900;
	$mon += 1;
	return "$year-$mon-$day";
}
