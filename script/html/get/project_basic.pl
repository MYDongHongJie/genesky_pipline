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
die "[ERROR] 文件 $file 中 id 值未定义\n" if(not exists $hash{"id"});
die "[ERROR] 文件 $file 中 project 值未定义\n" if(not exists $hash{"project"});
die "[ERROR] 文件 $file 中 department 值未定义\n" if(not exists $hash{"department"});

print "
	<table border=0 style='margin:auto;font-size:24px;'>
		<tr height='40px'><td width='130px' style='text-align:left'>项目编号</td><td style='text-align:left'>".$hash{"id"}."</td></tr>
		<tr height='40px'><td style='text-align:left'>项目名称</td><td style='text-align:left'>".$hash{"project"}."</td></tr>
		<tr height='40px'><td style='text-align:left'>客户单位</td><td style='text-align:left'>".$hash{"department"}."</td></tr>
	</table>
";

