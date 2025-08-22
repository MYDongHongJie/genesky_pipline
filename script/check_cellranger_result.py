#!/home/donghj/snakemake/bin/python
#check_cellranger_result.py
#version: 0.1
#pipeline:10X_singlecell3
#author: dhj
#update: 2025/07/31
import os
import csv
import argparse
import sys
import argparse
from datetime import datetime

def term_check_log(sampledict, term, cmp, threshold):
    if term in sampledict.keys():
        term_value = sampledict[term]
    else:
        print("Term: " + term + " is not exists in metrics_summary.csv")
        exit(1)
        
    term_value = term_value.replace(',', '')
    ratio_fix = '%' if term_value.endswith('%') else ''
    term_value = term_value.rstrip("%")
    term_value = float(term_value)
    
    if cmp == 'is larger than':
        check_pass = term_value > threshold
    elif cmp == 'is smaller than':
        check_pass = term_value < threshold
    elif cmp == 'is equals to':
        check_pass = term_value == threshold
    else:
        print("cmp must be one of 'is larger than', 'is smaller than', 'is equals to'")
        exit(1)
        
    if check_pass:
        return "Error in " + sampledict['Category'] + ": " + term + ' ' + cmp + ' ' + str(threshold) + ratio_fix + "\n"
    return ""

def check_log(sampleterm,platform):
    log=''
    if platform=='10X':
        log = log + term_check_log(sampleterm,'Estimated Number of Cells','is smaller than',4000)
        log = log + term_check_log(sampleterm,'Estimated Number of Cells','is larger than',18000)
        log = log + term_check_log(sampleterm,'Median Genes per Cell','is smaller than',800)
        log = log + term_check_log(sampleterm,'Reads Mapped to Genome','is smaller than',80)
        log = log + term_check_log(sampleterm,'Q30 Bases in RNA Read','is smaller than',80)
        log = log + term_check_log(sampleterm,'Q30 Bases in Barcode','is smaller than',85)
        log = log + term_check_log(sampleterm,'Q30 Bases in UMI','is smaller than',85)
        log = log + term_check_log(sampleterm,'Fraction Reads in Cells','is smaller than',70)
    elif platform=='mobi':
        log = log + term_check_log(sampleterm,'Estimated Number of Cells','is smaller than',4000)
        log = log + term_check_log(sampleterm,'Estimated Number of Cells','is larger than',18000)
        log = log + term_check_log(sampleterm,'Median Genes per Cell','is smaller than',800)
        log = log + term_check_log(sampleterm,'Reads Mapped to Genome','is smaller than',80)
        log = log + term_check_log(sampleterm,'Q30 Bases in RNA Read','is smaller than',80)
        log = log + term_check_log(sampleterm,'Q30 Bases in Barcode+UMI','is smaller than',85)
        log = log + term_check_log(sampleterm,'Fraction of Unique Reads in Cells','is smaller than',70)
    elif platform=='huadaC4':
        log = log + term_check_log(sampleterm,'Estimated number of cell','is smaller than',4000)
        log = log + term_check_log(sampleterm,'Estimated number of cell','is larger than',18000)
        log = log + term_check_log(sampleterm,'Median genes per cell','is smaller than',800)
        log = log + term_check_log(sampleterm,'Reads mapped to genome','is smaller than',80)
        log = log + term_check_log(sampleterm,'cDNA Q30 bases in reads','is smaller than',80)
        log = log + term_check_log(sampleterm,'Fraction Reads in cell','is smaller than',70)
    else:
        print("Platform: "+platform+" not support yet")
        exit(1)
    return log

# 单位转换表
unit = {
    "K": 10**3,
    "M": 10**6,
    "G": 10**9,
}

def N2UN(unit_dict, n):
    import re

    # 判断是否为纯数字（整数或小数）
    if not re.match(r'^\d+(\.\d+)?$', str(n)):
        print(f"[WARNING] 无法转换为 KMG 格式, 不是纯数字 {n}")
        return n

    n = float(n)
    un = n
    u = ""

    # 按单位换算值升序排列
    for key in sorted(unit_dict, key=lambda k: unit_dict[k]):
        if n / unit_dict[key] < 1:
            break
        un = n / unit_dict[key]
        u = key

    return f"{int(un)}{u}"

def UN2N(un):
    if isinstance(un, (int, float)):
        return un
    try:
        for k in unit:
            if un.endswith(k):
                return float(un.rstrip(k)) * unit[k]
    except:
        pass
    print(f"[WARNING] 无法转换为纯数字, 不是 KMG 格式 {un}")
    return un

def read_status(file):
    result = {}
    with open(file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                result[parts[0]] = parts[1]
    return result

def read_metrics(file):
    with open(file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            return row
    return {}







def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input',required=True,type=str,nargs='+',help='input metrics_summary.csv file of samples')
    parser.add_argument('--platform',required=True,type=str,help='platform of samples')
    parser.add_argument('--output',default='',type=str,help='output directory')
    parser.add_argument("--sample", required=True,type=str,nargs='+',help='input sample name')
    parser.add_argument("--qcfile", required=True,type=str,nargs='+',help='input fataQC status.txt file of samples')
    #后面的参数
 #   parser.add_argument("--config", "-c", required=True)
    #parser.add_argument("--qc-dir", "-q", required=True)
    #parser.add_argument("--cellranger-dir", "-cd", required=True)
    #parser.add_argument("--aggregate-dir", "-a", required=True)
    parser.add_argument("--output_file", "-o", required=True)
    parser.add_argument("--cell_count", "-cc", type=int, default=2000)
    parser.add_argument("--base_count_gt", "-bcgt", default="120G")
    parser.add_argument("--base_count_le", "-bcle", default="30G")
    args = parser.parse_args()
    
    if(args.platform=='huadaC4'):
        delimiter="\t"
    else:
        delimiter=","
    check_file=open(args.output,'w')
    for summary in args.input:
        file=open(summary,'r')
        paras,line2=csv.reader(file, delimiter=delimiter)
        sampledict=dict(zip(paras, line2))
        file.close()
        sampledict['Category']=summary.split('/')[-2]
        check_file.write(check_log(sampledict,args.platform))
    check_file.close()
    #第二个文件夹
    failed = []
    status_dic = dict(zip(args.sample, args.qcfile))
    metrics_dic = dict(zip(args.sample, args.input))
    with open(args.output_file, 'w') as out:
        out.write("Sample\tCell_Count\tBase_Count\tQC\n")
        for sample in args.sample:
            metrics = read_metrics(metrics_dic[sample])
            cell_count = int(metrics.get("Estimated Number of Cells", "0").replace(",", ""))

            min_bc = args.base_count_gt if cell_count >= args.cell_count else args.base_count_le
            raw_bases = read_status(status_dic[sample]).get("raw_bases", 0)
            raw_bases = float(raw_bases)

            qc = "PASS"
            add = f"(>={min_bc})"
            if raw_bases < UN2N(min_bc):
                qc = "Low_Base_Count"
                add = f"(<{min_bc})"
                failed.append(sample)

            out.write(f"{sample}\t{cell_count}\t{N2UN(unit,raw_bases)}{add}\t{qc}\n")

        out.write("\n")
        if not failed:
            out.write(f"共计{len(args.sample)}个样本, 质控合格\nCondition:0\n")
        else:
            out.write(f"共计{len(args.sample)}个样本, {len(failed)}个样本质控不合格, 不合格样本: {' '.join(failed)}\nCondition:1\n")



if __name__ == '__main__':
    main() 
