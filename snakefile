#STDv3.smk.py
#version: 0.1
#pipeline: single cell 3end
#author: dhj
#update: 2025/07/18

#put all used tools into bin folder
#FastQC : v0.12.1
#cellranger : 版本控制

#python : 3.13.2
#python lib:
#treelib 1.7.0
#argparse
#json 

#R/Rscript : 4.4.2
#R lib:
#argparse
#patchwork
#ggsci
#circlize
#ggraph
#data.table
#future
#tibble
#dplyr
#Seurat
#cowplot
#ggplot2
#scales
#jsonlite
#RColorBrewer
#BiocGenerics
#Cairo
#ragg
#tidyverse
#ComplexHeatmap
#DoubletFinder
#harmony
#dittoSeq

#------------------------------------加载模块------------------------------------
import os
import re
import glob
import json
import time
from treelib import Tree, Node
import pandas as pd
# ---------------------------------- 配置信息 ----------------------------------
configfile: "config.yaml"

species = config['header']['species']
all_specise_load = "/home/donghj/scRNA_snakemake/extra_dir/ref.json"
all_specis = json.load(open(all_specise_load))
platforms = config['header']['platform']
ref = all_specis[species][platforms+'_dir']

sample_info = config['header']['sample_info']
pattle = config['header']['pattle']
cellranger = config['header']['cellranger_version']
cellranger_path = f"/home/genesky/software/cellranger/{cellranger}/cellranger"
used_dnbc4tools ="/home/pub/software/dnbc4tools2.1.3/dnbc4tools"
contract_number = config['header']['contract_number']
result_output = config['header']['output']
raw_data_dir = "/home/pub/project/research"
TMP_DIR = config["header"]["tmp_dir"]
os.makedirs(TMP_DIR, exist_ok=True)
person_list = {'donghj': "董宏杰", 'hangxh': "韩晓和", "zhouyh": "周元昊"}
#工作人员确定
user=os.getenv("USER")
person = person_list[user]
cwd = os.getcwd()
abs_path = os.path.join(cwd,result_output)
# ------------------------ 生成样本树（2层：group-sample） ------------------------
#检测是否有sample_info文件
if not os.path.exists(sample_info):
    rawdatapath = os.path.join("/home/pub/project/research", contract_number)
    files = os.listdir(rawdatapath)
    samples = set()
    for f in files:
        match = re.sub(r'_(R1|R2)\.fastq\.gz$', '', f)
        samples.add(match)
    with open(sample_info, "w") as out:
        for s in sorted(samples):
            out.write(f"{s}\t{s}\n")
def TreeSampleInfo_2level(file):
    sampletree = Tree()
    sampletree.create_node(tag='root', identifier='0_root')
    for line in open(file):
        if not line.startswith("#"):
            line = line.strip().split()
            if len(line) < 2:
                continue
            group, sample = line[0], line[1]
            group_id = f"1_{group}"
            if not sampletree.contains(group_id):
                sampletree.create_node(tag=group, identifier=group_id, parent='0_root')
            sample_id = f"2_{sample}"
            if not sampletree.contains(sample_id):
                sampletree.create_node(tag=sample, identifier=sample_id, parent=group_id)
    return sampletree


def log_simplify(log):
    if os.path.exists(log):
        file=open(log,'r')
        logs=file.readlines()
        file.close()
    else:
        logs=log.splitlines()
    errors="\n".join(filter(lambda l:l.startswith(('Error','Missing')), logs))
    return("\""+errors+"\"")


#定义一个函数，如果load_dir存在，则返回fastq_load，否则返回fastp的output
#load_dir = config['cellranger']['fastq_load']
def get_fastq_path(wildcards):
    base = config.get('cellranger', {}).get('fastq_load', None)
    if base and os.path.isdir(base):
        return base
    else:
        return  os.path.abspath(TMP_DIR)

def convert_to_r_format(subset_str):
    """
    将 "celltype:Bcell,Tcell;sample.name:sample1" 转换为R格式
    """
    conditions = []
    
    for condition in subset_str.split(';'):
        if ':' in condition:
            field, values_str = condition.split(':', 1)
            field = field.strip()
            values = [v.strip() for v in values_str.split(',') if v.strip()]
            
            if values:
                # 使用两个单引号表示R中的单引号
                values_r = ', '.join([f"\"{v}\"" for v in values])
                conditions.append(f"{field} %in% c({values_r})")
    # ext = 
    # ext = '\'\''+ext+'\'\''
    return ' & '.join(conditions)
def find_subfolder_names_with_bam(root_dir):
    """
    返回 root_dir 下所有包含 'possorted_genome_bam.bam' 文件的子文件夹名称列表
    """
    result = []
    target = "possorted_genome_bam.bam"

    for dirpath, _, filenames in os.walk(root_dir):
        if target in filenames:
            # 只取当前文件夹的名字
            result.append(os.path.basename(dirpath))
    return result


#报错处理
onerror:
    # shell("echo -e \""+task_number+"\terror\" >> /PERSONALBIO/work/singlecell/s07/auto_script/auto_lims_log.tsv && wechat_notify.py --AT 刘少芬 "+contact_person+" --task "+task_number+" --status error --log "+log_simplify(log))
     shell("echo -e 报错" )


# ------------------------ 规则定义 ------------------------
if config["header"].get("cellranger_Analysis",False):
    SampleTree = TreeSampleInfo_2level(sample_info)
    SAMPLES_NAME = [n.tag for n in SampleTree.filter_nodes(lambda x: SampleTree.depth(x) == 2)]
    SAMPLES_GROUP = [n.tag for n in SampleTree.filter_nodes(lambda x: SampleTree.depth(x) == 1)]
    print(SampleTree)
    include: "script/rules/cellranger.smk"
    CR_REQUIRED = [
    file 
    for sublist in [
        expand("{result_output}/cellranger/{sample}/web_summary.html", sample=SAMPLES_NAME, result_output=result_output),
        expand("{result_output}/cellranger/{sample}/metrics_summary.csv", sample=SAMPLES_NAME, result_output=result_output),
        expand("{result_output}/cellranger/{sample}/cloupe.cloupe", sample=SAMPLES_NAME, result_output=result_output),
        expand("{result_output}/cellranger/{sample}/filtered_feature_bc_matrix/{matrixfile}", 
               sample=SAMPLES_NAME, 
               matrixfile=['barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz'], 
               result_output=result_output),
        [f"{result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/CellRanger_Summary.xlsx","log/cellranger_aggr.log"]
    ] 
    for file in sublist
    ]
    HUADA_C4_REQUIRED = [
    file 
    for sublist in [
        expand("{result_output}/cellranger/{sample}/output/{sample}_scRNA_report.html", sample=SAMPLES_NAME, result_output=result_output),
        expand("{result_output}/cellranger/{sample}/output/metrics_summary.xls", sample=SAMPLES_NAME, result_output=result_output),
        expand("{result_output}/cellranger/{sample}/output/filtered_feature_bc_matrix/{matrixfile}", 
               sample=SAMPLES_NAME, 
               matrixfile=['barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz'], 
               result_output=result_output),
        [f"{result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/CellRanger_Summary.xlsx"]
    ] 
    for file in sublist]
    MAPPING_SUMMARY={
        '10X':CR_REQUIRED,
        #'mobi':"summary/02_Mobivision/{sample}/{sample}_summary.csv",
        'huadaC4':HUADA_C4_REQUIRED
    }

SUBSEQYENT_REQUIRED = [f"{result_output}/create/create.rds",
                       f"{result_output}/01_QC/after_QC.rds",
                       *expand("{result_output}/03_seurat_clusters_proportion_visualization/SummaryCluster_seurat_clusters_splitby_{label}_alluvial{ext}",label=["sample.name","group"],ext=[".pdf",".png"],result_output=result_output),
                       f"{result_output}/04_Marker/allmarkers.xlsx"]

if not config["singleR"]["ref"] is None:
        singleR_database = config["singleR"]["ref"]
        SUBSEQYENT_REQUIRED.append(f"{result_output}/06_singleR/singleR_celltype.rds")
else:
        singleR_database = "/home/pub/singlecell_ref/singleR_database"
        if not config['header']['organization'] is None:
            singleR_database = singleR_database + "/" + config["header"]["species"] + "/" + config["header"]["organization"]+".rds"
            if  os.path.exists(singleR_database):
                SUBSEQYENT_REQUIRED.append(f"{result_output}/06_singleR/singleR_celltype.rds")
            else:
                print("singleR database not exists: " + singleR_database +" Skip singleR analysis")    
        else:
            print("organization is None, skip singleR analysis")
#高级分析模块判断
def get_Advanced_analytics_input(config_file):
    subset = config_file["subcluster"]["subset"].lower()
    params_list = {}
    if subset == "all":
        params_list["input"] = f"{result_output}/02_Cluster/cluster.rds"
        params_list["output"] = f"{result_output}/Advanced_analytics/all"
    else:
        # 例如 "celltype:Bcell,Tcell;sample.name:sample1"
        sub_rds_load = config_file["subcluster"]["subset"].replace(";", "_").replace(":", "_").replace(",", "_")
        params_list["input"] = f"{result_output}/Advanced_analytics/{sub_rds_load}/recluster/sub.rds"
        params_list["output"] = f"{result_output}/Advanced_analytics/{sub_rds_load}"
    #修改分组和新的样本名
    if config_file["subcluster"]["re_sample_diff_file"] is not None and config_file["subcluster"]["re_sample_diff_file"] != "":
        re_sample_diff_file = config_file["subcluster"]["re_sample_diff_file"]
        if not os.path.exists(re_sample_diff_file):
            raise ValueError(f"分组和重命名样本名文件 {re_sample_diff_file} 不存在！")
        else:
            params_list["re_sample_diff_file"] = re_sample_diff_file
    else:
        params_list["re_sample_diff_file"] = ""
    return params_list



ADVANCED_RESULT=[]

if config["subcluster"]["subset"].lower() !="all":
    params_list = get_Advanced_analytics_input(config)
    if config["subcluster"]["subcluster_rds"] =='all':
        rule_subcluster_input = f"{result_output}/Advanced_analytics/all/celltype_annoted.rds"
    else:
        rule_subcluster_input = config["subcluster"]["subcluster_rds"]
#检查文件是否存在
    if not os.path.exists(rule_subcluster_input):
        raise ValueError(f"文件 {rule_subcluster_input} 不存在！")
    subset_ext=convert_to_r_format(config["subcluster"]["subset"])
    print(f"选择的细胞类型为{subset_ext}")
    include: "script/rules/subcluster.smk" 
    ADVANCED_RESULT.append(f"{params_list["input"]}")
    ADVANCED_RESULT.append(params_list["output"]+"/recluster/sub.cloupe")

if config['subcluster']['module'].get("celltype",False):
    params_list = get_Advanced_analytics_input(config)
    print(params_list)
    include: "script/rules/celltype.smk"
    ADVANCED_RESULT.append(f"{params_list["output"]}/celltype_annoted.rds")
    ADVANCED_RESULT.append(f"{params_list['output']}/03_Marker/allmarkers.xlsx")
    ADVANCED_RESULT = ADVANCED_RESULT +expand("{outdir}/02_celltype_proportion_visualization/SummaryCluster_celltype_splitby_{label}_alluvial{ext}",
            label=["sample.name", "group"], ext=[".pdf", ".png"], outdir=params_list["output"])

if config['subcluster']['module'].get('diffexp',False):
    params_list = get_Advanced_analytics_input(config)
    contract = config["subcluster"]["Diff_contract"]
    contract = contract.split(",")
    contract_splits =  [s.replace(":", "_", 1).replace(":", "-vs-", 1) for s in contract]
    include: "script/rules/diffexp.smk"
    #Diff_contract = config['subcluster']['diffexp']['Diff_contract']
    ADVANCED_RESULT.append("log/diffexp.log")
    ADVANCED_RESULT=ADVANCED_RESULT+ expand("log/diff_enrichment_{contract_split}.log", contract_split=contract_splits)
    ADVANCED_RESULT=ADVANCED_RESULT+ expand("log/diff_GSEA_{contract_split}.log", contract_split=contract_splits)
if config['subcluster']['module'].get('cellchat',False):
    if species not in ["mouse","human"]:
        raise ValueError("Cellchat analysis only support mouse and human")
    params_list = get_Advanced_analytics_input(config)
    if not config['subcluster']["cellchat_contract"] is None:
        cellchat_contract = config['subcluster']["cellchat_contract"]
        cellchat_contract = cellchat_contract.split(",")
        cellchat_groups = [s.split(":")[0] for s in cellchat_contract]
        cellchat_groups = list(set(cellchat_groups))
        cellchat_groups_dic = dict()
        for group in cellchat_groups:
            temp = [s.split(":", 1)[1] for s in cellchat_contract if s.split(":", 1)[0] == group]
            #list装成字符 
            cellchat_groups_dic[group] = "+".join(temp)
    else:
        cellchat_groups = ["AllCells_Unsupervised"]
    include: "script/rules/cellchat.smk"
    ADVANCED_RESULT=ADVANCED_RESULT+ expand("{outdir}/CellChat/{cellchat_group}/cellchat_results.rds",outdir = params_list["output"],cellchat_group = cellchat_groups)
    ADVANCED_RESULT=ADVANCED_RESULT+ expand("{outdir}/CellChat/{cellchat_group}/communication.xls", outdir = params_list["output"],cellchat_group = cellchat_groups)

if config['subcluster']['module'].get('monocle',False):
    params_list = get_Advanced_analytics_input(config)
    M2_REQUIRED=[  f"{params_list['output']}/Monocle/monocle2.rds",
          f"{params_list['output']}/Monocle/genes_for_order.csv",
          f"{params_list['output']}/Monocle/diff_test_Pseudotime.txt"]
    M2P_REQUIRED=[f"{params_list['output']}/Monocle/trajectory_order_tree/cell_trajectory_Pseudotime.pdf"]
    BEAM_REQUIRED=[  f"{params_list['output']}/Monocle/BEAM/BEAM_order.xls"]
    method = config['subcluster']['monocle_method']
    idents = config['subcluster']['monocle_idents']
    pcount = config['subcluster']['monocle_pcount'] if not config['subcluster']["monocle_pcount"] is None else ""
    maxcell = config['subcluster']['monocle_maxcell'] if not config['subcluster']["monocle_maxcell"] is None else ""
    root= config['subcluster']['monocle_root'] if not config['subcluster']["monocle_root"] is None else ""
    nselect = config['subcluster']['monocle_nselect'].split(",") if not config['subcluster']["monocle_nselect"] is None else ""
    ngroup = config['subcluster']['monocle_ngroup'].split(",") if not config['subcluster']["monocle_ngroup"] is None else ""
    Monocle2Plot_ext = f"--colorby {config['subcluster']['monocle_colorby']}"
    Monocle2Plot_ext = Monocle2Plot_ext + f" --merge {config['subcluster']['monocle_merge']}" if not config['subcluster']["monocle_merge"] is None else Monocle2Plot_ext
    ADVANCED_RESULT=ADVANCED_RESULT+M2_REQUIRED+M2P_REQUIRED
    if config['subcluster'].get('monocle_BEAM',False):
      ADVANCED_RESULT=ADVANCED_RESULT+BEAM_REQUIRED
    include: "script/rules/monocle.smk"

if config['subcluster']['module'].get('scenic',False):
    params_list = get_Advanced_analytics_input(config)
    tfdb = all_specis[species]["tf"]
    featuredb = all_specis[species]["feather"]
    motifdb = all_specis[species]["tbl"]

    if config['subcluster']['scenic_subset'] is None:
        subset = ""
    else:
        subset = f"--subset {config['subcluster']['scenic_subset']}"

    if config["subcluster"]["scenic_downsample"] is None:
        downsample = ""
    else:
        downsample = f"--downsample {config["subcluster"]["scenic_downsample"]}"
    path = os.path.join(params_list['output'],"SCENIC/temp")
    file_exists = os.path.isfile(os.path.join(path, 'sub_cells.rds'))
    if file_exists:
        rds = os.path.join(path, 'sub_cells.rds')
    else:
        rds = f"{params_list['output']}/celltype_annoted.rds"
    get_CSI_group_raw = config["subcluster"]["scenic_groupby"]
    get_CSI_group = get_CSI_group_raw.split(",")[0]
    scenic_result_plot = f"{params_list['output']}/SCENIC/plot"
    include: "script/rules/scenic.smk"
    ADVANCED_RESULT.append(f"{scenic_result_plot}/3.Cisplot/CSI_heatmap.csv")

if config['subcluster']['module'].get('scvelo',False):
    params_list = get_Advanced_analytics_input(config)
    CellRanger_Res = os.path.join(result_output,"cellranger")
    sample_cellrangers=find_subfolder_names_with_bam(CellRanger_Res)
    gtf_file = os.path.join(all_specis[species]["10X_dir"],"genes")
    if os.path.exists(os.path.join(gtf_file,"genes.gtf.gz")):
        os.makedirs(os.path.join(result_output,"loom"),exist_ok=True)
        subprocess.call([f"gunzip -c {gtf_file}/genes.gtf.gz >{result_output}/loom/genes.gtf"],shell=True)
        gtf = os.path.join(result_output,"loom/genes.gtf")
    else:
        gtf = os.path.join(gtf_file,"genes.gtf")
    loom_want = expand("{result_output}/loom/{sample}.loom",sample=sample_cellrangers,result_output=result_output)
    if config['subcluster']['scvelo_subset'] is None:
        scvelo_subset = ""
    else:
        scvelo_subset = f"--subset {config['subcluster']['scvelo_subset']}"
    ADVANCED_RESULT=ADVANCED_RESULT+ loom_want
    ADVANCED_RESULT.append(f"{params_list['output']}/SCVELO/adata_with_scvelo.h5ad")
    ADVANCED_RESULT.append(f"{params_list['output']}/SCVELO/velocity_data.xls")
    include: "script/rules/scvelo.smk"



def all_input(wildcards):
    wanted_input = []
    if config['header'].get('cellranger_Analysis', False):
        wanted_input.extend(
            expand("log/{sample}_fastqc.log", sample=SAMPLES_NAME) +
            expand("log/{sample}_fqchk.log", sample=SAMPLES_NAME) +
            expand("log/{sample}_fastp.log", sample=SAMPLES_NAME) +
            expand("log/{sample}_fastqc_clean.log", sample=SAMPLES_NAME) +
            expand("log/{sample}_fqchk_clean.log", sample=SAMPLES_NAME) +
            expand("log/{sample}_stat.log", sample=SAMPLES_NAME)
        )
        
        # 添加固定文件
        wanted_input.extend([
            "log/wechat_notice.log",
            "log/wechat_notice2.log"]
        )
        wanted_input.extend( MAPPING_SUMMARY[platforms])
    if config["header"].get("Standard_Analysis",False):
        wanted_input.extend(SUBSEQYENT_REQUIRED)
    if config["header"].get("Report_create",False):
        wanted_input.extend([f"{result_output}/report/数据分析结果报告.pdf"])
    wanted_input.extend(ADVANCED_RESULT)
    return wanted_input

rule all:
    input:
        all_input

#create
if config["create"].get("input_file") is None:
    Cellranger_check = "log/cellranger_finished.log"
    Create_input_file = f"{result_output}/cellranger"
else:
    Create_input_file = config["create"]['input_file']

rule Seurat_Create:
    input:
        lambda wildcards: Create_input_file
    output:
        rds=f"{result_output}/create/create.rds"
    params:
        script_file = config["create"]["script"],
        mincells = config["create"]["mincells"],
        minfeatures = config["create"]["minfeatures"],
        outdir =  os.path.join(result_output,"create"),
        report_log = "log/report_need.log"
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo 'Rscript {params.script_file} --prefix create -o {params.outdir} create -s {input} "
        "--genecolumn 2 --minCells {params.mincells} --minfeatures {params.minfeatures} "
        "-m {sample_info} --type genesky' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{"
        "  if Rscript {params.script_file} --prefix create -o {params.outdir} create -s {input} "
        "      --genecolumn 2 --minCells {params.mincells} --minfeatures {params.minfeatures} "
        "      -m {sample_info} --type genesky 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule Seurat_Create OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule Seurat_Create FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"
 


#QC
filter_parmams = "nFeature_RNA,nCount_RNA"
nFeature_RNA_extent = config["QC"]["nFeature_RNA"].split(",")
nCount_RNA_extent = config["QC"]["nCount_RNA"].split(",")
filter_max =f"{nFeature_RNA_extent[1]},{nCount_RNA_extent[1]}"
filter_min =f"{nFeature_RNA_extent[0]},{nCount_RNA_extent[0]}"

extra_file = config["QC"]["extra_file"]
if not extra_file is None:
    qc_ext = f"--gmt {extra_file}"
    if not config["QC"].get("extra_col") is None:
        print("提供了额外过滤文件，但是没有给过滤参数，可以选择过滤文件的列名作为过滤参数")
    else:
        extra_col = config["QC"]["extra_col"].split(",")
        extra_percent = config["QC"]["extra_percent"].split(",")
        if len(extra_col) != len(extra_percent):
            print("过滤参数和过滤文件列数不匹配，请检查")
            exit(1)
        else:
            filter_parmams = f"{filter_parmams},{','.join(extra_col)}" 
            filter_max = f"{filter_max},{','.join(extra_percent)}"
            extra_percent_min = ["0"]*len(extra_percent)
            filter_min = f"{filter_min},{','.join(extra_percent_min)}"
else:
    qc_ext = ""

if not config["QC"].get("mtpercent") is None:
    if "mito" in all_specis[species]:
        mt_load = all_specis[species]['mito']
        #检测mt_load是否存在
        if os.path.exists(mt_load):
            #species_mito = mt_load
            scFeature = "nCount_RNA,nFeature_RNA,percent.mt"
            filter_parmams = f"{filter_parmams},percent.mt"
            filter_max = f"{filter_max},{config['QC']['mtpercent']}"
            filter_min = f"{filter_min},0"
        else:
            print(f"没有在 {mt_load} 找到该物种的线粒体基因文件，请检查")
            exit(1)
    else:
        print(f"没有在 {ref} 找到该物种的线粒体基因文件，请检查，本次分析将不计算线粒体基因含量")
        scFeature = "nCount_RNA,nFeature_RNA"
else:
    scFeature = "nCount_RNA,nFeature_RNA"
scFeatures = scFeature.split(",")
if config["QC"]["doublet"]:
    qc_ext = qc_ext + " --doublet TRUE"
else:
    qc_ext = qc_ext + " --doublet FALSE"
chem_version = config["QC"]["chem_version"]
if chem_version == "v4":
    qc_ext = f"{qc_ext} --chem_version {chem_version}"
else:
    qc_ext = f"{qc_ext}  --chem_version other"



rule SeuratQC:
    input:
        rules.Seurat_Create.output.rds
    output:
        f"{result_output}/01_QC/QC_cell_counts_per_step.xlsx",
        multiext(f"{result_output}/01_QC/QC_metrics_violin_before_vs_after_filtering", ".pdf", ".png"),
        rds=f"{result_output}/01_QC/after_QC.rds"
    params:
        script_file = config["QC"]["script"],
        pattle = pattle,
        species = species,
        report_log = "log/report_need.log"
    resources:
        mpi="pmi2"
    shell:
        # 先记录运行的命令
        "echo '运行的代码是:' >> {params.report_log};"
        "echo 'Rscript {params.script_file} --input {input} --output {result_output}/01_QC "
        "--prefix after_QC QC -s {params.species} --filter {filter_parmams} --high {filter_max} "
        "--low {filter_min} --pattle {params.pattle} {qc_ext}' >> {params.report_log};"

        # 激活环境并执行命令，记录完整输出和错误
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "  if Rscript {params.script_file} --input {input} --output {result_output}/01_QC "
        "      --prefix after_QC QC -s {params.species} --filter {filter_parmams} --high {filter_max} "
        "      --low {filter_min} --pattle {params.pattle} {qc_ext} 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule SeuratQC OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule SeuratQC FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"


#解析降维聚类的config参数

Seurat_Cluster_ext = ""
dim = config["Seurat_Cluster"]["dim"]
if not config["Seurat_Cluster"]["dim"] is None:
    Seurat_Cluster_ext = f" --dim {dim}" 
if not config["Seurat_Cluster"]["regression"] is None:
    Seurat_Cluster_ext = f"{Seurat_Cluster_ext} --regression {config['Seurat_Cluster']['regression']}"

rule Seurat_Cluster:
    input:
        rules.SeuratQC.output.rds
    output:
        rds=f"{result_output}/02_Cluster/cluster.rds",
        ElbowPlot=multiext(f"{result_output}/02_Cluster/pre_clustering/PCA_ElbowPlot", ".pdf", ".png"),
        VariableFeatures_distribution=multiext(f"{result_output}/02_Cluster/pre_clustering/VariableFeatures_distribution", ".pdf", ".png"),
        scFeature_result=expand("{result_output}/02_Cluster/pre_clustering/{feature}_umap{ext}",
                                feature=scFeatures, ext=[".pdf", ".png"], result_output=result_output)
    params:
        script_file=config["Seurat_Cluster"]["script"],
        pattle=pattle,
        cycle=config["Seurat_Cluster"]["cycle"],
        gather=config["Seurat_Cluster"]["gather"],
        resolution=config["Seurat_Cluster"]["res"],
        rely=config["Seurat_Cluster"]["rely"],
        scFeature=scFeature,
        report_log="log/report_need.log"
    resources:
        mpi="pmi2"
    shell:
        "echo 运行的代码是: >> {params.report_log};"
        "echo 'Rscript {params.script_file} --input {input} --output {result_output}/02_Cluster --prefix cluster Clustering -g {params.gather} --cycle {params.cycle} --rely {params.rely} --res {params.resolution} --pattle {params.pattle} --scFeature {params.scFeature} {Seurat_Cluster_ext} ' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ if Rscript {params.script_file} --input {input} --output {result_output}/02_Cluster --prefix cluster Clustering "
        "-g {params.gather} --cycle {params.cycle} --rely {params.rely} --res {params.resolution} "
        "--pattle {params.pattle} --scFeature {params.scFeature} {Seurat_Cluster_ext} 2>> {params.report_log} ; then "
        "echo '[SUCCESS] rule Seurat_Cluster OK' >> {params.report_log}; "
        "else echo '[ERROR] rule Seurat_Cluster FAILED' >> {params.report_log}; exit 1; fi; }}"


Cluster_plot_ext = ""
if not config["visualize"]["pointsize"] is None:
    Cluster_plot_ext = f"  --pointsize {config["visualize"]["pointsize"]}"

rule Cluster_plot:
    input:
        rules.Seurat_Cluster.output.rds
    output:
        expand("{result_output}/02_Cluster/{reduction}/{reduction}_reduction_groupby_seurat_clusters_splitby_{label}{ext}",
               reduction=["umap", "tsne"], label=["sample.name", "group"], ext=[".pdf", ".png"], result_output=result_output),
        expand("{result_output}/02_Cluster/{reduction}/{reduction}_reduction_density_isobars_groupby_seurat_clusters{ext}",
               reduction=["umap", "tsne"], ext=[".pdf", ".png"], result_output=result_output),
        expand("{result_output}/03_seurat_clusters_proportion_visualization/SummaryCluster_seurat_clusters_splitby_{label}{ext}",
               label=["sample.name", "group"], ext=[".pdf", ".png", ".xls"], result_output=result_output),
        expand("{result_output}/03_seurat_clusters_proportion_visualization/SummaryCluster_seurat_clusters_splitby_{label}_alluvial{ext}",
               label=["sample.name", "group"], ext=[".pdf", ".png"], result_output=result_output)
    params:
        script_file=config["visualize"]['script'],
        pattle=pattle,
        report_log="log/report_need.log"
    priority: 100
    resources:
        mpi="pmi2"
    shell:
        "echo 运行的代码是: >> {params.report_log};"
        "echo 'Rscript {params.script_file} -i {input} -o {result_output} Visualize --pattle {params.pattle} {Cluster_plot_ext}' >> {params.report_log};"
        "echo 生成loupe的代码是 >> {params.report_log};"
        "echo 'Rscript {params.script_file} -i {input} -o {result_output}/02_Cluster loupeR -n dimensional_reduction' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ if Rscript {params.script_file} -i {input} -o {result_output} Visualize --pattle {params.pattle} {Cluster_plot_ext} 2>> {params.report_log}  "
        "&& Rscript {params.script_file} -i {input} -o {result_output}/02_Cluster loupeR -n dimensional_reduction 2>> {params.report_log} ; "
        "then echo '[SUCCESS] rule Cluster_plot OK' >> {params.report_log}; "
        "else echo '[ERROR] rule Cluster_plot FAILED' >> {params.report_log}; exit 1; fi; }}"



rule marker:
    input:
        rules.Seurat_Cluster.output.rds
    output:
        marker_result = f"{result_output}/04_Marker/allmarkers.xlsx",
        marker_other = expand("{result_output}/04_Marker/{ext}",ext=["cluster_top5_dotplot.pdf","marker_number.pdf"],result_output=result_output),
        report_log ="log/report_need.log"
    params:
        min_pct = config["marker"]["min_pct"],
        logfcthreshold = config["marker"]["logfcthreshold"],
        topn = config["marker"]["topn"],
        avetopn = config["marker"]["avetopn"],
        #pattle =pattle,
        script_file = config["marker"]['script']
    resources:
        mpi="pmi2"
    shell:
        "echo '运行的代码是:' >> {output.report_log};"
        "echo 'Rscript {params.script_file} -i {input} -o {result_output}/04_Marker marker "
        "--plot TRUE --min_pct {params.min_pct} --logfcthreshold {params.logfcthreshold} "
        "--topn {params.topn} --avetopn {params.avetopn} --pattle customecol2_light' >> {output.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "  if Rscript {params.script_file} -i {input} -o {result_output}/04_Marker marker "
        "      --plot TRUE --min_pct {params.min_pct} --logfcthreshold {params.logfcthreshold} "
        "      --topn {params.topn} --avetopn {params.avetopn} --pattle customecol2_light "
        "      2>> {output.report_log}; then "
        "    echo '[SUCCESS] rule marker OK' >> {output.report_log}; "
        "  else "
        "    echo '[ERROR] rule marker FAILED' >> {output.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"


rule enrichment:
    input:
        rules.marker.output.marker_result
    output:
        log="log/enrichment_finshed.log",
        report_log="log/report_need.log"
    params:
        topn=config["enrichment"]["topn"],
        script_file = config["enrichment"]['script'],
        rankby = config["enrichment"]["rankby"]
    priority: 100
    shell:
        "echo '运行的代码是:' >> {output.report_log};"
        "echo 'Rscript {params.script_file} -o {result_output}/05_Enrichment Enrichments "
        "--file {input} -s {species} --topn {params.topn} --rankby {params.rankby}' >> {output.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "  if Rscript {params.script_file} -o {result_output}/05_Enrichment enrichment "
        "--file {input} -s {species} --topn {params.topn} --rankby {params.rankby} "
        " 2>> {output.report_log} ; then "
        "    touch {output.log}; "
        "    echo '[SUCCESS] rule marker enrichment OK' >> {output.report_log}; "
        "  else "
        "    echo '[ERROR] rule marker enrichment FAILED' >> {output.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"


rule singleR:
    input:
        rules.Seurat_Cluster.output.rds
    output:
        "{result_output}/06_singleR/singleR_celltype.rds"
    params:
        script_file = config["singleR"]['script'],
        ref = singleR_database,
        singleR_groupby = config["singleR"]["groupby"],
        downsample=config["singleR"]["downsample"],
        pattle =pattle,
        report_log="log/report_need.log"
    priority: 100
    resources:
        mpi="pmi2"
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo 'Rscript {params.script_file} -i {input} -o {result_output}/06_singleR "
        "--prefix singleR_celltype singleR --ref {params.ref} --groupby {params.singleR_groupby} "
        "--downsample {params.downsample} --pattle {params.pattle}' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "  if Rscript {params.script_file} -i {input} -o {result_output}/06_singleR "
        "--prefix singleR_celltype singleR --ref {params.ref} --groupby {params.singleR_groupby} "
        "--downsample {params.downsample} --pattle {params.pattle} "
        "2>> {params.report_log}; then "
        "    echo '[SUCCESS] rule singleR OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule singleR FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"



#目前只做到标准分析的报告整理
localrules: report_arrange
rule report_arrange:
    input:
        rules.Seurat_Cluster.output.rds
    output:
        report_html = f"{result_output}/report/数据分析结果报告.html",
        report_pdf = f"{result_output}/report/数据分析结果报告.pdf"
    params:
        result_output = result_output,
        report_log = "log/report_need.log"
    priority: 0
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo 'source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "python ./script/report_arranger.py --result_dir {result_output} --out_dir {result_output}/report "
        "--config config.yaml && perl /home/genesky/pipeline/html_markdown/v6.2.0/pipeline.pl "
        "-t {cwd}/script/html/ -r {abs_path}/report -p {abs_path}/report/project.txt' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "  if python ./script/report_arranger.py --result_dir {result_output} --out_dir {result_output}/report "
        "--config config.yaml && perl /home/genesky/pipeline/html_markdown/v6.2.0/pipeline.pl "
        "-t {cwd}/script/html/ -r {abs_path}/report -p {abs_path}/report/project.txt "
        "2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule Report_Arranger OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule Report_Arranger FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"


