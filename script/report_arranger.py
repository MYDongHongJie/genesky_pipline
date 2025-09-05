#!/home/donghj/snakemake/bin/python
#report_arranger.py
#version: 0.1
#pipeline: none
#author: dhj
#update: 2025/07/31
#整理结果报告，用以生成最后的报告结构
import argparse
#import requests
import yaml
import json
import os
from pathlib import Path
import subprocess
def ckeck_dir(result_dir,out_dir):
  cellranger_dir = os.path.join(result_dir , "cellranger")
  if os.path.isdir(cellranger_dir):
    print("检测到cellranger文件夹")
    SUPPLEMENTARY_FILES = os.path.join(out_dir ,"../","Supplementary_Files","CellRanger")
    subfolders = [d.name for d in Path(cellranger_dir).iterdir() if d.is_dir() and d.name != "aggr"]
    os.makedirs(SUPPLEMENTARY_FILES, exist_ok=True)
    for subfolder in subfolders:
      os.makedirs(os.path.join(SUPPLEMENTARY_FILES,subfolder), exist_ok=True)
      #subprocess.call(f'mkdir  {SUPPLEMENTARY_FILES}/{subfolder}', shell=True)
      subprocess.call(f'cp -r {cellranger_dir}/{subfolder}/{{filtered_feature_bc_matrix,molecule_info.h5,web_summary.html}} {SUPPLEMENTARY_FILES}/{subfolder}', shell=True)
  else:
    print("未检测到cellranger文件夹")
  qc_dir = os.path.join(result_dir , "01_QC") 
  if os.path.isdir(qc_dir):
    print("检测到01_QC文件夹")
    QC_report_load = os.path.join(out_dir ,"03.Cell_QC" ,"Cell_Quality_Images")
    os.makedirs(QC_report_load, exist_ok=True)
    subprocess.call(f'cp {qc_dir}/*.p* {QC_report_load}', shell=True)
    subprocess.call(f'cp {qc_dir}/QC_cell_counts_per_step.xlsx {QC_report_load}/../', shell=True)
  else:
    print("未检测到01_QC文件夹")
  
  Cluster_dir = os.path.join(result_dir , "02_Cluster")
  if os.path.isdir(Cluster_dir):
    print("检测到02_Cluster文件夹")
    CLUSTER_report_load = os.path.join(out_dir ,"04.Dimensional_Reduction")
    os.makedirs(CLUSTER_report_load, exist_ok=True)
    subprocess.call(f'cp -r {Cluster_dir}/umap {Cluster_dir}/tsne {CLUSTER_report_load} 2>/dev/null', shell=True)
    subprocess.call(f'cp -r {Cluster_dir}/pre_clustering/PCA_ElbowPlot.p* {CLUSTER_report_load}', shell=True)
    subprocess.call(f'cp -r {Cluster_dir}/dimensional_reduction.cloupe {CLUSTER_report_load}', shell=True)
    subprocess.call(f'cp -r {Cluster_dir}/clustering_results.csv {CLUSTER_report_load}', shell=True)
    if os.path.isdir(qc_dir):
      subprocess.call(f'cp -r {Cluster_dir}/pre_clustering/VariableFeatures_distribution.p* {QC_report_load}', shell=True)
      subprocess.call(f'cp -r {Cluster_dir}/pre_clustering/{{nCount_RNA_umap,nFeature_RNA_umap}}.p* {QC_report_load}', shell=True)
  else:
    print("未检测到02_Cluster文件夹")
  Vis_dir = os.path.join(result_dir , "03_seurat_clusters_proportion_visualization")
  if os.path.isdir(Vis_dir):
    print("检测到03_seurat_clusters_proportion_visualization文件夹")
    CLUSTER_report_load = os.path.join(out_dir ,"04.Dimensional_Reduction")
    CLUSTER_VIS_report_load =os.path.join(CLUSTER_report_load,"Cluster_Visualization")
    os.makedirs(CLUSTER_VIS_report_load, exist_ok=True)
    subprocess.call(f'cp -r {Vis_dir}/* {CLUSTER_VIS_report_load}', shell=True)
  else:
    print("未检测到03_seurat_clusters_proportion_visualization文件夹")
  Marker_dir = os.path.join(result_dir , "04_Marker")
  if os.path.isdir(Marker_dir):
    print("检测到04_Marker文件夹")
    MAKER_report_dir = os.path.join(out_dir , "05.Cluster_Marker_Gene")
    os.makedirs(MAKER_report_dir, exist_ok=True)
    subprocess.call(f'cp -r {Marker_dir}/* {MAKER_report_dir}', shell=True)
  else:
    print("未检测到04_Marker文件夹")
  Enrich_dir = os.path.join(result_dir , "05_Enrichment")
  if os.path.isdir(Enrich_dir):
    print("检测到05_Enrichment文件夹")
    ENRICH_report_dir = os.path.join(out_dir , "06.Enrichment")
    os.makedirs(ENRICH_report_dir, exist_ok=True)
    Enrich_dir_temp = os.path.join(Enrich_dir,"temp")
    if os.path.isdir(Enrich_dir_temp):
      subprocess.call(f'rm -rf  {Enrich_dir_temp}', shell=True)
    subprocess.call(f'cp -r {Enrich_dir}/* {ENRICH_report_dir}', shell=True)
    subprocess.call(f'rename ".txt" ".xls"  {ENRICH_report_dir}/*/*', shell=True)
  else:
    print("未检测到05_Enrichment文件夹")
  singleR_dir = os.path.join(result_dir , "06_singleR")
  if os.path.isdir(singleR_dir):
    print("检测到06_singleR文件夹")
    SINGLER_report_dir = os.path.join(out_dir , "07.SingleR")
    os.makedirs(SINGLER_report_dir, exist_ok=True)
    subprocess.call(f'cp -r {singleR_dir}/*{{png,pdf,xls,xlsx}} {SINGLER_report_dir}', shell=True)
  else:
    print("未检测到06_singleR文件夹")
  dir_name = ["01_Cluster","02_celltype_proportion_visualization","03_Marker","04_diffexp","05_Enrichment","Ident_CellType_Markers","06_Gene_GSEA"]
  Advanced_result = []
  for root, dirs, files in os.walk(os.path.join(result_dir ,"Advanced_analytics")):
    if any (dir in dirs for dir in dir_name):
      Advanced_result.append(os.path.basename(root))
  if Advanced_result:
    for group_dir in Advanced_result:
      all_celltype_dir = os.path.join(result_dir ,"Advanced_analytics", group_dir)
      if os.path.isdir(all_celltype_dir):
        print(f"检测到Advanced_analytics/{group_dir}文件夹")
        ALL_CELLTYPE_report_dir = os.path.join(out_dir , "Advanced_analytics", group_dir)
        os.makedirs(ALL_CELLTYPE_report_dir, exist_ok=True)
        if os.path.isfile(os.path.join(all_celltype_dir,"celltype_annoted.rds")):
          os.makedirs(os.path.join(out_dir ,"../","Supplementary_Backup",group_dir),exist_ok=True)
          subprocess.call(f'cp -r {os.path.join(all_celltype_dir,"celltype_annoted.rds")} {os.path.join(out_dir ,"../","Supplementary_Backup",group_dir)}', shell=True)
          
      for dir in dir_name:
        nedd_dif = os.path.join(all_celltype_dir, dir)
        if os.path.isdir(nedd_dif):
          if os.path.isdir(os.path.join(ALL_CELLTYPE_report_dir,dir)):
            subprocess.call(f'cp -r {nedd_dif}/* {ALL_CELLTYPE_report_dir}/{dir}', shell=True)
          else:
            subprocess.call(f'cp -r {nedd_dif} {ALL_CELLTYPE_report_dir}/{dir}', shell=True)
      if os.path.isdir(os.path.join(all_celltype_dir,"CellChat")):
        print(f"检测到Advanced_analytics/{group_dir}/CellChat文件夹")
        cellchat_group_result =[]
        for root, dirs, files in os.walk(os.path.join(all_celltype_dir,"CellChat")):
          if "cellchat_results.rds" in files:
            cellchat_group_result.append(os.path.basename(root))
        for cellchat_group in cellchat_group_result:
          src = os.path.join(all_celltype_dir, "CellChat", cellchat_group) + "/"
          dst = os.path.join(ALL_CELLTYPE_report_dir, "07_CellChat", cellchat_group)
          os.makedirs(dst, exist_ok=True)
          subprocess.call(f'rsync -av --exclude="*.rds" {src} {dst}', shell=True)
  else:
    print("未在Advanced_analytics检测到任何文件夹")

def project_make(id, project, department, sequencing, Ref, out_dir):
    """
    将项目信息写入格式化的文本文件（project.txt）
    
    参数:
        id (str): 项目ID
        project (str): 项目名称
        department (str): 部门名称
        sequencing (str): 测序信息
        Ref (str): 参考基因组
        out_dir (str): 输出目录路径
    """
    # 创建格式化内容
    content = f"""\
		id\t\t\t= {id}
		project\t\t= {project}
		department\t= {department}
		sequencing\t= {sequencing}
		Ref\t\t= {Ref}
		"""
    # 确保输出目录存在
    os.makedirs(out_dir, exist_ok=True)
    
    # 写入文件
    output_file = os.path.join(out_dir, "project.txt")
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f"项目信息已保存至: {output_file}")



  
def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--result_dir',required=True,type=str, help='代码的结果输出路径')
    parser.add_argument('--out_dir',required=True,type=str, help='report的结果输出路径')
    parser.add_argument('--config',default='no_id',type=str, help='config路径')
    #parser.add_argument('--log',default='',type=str, help='log file of message')
    args = parser.parse_args()
    
    config = args.config
    with open(config) as f:
        config = yaml.safe_load(f)
    
    result_load = os.path.join(args.out_dir , "SCGES_单细胞转录组测序分析报告")
    ckeck_dir(args.result_dir, result_load)
    #project_txt_load = os.path.join(args.out_dir , "")
    species_ref_dict = {"human":"GRCh38", "mouse":"GRCm39","pig":"susScr11","rabbit":"oryCun2","rat":"GRCr8"}
    if config['header']['species'] in species_ref_dict.keys():
        ref = species_ref_dict[config['header']['species']]
    else:
        ref = config['header']['species']
    project_make(config['header']['contract_number'], config['header']['project_name'], config['header']['department'], ref, config['header']['sequencing'] , args.out_dir)

if __name__ == '__main__':
    main() 