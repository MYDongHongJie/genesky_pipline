ref = all_specis[species][platforms+'_dir']
gtf_file_load = os.join.path(ref,"genes")
#检测gtf_file_load目录下是否存在infercnv_need.gtf文件
if not os.path.exists(os.join.path(gtf_file_load,"infercnv_need.gtf")):
    pattern = os.path.join(gtf_file_load, "genes.gtf*")
    gtf_files = glob.glob(pattern)[0]
    if not gtf_files:  # 列表为空
        raise FileNotFoundError(f"No files starting with 'genes.gtf' found in folder: {gtf_file_load}")
else:
    gtf_files = os.join.path(gtf_file_load,"infercnv_need.gtf")
    import 


if ref = all_specis[species][platforms+'_dir']


rule infercnv:
    input:
        rds=f"{params_list['output']}/celltype_annoted.rds",
    output:
        loom=f"{result_output}/loom/{{sample}}.loom",
        tmp  =f"{result_output}/cellranger/{{sample}}/cellsorted_possorted_genome_bam.bam"