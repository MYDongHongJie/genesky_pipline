anno_file = config["subcluster"]["celltype_file"]
#检查是否存在这个文件
if not os.path.exists(anno_file):
    raise ValueError(f"注释文件 {anno_file} 不存在！")


if params_list["re_sample_diff_file"] is not None and os.path.exists(params_list["re_sample_diff_file"]):
    celltype_ext = f"--re_sample_diff {params_list['re_sample_diff_file']}"
    with open(params_list["re_sample_diff_file"], "r") as f:
        header = f.readline().strip().split("\t")
    split = ["sample.name", "group","raw"]
    split = split + [col for col in header if col not in split]
    split.remove("raw")
    split.remove("sample")
    celltype_Vis_ext = f" -s {','.join(split)} "
else:
    celltype_ext = ""
    celltype_Vis_ext=""
rule celltype:
    input:
        params_list["input"]
    output:
        f"{params_list['output']}/celltype_annoted.rds"
    resources:
        mem_mb=50000,
        mpi="pmi2",
    threads:10
    params:
        anno_file = anno_file,
        outdir = params_list["output"],
        report_log = "log/report_need.log"
    shell:
        "echo 运行的代码是: >> {params.report_log};"
        "echo '/home/donghj/snakemake/bin/Rscript  ./script/genesky_singlcell_tools.r -i {input} -o {params.outdir}  --prefix celltype_annoted --report FALSE  Celltyping --cluster {params.anno_file} {celltype_ext}' >> {params.report_log};"
         "{{ "
        "  if /home/donghj/snakemake/bin/Rscript  ./script/genesky_singlcell_tools.r -i {input} -o {params.outdir} --prefix celltype_annoted --report FALSE  Celltyping --cluster {params.anno_file} {celltype_ext} "
        "      2>> {params.report_log}; then "
        "    echo '[SUCCESS] rule celltype OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule celltype FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

pointsize = config["subcluster"]["pointsize"]
if pointsize is None:
  celltype_Vis_ext = f"{celltype_Vis_ext} --groupby celltype"
else:
  celltype_Vis_ext = f"{celltype_Vis_ext} --groupby celltype --pointsize {pointsize}"

rule celltype_Vis:#可视化
    input:
        rules.celltype.output
    output:
        expand("{outdir}/01_Cluster/{reduction}/{reduction}_reduction_groupby_celltype_splitby_{label}{ext}",outdir=params_list["output"],reduction=["umap", "tsne"], label=["sample.name", "group"], ext=[".pdf", ".png"]),
        expand("{outdir}/01_Cluster/{reduction}/{reduction}_reduction_density_isobars_groupby_celltype{ext}",
               reduction=["umap", "tsne"], ext=[".pdf", ".png"], outdir=params_list["output"]),
        expand("{outdir}/02_celltype_proportion_visualization/SummaryCluster_celltype_splitby_{label}{ext}",
               label=["sample.name", "group"], ext=[".pdf", ".png", ".xls"],outdir=params_list["output"]),
        expand("{outdir}/02_celltype_proportion_visualization/SummaryCluster_celltype_splitby_{label}_alluvial{ext}",
               label=["sample.name", "group"], ext=[".pdf", ".png"], outdir=params_list["output"])
    resources:
        mem_mb=50000,
        mpi="pmi2",
    threads:10
    params:
        pattle=pattle,
        report_log="log/report_need.log",
        outdir = params_list["output"]
    priority: 100
    shell:
        "echo 运行的代码是: >> {params.report_log};"
        "echo '/home/donghj/snakemake/bin/Rscript  ./script/genesky_singlcell_tools.r -i {input} -o {params.outdir} --report FALSE  Visualize --pattle {params.pattle} {celltype_Vis_ext}' >> {params.report_log};"
        "{{ if /home/donghj/snakemake/bin/Rscript  ./script/genesky_singlcell_tools.r -i {input} -o {params.outdir} --report FALSE  Visualize --pattle {params.pattle} {celltype_Vis_ext} 2>> {params.report_log} ;"
        "then echo '[SUCCESS] rule celltype_Vis OK' >> {params.report_log}; "
        "else echo '[ERROR] rule celltype_Vis FAILED' >> {params.report_log}; exit 1; fi; }}"

rule marker_celltype:
    input:
        rules.celltype.output
    output:
        marker_result = f"{params_list['output']}/03_Marker/allmarkers.xlsx",
        marker_other = expand("{outdir}/03_Marker/{ext}",ext=["cluster_top5_dotplot.pdf","marker_number.pdf"],outdir=params_list["output"])
    resources:
        mem_mb=50000,
        mpi="pmi2",
    threads:10
    params:
        outdir = params_list["output"],
        min_pct = config["marker"]["min_pct"],
        logfcthreshold = config["marker"]["logfcthreshold"],
        topn = config["marker"]["topn"],
        avetopn = config["marker"]["avetopn"],
        #pattle =pattle,
        report_log ="log/report_need.log",
        script_file = config["marker"]['script']
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo '/home/donghj/snakemake/bin/Rscript {params.script_file} -i {input} -o {params.outdir}/03_Marker marker "
        "--plot TRUE --min_pct {params.min_pct} --logfcthreshold {params.logfcthreshold} "
        "--topn {params.topn} --avetopn {params.avetopn} --pattle customecol2_light -g celltype' >> {params.report_log};"
       # "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "  if /home/donghj/snakemake/bin/Rscript {params.script_file} -i {input} -o {params.outdir}/03_Marker marker "
        "      --plot TRUE --min_pct {params.min_pct} --logfcthreshold {params.logfcthreshold} "
        "      --topn {params.topn} --avetopn {params.avetopn} --pattle customecol2_light -g celltype "
        "      2>> {params.report_log}; then "
        "    echo '[SUCCESS] rule celltype marker vis OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule celltype marker vis FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

