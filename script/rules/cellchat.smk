cellchat_ext = {}
if cellchat_groups == ["AllCells_Unsupervised"]:
    cellchat_ext["AllCells_Unsupervised"] = ""
else:
    for group_use in cellchat_groups:
        cellchat_ext[group_use] = f" --group {group_use} --contrast {cellchat_groups_dic[group_use]} "


rule cellchat:
    input:
        f"{params_list['output']}/celltype_annoted.rds"
    output:
        params_list['output']+"/CellChat/{cellchat_group}/communication.xls",
        params_list['output']+"/CellChat/{cellchat_group}/cellchat_results.rds",
    params:
        pattle=pattle,
        outdir = params_list['output']+"/CellChat/{cellchat_group}",
        ext    = lambda wildcards: cellchat_ext[wildcards.cellchat_group],
        report_log = "log/report_need.log"
    resources:
        mpi="pmi2",
        mem_mb=40000
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo ' source /home/genesky/software/conda/23.3.1/bin/activate /home/genesky/software/cellchat/1.6.1 && Rscript /home/donghj/cellchat/cellchat_main.r --input {input} --species {species} --palette {params.pattle} --column4cell celltype -o {params.outdir} {params.ext} ' >> {params.report_log};"
        "source /home/genesky/software/conda/23.3.1/bin/activate /home/genesky/software/cellchat/1.6.1 && "
        "{{ "
        "if Rscript /home/donghj/cellchat/cellchat_main.r --input {input} --species {species} --palette {params.pattle} --column4cell celltype -o {params.outdir} {params.ext} " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule CellChat {wildcards.cellchat_group} OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule CellChat {wildcards.cellchat_group}  FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"