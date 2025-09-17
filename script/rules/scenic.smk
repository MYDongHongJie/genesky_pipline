rule scenic_loom_prepare:
    input:
        f"{params_list['output']}/celltype_annoted.rds"
    output:
        params_list['output']+"/SCENIC/temp/count.csv",
    params:
        report_log = "log/report_need.log",
        subset = subset,
        downsample= downsample,
        outdir = params_list['output']+"/SCENIC/temp/",
    resources:
        mpi="pmi2",
        mem_mb=10000
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo '   source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && Rscript /home/donghj/scenic/1.loom_prepare.r -i {input} -n count.csv  -o {params.outdir} {params.downsample} {params.subset} ' >> {params.report_log};"
        "  source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "if Rscript /home/donghj/scenic/1.loom_prepare.r -i {input} -n count.csv  -o {params.outdir} {params.downsample} {params.subset} " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule scenic_loom_prepare  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule scenic_loom_prepare   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule loom_create:
    input:
        check = rules.scenic_loom_prepare.output
    output:
        loom_file=params_list['output']+"/SCENIC/temp/marix.loom"
    params:
        report_log = "log/report_need.log",
        outdir = params_list['output']+"/SCENIC/temp/"
    resources:
        mem_mb=20000,
        mpi="pmi2",
    threads: 1
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo '   source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && python /home/donghj/scenic/2.loom_create.py -c {input.check} -n marix.loom -o {params.outdir} ' >> {params.report_log};"
        "  source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "if python /home/donghj/scenic/2.loom_create.py -c {input.check} -n marix.loom -o {params.outdir} " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule loom_create  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule loom_create   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule pyscenic_grn:
    input:
        loom_file=rules.loom_create.output.loom_file
    output:
        grn_file=params_list['output']+"/SCENIC/grn.csv"
    params:
        tfs_path=tfdb,
        threads=10,
        report_log="log/report_need.log"
    threads: 10
    resources:
        mem_mb=100000,
        mpi="pmi2",
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo '   /home/genesky/software/pyscenic/0.12.1/bin/pyscenic grn --output {output.grn_file} {input.loom_file} {params.tfs_path} --num_workers {params.threads}' >> {params.report_log};"
        "{{ "
        "if /home/genesky/software/pyscenic/0.12.1/bin/pyscenic grn --output {output.grn_file} {input.loom_file} {params.tfs_path} --num_workers {params.threads} " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule pyscenic_grn  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule pyscenic_grn   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule pyscenic_ctx:
    input:
        grn_file=rules.pyscenic_grn.output.grn_file,
        loom_file=rules.loom_create.output.loom_file
    output:
        ctx_file=params_list['output']+"/SCENIC/ctx.csv"
    params:
        feather_paths=featuredb,
        table_path=motifdb,
        threads=10,
        report_log="log/report_need.log"
    resources:
        mem_mb=100000,
        mpi="pmi2",
    threads: 10
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo ' /home/genesky/software/pyscenic/0.12.1/bin/pyscenic ctx {input.grn_file} {params.feather_paths}  --annotations_fname {params.table_path} --expression_mtx_fname {input.loom_file} --output {output.ctx_file} --num_workers {params.threads}  --mask_dropouts' >> {params.report_log};"
        "{{ "
        "if /home/genesky/software/pyscenic/0.12.1/bin/pyscenic ctx {input.grn_file} {params.feather_paths}  --annotations_fname {params.table_path} --expression_mtx_fname {input.loom_file} --output {output.ctx_file} --num_workers {params.threads}  --mask_dropouts " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule pyscenic_ctx  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule pyscenic_ctx   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule pyscenic_aucell:
    input:
        ctx_file=rules.pyscenic_ctx.output.ctx_file,
        loom_file=rules.loom_create.output.loom_file
    output:
        aucell_loom=params_list['output']+"/SCENIC/aucell.loom"
    resources:
        mem_mb=50000,
        mpi="pmi2",
    params:
        threads=10,
        report_log="log/report_need.log"
    threads:10
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo '   /home/genesky/software/pyscenic/0.12.1/bin/pyscenic aucell {input.loom_file} {input.ctx_file} --output {output.aucell_loom} --num_workers {params.threads} ' >> {params.report_log};"
        "{{ "
        "if /home/genesky/software/pyscenic/0.12.1/bin/pyscenic aucell {input.loom_file} {input.ctx_file} --output {output.aucell_loom} --num_workers {params.threads} " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule pyscenic_aucell  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule pyscenic_aucell   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

#------------------------------绘图部分--------------------------------------



rule pyscenic_calcRSS_plot:
    input:
        auc_loom=rules.pyscenic_aucell.output.aucell_loom
    output:
        rss_data=f"{scenic_result_plot}/1.rssplot/{get_CSI_group}_scenic_rss_data.xls",
        regulon_gene_file = f"{scenic_result_plot}/2.AUCplot/regulon_gene_Relationship.xls"
    resources:
        mem_mb=50000,
        mpi="pmi2",
    threads:10
    params:
        rds_calcRSS = rds,
        groupby = config["subcluster"]["scenic_groupby"],
        heatpattle = config['header']["pattle"],
        report_log="log/report_need.log"
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo ' /home/genesky/software/rscenic/1.1.2/bin/Rscript /home/donghj/scenic/3.calcRSS_by_scenic.R -l {input.auc_loom} -i {params.rds_calcRSS} --groupby {params.groupby} -o {scenic_result_plot} --heatpattle {params.heatpattle} ' >> {params.report_log};"
        "{{ "
        "if /home/genesky/software/rscenic/1.1.2/bin/Rscript /home/donghj/scenic/3.calcRSS_by_scenic.R -l {input.auc_loom} -i {params.rds_calcRSS} --groupby {params.groupby} -o {scenic_result_plot} --heatpattle {params.heatpattle} " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule pyscenic_calcRSS_plot  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule pyscenic_calcRSS_plot   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule get_CSI_and_regulon_gene_relationship:
    input:
        get_CSI_check = f"{scenic_result_plot}/2.AUCplot/regulon_gene_Relationship.xls"
    output:
        cis_heatmap= f"{scenic_result_plot}/3.Cisplot/CSI_heatmap.csv",
        filtered_grn_network = f"{scenic_result_plot}/4.regulon_gene/filtered_grn_network.csv"
    resources:
        mem_mb=10000,
        mpi="pmi2",
    threads:8
    params:
        topn = 5,
        rss_data = f"{scenic_result_plot}/1.rssplot/{get_CSI_group}_scenic_rss_data.xls",
        auc_loom = f"{params_list['output']}/SCENIC/aucell.loom",
        regulon_gene_file = f"{scenic_result_plot}/2.AUCplot/regulon_gene_Relationship.xls",
        ingrn = f"{params_list['output']}/SCENIC/grn.csv",
        report_log="log/report_need.log"
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo '   /home/donghj/snakemake/bin/python3.13 /home/donghj/scenic/4.get_CSI_and_regulon_gene_relationship.py --loom_dir {params.auc_loom} -i {params.rss_data} --regulon_gene_file {params.regulon_gene_file} -o {scenic_result_plot} --ingrn {params.ingrn} --topn {params.topn} ' >> {params.report_log};"
        "{{ "
        "if /home/donghj/snakemake/bin/python3.13 /home/donghj/scenic/4.get_CSI_and_regulon_gene_relationship.py --loom_dir {params.auc_loom} -i {params.rss_data} --regulon_gene_file {params.regulon_gene_file} -o {scenic_result_plot} --ingrn {params.ingrn} --topn {params.topn} " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule pyscenic_calcRSS_plot  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule pyscenic_calcRSS_plot   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"
        
