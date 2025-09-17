rule loom_to_bam:
    input:
        f"{result_output}/cellranger",
    output:
        loom=f"{result_output}/loom/{{sample}}.loom",
        tmp  =f"{result_output}/cellranger/{{sample}}/cellsorted_possorted_genome_bam.bam"
    params:
        report_log = "log/report_need.log",
        gtf = gtf,
        outdir = f"{result_output}/loom/"
    resources:
        mpi="pmi2",
        mem_mb=100000
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo ' source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && Rscript /home/donghj/scvelo/bam2loom.r -i {input}/{wildcards.sample} -n {wildcards.sample} -o {params.outdir} -g {params.gtf} && rm -rf {output.tmp} ' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "if Rscript /home/donghj/scvelo/bam2loom.r -i {input}/{wildcards.sample} -n {wildcards.sample} -o {params.outdir} -g {params.gtf} && rm -rf {output.tmp}" 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule loom_to_bam  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule loom_to_bam   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule scvelo:
    input:
        loom=expand("{result_output}/loom/{sample}.loom", sample=sample_cellrangers,result_output= result_output),
        rds=f"{params_list['output']}/celltype_annoted.rds"
    output:
        params_list['output']+"/SCVELO/velocity_data.xls",
        params_list['output']+"/SCVELO/adata_with_scvelo.h5ad"
    params:
        report_log = "log/report_need.log",
        loom_dir = f"{result_output}/loom/",
        groupby = config["subcluster"]["scvelo_groupby"],
        out_dir =  params_list['output']+"/SCVELO/",
        subset = scvelo_subset
    resources:
        mpi="pmi2",
        mem_mb=50000
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "Rscript /home/donghj/scvelo/scvelo.r -i {input.rds} -o {params.out_dir} -g {params.groupby} -l {params.loom_dir}  2>> {params.report_log};"