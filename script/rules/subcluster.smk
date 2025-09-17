rule subcluster:
    input:
        rule_subcluster_input
    output:
        params_list["input"]
    resources:
        mem_mb=50000,
        mpi="pmi2",
    threads:10
    params:
        res = config['subcluster']['subcluster_res'],
        batch = config['subcluster']['subcluster_batch'],
        outdir = params_list["output"]+"/recluster",
        report_log="log/report_need.log",
    shell:
        "echo '/home/donghj/snakemake/bin/Rscript  ./script/genesky_singlcell_tools.r -i {input} -o {params.outdir} --subset {subset_ext} --prefix sub Subcluster -b {params.batch} -s {params.res}'  >> {params.report_log} ;"
        "{{ "
        "if /home/donghj/snakemake/bin/Rscript  ./script/genesky_singlcell_tools.r -i {input} -o {params.outdir}  --subset '{subset_ext}'  --prefix sub Subcluster -b {params.batch} -s {params.res} " 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule subcluster  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule subcluster   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule subcluster_loupe:
    input:
        rules.subcluster.output
    output:
        params_list["output"]+"/recluster/sub.cloupe"
    resources:
        mem_mb=50000,
        mpi="pmi2",
    threads:10
    params:
        report_log="log/report_need.log",
        outdir = params_list["output"]+"/recluster",
    shell:
      "echo '/home/donghj/snakemake/bin/Rscript ./script/genesky_singlcell_tools.r -i {input} -o {params.outdir} loupeR -n sub' >> {params.report_log} ;"
      "{{ "
        "if /home/donghj/snakemake/bin/Rscript  ./script/genesky_singlcell_tools.r  -i {input} -o {params.outdir} loupeR -n sub" 
        " 2>> {params.report_log} ; then "
        "    echo '[SUCCESS] rule subcluster_loupe  OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule subcluster_loupe   FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"