rule diffexp:
    input:
        f"{params_list['output']}/celltype_annoted.rds"
    output:
        "log/diffexp.log"
    params:
        output = params_list['output']+"/04_diffexp",
        logs="log/diffexp.log",
        report_log="log/report_need.log",
        pattle = pattle,
        pvalue = config["subcluster"]["Diff_pvalue"],
        split = config["subcluster"]["Diff_split"],
        logfc = config["subcluster"]["Diff_logfc"],
        fc = config["subcluster"]["Diff_FC"],
        contract = config["subcluster"]["Diff_contract"]
    shell:
        "echo 运行的代码是: >> {params.report_log};"
        "echo 'Rscript  ./script/genesky_singlcell_tools.r -i {input} -o {params.output}   diffexp --pattle {params.pattle} --pvalue {params.pvalue} --splitby {params.split} --logfc {params.logfc} --FC {params.fc} --contrasts {params.contract}  && echo rule diffexp OK > {params.logs}' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ if Rscript  ./script/genesky_singlcell_tools.r -i {input} -o {params.output}   diffexp --pattle {params.pattle} --pvalue {params.pvalue} --splitby {params.split} --logfc {params.logfc} --FC {params.fc} --contrasts {params.contract} 2>> {params.report_log} ;"
        "then echo '[SUCCESS] rule diffexp OK' >> {params.report_log} && echo '[SUCCESS] rule diffexp OK' > {params.logs} ; "
        "else echo '[ERROR] rule diffexp FAILED' >> {params.report_log}; exit 1; fi; }}"





rule diffexp_result_arranger:
    input:
        log = "log/diffexp.log" 
    output:
        allgene=expand("{outdir}/05_Enrichment/{contract_split}/all_gene_merged.xls",outdir=params_list['output'],contract_split=contract_splits),
        up =expand("{outdir}/05_Enrichment/{contract_split}/Up_gene_merged.xls",outdir=params_list['output'],contract_split=contract_splits),
        down=expand("{outdir}/05_Enrichment/{contract_split}/Down_gene_merged.xls",outdir=params_list['output'],contract_split=contract_splits),
        log = "log/diffexp_result_arranger.log"
    params:
        report_log="log/report_need.log",
        diff_out = directory(f"{params_list['output']}/04_diffexp/"),
        output = params_list['output']+"/05_Enrichment"
    shell:
        "echo 运行的代码是: >> {params.report_log};"
        "echo 'Rscript ./script/Diff_file_arranger.r --input {params.diff_out} --output {params.output}'  >> {params.report_log} ;"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ if Rscript  ./script/Diff_file_arranger.r --input {params.diff_out} --output {params.output} 2>> {params.report_log} ;"
        "then echo '[SUCCESS] rule diffexp_result_arranger OK' >> {params.report_log} && touch {output.log}; "
        "else echo '[ERROR] rule diffexp_result_arranger FAILED' >> {params.report_log}; exit 1; fi; }}"

rule diff_enrichment:
    input:
        allgene=f"{params_list['output']}/05_Enrichment/{{contract_split}}/all_gene_merged.xls",
        up =f"{params_list['output']}/05_Enrichment/{{contract_split}}/Up_gene_merged.xls",
        down=f"{params_list['output']}/05_Enrichment/{{contract_split}}/Down_gene_merged.xls"
    output:
        # expand(directory("{outdir}/05_Enrichment/{contract_split}/{Reg}"),outdir=params_list['output'],contract_split=contract_splits,Reg =["ALL", "Up", "Down"]),
        log="log/diff_enrichment_{contract_split}.log"
    params:
        report_log="log/report_need.log",
        result_output = params_list['output']+"/05_Enrichment",
        topn=config["enrichment"]["topn"],
        rankby = config["enrichment"]["rankby"]
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo 'Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/ALL   Enrichments --file {params.result_output}/{wildcards.contract_split}/all_gene_merged.xls -s {species} --topn {params.topn} --rankby {params.rankby} && Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/Up   Enrichments --file {params.result_output}/{wildcards.contract_split}/Up_gene_merged.xls -s {species} --topn {params.topn} --rankby {params.rankby} && Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/Down   Enrichments --file {params.result_output}/{wildcards.contract_split}/Down_gene_merged.xls -s {species} --topn {params.topn} --rankby {params.rankby}' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "  if Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/ALL   Enrichments --file {params.result_output}/{wildcards.contract_split}/all_gene_merged.xls -s {species} --topn {params.topn} --rankby {params.rankby} && Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/Up   Enrichments --file {params.result_output}/{wildcards.contract_split}/Up_gene_merged.xls -s {species} --topn {params.topn} --rankby {params.rankby} && Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/Down   Enrichments --file {params.result_output}/{wildcards.contract_split}/Down_gene_merged.xls -s {species} --topn {params.topn} --rankby {params.rankby} "
        " 2>> {params.report_log} ; then "
        "    touch {output.log}; "
        "    echo '[SUCCESS] rule diffexp enrichment OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule diffexp enrichment FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule diff_GSEA:
    input:
        allgene=f"{params_list['output']}/05_Enrichment/{{contract_split}}/all_gene_merged.xls",
        up =f"{params_list['output']}/05_Enrichment/{{contract_split}}/Up_gene_merged.xls",
        down=f"{params_list['output']}/05_Enrichment/{{contract_split}}/Down_gene_merged.xls"
    output:
        # expand(directory("{outdir}/05_Enrichment/{contract_split}/{Reg}"),outdir=params_list['output'],contract_split=contract_splits,Reg =["ALL", "Up", "Down"]),
        log="log/diff_GSEA_{contract_split}.log"
    params:
        report_log="log/report_need.log",
        result_output = params_list['output']+"/06_Gene_GSEA",
        enrichment = params_list['output']+"/05_Enrichment"
    shell:
        "echo '运行的代码是:' >> {params.report_log};"
        "echo 'Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/ALL   GSEA --file {params.enrichment}/{wildcards.contract_split}/all_gene_merged.xls -s {species}  && Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/Up   GSEA --file {params.enrichment}/{wildcards.contract_split}/Up_gene_merged.xls -s {species}  && Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/Down   GSEA --file {params.enrichment}/{wildcards.contract_split}/Down_gene_merged.xls -s {species} ' >> {params.report_log};"
        "source /home/genesky/software/conda/4.9.2/bin/activate /home/donghj/snakemake && "
        "{{ "
        "  if Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/ALL   GSEA --file {params.enrichment}/{wildcards.contract_split}/all_gene_merged.xls -s {species}  && Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/Up   GSEA --file {params.enrichment}/{wildcards.contract_split}/Up_gene_merged.xls -s {species}  && Rscript ./script/genesky_singlcell_tools.r -o {params.result_output}/{wildcards.contract_split}/Down   GSEA --file {params.enrichment}/{wildcards.contract_split}/Down_gene_merged.xls -s {species}  "
        " 2>> {params.report_log} ; then "
        "    touch {output.log}; "
        "    echo '[SUCCESS] rule diffexp enrichment OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule diffexp enrichment FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"