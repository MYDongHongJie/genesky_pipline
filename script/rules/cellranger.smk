rule fastqc_raw:
    input:
        r1=lambda wildcards: f"/home/pub/project/research/{contract_number}/{wildcards.sample}_R1.fastq.gz"
            if platforms == "10X"
            else f"/home/pub/project/research/{contract_number}/{wildcards.sample}_cDNA_R1.fastq.gz",
        r2=lambda wildcards: f"/home/pub/project/research/{contract_number}/{wildcards.sample}_R2.fastq.gz"
            if platforms == "10X"
            else f"/home/pub/project/research/{contract_number}/{wildcards.sample}_cDNA_R2.fastq.gz"
    output:
        r1_dir = directory(f"{result_output}/fataQC/{{sample}}/R1_fastqc"),
        r2_dir = directory(f"{result_output}/fataQC/{{sample}}/R2_fastqc"),
        log = f"log/{{sample}}_fastqc.log"
    threads: 10
    priority: 100
    params:
        use_fastqc = config["software"]["fastqc"]
    shell:
        """
        {params.use_fastqc} -t {threads} --extract -o {TMP_DIR} {input.r1} {input.r2}

        # 复制 R1 文件夹内容
        mkdir -p {output.r1_dir}
        for d in {TMP_DIR}/{wildcards.sample}_R1_fastqc {TMP_DIR}/{wildcards.sample}_cDNA_R1_fastqc; do
            if [ -d "$d" ]; then
                cp -r "$d"/* {output.r1_dir}/
            fi
        done

        # 复制 R2 文件夹内容
        mkdir -p {output.r2_dir}
        for d in {TMP_DIR}/{wildcards.sample}_R2_fastqc {TMP_DIR}/{wildcards.sample}_cDNA_R2_fastqc; do
            if [ -d "$d" ]; then
                cp -r "$d"/* {output.r2_dir}/
            fi
        done

        echo 'fastqc done' > {output.log}
        """




rule fqchk_raw:
    input:
        r1=lambda wildcards: f"/home/pub/project/research/{contract_number}/{wildcards.sample}_R1.fastq.gz"
            if platforms == "10X"
            else f"/home/pub/project/research/{contract_number}/{wildcards.sample}_cDNA_R1.fastq.gz",
        r2=lambda wildcards: f"/home/pub/project/research/{contract_number}/{wildcards.sample}_R2.fastq.gz"
            if platforms == "10X"
            else f"/home/pub/project/research/{contract_number}/{wildcards.sample}_cDNA_R2.fastq.gz"
    output:
        r1_quality = f"{result_output}/fataQC/{{sample}}/R1.quality.txt",
        r2_quality = f"{result_output}/fataQC/{{sample}}/R2.quality.txt",
        log= f"log/{{sample}}_fqchk.log"
    params:
        use_seqtk = config["software"]["Seqtk"]
    priority: 100
    shell:
        """
        {params.use_seqtk} fqchk {input.r1} > {output.r1_quality} 
        {params.use_seqtk} fqchk {input.r2} > {output.r2_quality}
        echo 'fqchk done' > {output.log}
        """

rule fastp:
    input:
        r1=lambda wildcards: f"/home/pub/project/research/{contract_number}/{wildcards.sample}_R1.fastq.gz"
            if platforms == "10X"
            else f"/home/pub/project/research/{contract_number}/{wildcards.sample}_cDNA_R1.fastq.gz",
        r2=lambda wildcards: f"/home/pub/project/research/{contract_number}/{wildcards.sample}_R2.fastq.gz"
            if platforms == "10X"
            else f"/home/pub/project/research/{contract_number}/{wildcards.sample}_cDNA_R2.fastq.gz"
    output:
        r1_clean = f"{TMP_DIR}/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2_clean = f"{TMP_DIR}/{{sample}}_S1_L001_R2_001.fastq.gz",
        html = f"{result_output}/fataQC/{{sample}}/fastp.html",
        json = f"{result_output}/fataQC/{{sample}}/fastp.json",
        log = f"log/{{sample}}_fastp.log"
    threads: 10
    priority: 100
    params:
        use_fastp = config["software"]["fastp"],
        length_required = 28 if platforms == "10X" else 150,
        fp_exp = "" if platforms == "10X" else " --disable_adapter_trimming "
    shell:
        "{params.use_fastp} --thread {threads} --cut_right --length_required {params.length_required} --n_base_limit 1 --qualified_quality_phred 20 --unqualified_percent_limit 40 --dont_eval_duplication --in1 {input.r1} --in2 {input.r2} --out1 {output.r1_clean} --out2 {output.r2_clean} --html {output.html} --json {output.json} {params.fp_exp} && echo 'fastp done' > {output.log}"

rule fastqc_clean:
    input:
        r1 = rules.fastp.output.r1_clean,
        r2 = rules.fastp.output.r2_clean
    output:
        r1_dir = directory(f"{result_output}/fataQC/{{sample}}/final_R1_fastqc"),
        r2_dir = directory(f"{result_output}/fataQC/{{sample}}/final_R2_fastqc"),
        log = f"log/{{sample}}_fastqc_clean.log"
    threads: 10
    priority: 100
    params:
        use_fastqc = config["software"]["fastqc"]
    shell:
        "{params.use_fastqc} -t {threads} --extract -o {TMP_DIR} {input.r1} {input.r2} && mkdir -p {output.r1_dir} && cp -r {TMP_DIR}/{wildcards.sample}_S1_L001_R1_001_fastqc/* {output.r1_dir}/ && mkdir -p {output.r2_dir} && cp -r {TMP_DIR}/{wildcards.sample}_S1_L001_R2_001_fastqc/* {output.r2_dir}/ && echo 'fastqc_clean done' > {output.log}"

rule fqchk_clean:
    input:
        r1 = rules.fastp.output.r1_clean,
        r2 = rules.fastp.output.r2_clean
    output:
        r1_quality = f"{result_output}/fataQC/{{sample}}/final_R1.quality.txt",
        r2_quality = f"{result_output}/fataQC/{{sample}}/final_R2.quality.txt",
        log = f"log/{{sample}}_fqchk_clean.log"
    priority: 100
    params:
        use_seqtk = config["software"]["Seqtk"]
    shell:
        "{params.use_seqtk} fqchk {input.r1} > {output.r1_quality} && {params.use_seqtk} fqchk {input.r2} > {output.r2_quality} && echo 'fqchk_clean done' > {output.log}"

rule stat_txt:
    input:
        json = rules.fastp.output.json
    output:
        txt = f"{result_output}/fataQC/{{sample}}/stat.txt",
        log = f"log/{{sample}}_stat.log"
    priority: 100
    shell:
        """
        python script/parse_fastp_json.py {input.json} {output.txt}
        echo -e 'stat_txt done' > {output.log}
        echo -e '序列质控 ... OK' > log/report_need.log
        """





if platforms == "10X":
    nobam = config['cellranger']['nobam']
    include_introns= config['cellranger']['include_introns']
    if not config['cellranger']['forcell'] is None :
        forcecells=config['cellranger']['forcell']
        cellranger_ext = f" --force-cells {forcecells}"
    else :
        cellranger_ext=""
    rule CellRangerCount:
        input:
            fastq_dir = lambda wildcards: get_fastq_path(wildcards)
        output:
            web_summary = f"{result_output}/cellranger/{{sample}}/web_summary.html",
            metrics_csv = f"{result_output}/cellranger/{{sample}}/metrics_summary.csv",
            cloupe_file = f"{result_output}/cellranger/{{sample}}/cloupe.cloupe",
            barcodes = f"{result_output}/cellranger/{{sample}}/filtered_feature_bc_matrix/barcodes.tsv.gz",
            features = f"{result_output}/cellranger/{{sample}}/filtered_feature_bc_matrix/features.tsv.gz",
            matrix = f"{result_output}/cellranger/{{sample}}/filtered_feature_bc_matrix/matrix.mtx.gz"
        threads: 12
        priority: 0
        resources:
            mem_mb=90000
            #slurm_partition=used_batch
        params: 
            no_create_bam = '--no-bam' if nobam else '',
            include_introns = "true" if include_introns else "false",
            report_log2 = "../../log/report_need.log",
            report_log = "log/report_need.log"
        shell:
            "echo 运行的代码是: >> {params.report_log};"
            "echo 'cd {result_output}/cellranger && rm -rf {wildcards.sample} && "
            "{cellranger_path} count {params.no_create_bam} --id={wildcards.sample} "
            "--localcores {threads} --localmem 80 --fastqs={input.fastq_dir} --sample={wildcards.sample} "
            "--transcriptome={ref} --include-introns={params.include_introns} {cellranger_ext} "
            "&& mv {wildcards.sample}/outs/* {wildcards.sample}/ && rm -rf {wildcards.sample}/outs' "
            ">> {params.report_log};"
            "cd {result_output}/cellranger && rm -rf {wildcards.sample}; "
            "if {cellranger_path} count {params.no_create_bam} --id={wildcards.sample} "
            "--localcores {threads} --localmem 80 --fastqs={input.fastq_dir} --sample={wildcards.sample} "
            "--transcriptome={ref} --include-introns={params.include_introns} {cellranger_ext} "
            "&& mv {wildcards.sample}/outs/* {wildcards.sample}/ && rm -rf {wildcards.sample}/outs "
            "2>> {params.report_log2}; then "
            "echo '[SUCCESS] cellranger {wildcards.sample} OK' >> {params.report_log2}; "
            "else echo '[ERROR] cellranger {wildcards.sample} FAILED' >> {params.report_log2}; exit 1; fi"

    #AGGR,只在10X中用
    localrules: cellranger_aggr_prep
    rule cellranger_aggr_prep:
        """
        getting cellranger aggr libraries
        """ 
        input:
            metrics = expand(
                f"{result_output}/cellranger/{{sample}}/metrics_summary.csv",
                sample=SAMPLES_NAME
            ),
            samples_file = config['header']['sample_info']
        params:
            cellranger_version = config['header']['cellranger_version']
        output:
            libraries_file = os.path.abspath("libraries.csv")
            #libraries_file = config['report_params']['libraries_file']
        run:
            samples = pd.read_csv(input.samples_file,sep ='\t',names=["group","sample"], dtype=str).set_index("sample", drop=False)
            if params.cellranger_version == "7.1.0":
                libraries = pd.DataFrame(columns=['sample_id', 'molecule_h5'])
                for sample in list(samples.index):
                    libraries.loc[sample, 'sample_id'] = sample
                    libraries.loc[sample, 'molecule_h5'] = f"{sample}/molecule_info.h5"
                    libraries.to_csv(output.libraries_file, encoding='utf-8', index=False)
            else:
                libraries = pd.DataFrame(columns=['library_id', 'molecule_h5'])
                for sample in list(samples.index):
                    libraries.loc[sample, 'library_id'] = sample
                    libraries.loc[sample, 'molecule_h5'] = f"{sample}/molecule_info.h5"
                    libraries.to_csv(output.libraries_file, encoding='utf-8', index=False)

    rule cellranger_aggr:
        input:
            libraries_file = rules.cellranger_aggr_prep.output.libraries_file
        output:
            log = "log/cellranger_aggr.log",
            aggr_web = f"{result_output}/cellranger/aggr/web_summary.html",
            aggr_cloup = f"{result_output}/cellranger/aggr/cloupe.cloupe",
            tmpfiles = temp([f"{result_output}/cellranger/aggr/_cmdline",
                        f"{result_output}/cellranger/aggr/_filelist",
                        f"{result_output}/cellranger/aggr/_finalstate",
                        f"{result_output}/cellranger/aggr/_invocation",
                        f"{result_output}/cellranger/aggr/_jobmode",
                        f"{result_output}/cellranger/aggr/_log",
                        f"{result_output}/cellranger/aggr/_mrosource",
                        f"{result_output}/cellranger/aggr/_perf",
                        f"{result_output}/cellranger/aggr/_sitecheck",
                        f"{result_output}/cellranger/aggr/_tags",
                        f"{result_output}/cellranger/aggr/_timestamp",
                        f"{result_output}/cellranger/aggr/_uuid",
                        f"{result_output}/cellranger/aggr/_vdrkill",
                        f"{result_output}/cellranger/aggr/_versions",
                        directory(f"{result_output}/cellranger/aggr/SC_RNA_AGGREGATOR_CS"),
                        f"{result_output}/cellranger/aggr/aggr.mri.tgz"])
        shell:
            """
            # 1. 清理并创建目录
            cd {result_output}/cellranger/ && \
            rm -rf aggr && \
                
            # 2. 运行cellranger
            {cellranger_path} aggr \
                --id=aggr \
                --csv={input.libraries_file} \
                --normalize=none && \
                
            # 3. 移动输出文件
            mv aggr/outs/count/*  ./aggr/outs
            mv aggr/outs/* ./aggr/ && \
            rm -rf aggr/outs && \
            rm -rf aggr/count && \
            # 4. 记录日志
            cd - && \
            echo "cellranger_aggr done" > {output.log}
            """
    #TODO: 这里的input里面有sample变量，识别不出来，要成显示的方式
    rule Cellranger_report_arrange_10X:
        input:
            r1_quality_raw=expand(f"{result_output}/fataQC/{{sample}}/R1.quality.txt",sample=SAMPLES_NAME),
            r2_quality_raw=expand(f"{result_output}/fataQC/{{sample}}/R2.quality.txt",sample=SAMPLES_NAME),
            r1_quality_clean = expand(f"{result_output}/fataQC/{{sample}}/final_R1.quality.txt",sample=SAMPLES_NAME),
            r2_quality_clean = expand(f"{result_output}/fataQC/{{sample}}/final_R2.quality.txt",sample=SAMPLES_NAME),
            r1_dir_raw = expand(f"{result_output}/fataQC/{{sample}}/R1_fastqc", sample=SAMPLES_NAME),
            r2_dir_raw = expand(f"{result_output}/fataQC/{{sample}}/R2_fastqc", sample=SAMPLES_NAME),
            r1_dir_clean = expand(f"{result_output}/fataQC/{{sample}}/final_R1_fastqc", sample=SAMPLES_NAME),
            r2_dir_clean = expand(f"{result_output}/fataQC/{{sample}}/final_R2_fastqc", sample=SAMPLES_NAME),
            web_summary = expand(f"{result_output}/cellranger/{{sample}}/web_summary.html",sample=SAMPLES_NAME),
            metrics_csv = expand(f"{result_output}/cellranger/{{sample}}/metrics_summary.csv",sample=SAMPLES_NAME),
            satat = expand(f"{result_output}/fataQC/{{sample}}/stat.txt",sample=SAMPLES_NAME),
            aggr_web = f"{result_output}/cellranger/aggr/web_summary.html"
        output:
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images/{sample}_base_content{ext}",sample=SAMPLES_NAME,ext=[".png",".pdf"],result_output=result_output),
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images/{sample}_final_base_content{ext}",sample=SAMPLES_NAME,ext=[".png",".pdf"],result_output=result_output),
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images/{sample}_error_rate{ext}",sample=SAMPLES_NAME,ext=[".png",".pdf"],result_output=result_output),
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images/{sample}_final_error_rate{ext}",sample=SAMPLES_NAME,ext=[".png",".pdf"],result_output=result_output),
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/CellRanger_Summary.xlsx",result_output=result_output)
        shell:
            '''
            Rscript ./script/Cellranger_report_arranger_10X.r  --read1 '{input.r1_quality_raw}'  --read2 '{input.r2_quality_raw}' -s '{SAMPLES_NAME}'  -f false  -o {result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images &&
            Rscript ./script/Cellranger_report_arranger_10X.r --read1 '{input.r1_quality_clean}' --read2 '{input.r2_quality_clean}' -s '{SAMPLES_NAME}' -f true -o {result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images &&
            Rscript ./script/Cellranger_report_arranger_10_web_mer.r --webs '{input.web_summary}' --metrics '{input.metrics_csv}' --samples '{SAMPLES_NAME}' -o {result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/  &&
            Rscript ./script/Cellranger_report_arranger_10X_fastqc_stat.r --r1_dir_raw '{input.r1_dir_raw}' --r2_dir_raw '{input.r2_dir_raw}' --r1_dir_final '{input.r1_dir_clean}' --r2_dir_final '{input.r2_dir_clean}' --stat '{input.satat}' --sample '{SAMPLES_NAME}' -o {result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/ && 
            cp -r {result_output}/cellranger/aggr/filtered_feature_bc_matrix {result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/ && \
            cp -r {result_output}/cellranger/aggr/web_summary.html {result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/filtered_feature_bc_matrix && \
            Rscript ./script/Cellranger_report_aggr.r --barcodes {result_output}/cellranger/aggr/filtered_feature_bc_matrix/barcodes.tsv.gz --csv {result_output}/cellranger/aggr/aggregation.csv --outfile {result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/filtered_feature_bc_matrix/cell_sample.txt
        '''




if platforms == "huadaC4":
    rule Dnbc4toolsCount:
        input:
            f1=f"{TMP_DIR}/{{sample}}_S1_L001_R1_001.fastq.gz",
            f2=f"{TMP_DIR}/{{sample}}_S1_L001_R2_001.fastq.gz",
            o1=f"/home/pub/project/research/{contract_number}/{{sample}}_Oligo_R1.fastq.gz",
            o2=f"/home/pub/project/research/{contract_number}/{{sample}}_Oligo_R2.fastq.gz"
        output:
            web_summary = f"{result_output}/cellranger/{{sample}}/output/{{sample}}_scRNA_report.html",
            metrics_csv = f"{result_output}/cellranger/{{sample}}/output/metrics_summary.xls",
            barcodes = f"{result_output}/cellranger/{{sample}}/output/filtered_feature_bc_matrix/barcodes.tsv.gz",
            features = f"{result_output}/cellranger/{{sample}}/output/filtered_feature_bc_matrix/features.tsv.gz",
            matrix = f"{result_output}/cellranger/{{sample}}/output/filtered_feature_bc_matrix/matrix.mtx.gz"
        threads: 12
        priority: 0
        resources:
            mem_mb=180000
            #slurm_partition=used_batch
        shell:
            "{used_dnbc4tools} rna run --threads 12 --expectcells 10000 --outdir {result_output}/cellranger/ --cDNAfastq1 {input.f1} --cDNAfastq2 {input.f2} --oligofastq1 {input.o1} --oligofastq2 {input.o2} --name {wildcards.sample} --genomeDir {ref} && rename filter_matrix filtered_feature_bc_matrix {result_output}/cellranger/{wildcards.sample}/output/*"

    #TODO: 这里的input里面有sample变量，识别不出来，要成显示的方式
    rule Cellranger_report_arrange_huadaC4:
        input:
            r1_quality_raw=expand(f"{result_output}/fataQC/{{sample}}/R1.quality.txt",sample=SAMPLES_NAME),
            r2_quality_raw=expand(f"{result_output}/fataQC/{{sample}}/R2.quality.txt",sample=SAMPLES_NAME),
            r1_quality_clean = expand(f"{result_output}/fataQC/{{sample}}/final_R1.quality.txt",sample=SAMPLES_NAME),
            r2_quality_clean = expand(f"{result_output}/fataQC/{{sample}}/final_R2.quality.txt",sample=SAMPLES_NAME),
            r1_dir_raw = expand(f"{result_output}/fataQC/{{sample}}/R1_fastqc", sample=SAMPLES_NAME),
            r2_dir_raw = expand(f"{result_output}/fataQC/{{sample}}/R2_fastqc", sample=SAMPLES_NAME),
            r1_dir_clean = expand(f"{result_output}/fataQC/{{sample}}/final_R1_fastqc", sample=SAMPLES_NAME),
            r2_dir_clean = expand(f"{result_output}/fataQC/{{sample}}/final_R2_fastqc", sample=SAMPLES_NAME),
            web_summary = expand(f"{result_output}/cellranger/{{sample}}/output/{{sample}}_scRNA_report.html",sample=SAMPLES_NAME),
            metrics_csv = expand(f"{result_output}/cellranger/{{sample}}/output/metrics_summary.xls",sample=SAMPLES_NAME),
            satat = expand(f"{result_output}/fataQC/{{sample}}/stat.txt",sample=SAMPLES_NAME)
        output:
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images/{sample}_base_content{ext}",sample=SAMPLES_NAME,ext=[".png",".pdf"],result_output=result_output),
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images/{sample}_final_base_content{ext}",sample=SAMPLES_NAME,ext=[".png",".pdf"],result_output=result_output),
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images/{sample}_error_rate{ext}",sample=SAMPLES_NAME,ext=[".png",".pdf"],result_output=result_output),
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images/{sample}_final_error_rate{ext}",sample=SAMPLES_NAME,ext=[".png",".pdf"],result_output=result_output),
            expand("{result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/CellRanger_Summary.xlsx",result_output=result_output)
        shell:
            '''
            Rscript ./script/Cellranger_report_arranger_10X.r  --read1 '{input.r1_quality_raw}'  --read2 '{input.r2_quality_raw}' -s '{SAMPLES_NAME}'  -f false  -o {result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images &&
            Rscript ./script/Cellranger_report_arranger_10X.r --read1 '{input.r1_quality_clean}' --read2 '{input.r2_quality_clean}' -s '{SAMPLES_NAME}' -f true -o {result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/Quality_Images &&
            Rscript ./script/Cellranger_report_arranger_10_web_mer.r --platforms huadaC4 --webs '{input.web_summary}' --metrics '{input.metrics_csv}' --samples '{SAMPLES_NAME}' -o {result_output}/report/SCGES_单细胞转录组测序分析报告/02.CellRanger/  &&
            Rscript ./script/Cellranger_report_arranger_10X_fastqc_stat.r --r1_dir_raw '{input.r1_dir_raw}' --r2_dir_raw '{input.r2_dir_raw}' --r1_dir_final '{input.r1_dir_clean}' --r2_dir_final '{input.r2_dir_clean}' --stat '{input.satat}' --sample '{SAMPLES_NAME}' -o {result_output}/report/SCGES_单细胞转录组测序分析报告/01.Quality_Statistics/ 
        '''


rule Cellranger_report:
    input: 
        metrics = expand(
            f"{result_output}/cellranger/{{sample}}/metrics_summary.csv",
            sample=SAMPLES_NAME
        ) if platforms == "10X" else expand(
            f"{result_output}/cellranger/{{sample}}/output/metrics_summary.xls",
            sample=SAMPLES_NAME
        ),
        status = expand(f"{result_output}/fataQC/{{sample}}/stat.txt",sample=SAMPLES_NAME)
    output:
        log = "log/cellranger_report_check.log",
        qc_log = "log/qc.txt"
    params:
        cell_count = config['cellranger']['cell_count'],
        base_count_gt = config['cellranger']['base_count_gt'],
        base_count_le = config['cellranger']['base_count_le']
    shell:
         """
            script/check_cellranger_result.py \
                --input {input.metrics} \
                --output {output.log} \
                --platform {platforms} \
                --sample {SAMPLES_NAME} \
                --output_file {output.qc_log} \
                --qcfile {input.status} \
                --base_count_gt {params.base_count_gt} \
                --base_count_le {params.base_count_le} \
                --cell_count {params.cell_count} 
        """








rule wechat_notice:
    input:
        log = "log/cellranger_report_check.log"
    output:
        log = "log/wechat_notice.log"
    params:
        log_simple = lambda wildcards, input, output: log_simplify(input.log)
    shell:
        "python script/wechat_notice.py \
                --AT {person} \
                --task {contract_number} \
                --status success \
                --log {params.log_simple} && \
        echo 'wechat_notice done' > {output.log}"

rule wechat_notice2:
    input:
        rules.Cellranger_report.output.qc_log
    output:
        log = "log/wechat_notice2.log"
    shell:
        "python script/wechat_notice2.py \
                --qc_files {input} \
                --platforms {platforms} \
                --result_dir {abs_path} \
                --project {contract_number} &&  \
        echo 'wechat_notice2 done' > {output.log}"