# Python 函数：拼接 Monocle2 命令
def build_monocle2_cmd(input_seurat_rds, output_dir, method,
                       idents=None, nselect=None, ngroup=None,
                       pcount=None, pratio=None, maxcell=None, root=None):
    """
    根据输入参数拼接 RunMonocle2.R 的命令行
    """
    cmd = [
        "Rscript /home/donghj/monocle/RunMonocle2.R",
        f"--rds {input_seurat_rds}",
        f"--method {method}",
        f"--output {output_dir}"
    ]
    if idents:
        cmd.append(f"--idents {idents}")
    if nselect:
        cmd.append(f"--nselect {nselect}")
    if ngroup:
        cmd.append(f"--ngroup {ngroup}")
    if pcount:
        cmd.append(f"--pcount {pcount}")
    if maxcell:
        cmd.append(f"--maxcell {maxcell}")
    if root:
        cmd.append(f"--root {root}")
    return " ".join(cmd)


# Snakemake Rule
if config['subcluster']['module'].get('monocle',False):
    monocle_cmd = build_monocle2_cmd(
            f"{params_list['output']}/celltype_annoted.rds",
            params_list['output'] + "/Monocle/",
            method,
            idents=idents,
            nselect=nselect,
            ngroup=ngroup,
            pcount=pcount,
            maxcell=maxcell,
            root=root
        )
    rule Monocle2Run:
        input:
            seurat_rds = f"{params_list['output']}/celltype_annoted.rds"
        output:
            rds = M2_REQUIRED
        resources:
            mem_mb=20000,
            mpi="pmi2"
        params:
            outdir = params_list['output'] + "/Monocle/",
            method = "VariableFeatures",
            idents = idents,
            nselect = nselect,
            ngroup = ngroup,
            pcount = pcount,
            maxcell = maxcell,
            root = root,
            report_log = "log/report_need.log",
            cmd = monocle_cmd
        shell:
            """
            echo "运行的代码是:" >> {params.report_log}
            echo "  source /home/genesky/software/conda/23.3.1/bin/activate /home/genesky/software/monocle/2.30.1 && {params.cmd}" >> {params.report_log}
              source /home/genesky/software/conda/23.3.1/bin/activate /home/genesky/software/monocle/2.30.1 && \
            {{ \
                if {params.cmd} 2>> {params.report_log}; then \
                    echo "[SUCCESS] rule Monocle2Run OK" >> {params.report_log}; \
                else \
                    echo "[ERROR] rule Monocle2Run FAILED" >> {params.report_log}; \
                    exit 1; \
                fi \
            }}
            """



rule Monocle2Plot:
    input:
        M2_REQUIRED
    output:
        M2P_REQUIRED
    params:
        report_log = "log/report_need.log",
        outdir = params_list['output']+"/Monocle/",  # 你之前定义的目录变量
    resources:
        mem_mb=20000,
        mpi="pmi2"
        #slurm_partition=used_batch
    shell:
        "echo '运行的代码是:' >> {params.report_log}; "
        "echo   source /home/genesky/software/conda/23.3.1/bin/activate /home/genesky/software/monocle/2.30.1 &&"
        "echo ' Rscript   /home/donghj/monocle/PlotMonocle2.R --rds {input[0]} --diff {input[2]} --output {params.outdir} ' >> {params.report_log}; "
        "  source /home/genesky/software/conda/23.3.1/bin/activate /home/genesky/software/monocle/2.30.1 && "
        "{{ "
        "if Rscript   /home/donghj/monocle/PlotMonocle2.R --rds {input[0]} --diff {input[2]} --output {params.outdir} "
        " 2>> {params.report_log}; then "
        "    echo '[SUCCESS] rule Monocle2Plot OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule Monocle2Plot FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

rule Monocle2BEAM:
    input:
        M2_REQUIRED
    output:
        BEAM_REQUIRED
    resources:
        mem_mb=50000,
        mpi="pmi2"
    params:
        report_log = "log/report_need.log",
        outdir = params_list['output']+"/Monocle/BEAM"  # 你之前定义的目录变量
    shell:
        "echo '运行的代码是:' >> {params.report_log}; "
        "echo   source /home/genesky/software/conda/23.3.1/bin/activate /home/genesky/software/monocle/2.30.1 &&"
        "echo ' Rscript   /home/donghj/monocle/BEAMMonocle2.R --rds {input[0]} --species {species} --output {params.outdir} ' >> {params.report_log}; "
        "  source /home/genesky/software/conda/23.3.1/bin/activate /home/genesky/software/monocle/2.30.1 && "
        "{{ "
        "if Rscript   /home/donghj/monocle/BEAMMonocle2.R --rds {input[0]} --species {species} --output {params.outdir} "
        " 2>> {params.report_log}; then "
        "    echo '[SUCCESS] rule Monocle2BEAM OK' >> {params.report_log}; "
        "  else "
        "    echo '[ERROR] rule Monocle2BEAM FAILED' >> {params.report_log}; "
        "    exit 1; "
        "  fi "
        "}}"

