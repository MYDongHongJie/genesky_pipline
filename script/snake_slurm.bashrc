 cp -r /home/donghj/scRNA_snakemake/10X_singlecell3/genesky_scrna_pipline ./

smq() {
    if [[ $# -lt 1 ]]; then
        echo "Usage: smq <Snakefile> [more_options]"
        return 1
    fi

    # Snakefile 路径
    snakefile="$1"
    snakefile=$(readlink -e "$snakefile")
    shift

    # 工作目录固定为当前目录
    wd=$(pwd)

    # 额外参数
    more="$@"

    echo "运行 snakemake文件: $snakefile 在: $wd 路径下，运行开始时间是 $(date +%Y年%m月%d日%H小时%M分) 额外参数是: $more"

    # 调用 Snakemake
    /home/donghj/snakemake/bin/snakemake -s "$snakefile" \
              -d "$wd" \
              --executor slurm \
              -k \
              $more 
}

smk() {
    /home/donghj/snakemake/bin/snakemake -npr -s "$1"
}

smu(){
if [[ $# -eq 0 ]]; then snakemake --unlock
else snakemake --unlock -s $1
fi
}