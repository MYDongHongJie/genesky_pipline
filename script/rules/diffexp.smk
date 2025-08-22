
rule diffexp:
    input:
        f"{params_list['output']}/celltype_annoted.rds"
    output:
        
    shell:
        "Rscript {params_list['script']}/diffexp.R {input[0]} {input[1]} {output[0]} {params_list['output']}"

