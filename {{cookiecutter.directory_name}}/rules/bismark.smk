rule bismark_mapping_pe:
    input:
        align_pe_find_input
    output:
        bam="output/mapped/{sample}-{rep}-{unit}_bt2_PE.bam",
        report="output/mapped/{sample}-{rep}-{unit}_bt2_PE_PE_report.txt"
    log:
        "logs/bismark/{sample}-{rep}-{unit}.log"
    params:
        index=config["bismark"]["index"],
        extra=config["bismark"]["params"],
        basename="{sample}-{rep}-{unit}_bismark_bt2_PE"
    threads:
        config["threads"]
    resources:
        nodes=config["nodes"],
        mem=config["mem"]
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark {params.extra} -o output/mapped/ -B {params.basename} {params.index} -1 {input[0]} -2 {input[1]} > {log}
        """
    
rule bismark_mapping_se:
    input:
        align_se_find_input
    output:
        bam="output/mapped/{sample}-{rep}-{unit}_bt2_SE.bam",
        report="output/mapped/{sample}-{rep}-{unit}_bt2_SE_SE_report.txt"
    log:
        "logs/bismark/{sample}-{rep}-{unit}.log"
    params:
        index=config["bismark"]["index"],
        extra=config["bismark"]["params"],
        basename="{sample}-{rep}-{unit}_bt2_SE"
    threads:
        config["threads"]
    resources:
        nodes=config["nodes"],
        mem=config["mem"]
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark {params.extra} -o output/mapped/ -B {params.basename}  {params.index} {input} > {log}
        """

    

rule bismark_methylation_extractor:
    input:
        lambda wildcards: "output/mapped/{sample}-{rep}-{unit}_bt2_{type}.bam"
    output:
        CpG="output/bismark_methylation_extract/CpG_context_{sample}-{rep}-{unit}_bt2_{type}.txt",
        CHG="output/bismark_methylation_extract/CHG_context_{sample}-{rep}-{unit}_bt2_{type}.txt",
        CHH="output/bismark_methylation_extract/CHH_context_{sample}-{rep}-{unit}_bt2_{type}.txt"
    log:
        "logs/bismark_methylation_extract/{sample}-{rep}-{unit}_{type}.log"
    params:
        genome=config["bismark"]["index"],
        extra=config["bismark_methylation_extractor"]["params"]
    threads:
        config["threads"]
    resources:
        nodes=config["nodes"],
        mem=config["mem"]
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark_methylation_extractor {params.extra} --parallel 6 -o output/bismark_methylation_extract \
            --genome_folder {params.genome} {input} > {log}
        """

def bismark2bedGraph_find_input(wildcards):
    # keep rows with specified sample index
    slices = samples.loc[wildcards.sample, :]
    # get the bam files corresponding to the specified sample
    CpG_inputs = []
    CHG_inputs = []
    CHH_inputs = []
    for row in slices.itertuples():
        type = "SE" if is_single_end(row.sample, row.rep, row.unit) else "PE"
        CpG_inputs.append(f"output/bismark_methylation_extract/CpG_context_{row.sample}-{row.rep}-{row.unit}_bt2_{type}.txt") 
        CHG_inputs.append(f"output/bismark_methylation_extract/CHG_context_{row.sample}-{row.rep}-{row.unit}_bt2_{type}.txt") 
        CHH_inputs.append(f"output/bismark_methylation_extract/CHH_context_{row.sample}-{row.rep}-{row.unit}_bt2_{type}.txt") 

    if wildcards.context == "CpG":
        return CpG_inputs
    elif wildcards.context == "CHG":
        return CHG_inputs
    elif wildcards.context == "CHH":
        return CHH_inputs

rule bismark2bedGraph_context:
    input:
        bismark2bedGraph_find_input
    output:
#        directory("output/bismark2bedGraph/{sample}_{context}/"),
        "output/bismark2bedGraph/{sample}_{context}/{sample}_bismark_bt2.{context}.bedGraph.gz"
    log:
        "logs/bismark_methylation_extract/{sample}_{context}.bismark2bedGraph.log"
    params:
        basename="{sample}_bismark_bt2.{context}.bedGraph",
        extra=config["bismark2bedGraph"]["params"]
    threads:
        config["threads"]
    resources:
        nodes=1,
        mem=20,
        walltime="40:00:00"
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark2bedGraph {params.extra} --CX -o {params.basename} --dir $(dirname {output}) {input}
        """

rule bismark_bedGraph2bigwig:
    input:
        bg = "output/bismark2bedGraph/{sample}_{context}/{sample}_bismark_bt2.{context}.bedGraph.gz",
        chrom_size = config["bedGraphToBigWig"]["chrom"]
    output:
        "output/bismark2bedGraph/{sample}_{context}/{sample}_bismark_bt2.{context}.bw"
    conda:
        "../snakemake-pipeline-general/envs/ucsc.yaml"
    resources:
        mem = 40,
        walltime = "12:00:00"
    shell:
        """
        gunzip -c {input.bg} | tail -n +2 | sort -k1,1 -k2,2n > {input.bg}.tmp
        bedGraphToBigWig {input.bg}.tmp {input.chrom_size} {output}
        """
        
rule bismark2bedGraph_all:
    input:
        bismark2bedGraph_find_input
    output:
        directory("output/bismark2bedGraph/{sample}_all/")
    log:
        "logs/bismark_methylation_extract/{sample}.bismark2bedGraph.log"
    params:
        basename="{sample}_bismark_bt2.all.bedGraph",
        extra=config["bismark2bedGraph"]["params"]
    threads:
        config["threads"]
    resources:
        nodes=1,
        mem=80,
        walltime="120:00:00"
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark2bedGraph {params.extra} --CX -o {params.basename} --dir {output} {input}
        """
