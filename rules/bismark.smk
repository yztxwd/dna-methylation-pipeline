rule bismark_mapping_pe:
    input:
        align_pe_find_input(wildcards)
    output:
        bam="output/mapped/{sample}-{rep}-{unit}_bismark_bt2_PE.bam",
        report="output/mapped/{sample}-{rep}-{unit}_bismark_bt2_PE_PE_report.txt"
    cluster:
        log="logs/bismark/{sample}-{rep}-{unit}.log"
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
        bismark {params.extra} -o output/mapped/ -S {params.basename} {params.index} -1 {input[0]} -2 {input[1]}
        """
    
rule bismark_mapping_se:
    input:
        align_se_find_input(wildcards)
    output:
        bam="output/mapped/{sample}-{rep}-{unit}_bismark_bt2_SE.bam",
        report="output/mapped/{sample}-{rep}-{unit}_bismark_bt2_SE_SE_report.txt"
    cluster:
        log="logs/bismark/{sample}-{rep}-{unit}.log"
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
        bismark {params.extra} -o output/mapped/ -S {params.basename} {params.index} {input}
        """

rule bismark_methylation_extractor:
    input:
        lambda wildcards: "output/mapped/{sample}-{rep}-{unit}_bismark_bt2_" + ("SE.bam" if is_single_end(**wildcards) else "PE.bam")
    output:
        CpG="output/bismark_methylation_extract/CpG_context_{sample}-{rep}-{unit}_bismark_bt2.txt",
        CHG="output/bismark_methylation_extract/CHG_context_{sample}-{rep}-{unit}_bismark_bt2.txt",
        CHH="output/bismark_methylation_extract/CHH_context_{sample}-{rep}-{unit}_bismark_bt2.txt"
    cluster:
        log="logs/bismark_methylation_extract/{sample}-{rep}-{unit}.log"
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
            --genome_folder {params.genome} {input}
        """

def bismark2bedGraph_find_input(wildcards):
    # keep rows with specified sample index
    slices = samples.loc[sample, :].reset_index()
    # get the bam files corresponding to the specified sample
    CpG_inputs = list(slices.apply(lambda row: f"output/bismark_methylation_extract/CpG_context_{wildcards.sample}-{row["rep"]}-{row["unit"]}_bismark_bt2.txt"))
    CHG_inputs = list(slices.apply(lambda row: f"output/bismark_methylation_extract/CHG_context_{wildcards.sample}-{row["rep"]}-{row["unit"]}_bismark_bt2.txt"))
    CHH_inputs = list(slices.apply(lambda row: f"output/bismark_methylation_extract/CHH_context_{wildcards.sample}-{row["rep"]}-{row["unit"]}_bismark_bt2.txt"))

    return CpG_inputs + CHG_inputs + CHH_inputs 

rule bismark2bedGraph_CpG:
    input:
        bismark2bedGraph_find_input
    output:
        "output/bismark_methylation_extract/{sample}_bismark_bt2.CpG.bedGraph.gz"
    cluster:
        log="logs/bismark_methylation_extract/{sample}.bismark2bedGraph.log"
    params:
        basename="{sample}_bismark_bt2.CpG.bedGraph"
        extra=config["bismark2bedGraph"]["params"]
    threads:
        4
    resources:
        nodes=1,
        mem=20
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark2bedGraph {params.extra} -o {params.basename} --dir output/bismark_methylation_extract {input}
        """
        
rule bismark2bedGraph_all:
    input:
        bismark2bedGraph_find_input
    output:
        "output/bismark_methylation_extract/{sample}_bismark_bt2.all.bedGraph.gz"
    cluster:
        log="logs/bismark_methylation_extract/{sample}.bismark2bedGraph.log"
    params:
        basename="{sample}_bismark_bt2.all.bedGraph"
        extra=config["bismark2bedGraph"]["params"]
    threads:
        4
    resources:
        nodes=1,
        mem=20
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark2bedGraph {params.extra} --CX -o {params.basename} --dir output/bismark_methylation_extract {input}
        """
