import pandas as pd
from snakemake.utils import validate, min_version

#### Set minimum snakemake version ####
min_version("5.1.2")

#### Load config and sample sheets ####

configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index(["sample", "rep", "unit"], drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

#### target rules ####

# default only report CpG, uncomment to report all cytosine 
rule all:
    input:
        "output/qc/multiqc/multiqc.html",
        expand("output/bismark_methylation_extract/{sample}_bismark_bt2.CpG.bedGraph.gz", sample = samples["sample"])
#        expand("output/bismark_methylation_extract/{sample}_bismark_bt2.all.bedGraph.gz", sample = samples["sample"])

#### setup singularity ####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

#### setup report ####
report: "report/workflow.rst"

#### load rules ####
include: "snakemake-pipeline-general/rules/global.smk"
include: "snakemake-pipeline-general/rules/qc.smk"
include: "snakemake-pipeline-general/rules/trim.smk"
include: "snakemake-pipeline-general/rules/fastp.smk"
include: "snakemake-pipeline-general/rules/bowtie2.smk"
include: "snakemake-pipeline-general/rules/bwa.smk"
include: "snakemake-pipeline-general/rules/filter.smk"
include: "snakemake-pipeline-general/rules/merge.smk"
include: "snakemake-pipeline-general/rules/coverage.smk"
include: "rules/bismark.smk"
