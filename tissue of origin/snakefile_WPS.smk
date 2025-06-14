from snakemake.utils import validate
import pandas as pd
from os.path import exists

configfile: "config/cfdna_FDP.yml"


validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")
regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
validate(regions, schema="workflow/schemas/regions.schema.yaml")


def get_WPS_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_WPS.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_COV_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_COV.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_STARTS_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_STARTS.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_WPS_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{ref_SAMPLE}_WPS.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_COV_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{ref_SAMPLE}_COV.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_STARTS_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/background/{{GENOME}}/table/{{target_region}}--{ref_SAMPLE}_STARTS.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_length(input):
    if exists(input):
        df = pd.read_csv(input, sep="\t", header=None)
        length = df[2] - df[1]
        return length[0]
    else:
        length=2000
        #print(f"{input} does not exist: using length of {length}")
        return length


rule all:
    input:
        expand(
            expand(
                "results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
                zip,
                ID=samples["ID"],
                GENOME=samples["genome_build"],
                allow_missing=True,
            ),
            target_region=regions["target"],
        ),
        expand(
            expand(
                "results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_background_regions.bed.gz",
                zip,
                ID=samples["ID"],
                GENOME=samples["genome_build"],
                allow_missing=True,
            ),
            target_region=regions["target"],
        ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_COV.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_STARTS.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),
        expand(
            expand(
                "results/plots/overlays/{ID}/{target_region}--{SAMPLE}_overlays.pdf",
                zip,
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                allow_missing=True,
            ),
            target_region=regions["target"],
        ),

rule add_flanks:
    input:
        Region_file=lambda wildcards: regions["path"][wildcards.target_region],
    output:
        "results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_flanking.bed.gz",
    params:
        flank=1000,
    shell:
        """
        (fname={input.Region_file};
        if [[ $fname == *.gz ]]; 
        then zcat $fname; 
        else cat $fname; 
        fi) | 
        awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{start = $2 - {params.flank}; if (start < 0) start = 0; end = $3 + {params.flank}; print $1, start, end, $4, $5, $6}}' |
        gzip -c > {output}
        """


rule exclude_blacklist:
    input:
        Region_file_gz="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_flanking.bed.gz",
        blacklist=lambda wildcards: config[wildcards.GENOME]["universal_blacklist"]
    output:
        Region_file="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_flanking.bed",
        Region_file2="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz"
    conda:
        "workflow/envs/background.yml"
    shell:
        """
        #set -euo pipefail  # Ensure strict bash mode for error handling
        echo "Decompressing {input.Region_file_gz} to {output.Region_file}"
        gzip -d -c {input.Region_file_gz} > {output.Region_file}
        
        echo "Running bedtools intersect"           
        bedtools intersect -v -a {output.Region_file} -b {input.blacklist} | gzip -c > {output.Region_file2}
        
        echo "Output generated at {output.Region_file2}"
        """

rule generate_random_background:
    input:
        region="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
        genome=lambda wildcards: config[wildcards.GENOME]["genome_autosomes"],
        gap=lambda wildcards: config[wildcards.GENOME]["universal_blacklist"],
    output:
        "results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_background_regions.bed.gz",
    params:
        length=lambda wildcards, input: get_length(input.region),
    conda:
        "workflow/envs/background.yml"
    shell:
        """
        bedtools random -n 1000 -l {params.length} -g {input.genome} | \
        bedtools shuffle -i stdin -g {input.genome} -excl {input.gap} -noOverlapping | \
        gzip -c > {output}
        echo "End generate random background"
        """

rule extract_counts:
    input:
        target="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS.csv.gz",
        COV="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_COV.csv.gz",
        STARTS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_STARTS.csv.gz",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_%s.csv.gz",
    conda:
        "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/WPS/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.target} \
        -o {params.out_pre} {input.BAMFILE}
        """


rule extract_counts_background:
    input:
        background="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_background_regions.bed.gz",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.background.csv.gz",
        COV="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.background.csv.gz",
        STARTS="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.background.csv.gz",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_%s.background.csv.gz",
    conda:
        "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/WPS/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.background} \
        -o {params.out_pre} {input.BAMFILE}
        """


rule plot_overlays:
    input:
        WPS=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{{SAMPLE}}_WPS.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        WPS_ref=lambda wildcards: get_WPS_ref(wildcards.SAMPLE),
        COV=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{{SAMPLE}}_COV.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        COV_ref=lambda wildcards: get_COV_ref(wildcards.SAMPLE),
        WPS_back=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{{SAMPLE}}_WPS.background.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        WPS_back_ref=lambda wildcards: get_WPS_background_ref(wildcards.SAMPLE),
        COV_back=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{{SAMPLE}}_COV.background.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        COV_back_ref=lambda wildcards: get_COV_background_ref(wildcards.SAMPLE),
    output:
        "results/plots/overlays/{ID}/{target_region}--{SAMPLE}_overlays.pdf",
    params:
        target="{target_region}",
        sample="{SAMPLE}",
        ref_IDs=lambda wildcards: samples["ref_samples"][wildcards.SAMPLE].split(","),
        overlay_mode = config["plotting"]["overlay_mode"],
        smoothing = config["plotting"]["smoothing"],
        rolling = config["plotting"]["rolling"],
        background_norm = config["plotting"]["background_norm"]
    conda:
        "workflow/envs/overlays.yml"
    script:
        "workflow/scripts/WPS/overlays.py"
