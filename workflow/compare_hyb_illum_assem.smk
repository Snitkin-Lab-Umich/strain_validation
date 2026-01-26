# Author: Dhatri Badri

configfile: "config/config_compare.yaml"

import pandas as pd
import os
import shutil, glob

# Load configuration
samples_df = pd.read_csv(config["samples"])
PREFIX = config["prefix"]
ASSEMBLY_DIR = config["assembly"]
SHORT_READS_DIR = config["short_reads"]

SAMPLE = list(samples_df["Sample_name"])

sample_ref_genome_dict = dict(zip(samples_df["Sample_name"], samples_df["ref_genome_path_fasta"]))
# sample_ref_genome_dict_gbk = dict(zip(samples_df["Sample_name"], samples_df["ref_genome_path_gbk"]))

if not os.path.exists("results/" + PREFIX):
    os.makedirs("results/" + PREFIX, exist_ok=True)

rule all:
    input:
        # expand("results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta", prefix=PREFIX, sample=SAMPLE),
        # expand("results/{prefix}/cauris_dotplot_repo/{sample}/plots/{sample}.pdf", prefix=PREFIX, sample=SAMPLE),
        # expand("results/{prefix}/cauris_dotplot_repo/dotplot.py", prefix=PREFIX),
        expand("results/{prefix}/RagTag/{sample}/ragtag.scaffold.fasta", prefix=PREFIX, sample=SAMPLE),
        expand("results/{prefix}/final_results/{sample}.pdf", prefix=PREFIX, sample=SAMPLE),

rule trimmomatic_pe:
    input:    
        r1=lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R1.fastq.gz"),
        r2=lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R2.fastq.gz"),
    output:
        r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_unpaired.fastq.gz",
    params:
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
        threads = config["ncores"],
    log:
        "logs/{prefix}/trimmomatic/{sample}/{sample}.log"
    #conda:
    #    "envs/trimmomatic.yaml"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    #envmodules:
    #    "Bioinformatics",
    #    "trimmomatic"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {params.threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}"

rule downsample:
    input:
        r1 ="results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz",
    output:
        outr1 = "results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz",
        outr2 = "results/{prefix}/downsample/{sample}/{sample}_R2_trim_paired.fastq.gz",
    params:
        gsize = config["genome_size"],
    # run:
    #     downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})
    wrapper:
        "file:workflow/wrapper_functions/downsample"

rule assembly:
    input:
        r1 = "results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/downsample/{sample}/{sample}_R2_trim_paired.fastq.gz",
    output:
        spades_assembly = "results/{prefix}/spades/{sample}/contigs.fasta", 
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        threads= 10,
        prefix_dir = "results/{prefix}/",
        sample = "{sample}"
    #conda:
    #    "envs/spades.yaml"
    singularity:
        "docker://staphb/spades:4.0.0"
    #envmodules:
    #    "Bioinformatics",
    #    "spades/4.0.0"
    shell:
        "spades.py --isolate --pe1-1 {input.r1} --pe1-2 {input.r2} -o {params.out_dir} --threads {params.threads}"

rule bioawk:
    input:
        spades_assembly = "results/{prefix}/spades/{sample}/contigs.fasta"
    output:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta"
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        bioawk_params = config["bioawk"],
        prefix = "{sample}"
    conda:
        "envs/bioawk.yaml"
    # singularity:
    #     "docker://lbmc/bioawk:1.0"
    shell:
        """
        workflow/scripts/bioawk.sh {input.spades_assembly} {output.spades_l1000_assembly} {params.out_dir} {params.prefix}
        """

# Install cauris_dotplot GitHub repo
rule install_cauris_dotplot:
    output:
        "results/{prefix}/cauris_dotplot_repo/dotplot.py",
    params:
        dotplot_github_outdir = "results/{prefix}/cauris_dotplot_repo"
    shell:
        """
        git clone https://github.com/Snitkin-Lab-Umich/cauris_dotplot_pipeline.git {params.dotplot_github_outdir}
        """

rule ragtag:
    input:
        # ragtag_script = "results/{prefix}/RagTag/ragtag.py",
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta"
    output:
        ragtag_scaffold = "results/{prefix}/RagTag/{sample}/ragtag.scaffold.fasta"
    params:
        ref_genome=lambda wc: sample_ref_genome_dict[wc.sample],
        sample_out_dir = "results/{prefix}/RagTag/{sample}/",
        # ragtag_out_dir = "results/{prefix}/RagTag",
        threads= 10,
        # basedir=workflow.basedir
    # conda:
    #     "envs/ragtag.yaml"
    singularity:
        "docker://quay.io/biocontainers/ragtag:2.1.0--pyhb7b1952_0"
    shell:
        """
        ragtag.py scaffold {params.ref_genome} {input.spades_l1000_assembly} -t {params.threads} -o {params.sample_out_dir}
        """

# Run dotplot for each query
rule run_dotplot:
    input:
        ragtag_scaffold = "results/{prefix}/RagTag/{sample}/ragtag.scaffold.fasta",
        subject_path=lambda wc: sample_ref_genome_dict[wc.sample],
        script= "results/{prefix}/cauris_dotplot_repo/dotplot.py",
    output:
        "results/{prefix}/cauris_dotplot_repo/{sample}/plots/{sample}.pdf"
    params:
        sample="{sample}",
        # table_of_regions=config["table_of_regions"],
        basedir=workflow.basedir
    envmodules:
       "Bioinformatics",
       "mummer/4.0.0rc1",
       "R/4.4.0",
       "python3.10-anaconda/2023.03",
       "Rtidyverse/4.4.3"
    shell:
        """
        strainval_dir=$(dirname {params.basedir})
        python $strainval_dir/{input.script} \
            --name {params.sample} \
            --query $strainval_dir/{input.ragtag_scaffold} \
            --subject {input.subject_path} 
        """
        #strainval_dir=$(dirname {params.basedir})

rule organize_results_dir:
    input:
        dotplot_results = "results/{prefix}/cauris_dotplot_repo/{sample}/plots/{sample}.pdf",
    output:
        "results/{prefix}/final_results/{sample}.pdf"
    params:
        cauris_repo_outdir = "results/{prefix}/cauris_dotplot_repo",
        dotplot_outdir = "results/{prefix}/dotplots",
        final_results_outdir = "results/{prefix}/final_results",
        sample = "{sample}"
    shell:
        """
        mv {input.dotplot_results} {params.final_results_outdir}
        if [ ! -d "{params.dotplot_outdir}" ]; then
            mkdir -p "{params.dotplot_outdir}"
        fi
        mv {params.cauris_repo_outdir}/{params.sample} {params.dotplot_outdir}
        """
