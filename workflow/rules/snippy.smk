# Author: Dhatri Badri

# Run variant calling (snippy) to get variant information        
rule run_snippy:
    input:
        r1=lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R1.fastq.gz"),
        r2=lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R2.fastq.gz"),
        ref=lambda wc: sample_ref_genome_dict[wc.sample],
    output:
        "results/{prefix}/snippy/{sample}/{sample}.csv",
    params:
        cpus = config["cores"],
        sample ="{sample}",
        outdir = "results/{prefix}/snippy/{sample}"
    # singularity:
    #     "docker://staphb/snippy:4.6.0"
    envmodules:
       "Bioinformatics",
       "snippy"
    shell:
        """       
        cleaned_ref=$(mktemp --suffix=.fasta)
        sed 's/;.*//' {input.ref} > $cleaned_ref
        snippy --cpus {params.cpus} \
               --prefix {params.sample} \
               --outdir {params.outdir} \
               --ref $cleaned_ref \
               --R1 {input.r1} \
               --R2 {input.r2} \
               --force
        """

