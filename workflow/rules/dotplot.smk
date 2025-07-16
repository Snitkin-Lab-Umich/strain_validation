# Author: Dhatri Badri

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

# Run dotplot for each query
rule run_dotplot:
    input:
        query_path=lambda wc: os.path.join(ASSEMBLY_DIR, f"{wc.sample}.fna"),
        subject_path=lambda wc: sample_ref_genome_dict[wc.sample],
        script="results/{prefix}/cauris_dotplot_repo/dotplot.py",
    output:
        "results/{prefix}/cauris_dotplot_repo/{sample}/plots/{sample}.pdf"
    params:
        sample="{sample}"
    envmodules:
       "Bioinformatics",
       "mummer/4.0.0rc1",
       "R/4.4.0",
       "python3.10-anaconda/2023.03",
       "Rtidyverse/4.4.3"
    shell:
        """
        python {input.script} \
            --name {params.sample} \
            --query {input.query_path} \
            --subject {input.subject_path}
        """


