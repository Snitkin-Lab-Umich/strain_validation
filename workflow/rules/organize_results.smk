# Author: Dhatri Badri

# Move synteny plots and variant information to final results folder
rule organize_results_dir:
    input:
        dotplot_results = "results/{prefix}/cauris_dotplot_repo/{sample}/plots/{sample}.pdf",
        snippy_results = "results/{prefix}/snippy/{sample}/{sample}.csv",
    output:
        "results/{prefix}/final_results/{sample}.csv",
        "results/{prefix}/final_results/{sample}.pdf"
    params:
        cauris_repo_outdir = "results/{prefix}/cauris_dotplot_repo",
        dotplot_outdir = "results/{prefix}/dotplots",
        final_results_outdir = "results/{prefix}/final_results",
        sample = "{sample}"
    shell:
        """
        mv {input.dotplot_results} {params.final_results_outdir}
        mv {input.snippy_results} {params.final_results_outdir}
        if [ ! -d "{params.dotplot_outdir}" ]; then
            mkdir -p "{params.dotplot_outdir}"
        fi
        mv {params.cauris_repo_outdir}/{params.sample} {params.dotplot_outdir}
        """

rule indel_nucmer:
    input:
        query_path=lambda wc: os.path.join(ASSEMBLY_DIR, f"{wc.sample}.fna"),
        subject_path=lambda wc: sample_ref_genome_dict[wc.sample],
        snippy_report = "results/{prefix}/final_results/{sample}.csv",
    output:
        ref_based_diff = "results/{prefix}/final_results/{sample}.rdiff",
        dnadiff_summary = "results/{prefix}/final_results/{sample}.report",
    params:
        sample="{sample}",
        dotplot_outdir = "results/{prefix}/dotplots",
        final_results_outdir = "results/{prefix}/final_results",
    envmodules:
       "Bioinformatics",
       "mummer/4.0.0rc1",
    shell:
        """
        cd {params.dotplot_outdir}/{params.sample}/nucmer
        dnadiff -p {params.sample} {input.subject_path} {input.query_path}
        cd ../../../../..
        mv {params.dotplot_outdir}/{params.sample}/nucmer/{params.sample}.rdiff {params.final_results_outdir}
        mv {params.dotplot_outdir}/{params.sample}/nucmer/{params.sample}.report {params.final_results_outdir}
        """

