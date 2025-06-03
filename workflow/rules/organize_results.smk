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


