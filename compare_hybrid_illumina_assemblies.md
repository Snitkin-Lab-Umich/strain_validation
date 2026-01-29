# The pipeline format is the same as before

# Strain validation
Snakemake pipeline to find the additional or missing scaffols present in the cured isolates

**The pipeline assumes you have created the sample sheet found [here](README.md#samples).**

## Summary
This pipeline performs the following:

1. **Trimming & quality control** of reads
2. **Downsampling** to manage sequencing depth
3. **Assembly** (SPAdes) and filtering large contigs (>1kb)
4. **Genome scaffolding** using RagTag against reference genomes
5. **Dotplot analysis** (with `cauris_dotplot_pipeline`)

The workflow generates all the output in the `prefix` folder set in  `config/config_compare.yaml`. Each workflow steps gets its own individual folder as shown below. This structure provides a general view of how outputs are organized, with each tool or workflow step having its own directory. **_Note that this overview does not capture all possible outputs from each tool; it only highlights the primary directories and _SOME_ of their contents._** 

**The synteny plots can be found in the `final_results` folder.**

```
.
├── cauris_dotplot_repo
│   ├── dotplot.py
│   ├── highlight_data.tsv
│   ├── make_plots.R
│   ├── README.md
│   ├── run_dotplot.job
│   └── sloop.py
├── dotplots
│   └── KPC_653_A
│       ├── contig_data
│       │   ├── MI_KPC_653_flye_medaka_polypolish_contig_data.csv
│       │   └── ragtag.scaffold_contig_data.csv
│       ├── KPC_653_A_debug_log.txt
│       ├── nucmer
│       │   ├── KPC_653_A.coord
│       │   └── KPC_653_A.delta
│       └── plots
├── downsample
│   └── KPC_653_A
│       ├── KPC_653_A_R1_trim_paired.fastq.gz
│       └── KPC_653_A_R2_trim_paired.fastq.gz
├── final_results
│   └── KPC_653_A.pdf
├── RagTag
│   ├── KPC_30_A
│   │   ├── ragtag.scaffold.confidence.txt
│   │   ├── ragtag.scaffold.err
│   │   ├── ragtag.scaffold.fasta
│   │   └── ragtag.scaffold.stats
├── spades
│   └── KPC_653_A
│       ├── contigs.fasta
│       ├── KPC_653_A_contigs_l1000.fasta
│       ├── KPC_653_A_contigs_l1000.fasta.fai
│       ├── run_spades.sh
│       ├── run_spades.yaml
└── trimmomatic
    └── KPC_653_A
        ├── KPC_653_A_R1_trim_paired.fastq.gz
        ├── KPC_653_A_R1_trim_paired.fastq.gz_fastqchk.txt
        ├── KPC_653_A_R1_trim_unpaired.fastq.gz
        ├── KPC_653_A_R2_trim_paired.fastq.gz
        └── KPC_653_A_R2_trim_unpaired.fastq.gz
```

## Installation

>If you are using Great Lakes HPC, ensure you are cloning the repository in your scratch directory. Change `your_uniqname` to your uniqname.

```
cd /scratch/esnitkin_root/esnitkin1/your_uniqname/
```

> Clone the github directory onto your system.

```

git clone https://github.com/Snitkin-Lab-Umich/strain_validation.git

```

> Ensure you have successfully cloned `strain_validation`. Type `ls` and you should see the newly created directory **_strain_validation_**. Move to the newly created directory.

```

cd strain_validation

```
> To ensure a clean starting environment, deactivate and/or remove any modules/packages loaded.

```
module purge
conda deactivate
```

> Load bioinformatics, snakemake, singularity and R modules from Great Lakes modules.

```

module load Bioinformatics snakemake singularity R

```


### Cluster file
Change the walltime according to the number of samples you have in `config/cluster.json` to ensure the jobs are being submitted in a timely manner. Update `email` flag to your email.

## Quick start

### Run pipeline on a set of samples.

>Preview the steps in Strain validation by performing a dryrun of the pipeline.

```

snakemake -s workflow/compare_hyb_illum_assem.smk --dryrun

```

> Submit the pipeline as a batch job. Change these SBATCH commands: `--job-name` to a more descriptive name like `run_strainval`, `--mail-user` to your email address, `--time` depending on the number of samples you have (should be more than what you specified in `config/cluster.json`). Feel free to make changes to the other flags if you are comfortable doing so. The sbat script can be found in the current directory—it's called `StrainValidation_compare.sbat`. Don't forget to submit the script to Slurm! `sbatch StrainValidation_compare.sbat`.

```
#!/bin/bash

#SBATCH --job-name=compare_assemblies
#SBATCH --mail-user=youremail@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --export=ALL
#SBATCH --partition=standard
#SBATCH --account=esnitkin1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=3 --mem=10g --time=08:15:00

# Load necessary modules
module load Bioinformatics snakemake singularity R mummer/4.0.0rc1 Rtidyverse/4.4.3

# Run pipeline
snakemake -s workflow/compare_hyb_illum_assem.smk --use-conda --use-singularity --use-envmodules -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock
```

### Near Essential
* [Snakemake>=7.32.4](https://snakemake.readthedocs.io/en/stable/#)

### Tool stack used in workflow

- [Snakemake >=7.32.4](https://snakemake.readthedocs.io/en/stable/)
- [R](https://www.r-project.org/)
- [snippy](https://github.com/tseemann/snippy)
- [nucmer](https://github.com/mummer4/mummer)
- [Pandas](https://pandas.pydata.org/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [SPAdes](https://github.com/ablab/spades)
- [RagTag](https://github.com/malonge/RagTag)
- [Dotplot pipeline](https://github.com/Snitkin-Lab-Umich/cauris_dotplot_pipeline)
