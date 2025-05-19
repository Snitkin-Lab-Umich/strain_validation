# Strain validation
Snakemake pipeline to validate mutants made in the lab

## Summary
This pipeline performs two  steps:
- Firstly, it creates synteny dot plots using long reads based on an in house pipeline [cauris_dotplot](https://github.com/Snitkin-Lab-Umich/cauris_dotplot) 
- Secondly, it runs [snippy](https://github.com/tseemann/snippy), a variant calling tool, that takes paired end reads and a reference genome.

The workflow generates all the output in the `prefix` folder set in  `config/config.yaml`. Each workflow steps gets its own individual folder as shown below. This structure provides a general view of how outputs are organized, with each tool or workflow step having its own directory. **Note that this overview does not capture all possible outputs from each tool; it only highlights the primary directories and _SOME_ of their contents.** 

```
results/
└── 2025-05-12_Test_Strain_Validation_Pipeline
    ├── dotplots
    │   └── BXYRFN_3_sample_03
    │       ├── BXYRFN_3_sample_03_debug_log.txt
    │       ├── contig_data
    │       │   ├── BXYRFN_3_sample_03.fna_contig_data.csv
    │       │   └── MI_KPC_112_flye_medaka_polypolish_contig_data.csv
    │       ├── nucmer
    │       │   ├── BXYRFN_3_sample_03.fna_to_MI_KPC_112_flye_medaka_polypolish.coord
    │       │   └── BXYRFN_3_sample_03.fna_to_MI_KPC_112_flye_medaka_polypolish.delta
    │       └── plots
    │           └── BXYRFN_3_sample_03.fna_to_MI_KPC_112_flye_medaka_polypolish.pdf
    └── snippy
        └── BXYRFN_3_sample_03_illumina
            ├── BXYRFN_3_sample_03_illumina.aligned.fa
            ├── BXYRFN_3_sample_03_illumina.bam
            ├── BXYRFN_3_sample_03_illumina.bam.bai
            ├── BXYRFN_3_sample_03_illumina.bed
            ├── BXYRFN_3_sample_03_illumina.csv
            ├── BXYRFN_3_sample_03_illumina.filt.vcf
            ├── BXYRFN_3_sample_03_illumina.gff
            ├── BXYRFN_3_sample_03_illumina.html
            ├── BXYRFN_3_sample_03_illumina.log
            ├── BXYRFN_3_sample_03_illumina.raw.vcf
            ├── ref.fa -> reference/ref.fa
            └── ref.fa.fai -> reference/ref.fa.fai
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

> Load bioinformatics, snakemake, singularity and R modules from Great Lakes modules.

```

module load Bioinformatics snakemake singularity mummer/4.0.0rc1 R

```

Check to see if you have ggplot2 installed in your home directory 

> Open an R session by typing the below command on the terminal.

```
R
```
> Once you open an R session, load the ggplot2 library

```
library(ggplot2)
```

If you do not see any message(s) show up after loading the library, that means you have it installed in your home, quit the session `q()` and skip to this [section](#setup-config-and-samples-files). However if you see this message `Error in library(ggplot2) : There is no package called 'ggplot2'`, you have to install it in your home directory manually with this command. 

```
install.packages("ggplot2", repos = "https://repo.miserver.it.umich.edu/cran/")
```

>Check to see if the installation worked by running this command

```
library(ggplot2)
```

Quit the R session `q()` after successful installation of `ggplot2` and move onto the next section of the pipeline.  

## Setup config, sample and cluster files

**_If you are just testing this pipeline, the config and sample files are already loaded with test data, so you do not need to make any additional changes to them. However, it is a good idea to change the prefix (name of your output folder) in the config file to give you an idea of what variables need to be modified when running your own samples._**

### Customize config.yaml and set tool specific parameters
As an input, the snakemake file takes a config file where you can set the path to `sample_sheet.csv`, path to ONT long reads and illumina short reads, etc. Instructions on how to modify `config/config.yaml` is found in `config.yaml`. 

### Samples
This sample file should contain 6 columns: `Reference_genome`, `Reference_genome_path`, `Sample_name`,  `Illumina_F`, `Illumina_R` and `ONT_assembly`. An example of how the `sample_sheet.csv` should look like can be found here `/nfs/turbo/umms-esnitkin/Project_MIDGE_Bac/Analysis/Plasmid_curing/2025-04-16_EXAMPLE_CURING_EXP/sample_sheet.csv`.

The reference genome column should contain the name of the reference genome, the reference_genome_path should contain the path to the reference genome which can be found here `/nfs/turbo/umms-esnitkin/Project_MIDGE_Bac/Analysis/Plasmid_curing/MDHHS_hybrid_genome_assembly_paths.txt`. The Sample_name should be the name of the reference genome with a suffix `_C`. `Illumina_F` and `Illumina_R` are the names of the paired end reads with the suffix `_R1/R2.fastq.gz`. ONT_assembly is the name of the long read. 


## Quick start

### Run pipeline on a set of samples.

>Preview the steps in Strain validation by performing a dryrun of the pipeline.

```

snakemake -s workflow/Snakefile --dryrun

```

> Submit the pipeline as a batch job. Change these SBATCH commands: `--job-name` to a more descriptive name like `run_strainval`, `--mail-user` to your email address, `--time` depending on the number of samples you have (should be more than what you specified in `config/cluster.json`). Feel free to make changes to the other flags if you are comfortable doing so. The sbat script can be found in the current directory—it's called `StrainValidation.sbat`. Don't forget to submit the script to Slurm! `sbatch StrainValidation.sbat`.

```
#!/bin/bash

#SBATCH --job-name=strain_val
#SBATCH --mail-user=youremail@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --export=ALL
#SBATCH --partition=standard
#SBATCH --account=esnitkin1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=3 --mem=10g --time=08:15:00

# Load necessary modules
module load Bioinformatics snakemake singularity R mummer/4.0.0rc1

# Run pipeline
snakemake -s workflow/Snakefile --use-envmodules -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 100 --nolock
```

## DAG of pipeline
![Alt text](images/dag.svg)

## Dependencies

### Near Essential
* [Snakemake>=7.32.4](https://snakemake.readthedocs.io/en/stable/#)

### Tool stack used in workflow

* [R](https://www.r-project.org/)
* [snippy](https://github.com/tseemann/snippy)
* [nucmer](https://github.com/mummer4/mummer)
* [Pandas](https://pandas.pydata.org/)

