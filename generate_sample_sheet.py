import pandas as pd
import shutil
from pathlib import Path
import glob
from tqdm import tqdm
import os
import argparse

# ------------- Set up command-line argument parser ---------------------------------------
parser = argparse.ArgumentParser(description="Move and log files for assembly processing")
parser.add_argument('--dryrun', action='store_true', help='Simulate file operations without actually moving files')
args = parser.parse_args()
dryrun = args.dryrun

# ------------- Define paths ---------------------------------------------------------------
base_dir = Path.cwd()
plasmidsaurus_data = base_dir / "Plasmidsaurus_data"
ont_out = base_dir / "ONT_assemblies"
illumina_out = base_dir / "illumina_reads"

# Create output directories if they don't exist
ont_out.mkdir(parents=True, exist_ok=True)
illumina_out.mkdir(parents=True, exist_ok=True)

lookup = pd.read_csv(plasmidsaurus_data / "2025-05-20_curing_batch_1_sample_lookup.csv")

# Load and prepare GBK assemblies
gbk_assemblies = pd.read_csv(base_dir / "MDHHS_hybrid_gbf_paths.txt", sep="\t")
gbk_assemblies = gbk_assemblies.rename(columns={"Assembly_path": "Assembly_path_gbk"})

# Load and prepare FASTA assemblies
fasta_assemblies = pd.read_csv(base_dir / "MDHHS_hybrid_genome_assembly_paths.txt", sep="\t")
fasta_assemblies = fasta_assemblies.drop(columns=["Sample_ID"], errors="ignore")  # drop if exists
fasta_assemblies = fasta_assemblies.rename(columns={"Assembly_path": "Assembly_path_fasta"})

# Merge both dataframes on Illumina_genome_ID and ONT_genome_ID
assembly = pd.merge(
    gbk_assemblies,
    fasta_assemblies,
    on=["Illumina_genome_ID", "ONT_genome_ID"],
    how="outer"
)

sample_sheet = []

# -------------- Main processing loop -------------------------------------------------
for _, row in tqdm(lookup.iterrows(), total=lookup.shape[0], desc="📦 Processing samples"):
    sample_name = row['Sample_ID']
    ref_genome = row['Reference_Genome_ID']
    ont_batch = row['ONT_batch']
    ont_sub_id = row['ONT_Submission_ID']
    illumina_batch = row['Illumina_batch']
    illumina_sub_id = row['Illumina_Submission_ID']

    # -------------- 1. Move and rename ONT assembly -------------------------------------------------
    ont_dir_pattern = f"{ont_batch}_results/{ont_batch}_*_{ont_sub_id}/ONT-only/annotation"
    ont_dirs = list(plasmidsaurus_data.glob(ont_dir_pattern))

    if ont_dirs:
        annotation_dir = ont_dirs[0]
        fna_pattern = f"{ont_batch}_*_{ont_sub_id}.fna"
        fna_files = list(annotation_dir.glob(fna_pattern))
    else:
        tqdm.write(f"❌ No matching directory found for: {ont_dir_pattern}")
        fna_files = []

    if fna_files:
        fna_file = fna_files[0]
        dest_path = ont_out / f"{sample_name}.fna"
        if not dest_path.exists():
            if dryrun:
                tqdm.write(f"[DRYRUN] Would move {fna_file} -> {dest_path}")
            else:
                shutil.copy(fna_file, dest_path)
        else:
            tqdm.write(f"⚠️ File already exists: {dest_path}")
    else:
        tqdm.write(f"⚠️ ONT .fna file not found for {sample_name} ({ont_batch}, {ont_sub_id})")

    # ------------- 2. Move and rename Illumina fastqs ----------------------------------------------------
    illumina_dir = plasmidsaurus_data / f"{illumina_batch}_Illumina_fastq"
    r1 = list(illumina_dir.glob(f"{illumina_batch}_*_{illumina_sub_id}_illumina_R1.fastq.gz"))
    r2 = list(illumina_dir.glob(f"{illumina_batch}_*_{illumina_sub_id}_illumina_R2.fastq.gz"))

    illumina_r1_name, illumina_r2_name = "", ""
    if r1 and r2:
        r1_file = r1[0]
        r2_file = r2[0]
        dest_r1 = illumina_out / f"{sample_name}_R1.fastq.gz"
        dest_r2 = illumina_out / f"{sample_name}_R2.fastq.gz"

        if dryrun:
            tqdm.write(f"[DRYRUN] Would move {r1_file} -> {dest_r1}")
            tqdm.write(f"[DRYRUN] Would move {r2_file} -> {dest_r2}")
        else:
            shutil.copy(r1_file, dest_r1)
            shutil.copy(r2_file, dest_r2)

        illumina_r1_name = dest_r1.name
        illumina_r2_name = dest_r2.name
    else:
        missing = []
        if not r1:
            missing.append("R1")
        if not r2:
            missing.append("R2")
        tqdm.write(f"⚠️ Warning: ❓ Missing Illumina {' and '.join(missing)} fastq file(s) for sample {sample_name} in batch {illumina_batch}, submission ID {illumina_sub_id}")

    # --------- 3. Match reference genome ---------------------------------------
    ont_assembly_name = f"{sample_name}.fna" if fna_files else ""

    matched_assembly = assembly[assembly["Illumina_genome_ID"] == ref_genome]
    if matched_assembly.empty:
        tqdm.write(f"⚠️ No matching assembly found for reference genome {ref_genome}")
        continue

    matched_row = matched_assembly.iloc[0]

    row_data = {
        "Reference_genome": ref_genome,
        "Sample_name": sample_name,
        "Illumina_F": illumina_r1_name,
        "Illumina_R": illumina_r2_name,
        "ONT_assembly": ont_assembly_name,
    }

    # Add whichever paths are available
    if pd.notna(matched_row.get("Assembly_path_gbk")):
        row_data["ref_genome_path_gbk"] = matched_row["Assembly_path_gbk"]
        # tqdm.write(f"✅ GBK path used for {sample_name}")
    if pd.notna(matched_row.get("Assembly_path_fasta")):
        row_data["ref_genome_path_fasta"] = matched_row["Assembly_path_fasta"]
        # tqdm.write(f"✅ FASTA path used for {sample_name}")

    sample_sheet.append(row_data)

# ----- Save the sample sheet ----------------------------------
if not dryrun:
    df = pd.DataFrame(sample_sheet)
    df.to_csv("sample_sheet.csv", index=False)
