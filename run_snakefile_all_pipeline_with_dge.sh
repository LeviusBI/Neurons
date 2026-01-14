#!/bin/bash
#SBATCH --job-name=rnaseq_pipeline   
#SBATCH --cpus-per-task=10           
#SBATCH --mem=50G                     
#SBATCH --time=2-0                    
#SBATCH --output=rnaseq_%j.out        
#SBATCH --error=rnaseq_%j.err

# ---------- Salmon pipeline ----------
snakemake \
    --snakefile /mnt/tank/scratch/ezaitseva/snakefile_all_pipeline_with_dge.py \
    --configfile /mnt/tank/scratch/ezaitseva/config.yaml \
    --config quantification=salmon \
    --cores 10 \
    --use-conda \
    --rerun-incomplete \
    --latency-wait 30 \
    --printshellcmds

# ---------- STAR + HTSeq pipeline ----------
snakemake \
    --snakefile /mnt/tank/scratch/ezaitseva/snakefile_all_pipeline_with_dge.py \
    --configfile /mnt/tank/scratch/ezaitseva/config.yaml \
    --config quantification=star_htseq \
    --cores 10 \
    --use-conda \
    --rerun-incomplete \
    --latency-wait 30 \
    --printshellcmds



