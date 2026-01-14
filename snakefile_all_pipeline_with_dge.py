import pandas as pd
import os

# ---------- load config ----------
configfile: "config.yaml"

samples = pd.read_csv(
    os.path.join(config["work_dir"], config["sra_list"])
)["accession"].tolist()

samples.sort()

METHOD = config["quantification"]

# choose R script for DGE

if config["quantification"] == "salmon":
    rna_script = config["work_dir"] + "/rna-seq-analisis_salmon.R"
        
elif config["quantification"] == "star_htseq":
    rna_script = config["work_dir"] + "/rna-seq-analisis_htseq.R"
else:
    raise ValueError("quantification must be 'salmon' or 'star_htseq'")

# choose DGE inputs 

if config["quantification"] == "salmon":
    dge_input = expand(
        config["work_dir"] + "/salmon/{sample}/quant.sf",
        sample=samples
    )

elif config["quantification"] == "star_htseq":
    dge_input = expand(
        config["work_dir"] + "/htseq/{sample}.counts",
        sample=samples
    )

else:
    raise ValueError("quantification must be 'salmon' or 'star_htseq'")

# final target

rule all:
    input:
        up = os.path.join(
            config["work_dir"],
            "dif_expression_results",
            f"upregulation_{METHOD}.txt"
        ),
        down = os.path.join(
            config["work_dir"],
            "dif_expression_results",
            f"downregulation_{METHOD}.txt"
        )

# data download
# download fastq files using kingfisher

rule get_fastq:
    output: 
        config["work_dir"] + "/reads/raw/{sample}.fastq.gz"
    conda: "snake-files/envs/kingfisher.yaml"
    threads:
        config["kingfisher_threads"]
    shell: 
        """
        mkdir -p {config[work_dir]}/reads/raw
        kingfisher get \
            -r {wildcards.sample} \
            -m aws-http \
            -f fastq.gz \
            --download-threads {threads} \
            --output-directory {config[work_dir]}/reads/raw
        """

# quality control

rule fastp:
    input:
        R1 = config["work_dir"] + "/reads/raw/{sample}.fastq.gz"
    threads: 
        config["fastp_threads"]
    conda: "snake-files/envs/fastp.yaml"
    output:
        R1 = config["work_dir"] + "/reads/raw_filtered/{sample}.fastq.gz",
        html = config["work_dir"] + "/reports/fastp/{sample}.html"
    shell:
        """
        mkdir -p {config[work_dir]}/reads/raw_filtered
        mkdir -p {config[work_dir]}/reports/fastp
        fastp \
            -w {threads} \
            -i {input.R1} \
            -o {output.R1} \
            --length_required 25 \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 40 \
            --cut_front --cut_tail \
            --cut_right \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --n_base_limit 0 \
            --low_complexity_filter \
            --correction \
            --html {output.html} \
            2> {config[work_dir]}/reports/fastp/{wildcards.sample}.log || exit 1
        """

# build salmon index

rule salmon_index:
    input:
        transcripts=config["salmon"]["transcripts"]
    output:
        directory(config["salmon"]["index_dir"])
    threads:
        config["salmon"]["threads"]
    conda:
        "/mnt/tank/scratch/ezaitseva/snake-files/envs/salmon.yaml"
    shell:
        """
        mkdir -p {output}
        salmon index \
            -t {input.transcripts} \
            -i {output} \
            -k {config[salmon][kmer]}
        """

# salmon quantification

rule salmon:
    input: 
        index = config["salmon"]["index_dir"],
        R1 = config["work_dir"] + "/reads/raw_filtered/{sample}.fastq.gz"
    threads: 
        config["salmon"]["threads"]
    conda: "/mnt/tank/scratch/ezaitseva/snake-files/envs/salmon.yaml"
    params: 
        outdir=config["work_dir"] + "/salmon/{sample}"
    output: 
        config["work_dir"] + "/salmon/{sample}/quant.sf"
    shell: 
        """
        mkdir -p {params.outdir}
        salmon quant \
            -i {input.index} \
            -l U \
            -r {input.R1} \
            -o {params.outdir} \
            --validateMappings \
            --gcBias \
            --seqBias \
            --threads {threads} \
        """

# STAR

STAR_BIN = "/nfs/home/ezaitseva/STAR-2.7.11b/bin/Linux_x86_64/STAR"

# STAR index

rule star_index:
    input:
        genome=config["reference_star"]["genome_fasta"],
        gtf=config["reference_star"]["gtf"]
    output:
        directory(config["star"]["index"])
    threads:
        config["star"]["threads"]
    conda:
        "/mnt/tank/scratch/ezaitseva/snake-files/envs/star.yaml"
    shell:
        """
        mkdir -p {output}
        {STAR_BIN} \
          --runMode genomeGenerate \
          --runThreadN {threads} \
          --genomeDir {output} \
          --genomeFastaFiles {input.genome} \
          --sjdbGTFfile {input.gtf} \
          --sjdbOverhang {config[star][sjdbOverhang]}
        """

# STAR ALIGNMENT

rule star:
    input: 
        index = config["star"]["index"], 
        R1 = config["work_dir"] + "/reads/raw_filtered/{sample}.fastq.gz"
    threads:
        config["star"]["threads"]
    conda: "/mnt/tank/scratch/ezaitseva/snake-files/envs/star.yaml"
    params: 
        prefix = config["work_dir"] + "/star/{sample}_"
    output: 
        config["work_dir"] + "/star/{sample}_Aligned.sortedByCoord.out.bam"
    shell: 
        """
        mkdir -p {config[work_dir]}/star
        {STAR_BIN} \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --genomeDir {input.index} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn {input.R1} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 15000000000 \
            --outSAMunmapped Within \
            --outSAMattributes Standard
        """

# HTSeq COUNTING

rule htseq:
    input: 
        bam = config["work_dir"] + "/star/{sample}_Aligned.sortedByCoord.out.bam", 
        annotation = config["htseq"]["annot_dir"] + "/gencode.v49.chr_patch_hapl_scaff.annotation.gtf"
    output: 
        config["work_dir"] + "/htseq/{sample}.counts"
    conda: "/mnt/tank/scratch/ezaitseva/snake-files/envs/htseq.yaml"
    threads:
        config["htseq"]["threads"]
    shell: 
        """
        mkdir -p {config[work_dir]}/htseq
        htseq-count \
            -n {threads} \
            --idattr gene_id \
            -s no \
            --add-chromosome-info \
            --additional-attr=ID {input.bam} {input.annotation} > {output}
        """

# differential gene expression
rule dge:
    input:
        counts = dge_input,
        metadata = config["dge"]["metadata"]
    output:
        up = os.path.join(config["work_dir"], config["dge"]["output_dir"], f"upregulation_{METHOD}.txt"),
        down = os.path.join(config["work_dir"], config["dge"]["output_dir"], f"downregulation_{METHOD}.txt")
    conda: "rna_seq_r"  
    script:
        rna_script

