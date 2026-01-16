import pandas as pd
import os
from datetime import datetime

# Configure rule - update configuration via command line parameters
# IMPORTANT: configfile is loaded after this rule, so we use
# parameters directly via {config[...]} in shell command
rule configure:
    """
    Rule for dynamic configuration setup.
    Accepts parameters via --config and updates config.yaml
    
    Usage example:
    snakemake configure --config base_work_dir=/path/to/work reference_star.gtf=/path/to/gtf.gtf
    
    Or via separate script (recommended):
    python3 scripts/configure_pipeline.py --base-work-dir /path/to/work --gtf /path/to/gtf.gtf
    """
    output:
        config_file = "config.yaml"
    params:
        script = os.path.join(os.path.dirname(__file__), "scripts", "configure_pipeline.py")
    shell:
        """
        # Use Python script to update configuration
        # Parameters are passed via --config and available through {config[...]}
        python3 {params.script} \
            {config[base_work_dir]:+--base-work-dir } {config[base_work_dir]} \
            {config[reference_star][gtf]:+--gtf } {config[reference_star][gtf]} \
            {config[quantification]:+--quantification } {config[quantification]} \
            {config[sra_list]:+--sra-list } {config[sra_list]} \
            {config[dge][metadata]:+--metadata } {config[dge][metadata]} \
            {config[docker][image]:+--docker-image } {config[docker][image]}
        """

# Load and validate configuration
# Load configfile (may be updated by configure rule)
configfile: "config.yaml"

# Function to validate required parameters
def validate_config(config):
    """Validates that all required parameters are set"""
    required_params = {
        "base_work_dir": "Base working directory",
        "sra_list": "Path to SRA ID list file",
        "quantification": "Quantification method (salmon or star_htseq)",
        "salmon": {
            "transcripts": "Path to transcripts file",
            "index_dir": "Path to Salmon index",
            "kmer": "k-mer size",
            "threads": "Number of threads"
        },
        "reference_star": {
            "genome_fasta": "Path to genome file",
            "gtf": "Path to GTF file"
        },
        "star": {
            "index": "Path to STAR index",
            "threads": "Number of threads",
            "sjdbOverhang": "sjdbOverhang parameter",
            "star_bin": "Path to STAR executable"
        },
        "htseq": {
            "threads": "Number of threads"
        },
        "fastp_threads": "Number of threads for fastp",
        "kingfisher_threads": "Number of threads for kingfisher",
        "dge": {
            "metadata": "Path to metadata file",
            "output_dir": "Output directory name for DGE results"
        },
        "docker": {
            "image": "Docker image"
        },
        "conda_envs": {
            "kingfisher": "Path to kingfisher.yaml",
            "fastp": "Path to fastp.yaml",
            "salmon": "Path to salmon.yaml",
            "star": "Path to star.yaml",
            "htseq": "Path to htseq.yaml"
        }
    }
    
    errors = []
    
    # Check top-level parameters
    for key, desc in required_params.items():
        if isinstance(desc, dict):
            # Nested structure
            if key not in config:
                errors.append(f"Missing section: {key}")
            else:
                for subkey, subdesc in desc.items():
                    if subkey not in config[key] or config[key][subkey] is None:
                        errors.append(f"Missing parameter: {key}.{subkey} ({subdesc})")
        else:
            # Simple parameter
            if key not in config or config[key] is None:
                errors.append(f"Missing parameter: {key} ({desc})")
    
    if errors:
        raise ValueError("Configuration errors:\n" + "\n".join(f"  - {e}" for e in errors))
    
    # Check quantification value
    if config["quantification"] not in ["salmon", "star_htseq"]:
        raise ValueError(f"quantification must be 'salmon' or 'star_htseq', got: {config['quantification']}")
    
    return True

# Validate configuration
validate_config(config)

# Create timestamped directory for current run
timestamp = datetime.now().strftime("%H_%M_%S_%d_%m_%Y")
work_dir = os.path.join(config["base_work_dir"], timestamp)
os.makedirs(work_dir, exist_ok=True)

print(f"Work directory for this run: {work_dir}")

# Load sample list
sra_list_path = config["sra_list"]
if not os.path.isabs(sra_list_path):
    sra_list_path = os.path.join(work_dir, sra_list_path)

samples = pd.read_csv(sra_list_path)["accession"].tolist()
samples.sort()

METHOD = config["quantification"]

# Select R script for DGE
if METHOD == "salmon":
    rna_script = "rna-seq-analisis_salmon.R"
elif METHOD == "star_htseq":
    rna_script = "rna-seq-analisis_htseq.R"
else:
    raise ValueError("quantification must be 'salmon' or 'star_htseq'")

# Select input data for DGE
if METHOD == "salmon":
    dge_input = expand(
        os.path.join(work_dir, "salmon", "{sample}", "quant.sf"),
        sample=samples
    )
elif METHOD == "star_htseq":
    dge_input = expand(
        os.path.join(work_dir, "htseq", "{sample}.counts"),
        sample=samples
    )

# Final target
rule all:
    input:
        up = os.path.join(work_dir, config["dge"]["output_dir"], f"upregulation_{METHOD}.txt"),
        down = os.path.join(work_dir, config["dge"]["output_dir"], f"downregulation_{METHOD}.txt")

# Download FASTQ files using kingfisher
rule get_fastq:
    output: 
        fastq = os.path.join(work_dir, "reads", "raw", "{sample}.fastq.gz")
    conda: config["conda_envs"]["kingfisher"]
    threads: config["kingfisher_threads"]
    shell: 
        """
        mkdir -p $(dirname {output.fastq})
        kingfisher get \
            -r {wildcards.sample} \
            -m aws-http \
            -f fastq.gz \
            --download-threads {threads} \
            --output-directory $(dirname {output.fastq})
        
        # Rename file to expected name
        # Kingfisher may create {sample}_1.fastq.gz (paired-end) or {sample}.fastq.gz (single-end)
        if [ -f $(dirname {output.fastq})/{wildcards.sample}_1.fastq.gz ]; then
            mv $(dirname {output.fastq})/{wildcards.sample}_1.fastq.gz {output.fastq}
        elif [ ! -f {output.fastq} ]; then
            if [ -f $(dirname {output.fastq})/{wildcards.sample}.fastq.gz ]; then
                mv $(dirname {output.fastq})/{wildcards.sample}.fastq.gz {output.fastq}
            fi
        fi
        """

# Quality control with fastp
rule fastp:
    input:
        R1 = os.path.join(work_dir, "reads", "raw", "{sample}.fastq.gz")
    threads: config["fastp_threads"]
    conda: config["conda_envs"]["fastp"]
    output:
        R1 = os.path.join(work_dir, "reads", "raw_filtered", "{sample}.fastq.gz"),
        html = os.path.join(work_dir, "reports", "fastp", "{sample}.html")
    shell:
        """
        mkdir -p $(dirname {output.R1})
        mkdir -p $(dirname {output.html})
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
            --correction \
            --html {output.html} \
            2> {work_dir}/reports/fastp/{wildcards.sample}.log || exit 1
        """

# Build Salmon index
rule salmon_index:
    input:
        transcripts = config["salmon"]["transcripts"]
    output:
        index = directory(config["salmon"]["index_dir"])
    threads: config["salmon"]["threads"]
    conda: config["conda_envs"]["salmon"]
    shell:
        """
        mkdir -p {output}
        salmon index \
            -t {input.transcripts} \
            -i {output} \
            -k {config[salmon][kmer]}
        """

# Salmon quantification
rule salmon:
    input: 
        index = config["salmon"]["index_dir"],
        R1 = os.path.join(work_dir, "reads", "raw_filtered", "{sample}.fastq.gz")
    threads: config["salmon"]["threads"]
    conda: config["conda_envs"]["salmon"]
    params: 
        outdir = os.path.join(work_dir, "salmon", "{sample}")
    output: 
        quant = os.path.join(work_dir, "salmon", "{sample}", "quant.sf")
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
            --threads {threads}
        """

# STAR index
rule star_index:
    input:
        genome = config["reference_star"]["genome_fasta"],
        gtf = config["reference_star"]["gtf"]
    output:
        index = directory(config["star"]["index"])
    threads: config["star"]["threads"]
    conda: config["conda_envs"]["star"]
    params:
        star_bin = config["star"]["star_bin"]
    shell:
        """
        mkdir -p {output}
        {params.star_bin} \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {config[star][sjdbOverhang]}
        """

# STAR alignment
# HTSeq requires Unsorted BAM, so we use different settings
rule star:
    input: 
        index = config["star"]["index"], 
        R1 = os.path.join(work_dir, "reads", "raw_filtered", "{sample}.fastq.gz")
    threads: config["star"]["threads"]
    conda: config["conda_envs"]["star"]
    params: 
        prefix = os.path.join(work_dir, "star", "{sample}_"),
        star_bin = config["star"]["star_bin"],
        # HTSeq requires Unsorted BAM
        samtools = "samtools"  # assume samtools is available
    output: 
        bam = os.path.join(work_dir, "star", "{sample}_Aligned.out.bam")
    shell: 
        """
        mkdir -p $(dirname {output.bam})
        {params.star_bin} \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --genomeDir {input.index} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn {input.R1} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outSAMattributes Standard
        
        # Rename output file to expected name
        if [ -f {params.prefix}Aligned.out.bam ]; then
            mv {params.prefix}Aligned.out.bam {output.bam}
        fi
        """

# HTSeq counting
rule htseq:
    input: 
        bam = os.path.join(work_dir, "star", "{sample}_Aligned.out.bam"), 
        annotation = config["reference_star"]["gtf"]
    output: 
        counts = os.path.join(work_dir, "htseq", "{sample}.counts")
    conda: config["conda_envs"]["htseq"]
    threads: config["htseq"]["threads"]
    shell: 
        """
        mkdir -p $(dirname {output.counts})
        htseq-count \
            -f bam \
            -r pos \
            -n {threads} \
            --idattr gene_id \
            -s no \
            --add-chromosome-info \
            --additional-attr=ID {input.bam} {input.annotation} > {output.counts}
        """

# Differential gene expression
rule dge:
    input:
        counts = dge_input,
        metadata = config["dge"]["metadata"]
    output:
        up = os.path.join(work_dir, config["dge"]["output_dir"], f"upregulation_{METHOD}.txt"),
        down = os.path.join(work_dir, config["dge"]["output_dir"], f"downregulation_{METHOD}.txt")
    params:
        work_dir = work_dir,
        gtf_file = config["reference_star"]["gtf"],
        docker_image = config["docker"]["image"],
        # R script is located in the container
        r_script_path = os.path.join("/home/scripts", rna_script)
    shell:
        """
        # Mount work directory, GTF file and metadata directory
        docker run --rm \
            -v {params.work_dir}:{params.work_dir} \
            -v {params.gtf_file}:{params.gtf_file} \
            -v $(dirname {input.metadata}):$(dirname {input.metadata}) \
            -w {params.work_dir} \
            {params.docker_image} \
            /usr/bin/Rscript {params.r_script_path}
        """
