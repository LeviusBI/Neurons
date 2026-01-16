# RNA-seq Analysis Pipeline

A flexible Snakemake-based pipeline for RNA-seq data analysis with support for two quantification methods: Salmon (pseudoalignment) and STAR+HTSeq (alignment-based).

## Features

- **Flexible configuration**: All parameters are mandatory (no defaults) for explicit control
- **Automatic directory management**: Each run creates a timestamped directory (`HH_MM_SS_DD_MM_YYYY`)
- **Two quantification methods**: 
  - Salmon (fast pseudoalignment)
  - STAR + HTSeq (alignment-based counting)
- **Differential gene expression**: Automated DGE analysis using DESeq2 in Docker container
- **Dynamic configuration**: Update configuration via `configure` rule or Python script

## Pipeline Workflow

1. **Data download**: Download FASTQ files from SRA using Kingfisher
2. **Quality control**: Filter and trim reads using fastp
3. **Quantification**: 
   - Salmon: Pseudoalignment and transcript quantification
   - STAR+HTSeq: Genome alignment and gene counting
4. **Differential expression**: DESeq2 analysis in Docker container

## Requirements

- Snakemake
- Conda (for environment management)
- Docker (for DGE analysis)
- Python 3.6+
- R (in Docker container)

## Quick Start

### 1. Configure the pipeline

**Option A: Using the `configure` rule (recommended)**

```bash
snakemake configure \
    --config \
        base_work_dir=/path/to/work \
        reference_star.gtf=/path/to/gtf.gtf \
        quantification=salmon
```

**Option B: Using the configuration script**

```bash
python3 scripts/configure_pipeline.py \
    --base-work-dir /path/to/work \
    --gtf /path/to/gtf.gtf \
    --quantification salmon \
    --sra-list /path/to/sra_list.txt \
    --metadata /path/to/metadata.txt
```

**Option C: Manual configuration**

Edit `config.yaml` and fill in all required parameters (see `config.yaml` template).

### 2. Prepare input files

**SRA list file** (`sra_list.txt`):
```tsv
accession
SRR123456
SRR123457
SRR123458
```

**Metadata file** (`meta_data.txt`):
```tsv
sampleid    individual    diff_time
SRR123456   1             NPC
SRR123457   1             2 week neuron
SRR123458   2             NPC
SRR123459   2             2 week neuron
```

Required columns:
- `sampleid`: Sample identifier (must match SRA accession)
- `individual`: Replicate/donor factor
- `diff_time`: Experimental condition

### 3. Run the pipeline

```bash
snakemake \
    --snakefile snakefile_all_pipeline_with_dge.py \
    --configfile config.yaml \
    --cores 10 \
    --use-conda
```

## Configuration

All parameters in `config.yaml` are **mandatory** (no default values). The pipeline validates all parameters at startup.

### Key Parameters

- `base_work_dir`: Base directory where timestamped run directories will be created
- `sra_list`: Path to file with SRA accession IDs
- `quantification`: Method (`salmon` or `star_htseq`)
- `reference_star.gtf`: Path to GTF annotation file
- `dge.metadata`: Path to experimental metadata file
- `docker.image`: Docker image with R scripts

See `config.yaml` for the complete list of required parameters.

## Output Structure

Each run creates a timestamped directory:

```
base_work_dir/
└── HH_MM_SS_DD_MM_YYYY/
    ├── reads/
    │   ├── raw/              # Downloaded FASTQ files
    │   └── raw_filtered/     # Quality-filtered reads
    ├── salmon/               # Salmon quantification results (if method=salmon)
    │   └── {sample}/
    │       └── quant.sf
    ├── star/                 # STAR alignment results (if method=star_htseq)
    │   └── {sample}_Aligned.out.bam
    ├── htseq/                # HTSeq count files (if method=star_htseq)
    │   └── {sample}.counts
    ├── reports/
    │   └── fastp/           # Quality control reports
    └── dif_expression_results/
        ├── upregulation_{method}.txt
        └── downregulation_{method}.txt
```

## Configuration Methods

### Method 1: Snakemake `configure` rule

Update configuration via command line:

```bash
snakemake configure \
    --config \
        base_work_dir=/new/path \
        reference_star.gtf=/new/path/to/gtf.gtf \
        quantification=star_htseq
```

### Method 2: Python script

```bash
python3 scripts/configure_pipeline.py \
    --base-work-dir /path/to/work \
    --gtf /path/to/gtf.gtf \
    --quantification salmon
```

### Method 3: Manual editing

Edit `config.yaml` directly.

## Pipeline Rules

### Data Processing

- `get_fastq`: Download FASTQ files from SRA using Kingfisher
- `fastp`: Quality control and filtering
- `salmon_index`: Build Salmon index (if using Salmon)
- `salmon`: Salmon quantification
- `star_index`: Build STAR index (if using STAR+HTSeq)
- `star`: STAR alignment (outputs unsorted BAM for HTSeq)
- `htseq`: HTSeq counting (if using STAR+HTSeq)
- `dge`: Differential gene expression analysis in Docker container

### Configuration

- `configure`: Update configuration file dynamically

## Technical Details

### Docker Integration

The DGE analysis runs in a Docker container with automatic volume mounting:
- Work directory
- GTF annotation file
- Metadata directory

R scripts are located at `/home/scripts/` inside the container.

## Validation

The pipeline automatically validates:
- All required parameters are set
- `quantification` has a valid value (`salmon` or `star_htseq`)
- All file paths are specified

If any parameter is missing, the pipeline will report an error with a list of missing parameters.

## Troubleshooting

### Configuration errors

1. All parameters in `config.yaml` are filled (no `null` values)
2. All file paths are correct and accessible
3. Conda environment paths are valid

### File not found errors

Ensure that:
- SRA list file exists and is readable
- Metadata file exists and has correct format
- Reference files (GTF, genome, transcripts) are accessible
- Conda environment YAML files exist

## Examples

### Full workflow: Configure and run

```bash
# 1. Configure
snakemake configure \
    --config \
        base_work_dir=/data/rnaseq \
        reference_star.gtf=/data/annotations/gencode.v49.gtf \
        quantification=salmon \
        sra_list=/data/sra_list.txt \
        dge.metadata=/data/meta_data.txt

# 2. Run pipeline
snakemake \
    --snakefile snakefile_all_pipeline_with_dge.py \
    --configfile config.yaml \
    --cores 10 \
    --use-conda \
    --rerun-incomplete
```

### Switch quantification method

```bash
snakemake configure --config quantification=star_htseq
snakemake --configfile config.yaml --cores 10 --use-conda
```

## Authors

- Lev
- Katya
- Ekaterina


