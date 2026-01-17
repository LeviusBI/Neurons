# RNA-seq Analysis Pipeline

Snakemake-based pipeline for RNA-seq data analysis with support for two quantification methods: Salmon (pseudoalignment) and STAR+HTSeq (alignment-based).

## Features

- Two quantification methods: Salmon or STAR+HTSeq
- Automatic timestamped output directories
- Differential gene expression analysis using DESeq2
- Quality control with fastp
- PCA visualization

## Requirements

- Snakemake
- Conda (for environment management)
- Docker (for DGE analysis)
- Python 3.6+

## Installation

1. Clone the repository
2. Install Snakemake and Conda
3. Prepare Docker image: `rna_seq_v1.0:latest` (must contain R with DESeq2 and required libraries)

## Configuration

Edit `config.yaml` and fill in all required parameters. All paths must be absolute.

Key parameters:
- `base_work_dir`: Base directory for timestamped run directories
- `sra_list`: Path to file with SRA accession IDs (TSV with `accession` column)
- `quantification`: Method (`salmon` or `star_htseq`)
- `reference_star.gtf`: Path to GTF annotation file
- `dge.metadata`: Path to experimental metadata file
- `docker.image`: Docker image name

You can also update configuration via command line:

```bash
snakemake configure --config base_work_dir=/path/to/work reference_star.gtf=/path/to/gtf.gtf
```

Or using Python script:

```bash
python3 scripts/configure_pipeline.py --base-work-dir /path/to/work --gtf /path/to/gtf.gtf
```

## Input Files

**SRA list file** (`sra_list.txt`):
```tsv
accession
SRR123456
SRR123457
```

**Metadata file** (`meta_data.txt`):
```tsv
sampleid    individual    diff_time
SRR123456   1             NPC
SRR123457   1             2 week neuron
```

Required columns:
- `sampleid`: Must match SRA accession
- `individual`: Replicate/donor factor
- `diff_time`: Experimental condition

## Usage

```bash
snakemake \
    --snakefile snakefile_all_pipeline_with_dge.py \
    --configfile config.yaml \
    --cores 10 \
    --use-conda
```

## Pipeline Workflow

1. Download FASTQ files from SRA (Kingfisher)
2. Quality control and filtering (fastp)
3. Quantification:
   - **Salmon**: Pseudoalignment and transcript quantification
   - **STAR+HTSeq**: Genome alignment and gene counting
4. Differential expression analysis (DESeq2 in Docker container)

## Output Structure

Each run creates a timestamped directory (`HH_MM_SS_DD_MM_YYYY`):

```
base_work_dir/
└── HH_MM_SS_DD_MM_YYYY/
    ├── reads/
    │   ├── raw/              # Downloaded FASTQ files
    │   └── raw_filtered/     # Quality-filtered reads
    ├── salmon/               # Salmon results (if method=salmon)
    │   └── {sample}/
    │       └── quant.sf
    ├── star/                 # STAR alignment (if method=star_htseq)
    │   └── {sample}_Aligned.out.bam
    ├── htseq/                # HTSeq counts (if method=star_htseq)
    │   └── {sample}.counts
    ├── reports/
    │   └── fastp/           # Quality control reports
    └── dif_expression_results/
        ├── upregulation_{method}.txt
        ├── downregulation_{method}.txt
        ├── PCA_plot_{method}.png
        └── PCA_variance_{method}.png
```

## Docker Integration

DGE analysis runs in Docker container. R scripts are located in `R-scripts/` directory and are mounted into the container at runtime.

The container must:
- Have R installed with DESeq2 and required libraries
- Be accessible via `docker.image` name from config

## Project Structure

```
.
├── snakefile_all_pipeline_with_dge.py  # Main pipeline
├── config.yaml                         # Configuration file
├── R-scripts/                          # R scripts for DGE analysis
│   ├── rna-seq-analisis_salmon.R
│   └── rna-seq-analisis_htseq.R
├── scripts/
│   └── configure_pipeline.py           # Configuration helper
└── README.md
```

## Troubleshooting

- Ensure all paths in `config.yaml` are absolute
- Verify Docker image exists and is accessible
- Check that Conda environment YAML files exist at specified paths
- Ensure SRA list and metadata files are correctly formatted

## Authors

- Lev
- Katya
- Ekaterina
