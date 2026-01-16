#!/usr/bin/env python3
"""
Script for dynamic pipeline configuration setup.
Updates config.yaml based on provided parameters.
"""

import yaml
import sys
import os
import argparse


def load_config(config_path):
    """Loads existing config.yaml or creates base structure"""
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f) or {}
    else:
        # Base structure
        config = {
            'base_work_dir': None,
            'sra_list': None,
            'quantification': None,
            'salmon': {
                'transcripts': None,
                'index_dir': None,
                'kmer': None,
                'threads': None
            },
            'reference_star': {
                'genome_fasta': None,
                'gtf': None
            },
            'star': {
                'index': None,
                'threads': None,
                'sjdbOverhang': None,
                'star_bin': None
            },
            'htseq': {
                'threads': None
            },
            'fastp_threads': None,
            'kingfisher_threads': None,
            'dge': {
                'metadata': None,
                'output_dir': None
            },
            'docker': {
                'image': None
            },
            'conda_envs': {
                'kingfisher': None,
                'fastp': None,
                'salmon': None,
                'star': None,
                'htseq': None
            }
        }
    
    # Initialize nested structures if they don't exist
    if 'salmon' not in config:
        config['salmon'] = {}
    if 'reference_star' not in config:
        config['reference_star'] = {}
    if 'star' not in config:
        config['star'] = {}
    if 'htseq' not in config:
        config['htseq'] = {}
    if 'dge' not in config:
        config['dge'] = {}
    if 'docker' not in config:
        config['docker'] = {}
    if 'conda_envs' not in config:
        config['conda_envs'] = {}
    
    return config


def update_config(config, updates):
    """Updates configuration based on provided parameters"""
    for key, value in updates.items():
        if value is not None:  # Update only if value is provided
            if '.' in key:
                # Nested parameter (e.g., reference_star.gtf)
                parts = key.split('.')
                current = config
                for part in parts[:-1]:
                    if part not in current:
                        current[part] = {}
                    current = current[part]
                current[parts[-1]] = value
            else:
                # Simple parameter
                config[key] = value


def save_config(config, config_path):
    """Saves configuration to file"""
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False, 
                  allow_unicode=True, width=1000)
    print(f"✓ Configuration saved: {config_path}")


def main():
    parser = argparse.ArgumentParser(
        description='RNA-seq pipeline configuration setup'
    )
    
    # Frequently used parameters
    parser.add_argument('--base-work-dir', dest='base_work_dir',
                       help='Base working directory')
    parser.add_argument('--gtf', dest='reference_star.gtf',
                       help='Path to GTF annotation file')
    parser.add_argument('--quantification', 
                       choices=['salmon', 'star_htseq'],
                       help='Quantification method')
    parser.add_argument('--sra-list',
                       help='Path to SRA ID list file')
    parser.add_argument('--metadata', dest='dge.metadata',
                       help='Path to experimental metadata file')
    parser.add_argument('--docker-image', dest='docker.image',
                       help='Docker image with R scripts')
    
    # Additional parameters via --set
    parser.add_argument('--set', action='append', nargs=2, metavar=('KEY', 'VALUE'),
                       help='Set arbitrary parameter (can be used multiple times)')
    
    parser.add_argument('--config-file', default='config.yaml',
                       help='Path to configuration file (default: config.yaml)')
    
    args = parser.parse_args()
    
    # Load existing configuration
    config = load_config(args.config_file)
    
    # Prepare updates
    updates = {}
    
    # Update parameters from arguments
    if args.base_work_dir:
        updates['base_work_dir'] = args.base_work_dir
    if getattr(args, 'reference_star.gtf'):
        updates['reference_star.gtf'] = getattr(args, 'reference_star.gtf')
    if args.quantification:
        updates['quantification'] = args.quantification
    if args.sra_list:
        updates['sra_list'] = args.sra_list
    if getattr(args, 'dge.metadata'):
        updates['dge.metadata'] = getattr(args, 'dge.metadata')
    if getattr(args, 'docker.image'):
        updates['docker.image'] = getattr(args, 'docker.image')
    
    # Updates via --set
    if args.set:
        for key, value in args.set:
            updates[key] = value
    
    # Apply updates
    if updates:
        update_config(config, updates)
        save_config(config, args.config_file)
        print("\nUpdated parameters:")
        for key, value in updates.items():
            print(f"  {key}: {value}")
    else:
        print("No parameters specified for update.")
        print("Use --help for parameter help.")


if __name__ == '__main__':
    main()
