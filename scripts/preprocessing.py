#!/usr/bin/env python3
"""
Preprocessing script for single-cell RNA-seq data.

This script performs quality control, filtering, normalization, and
dimensionality reduction on raw scRNA-seq data.
"""

import argparse
import sys
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd
import yaml


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def load_data(input_path):
    """
    Load single-cell data from various formats.
    
    Parameters
    ----------
    input_path : str or Path
        Path to input data file or directory
        
    Returns
    -------
    adata : AnnData
        Annotated data matrix
    """
    input_path = Path(input_path)
    
    if input_path.suffix == '.h5ad':
        print(f"Loading h5ad file: {input_path}")
        adata = sc.read_h5ad(input_path)
    elif input_path.suffix == '.loom':
        print(f"Loading loom file: {input_path}")
        adata = sc.read_loom(input_path)
    elif input_path.is_dir():
        print(f"Loading 10X directory: {input_path}")
        adata = sc.read_10x_mtx(input_path, var_names='gene_symbols')
    else:
        raise ValueError(f"Unsupported file format: {input_path}")
    
    print(f"Loaded data: {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata


def calculate_qc_metrics(adata):
    """
    Calculate quality control metrics.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    """
    print("\nCalculating QC metrics...")
    
    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('mt-') | adata.var_names.str.startswith('MT-')
    
    # Identify ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl', 'RPS', 'RPL'))
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    print(f"Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.0f}")
    print(f"Median genes per cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"Mean counts per cell: {adata.obs['total_counts'].mean():.0f}")
    print(f"Mean mitochondrial %: {adata.obs['pct_counts_mt'].mean():.2f}%")


def filter_cells_and_genes(adata, config):
    """
    Filter cells and genes based on QC metrics.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    config : dict
        Configuration parameters
        
    Returns
    -------
    adata : AnnData
        Filtered data
    """
    print("\nFiltering cells and genes...")
    
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Get parameters from config
    qc_params = config['preprocessing']
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=qc_params['min_genes_per_cell'])
    sc.pp.filter_cells(adata, max_genes=qc_params['max_genes_per_cell'])
    
    # Filter based on mitochondrial content
    adata = adata[adata.obs['pct_counts_mt'] < qc_params['max_mito_percent'], :].copy()
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=qc_params['min_cells_per_gene'])
    
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    print(f"Cells removed: {n_cells_before - n_cells_after} ({(n_cells_before - n_cells_after) / n_cells_before * 100:.1f}%)")
    print(f"Genes removed: {n_genes_before - n_genes_after} ({(n_genes_before - n_genes_after) / n_genes_before * 100:.1f}%)")
    print(f"Remaining: {n_cells_after} cells × {n_genes_after} genes")
    
    return adata


def normalize_data(adata, config):
    """
    Normalize and log-transform data.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    config : dict
        Configuration parameters
    """
    print("\nNormalizing data...")
    
    norm_params = config['normalization']
    
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()
    
    # Normalize to target sum
    sc.pp.normalize_total(adata, target_sum=config['preprocessing']['target_sum'])
    
    # Log transform
    if norm_params['method'] == 'log1p':
        sc.pp.log1p(adata)
    
    # Store normalized data
    adata.layers['log1p_norm'] = adata.X.copy()
    
    print("Normalization complete")


def identify_highly_variable_genes(adata, config):
    """
    Identify highly variable genes.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    config : dict
        Configuration parameters
    """
    print("\nIdentifying highly variable genes...")
    
    n_top_genes = config['normalization']['highly_variable_genes']
    
    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=n_top_genes,
        subset=False,
        flavor='seurat'
    )
    
    n_hvg = adata.var['highly_variable'].sum()
    print(f"Identified {n_hvg} highly variable genes")


def scale_data(adata, config):
    """
    Scale data to unit variance and zero mean.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    config : dict
        Configuration parameters
    """
    print("\nScaling data...")
    
    max_value = config['normalization']['scale_max']
    
    sc.pp.scale(adata, max_value=max_value)
    
    print(f"Data scaled (clipped to max value: {max_value})")


def perform_pca(adata, config):
    """
    Perform PCA for dimensionality reduction.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    config : dict
        Configuration parameters
    """
    print("\nPerforming PCA...")
    
    n_pcs = config['cell_annotation']['n_pcs']
    
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
    
    print(f"PCA complete ({n_pcs} components)")


def main():
    parser = argparse.ArgumentParser(
        description='Preprocess single-cell RNA-seq data'
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input data file (.h5ad, .loom) or 10X directory'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output file path (.h5ad)'
    )
    parser.add_argument(
        '--config', '-c',
        default='config/analysis_config.yaml',
        help='Configuration file (default: config/analysis_config.yaml)'
    )
    parser.add_argument(
        '--skip-filtering',
        action='store_true',
        help='Skip cell and gene filtering'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    print(f"Loading configuration from {args.config}")
    config = load_config(args.config)
    
    # Set random seed for reproducibility
    np.random.seed(config['compute']['random_state'])
    
    # Configure scanpy
    sc.settings.verbosity = 1
    sc.settings.n_jobs = config['compute']['n_jobs']
    
    # Load data
    adata = load_data(args.input)
    
    # Calculate QC metrics
    calculate_qc_metrics(adata)
    
    # Filter cells and genes
    if not args.skip_filtering:
        adata = filter_cells_and_genes(adata, config)
        
    # Normalize data
    normalize_data(adata, config)
    
    # Identify highly variable genes
    identify_highly_variable_genes(adata, config)
    
    # Scale data (only on HVGs)
    scale_data(adata, config)
    
    # Perform PCA
    perform_pca(adata, config)
    
    # Save processed data
    print(f"\nSaving processed data to {args.output}")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(args.output, compression='gzip')
    
    print("\nPreprocessing complete!")
    print(f"Final dataset: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
