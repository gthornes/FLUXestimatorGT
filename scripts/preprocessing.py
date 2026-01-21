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
import scrublet as scr
import matplotlib.pyplot as plt


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
    
    # Identify mitochondrial genes - check multiple prefixes
    adata.var['mt'] = adata.var_names.str.startswith('mt-') | adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('Mt-')
    
    # Identify ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl', 'RPS', 'RPL'))
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    n_mt_genes = adata.var['mt'].sum()
    print(f"Found {n_mt_genes} mitochondrial genes")
    print(f"Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.0f}")
    print(f"Median genes per cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"Mean counts per cell: {adata.obs['total_counts'].mean():.0f}")
    print(f"Mean mitochondrial %: {adata.obs['pct_counts_mt'].mean():.2f}%")


def detect_doublets(adata, config, save_plot=None):
    """
    Detect doublets using Scrublet.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    config : dict
        Configuration parameters
    save_plot : str, optional
        Path to save doublet detection plot
        
    Returns
    -------
    adata : AnnData
        Data with doublet annotations added
    """
    print("\nDetecting doublets with Scrublet...")
    
    # Get expected doublet rate from config
    expected_doublet_rate = config['quality_control'].get('doublet_rate', 0.06)
    
    # Run Scrublet
    scrub = scr.Scrublet(
        adata.X, 
        expected_doublet_rate=expected_doublet_rate
    )
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=50
    )
    
    # Add results to adata
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets
    
    # Print detailed diagnostics
    print("=" * 60)
    print("DOUBLET DETECTION SUMMARY")
    print("=" * 60)
    print(f"Total cells analyzed: {len(doublet_scores)}")
    print(f"Expected doublet rate: {expected_doublet_rate*100:.1f}%")
    print(f"Expected doublets: ~{int(len(doublet_scores) * expected_doublet_rate)}")
    print(f"\nActual detected: {predicted_doublets.sum()} ({100 * predicted_doublets.sum() / len(predicted_doublets):.2f}%)")
    print(f"Automatic threshold: {scrub.threshold_:.3f}")
    print(f"\nDoublet score statistics:")
    print(f"  Min: {doublet_scores.min():.3f}")
    print(f"  25th percentile: {np.percentile(doublet_scores, 25):.3f}")
    print(f"  Median: {np.median(doublet_scores):.3f}")
    print(f"  75th percentile: {np.percentile(doublet_scores, 75):.3f}")
    print(f"  95th percentile: {np.percentile(doublet_scores, 95):.3f}")
    print(f"  99th percentile: {np.percentile(doublet_scores, 99):.3f}")
    print(f"  Max: {doublet_scores.max():.3f}")
    
    # Check how many cells would be flagged at different thresholds
    for threshold in [0.2, 0.3, 0.4, 0.5]:
        n_doublets = (doublet_scores > threshold).sum()
        pct = 100 * n_doublets / len(doublet_scores)
        print(f"\nAt threshold {threshold:.2f}: {n_doublets} cells ({pct:.2f}%)")
    
    print("\n" + "=" * 60)
    print("INTERPRETATION:")
    if predicted_doublets.sum() / len(predicted_doublets) < 0.01:
        print("⚠️  Very low detection suggests homotypic doublets dominate.")
        print("   These are two similar cells that look like singlets.")
        print("   Consider proceeding without aggressive doublet filtering.")
    else:
        print("✓ Detection rate looks reasonable.")
    print("=" * 60)
    
    # Optional: save plot
    if save_plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        
        # Histogram of doublet scores
        axes[0].hist(doublet_scores, bins=50, edgecolor='black')
        axes[0].set_xlabel('Doublet score')
        axes[0].set_ylabel('Number of cells')
        axes[0].set_title('Doublet score distribution')
        axes[0].axvline(scrub.threshold_, color='r', linestyle='--', 
                       label=f'Threshold ({scrub.threshold_:.3f})')
        axes[0].legend()
        
        # Doublet scores vs gene count
        scatter = axes[1].scatter(
            adata.obs['n_genes_by_counts'], 
            doublet_scores, 
            c=predicted_doublets, 
            cmap='coolwarm', 
            s=5, 
            alpha=0.5
        )
        axes[1].set_xlabel('Number of genes')
        axes[1].set_ylabel('Doublet score')
        axes[1].set_title('Genes vs Doublet Score')
        plt.colorbar(scatter, ax=axes[1], label='Predicted doublet')
        
        plt.tight_layout()
        plt.savefig(save_plot, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\nDoublet plot saved to {save_plot}")
    
    return adata


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
    adata = adata[adata.obs['pct_counts_mt'] < qc_params['max_mito_percent'], :]
    
    # Filter doublets if detected
    if 'predicted_doublet' in adata.obs.columns:
        n_doublets = adata.obs['predicted_doublet'].sum()
        adata = adata[~adata.obs['predicted_doublet'], :]
        print(f"Removed {n_doublets} predicted doublets")
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=qc_params['min_cells_per_gene'])
    
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    print(f"Cells removed: {n_cells_before - n_cells_after} ({(n_cells_before - n_cells_after) / n_cells_before * 100:.1f}%)")
    print(f"Genes removed: {n_genes_before - n_genes_after} ({(n_genes_before - n_genes_after) / n_genes_before * 100:.1f}%)")
    print(f"Remaining: {n_cells_after} cells × {n_genes_after} genes")
    
    return adata


def normalise_data(adata, config):
    """
    Normalise and log-transform data.
    
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
    
    # Normalise to target sum
    sc.pp.normalise_total(adata, target_sum=config['preprocessing']['target_sum'])
    
    # Log transform
    if norm_params['method'] == 'log1p':
        sc.pp.log1p(adata)
    
    # Store normalised data
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
    parser.add_argument(
        '--skip-doublet-detection',
        action='store_true',
        help='Skip doublet detection with Scrublet'
    )
    parser.add_argument(
        '--doublet-plot',
        type=str,
        default=None,
        help='Path to save doublet detection plot (optional)'
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
    
    # Detect doublets (before filtering)
    if not args.skip_doublet_detection:
        adata = detect_doublets(adata, config, save_plot=args.doublet_plot)
    
    # Filter cells and genes
    if not args.skip_filtering:
        adata = filter_cells_and_genes(adata, config)
        
    # Normalise data
    normalise_data(adata, config)
    
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
