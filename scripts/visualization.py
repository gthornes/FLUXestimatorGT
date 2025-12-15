#!/usr/bin/env python3
"""
Visualization script for flux analysis results.

This script generates publication-quality figures for metabolic flux analysis.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import yaml


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def setup_plotting(config):
    """Configure matplotlib settings."""
    vis_config = config['visualization']
    
    plt.rcParams['figure.dpi'] = vis_config['dpi']
    plt.rcParams['savefig.dpi'] = vis_config['dpi']
    plt.rcParams['figure.figsize'] = (10, 6)
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.labelsize'] = 11
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    
    sns.set_style("whitegrid")
    sns.set_palette(vis_config['color_palette'])


def plot_flux_heatmap(flux_df, output_path, top_n=50):
    """
    Create heatmap of top flux values across cell types.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results with columns: reaction_id, cell_type, flux
    output_path : Path
        Output file path
    top_n : int
        Number of top reactions to show
    """
    print(f"Creating flux heatmap...")
    
    # Pivot data for heatmap
    flux_pivot = flux_df.pivot(
        index='reaction_id', 
        columns='cell_type', 
        values='flux'
    ).fillna(0)
    
    # Select top N reactions by variance
    reaction_var = flux_pivot.var(axis=1)
    top_reactions = reaction_var.nlargest(top_n).index
    flux_pivot_top = flux_pivot.loc[top_reactions]
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(
        flux_pivot_top, 
        cmap='RdBu_r', 
        center=0,
        cbar_kws={'label': 'Flux (mmol/gDW/h)'},
        ax=ax
    )
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Reaction')
    ax.set_title(f'Top {top_n} Variable Reactions Across Cell Types')
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    
    print(f"Saved heatmap to {output_path}")


def plot_flux_distribution(flux_df, output_path):
    """
    Create violin plot of flux distributions per cell type.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results
    output_path : Path
        Output file path
    """
    print(f"Creating flux distribution plot...")
    
    # Filter for non-zero fluxes
    flux_df_nonzero = flux_df[flux_df['flux'].abs() > 1e-6].copy()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.violinplot(
        data=flux_df_nonzero,
        x='cell_type',
        y='flux',
        ax=ax
    )
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Flux (mmol/gDW/h)')
    ax.set_title('Distribution of Non-Zero Fluxes Across Cell Types')
    ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    
    print(f"Saved distribution plot to {output_path}")


def plot_pathway_comparison(flux_df, output_path, pathway_reactions=None):
    """
    Create bar plot comparing pathway activity across cell types.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results
    output_path : Path
        Output file path
    pathway_reactions : dict
        Mapping of pathway names to reaction IDs
    """
    print(f"Creating pathway comparison plot...")
    
    if pathway_reactions is None:
        # Default: show top reactions
        print("No pathway mapping provided, showing top reactions by absolute flux")
        
        # Calculate mean absolute flux per reaction
        reaction_mean = flux_df.groupby('reaction_id')['flux'].apply(
            lambda x: np.abs(x).mean()
        ).nlargest(20)
        
        # Filter to top reactions
        flux_df_top = flux_df[flux_df['reaction_id'].isin(reaction_mean.index)]
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Create grouped bar plot
        flux_pivot = flux_df_top.pivot(
            index='reaction_id',
            columns='cell_type',
            values='flux'
        )
        
        flux_pivot.plot(kind='bar', ax=ax)
        ax.set_xlabel('Reaction')
        ax.set_ylabel('Flux (mmol/gDW/h)')
        ax.set_title('Top 20 Reactions by Mean Absolute Flux')
        ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.tick_params(axis='x', rotation=45, labelsize=8)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    
    print(f"Saved pathway comparison to {output_path}")


def plot_metabolic_summary(flux_df, output_path):
    """
    Create summary statistics plot.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results
    output_path : Path
        Output file path
    """
    print(f"Creating metabolic summary plot...")
    
    # Calculate summary statistics per cell type
    summary_stats = flux_df.groupby('cell_type')['flux'].agg([
        ('Mean Flux', lambda x: np.abs(x).mean()),
        ('Active Reactions', lambda x: (np.abs(x) > 1e-6).sum()),
        ('Max Flux', lambda x: np.abs(x).max())
    ]).reset_index()
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot 1: Mean absolute flux
    axes[0].bar(summary_stats['cell_type'], summary_stats['Mean Flux'])
    axes[0].set_xlabel('Cell Type')
    axes[0].set_ylabel('Mean Absolute Flux')
    axes[0].set_title('Average Metabolic Activity')
    axes[0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Number of active reactions
    axes[1].bar(summary_stats['cell_type'], summary_stats['Active Reactions'])
    axes[1].set_xlabel('Cell Type')
    axes[1].set_ylabel('Number of Active Reactions')
    axes[1].set_title('Metabolic Pathway Breadth')
    axes[1].tick_params(axis='x', rotation=45)
    
    # Plot 3: Maximum flux
    axes[2].bar(summary_stats['cell_type'], summary_stats['Max Flux'])
    axes[2].set_xlabel('Cell Type')
    axes[2].set_ylabel('Maximum Absolute Flux')
    axes[2].set_title('Peak Metabolic Activity')
    axes[2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    
    print(f"Saved summary plot to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize flux analysis results'
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input flux results CSV file'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory for figures'
    )
    parser.add_argument(
        '--config', '-c',
        default='config/analysis_config.yaml',
        help='Configuration file (default: config/analysis_config.yaml)'
    )
    parser.add_argument(
        '--format',
        default='png',
        choices=['png', 'pdf', 'svg'],
        help='Output figure format (default: png)'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    print(f"Loading configuration from {args.config}")
    config = load_config(args.config)
    
    # Override format if specified
    if args.format:
        config['visualization']['figure_format'] = args.format
    
    # Setup plotting
    setup_plotting(config)
    
    # Load flux results
    print(f"Loading flux results from {args.input}")
    flux_df = pd.read_csv(args.input)
    print(f"Loaded {len(flux_df)} flux predictions for {flux_df['cell_type'].nunique()} cell types")
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    fmt = config['visualization']['figure_format']
    
    plot_flux_heatmap(
        flux_df, 
        output_dir / f'flux_heatmap.{fmt}'
    )
    
    plot_flux_distribution(
        flux_df,
        output_dir / f'flux_distribution.{fmt}'
    )
    
    plot_pathway_comparison(
        flux_df,
        output_dir / f'pathway_comparison.{fmt}'
    )
    
    plot_metabolic_summary(
        flux_df,
        output_dir / f'metabolic_summary.{fmt}'
    )
    
    print("\nVisualization complete!")
    print(f"Figures saved to {output_dir}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
