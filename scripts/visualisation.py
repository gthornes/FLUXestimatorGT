#!/usr/bin/env python3
"""
visualisation script for flux analysis results.

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


def load_module_annotations(annotations_path=None):
    """
    Load module pathway annotations.
    
    Parameters
    ----------
    annotations_path : str or Path, optional
        Path to annotations CSV file
        
    Returns
    -------
    pd.DataFrame or None
        Module annotations with pathway information
    """
    if annotations_path is None:
        annotations_path = Path('results/tables/module_annotations.csv')
    else:
        annotations_path = Path(annotations_path)
    
    if not annotations_path.exists():
        print(f"Warning: Annotations file not found at {annotations_path}")
        return None
    
    annotations = pd.read_csv(annotations_path)
    annotations.set_index('module_id', inplace=True)
    print(f"Loaded annotations for {len(annotations)} modules")
    return annotations

def setup_plotting(config):
    """Configure matplotlib settings."""
    vis_config = config['visualisation']
    
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


def plot_flux_heatmap(flux_df, output_path, annotations=None, top_n=30):
    """
    Create heatmap of top flux values across cell types.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results with columns: reaction_id, cell_type, flux
    output_path : Path
        Output file path
    annotations : pd.DataFrame, optional
        Module annotations
    top_n : int
        Number of top reactions to show
    """
    print(f"Creating flux heatmap...")
    
    # Aggregate by cell type and reaction (take mean flux across cells)
    flux_agg = flux_df.groupby(['reaction_id', 'cell_type'])['flux'].mean().reset_index()
    
    # Pivot data for heatmap
    flux_pivot = flux_agg.pivot(
        index='reaction_id', 
        columns='cell_type', 
        values='flux'
    ).fillna(0)
    
    # Select top N reactions by variance
    reaction_var = flux_pivot.var(axis=1)
    top_reactions = reaction_var.nlargest(top_n).index
    flux_pivot_top = flux_pivot.loc[top_reactions]
    
    # Add pathway labels if annotations available
    if annotations is not None:
        y_labels = []
        for module_id in flux_pivot_top.index:
            if module_id in annotations.index:
                pathway = annotations.loc[module_id, 'pathway']
                desc = annotations.loc[module_id, 'description']
                y_labels.append(desc)
            else:
                y_labels.append(module_id)
    else:
        y_labels = flux_pivot_top.index
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(14, 12))
    sns.heatmap(
        flux_pivot_top, 
        cmap='RdBu_r', 
        center=0,
        cbar_kws={'label': 'Flux (mmol/gDW/h)'},
        yticklabels=y_labels,
        ax=ax
    )
    ax.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax.set_ylabel('Metabolic Module', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {top_n} Variable Metabolic Modules Across Cell Types', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.tick_params(axis='y', labelsize=7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"Saved heatmap to {output_path}")


def plot_flux_distribution(flux_df, output_path):
    """
    Create distribution plot showing mean flux per module across cell types.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results
    output_path : Path
        Output file path
    """
    print(f"Creating flux distribution plot...")
    
    # Calculate mean absolute flux per module per cell for aggregation
    flux_summary = flux_df.groupby(['cell_type', 'reaction_id'])['flux'].apply(
        lambda x: np.abs(x).mean()
    ).reset_index(name='mean_abs_flux')
    
    # Filter to keep only meaningful fluxes (above 75th percentile)
    threshold = flux_summary['mean_abs_flux'].quantile(0.75)
    flux_filtered = flux_summary[flux_summary['mean_abs_flux'] > threshold]
    
    # Sort cell types by median of high-flux reactions
    cell_type_order = flux_filtered.groupby('cell_type')['mean_abs_flux'].median().sort_values(ascending=False).index
    
    # Calculate median flux for each cell type for color mapping
    cell_type_medians = flux_filtered.groupby('cell_type')['mean_abs_flux'].median()
    median_values = [cell_type_medians[ct] for ct in cell_type_order]
    
    # Normalize median values to [0, 1] for colourmap
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=min(median_values), vmax=max(median_values))
    cmap = plt.cm.Reds  # Red gradient colourmap (white to red)
    
    fig, ax = plt.subplots(figsize=(18, 8))
    
    # Create violin plot to show distribution shape
    parts = ax.violinplot(
        [flux_filtered[flux_filtered['cell_type'] == ct]['mean_abs_flux'].values 
         for ct in cell_type_order],
        positions=range(len(cell_type_order)),
        widths=0.7,
        showmeans=True,
        showmedians=True
    )
    
    # Colour the violins based on median flux (red gradient)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(cmap(norm(median_values[i])))
        pc.set_alpha(0.8)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
    
    # Style the statistical lines
    parts['cmeans'].set_color('red')
    parts['cmeans'].set_linewidth(2)
    parts['cmedians'].set_color('black')
    parts['cmedians'].set_linewidth(2),
    parts['cbars'].set_color('black'),
    parts['cmins'].set_color('black'),
    parts['cmaxes'].set_color('black')
    
    ax.set_xticks(range(len(cell_type_order)))
    ax.set_xticklabels(cell_type_order, rotation=45, ha='right')
    ax.set_ylabel('Mean Module Flux (mmol/gDW/h)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of High-Activity Metabolic Modules Across Cell Types\n(Top 25% most active modules)', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.grid(axis='y', alpha=0.3)
    
    # Add colorbar legend showing the flux gradient
    from matplotlib.cm import ScalarMappable
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label('Median Metabolic Flux\n(mmol/gDW/h)', fontsize=11, fontweight='bold')
    
    # Add legend for statistical lines
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='red', linewidth=2, label='Mean'),
        Line2D([0], [0], color='black', linewidth=2, label='Median')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"Saved distribution plot to {output_path}")


def plot_pathway_comparison(flux_df, output_path, annotations=None, pathway_reactions=None):
    """
    Create bar plot comparing pathway activity across cell types.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results
    output_path : Path
        Output file path
    annotations : pd.DataFrame, optional
        Module annotations
    pathway_reactions : dict
        Mapping of pathway names to reaction IDs
    """
    print(f"Creating pathway comparison plot...")
    
    if pathway_reactions is None:
        # Default: show top reactions for top cell types
        print("No pathway mapping provided, showing top reactions for most active cell types")
        
        # Calculate mean absolute flux per cell type and reaction
        flux_agg = flux_df.groupby(['cell_type', 'reaction_id'])['flux'].apply(
            lambda x: np.abs(x).mean()
        ).reset_index(name='mean_abs_flux')
        
        # Get top 5 cell types by overall metabolic activity
        cell_type_activity = flux_agg.groupby('cell_type')['mean_abs_flux'].mean().nlargest(5)
        top_cell_types = cell_type_activity.index.tolist()
        
        # Filter to top cell types
        flux_top_cells = flux_agg[flux_agg['cell_type'].isin(top_cell_types)]
        
        # Get top 10 reactions by variance across these cell types
        flux_pivot_temp = flux_top_cells.pivot(
            index='reaction_id',
            columns='cell_type',
            values='mean_abs_flux'
        ).fillna(0)
        
        reaction_var = flux_pivot_temp.var(axis=1).nlargest(10)
        top_reactions = reaction_var.index.tolist()
        
        # Filter to top reactions
        flux_final = flux_top_cells[flux_top_cells['reaction_id'].isin(top_reactions)]
        
        # Create final pivot for plotting
        flux_pivot = flux_final.pivot(
            index='reaction_id',
            columns='cell_type',
            values='mean_abs_flux'
        ).fillna(0)
        
        # Reorder by total flux
        flux_pivot = flux_pivot.loc[flux_pivot.sum(axis=1).sort_values(ascending=False).index]
        
        # Add pathway labels if annotations available
        if annotations is not None:
            y_labels = []
            for module_id in flux_pivot.index:
                if module_id in annotations.index:
                    pathway = annotations.loc[module_id, 'pathway']
                    desc = annotations.loc[module_id, 'description']
                    y_labels.append(desc)
                else:
                    y_labels.append(module_id)
        else:
            y_labels = flux_pivot.index
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        flux_pivot.plot(kind='barh', ax=ax, width=0.8)
        ax.set_yticklabels(y_labels, fontsize=9)
        ax.set_xlabel('Mean Absolute Flux (mmol/gDW/h)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Metabolic Module', fontsize=12, fontweight='bold')
        ax.set_title(f'Top 10 Variable Modules Across {len(top_cell_types)} Most Active Cell Types',
                     fontsize=13, fontweight='bold')
        ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
        ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"Saved pathway comparison to {output_path}")


def plot_metabolic_summary(flux_df, output_dir):
    """
    Create three separate summary statistics plots.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results
    output_dir : Path
        Output directory path
    """
    print(f"Creating metabolic summary plots...")
    
    # Calculate summary statistics per cell type
    summary_stats = flux_df.groupby('cell_type')['flux'].agg([
        ('Mean Flux', lambda x: np.abs(x).mean()),
        ('Active Reactions', lambda x: (np.abs(x) > 1e-6).sum()),
        ('Max Flux', lambda x: np.abs(x).max())
    ]).reset_index()
    
    # Sort by mean flux for consistent ordering
    summary_stats = summary_stats.sort_values('Mean Flux', ascending=False)
    
    # Get file extension from output_dir parent
    fmt = 'png'  # default
    
    # Plot 1: Mean absolute flux
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(summary_stats['cell_type'], summary_stats['Mean Flux'])
    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel('Mean Absolute Flux', fontsize=12)
    ax.set_title('Average Metabolic Activity', fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f'metabolic_summary_mean.{fmt}', dpi=300)
    plt.close()
    
    # Plot 2: Number of active reactions
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(summary_stats['cell_type'], summary_stats['Active Reactions'])
    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel('Number of Active Reactions', fontsize=12)
    ax.set_title('Metabolic Pathway Breadth', fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f'metabolic_summary_breadth.{fmt}', dpi=300)
    plt.close()
    
    # Plot 3: Maximum flux
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(summary_stats['cell_type'], summary_stats['Max Flux'])
    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel('Maximum Absolute Flux', fontsize=12)
    ax.set_title('Peak Metabolic Activity', fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f'metabolic_summary_peak.{fmt}', dpi=300)
    plt.close()
    
    print(f"Saved summary plots to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='Visualise flux analysis results'
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
        config['visualisation']['figure_format'] = args.format
    
    # Setup plotting
    setup_plotting(config)
    
    # Load module annotations
    annotations = load_module_annotations()
    
    # Load flux results
    print(f"Loading flux results from {args.input}")
    flux_df = pd.read_csv(args.input, index_col=0)
    
    # Load cell metadata
    metadata_path = Path(args.input).parent / 'cell_metadata.csv'
    if not metadata_path.exists():
        print(f"Error: Cell metadata file not found at {metadata_path}")
        print("Looking for cell_metadata.csv in the same directory as flux results")
        return 1
    
    cell_metadata = pd.read_csv(metadata_path, index_col=0)
    
    # Merge flux results with cell types
    flux_df = flux_df.join(cell_metadata['cell_type'], how='left')
    
    # Filter out 'Ignore' cell type
    flux_df = flux_df[flux_df['cell_type'] != 'Ignore']
    
    # Reshape from wide format (modules as columns) to long format (reaction_id, cell_type, flux)
    # Get module columns (M_1, M_2, etc.)
    module_cols = [col for col in flux_df.columns if col.startswith('M_')]
    
    # Melt the dataframe to long format
    flux_df = flux_df[module_cols + ['cell_type']].melt(
        id_vars=['cell_type'],
        value_vars=module_cols,
        var_name='reaction_id',
        value_name='flux'
    )
    
    print(f"Loaded {len(flux_df)} flux predictions for {flux_df['cell_type'].nunique()} cell types")
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    fmt = config['visualisation']['figure_format']
    
    plot_flux_heatmap(
        flux_df, 
        output_dir / f'flux_heatmap.{fmt}',
        annotations=annotations
    )
    
    plot_flux_distribution(
        flux_df,
        output_dir / f'flux_distribution.{fmt}'
    )
    
    plot_pathway_comparison(
        flux_df,
        output_dir / f'pathway_comparison.{fmt}',
        annotations=annotations
    )
    
    plot_metabolic_summary(
        flux_df,
        output_dir
    )
    
    print("\nvisualisation complete!")
    print(f"Figures saved to {output_dir}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
