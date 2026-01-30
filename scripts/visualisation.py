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
    
    # Calculate mean flux for each cell type for color mapping
    cell_type_means = flux_filtered.groupby('cell_type')['mean_abs_flux'].mean()
    mean_values = [cell_type_means[ct] for ct in cell_type_order]
    
    # Normalize mean values to [0, 1] for colourmap
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=min(mean_values), vmax=max(mean_values))
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
    
    # Colour the violins based on mean flux (red gradient)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(cmap(norm(mean_values[i])))
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
    ax.set_ylabel('Module Flux (mmol/gDW/h)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Cell Type (ordered by median flux)', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of High-Activity Metabolic Modules Across Cell Types\n(Top 25% most active modules)', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.grid(axis='y', alpha=0.3)
    
    # Add colorbar legend showing the flux gradient
    from matplotlib.cm import ScalarMappable
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label('Mean Metabolic Flux\n(mmol/gDW/h)', fontsize=11, fontweight='bold')
    
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


def plot_flux_distribution_by_stage(flux_df, output_path):
    """
    Create distribution plot showing mean flux per module across cell types, separated by stage.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results with 'stage' column
    output_path : Path
        Output file path
    """
    print(f"Creating flux distribution plot by stage...")
    
    # Check if stage column exists
    if 'stage' not in flux_df.columns:
        print("Warning: No 'stage' column found in flux data. Skipping stage-specific distribution plot.")
        return
    
    # Calculate mean absolute flux per module per cell for aggregation
    flux_summary = flux_df.groupby(['cell_type', 'stage', 'reaction_id'])['flux'].apply(
        lambda x: np.abs(x).mean()
    ).reset_index(name='mean_abs_flux')
    
    # Filter to keep only meaningful fluxes (above 75th percentile)
    threshold = flux_summary['mean_abs_flux'].quantile(0.75)
    flux_filtered = flux_summary[flux_summary['mean_abs_flux'] > threshold]
    
    # Get unique stages and ensure estrus is first
    stages = sorted(flux_filtered['stage'].unique())
    if 'estrus' in stages and 'diestrus' in stages:
        stages = ['estrus', 'diestrus']
    
    # Sort cell types by overall median flux
    cell_type_order = flux_filtered.groupby('cell_type')['mean_abs_flux'].median().sort_values(ascending=False).index
    
    # Create single figure with paired violins
    fig, ax = plt.subplots(figsize=(20, 8))
    
    # Color scheme: estrus = blue, diestrus = red
    stage_colors = {'estrus': 'cyan', 'diestrus': 'royalblue'}
    
    # Calculate positions for paired violins
    # Each cell type gets 2 positions (estrus and diestrus) with a gap between cell types
    spacing = 3  # Space between cell type pairs
    violin_width = 1.0
    
    all_parts = []
    all_positions = []
    all_colors = []
    
    for i, cell_type in enumerate(cell_type_order):
        base_pos = i * spacing
        
        for j, stage in enumerate(stages):
            stage_data = flux_filtered[(flux_filtered['stage'] == stage) & 
                                       (flux_filtered['cell_type'] == cell_type)]
            
            if len(stage_data) > 0:
                pos = base_pos + j * violin_width
                all_positions.append(pos)
                
                # Create individual violin plot
                parts = ax.violinplot(
                    [stage_data['mean_abs_flux'].values],
                    positions=[pos],
                    widths=violin_width * 0.9,
                    showmeans=True,
                    showmedians=True
                )
                
                # Color the violin
                for pc in parts['bodies']:
                    pc.set_facecolor(stage_colors[stage])
                    pc.set_alpha(0.7)
                    pc.set_edgecolor('black')
                    pc.set_linewidth(1)
                
                # Style the statistical lines
                parts['cmeans'].set_color('darkred')
                parts['cmeans'].set_linewidth(1.5)
                parts['cmedians'].set_color('black')
                parts['cmedians'].set_linewidth(1.5)
                parts['cbars'].set_color('black')
                parts['cmins'].set_color('black')
                parts['cmaxes'].set_color('black')
    
    # Set x-axis ticks at the center of each cell type pair
    cell_type_centers = [i * spacing + violin_width / 2 for i in range(len(cell_type_order))]
    ax.set_xticks(cell_type_centers)
    ax.set_xticklabels(cell_type_order, rotation=45, ha='right', fontsize=9)
    
    ax.set_ylabel('Mean Module Flux (mmol/gDW/h)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of High-Activity Metabolic Modules by Stage\n(Top 25% most active modules)', 
                 fontsize=16, fontweight='bold', pad=20)
    ax.grid(axis='y', alpha=0.3)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='cyan', alpha=0.7, edgecolor='black', label='Estrus'),
        Patch(facecolor='royalblue', alpha=0.7, edgecolor='black', label='Diestrus'),
        plt.Line2D([0], [0], color='darkred', linewidth=1.5, label='Mean'),
        plt.Line2D([0], [0], color='black', linewidth=1.5, label='Median')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.95)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved stage-specific distribution plot to {output_path}")


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


def plot_flux_by_category_and_stage(flux_df, output_path, annotations=None):
    """
    Create bar plot comparing metabolic category activity between cycle stages.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results with 'reaction_id', 'flux', 'cell_type', and 'stage' columns
    output_path : Path
        Output file path
    annotations : pd.DataFrame, optional
        Module annotations with 'pathway' column
    """
    print(f"\nCreating metabolic category comparison by stage...")
    
    # Check if stage column exists
    if 'stage' not in flux_df.columns:
        print("Warning: No 'stage' column found in flux data. Skipping category-by-stage plot.")
        return
    
    # Check if annotations available
    if annotations is None:
        print("Warning: No module annotations available. Cannot group by metabolic category.")
        return
    
    # Extract module ID from reaction_id (in case it's annotated with descriptions)
    import re
    def extract_module_id(name):
        match = re.search(r'(M_\d+)', str(name))
        return match.group(1) if match else name
    
    flux_df = flux_df.copy()
    flux_df['module_id'] = flux_df['reaction_id'].apply(extract_module_id)
    
    # Merge with annotations to get pathway category
    flux_with_category = flux_df.merge(
        annotations[['pathway']], 
        left_on='module_id', 
        right_index=True, 
        how='left'
    )
    
    # Remove rows without pathway annotation
    flux_with_category = flux_with_category.dropna(subset=['pathway'])
    
    if len(flux_with_category) == 0:
        print("Error: No modules matched with pathway annotations")
        return
    
    # Calculate mean absolute flux per category per stage
    category_flux = flux_with_category.groupby(['pathway', 'stage'])['flux'].apply(
        lambda x: np.abs(x).mean()
    ).reset_index(name='mean_flux')
    
    # Get unique stages
    stages = sorted(category_flux['stage'].unique())
    if len(stages) < 2:
        print(f"Warning: Only {len(stages)} stage(s) found. Need at least 2 for comparison.")
        return
    
    # Pivot for easier plotting
    category_pivot = category_flux.pivot(
        index='pathway',
        columns='stage',
        values='mean_flux'
    ).fillna(0)
    
    # Sort by total activity across both stages
    category_pivot['total'] = category_pivot.sum(axis=1)
    category_pivot = category_pivot.sort_values('total', ascending=False).drop('total', axis=1)
    
    # Create grouped bar plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = np.arange(len(category_pivot.index))
    width = 0.35
    
    # Color scheme matching stage distribution plot
    colors = {'estrus': 'steelblue', 'diestrus': 'coral'}
    
    # Plot bars for each stage
    for i, stage in enumerate(stages):
        offset = width * (i - 0.5)
        color = colors.get(stage, f'C{i}')
        ax.bar(x + offset, category_pivot[stage], width, 
               label=stage.capitalize(), color=color, alpha=0.8, edgecolor='black')
    
    ax.set_xlabel('Metabolic Category', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean Metabolic Flux (mmol/gDW/h)', fontsize=12, fontweight='bold')
    ax.set_title('Metabolic Activity by Category Across Cycle Stages', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(category_pivot.index, rotation=45, ha='right', fontsize=10)
    ax.legend(title='Cycle Stage', fontsize=11, title_fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print summary statistics
    print("\nMetabolic Category Activity Summary:")
    print(category_pivot.to_string())
    
    # Calculate fold changes
    if 'estrus' in stages and 'diestrus' in stages:
        fold_changes = (category_pivot['diestrus'] / (category_pivot['estrus'] + 1e-10))
        print("\nDiestrus/Estrus Fold Changes:")
        for pathway, fc in fold_changes.sort_values(ascending=False).items():
            print(f"  {pathway}: {fc:.2f}x")
    
    print(f"Saved category-by-stage plot to {output_path}")


def plot_flux_vs_gene_count(flux_df, output_path, module_gene_file=None):
    """
    Test correlation between module flux and number of genes per module.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux results in long format with 'reaction_id' and 'flux' columns
    output_path : Path
        Output file path
    module_gene_file : str or Path, optional
        Path to module-gene mapping file from scFEA
    """
    print(f"\nCreating flux vs. gene count analysis...")
    
    # Try to find the scFEA module gene file (m168 format)
    possible_paths = [
        Path('../../scFEA/data/module_gene_m168.csv'),
        Path('../../scFEA/data/module_gene_complete_mouse_m168.csv'),
        Path('../scFEA/data/module_gene_m168.csv'),
    ]
    
    if module_gene_file:
        possible_paths.insert(0, Path(module_gene_file))
    
    module_gene_path = None
    for path in possible_paths:
        if path.exists():
            module_gene_path = path
            print(f"Found module-gene file at: {path}")
            break
    
    if module_gene_path is None:
        print(f"Warning: Module-gene file not found. Tried:")
        for path in possible_paths:
            print(f"  - {path.absolute()}")
        print("Skipping flux vs. gene count analysis")
        return
    
    # Parse module-gene CSV to count genes per module
    # Format: M_1,gene1,gene2,gene3,...
    try:
        module_gene_counts = {}
        with open(module_gene_path, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) < 2:
                    continue
                module_id = parts[0]
                # Skip header rows or non-module rows
                if not module_id.startswith('M_'):
                    continue
                genes = parts[1:]  # All remaining columns are genes
                module_gene_counts[module_id] = len(genes)
        
        print(f"Counted genes for {len(module_gene_counts)} modules")
        print(f"Sample: {list(module_gene_counts.items())[:5]}")
        
    except Exception as e:
        print(f"Error parsing module-gene file: {e}")
        return
    
    # Calculate mean absolute flux per module across all cells
    mean_flux_per_module = flux_df.groupby('reaction_id')['flux'].apply(
        lambda x: np.abs(x).mean()
    ).reset_index()
    mean_flux_per_module.columns = ['module', 'mean_flux']
    
    print(f"Sample module names in flux data: {mean_flux_per_module['module'].head(10).tolist()}")
    
    # Extract M_XX pattern from annotated module names if needed
    import re
    def extract_module_id(name):
        match = re.search(r'(M_\d+)', str(name))
        return match.group(1) if match else name
    
    mean_flux_per_module['module_id'] = mean_flux_per_module['module'].apply(extract_module_id)
    
    # Add gene counts using extracted module IDs
    mean_flux_per_module['gene_count'] = mean_flux_per_module['module_id'].map(module_gene_counts)
    
    # Remove modules without gene count data
    analysis_df = mean_flux_per_module.dropna(subset=['gene_count'])
    
    if len(analysis_df) < 2:
        print(f"Error: Only {len(analysis_df)} modules matched")
        print(f"Tried to match: {mean_flux_per_module['module_id'].unique()[:10]}")
        print("Skipping analysis")
        return
    
    print(f"Analyzing {len(analysis_df)} modules")
    
    # Calculate correlations
    from scipy.stats import pearsonr, spearmanr, linregress
    pearson_r, pearson_p = pearsonr(analysis_df['gene_count'], analysis_df['mean_flux'])
    spearman_r, spearman_p = spearmanr(analysis_df['gene_count'], analysis_df['mean_flux'])
    
    print(f"Pearson correlation: r={pearson_r:.3f}, p={pearson_p:.2e}")
    print(f"Spearman correlation: r={spearman_r:.3f}, p={spearman_p:.2e}")
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Linear scale scatter plot
    axes[0].scatter(analysis_df['gene_count'], analysis_df['mean_flux'], 
                    alpha=0.6, s=50, color='steelblue')
    axes[0].set_xlabel('Number of Genes in Module', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Mean Metabolic Flux', fontsize=11, fontweight='bold')
    axes[0].set_title(f'Flux vs. Gene Count per Module\nPearson r={pearson_r:.3f}, p={pearson_p:.2e}', 
                      fontsize=12, fontweight='bold')
    axes[0].grid(alpha=0.3)
    
    # Add regression line
    slope, intercept, r_value, p_value, std_err = linregress(analysis_df['gene_count'], 
                                                               analysis_df['mean_flux'])
    x_line = np.linspace(analysis_df['gene_count'].min(), analysis_df['gene_count'].max(), 100)
    axes[0].plot(x_line, slope * x_line + intercept, 'r--', alpha=0.8, linewidth=2, 
                 label=f'y = {slope:.4f}x + {intercept:.4f}')
    axes[0].legend(loc='upper left', fontsize=9)
    
    # Log-scale version
    axes[1].scatter(analysis_df['gene_count'], analysis_df['mean_flux'], 
                    alpha=0.6, s=50, color='coral')
    axes[1].set_xlabel('Number of Genes in Module', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Mean Metabolic Flux (log scale)', fontsize=11, fontweight='bold')
    axes[1].set_yscale('log')
    axes[1].set_title(f'Flux vs. Gene Count (Log Scale)\nSpearman r={spearman_r:.3f}, p={spearman_p:.2e}', 
                      fontsize=12, fontweight='bold')
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print modules with highest gene counts
    print("\nModules with most genes:")
    top_gene_counts = analysis_df.nlargest(10, 'gene_count')[['module', 'gene_count', 'mean_flux']]
    print(top_gene_counts.to_string(index=False))
    
    print(f"Saved flux vs. gene count plot to {output_path}")
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
    
    # Use refined_cell_type if available (from subclustering), otherwise use cell_type
    if 'refined_cell_type' in cell_metadata.columns:
        celltype_col = 'refined_cell_type'
        print("✓ Using refined_cell_type from subclustering")
    else:
        celltype_col = 'cell_type'
        print("⚠ Using basic cell_type (no subclustering available)")
    
    # Merge flux results with cell types and stage (if available)
    metadata_cols_to_join = [celltype_col]
    if 'stage' in cell_metadata.columns:
        metadata_cols_to_join.append('stage')
        print("✓ Stage information available for stage-specific plots")
    else:
        print("⚠ No stage information - stage-specific plots will be skipped")
    
    flux_df = flux_df.join(cell_metadata[metadata_cols_to_join], how='left')
    flux_df.rename(columns={celltype_col: 'cell_type'}, inplace=True)
    
    # Filter out 'Ignore' and 'Inconclusive' cell types
    flux_df = flux_df[~flux_df['cell_type'].isin(['Ignore', 'Inconclusive'])]
    
    # Reshape from wide format (modules as columns) to long format (reaction_id, cell_type, flux)
    # Get module columns (M_1, M_2, etc.)
    module_cols = [col for col in flux_df.columns if col.startswith('M_')]
    
    # Determine which metadata columns to preserve during melt
    id_vars = ['cell_type']
    cols_to_keep = module_cols + ['cell_type']
    if 'stage' in flux_df.columns:
        id_vars.append('stage')
        cols_to_keep.append('stage')
    
    # Melt the dataframe to long format
    flux_df = flux_df[cols_to_keep].melt(
        id_vars=id_vars,
        value_vars=module_cols,
        var_name='reaction_id',
        value_name='flux'
    )
    
    print(f"Loaded {len(flux_df)} flux predictions for {flux_df['cell_type'].nunique()} cell types (with subclustering)")
    
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
    
    plot_flux_distribution_by_stage(
        flux_df,
        output_dir / f'flux_distribution_by_stage.{fmt}'
    )
    
    plot_flux_by_category_and_stage(
        flux_df,
        output_dir / f'flux_by_category_and_stage.{fmt}',
        annotations=annotations
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
    
    plot_flux_vs_gene_count(
        flux_df,
        output_dir / f'flux_vs_gene_count.{fmt}'
    )
    
    print("\nvisualisation complete!")
    print(f"Figures saved to {output_dir}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
