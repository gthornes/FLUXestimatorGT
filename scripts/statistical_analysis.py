#!/usr/bin/env python3
"""
Statistical analysis script for flux estimation results.

This script performs differential analysis comparing metabolic flux between:
- Estrus vs. Diestrus stages (within each cell type)
- Cell types (within each stage)
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sympy import re
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
        annotations_path = Path('data/reference/module_annotations.csv')
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


def calculate_effect_size(group1, group2):
    """
    Calculate Cohen's d effect size.
    
    Parameters
    ----------
    group1, group2 : array-like
        Two groups to compare
        
    Returns
    -------
    float
        Cohen's d effect size
    """
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    if pooled_std == 0:
        return 0.0
    
    return (np.mean(group1) - np.mean(group2)) / pooled_std


def differential_flux_analysis_by_stage(flux_df, cell_type, annotations=None, fdr_threshold=0.05):
    """
    Perform differential analysis comparing estrus vs. diestrus for a specific cell type.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux data in long format with 'reaction_id', 'cell_type', 'stage', 'flux' columns
    cell_type : str
        Cell type to analyze
    annotations : pd.DataFrame, optional
        Module annotations
    fdr_threshold : float
        FDR threshold for significance
        
    Returns
    -------
    pd.DataFrame
        Results table with statistics for each module
    """
    print(f"\nAnalyzing {cell_type}...")
    
    # Filter to this cell type
    ct_data = flux_df[flux_df['cell_type'] == cell_type].copy()
    
    # Check if both stages present
    stages = ct_data['stage'].unique()
    if len(stages) < 2:
        print(f"  Warning: Only {len(stages)} stage(s) found for {cell_type}")
        return None
    
    # Get estrus and diestrus labels
    estrus_stage = 'estrus' if 'estrus' in stages else stages[0]
    diestrus_stage = 'diestrus' if 'diestrus' in stages else stages[1]
    
    # Get all unique modules
    modules = ct_data['reaction_id'].unique()
    
    results = []
    
    for module in modules:
        module_data = ct_data[ct_data['reaction_id'] == module]
        
        estrus_flux = module_data[module_data['stage'] == estrus_stage]['flux'].values
        diestrus_flux = module_data[module_data['stage'] == diestrus_stage]['flux'].values
        
        # Skip if insufficient data
        if len(estrus_flux) < 3 or len(diestrus_flux) < 3:
            continue
        
        # Calculate statistics
        mean_estrus = np.mean(estrus_flux)
        mean_diestrus = np.mean(diestrus_flux)
        
        # Fold change (add small constant to avoid division by zero)
        epsilon = 1e-10
        fold_change = (mean_diestrus + epsilon) / (mean_estrus + epsilon)
        # Clip fold change to avoid log2(0) or log2(negative)
        fold_change = np.clip(fold_change, epsilon, None)
        log2_fc = np.log2(fold_change)
        
        # Wilcoxon rank-sum test (Mann-Whitney U)
        statistic, pval = stats.mannwhitneyu(
            estrus_flux, diestrus_flux, 
            alternative='two-sided'
        )
        
        # Effect size (Cohen's d)
        effect_size = calculate_effect_size(diestrus_flux, estrus_flux)
        
        results.append({
            'module': module,
            'cell_type': cell_type,
            'mean_estrus': mean_estrus,
            'mean_diestrus': mean_diestrus,
            'fold_change': fold_change,
            'log2_fold_change': log2_fc,
            'p_value': pval,
            'effect_size': effect_size,
            'n_estrus': len(estrus_flux),
            'n_diestrus': len(diestrus_flux)
        })
    
    if len(results) == 0:
        print(f"  No results for {cell_type}")
        return None
    
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction (Benjamini-Hochberg FDR)
    reject, pvals_corrected, _, _ = multipletests(
        results_df['p_value'], 
        alpha=fdr_threshold, 
        method='fdr_bh'
    )
    
    results_df['p_adjusted'] = pvals_corrected
    results_df['significant'] = reject
    
    # Add annotations if available
    if annotations is not None:
        # Extract module ID from potentially annotated names
        import re
        def extract_module_id(name):
            match = re.search(r'(M_\d+)', str(name))
            return match.group(1) if match else name
        
        results_df['module_id'] = results_df['module'].apply(extract_module_id)
        
        results_df = results_df.merge(
            annotations[['pathway', 'description']], 
            left_on='module_id', 
            right_index=True, 
            how='left'
        )
    
    # Sort by p-value
    results_df = results_df.sort_values('p_value')
    
    n_sig = results_df['significant'].sum()
    print(f"  Found {n_sig} significant modules (FDR < {fdr_threshold})")
    
    return results_df


def plot_volcano(results_df, cell_type, output_path, fdr_threshold=0.05, fc_threshold=1.5, xlim=(-10,10)):
    """
    Create volcano plot for differential analysis results.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Results from differential_flux_analysis_by_stage
    cell_type : str
        Cell type name for title
    output_path : Path
        Output file path
    fdr_threshold : float
        FDR significance threshold
    fc_threshold : float
        Fold-change threshold for highlighting
    xlim : tuple, optional
        X-axis limits as (min, max). If None, uses data range.
    """
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Calculate -log10(p-value)
    results_df = results_df.copy()
    results_df['-log10_pval'] = -np.log10(results_df['p_value'] + 1e-300)
    
    # Categorize points
    results_df['category'] = 'Not significant'
    
    # Significant and upregulated in diestrus
    results_df.loc[
        (results_df['significant']) & (results_df['fold_change'] >= fc_threshold),
        'category'
    ] = f'Diestrus enriched (FC≥{fc_threshold})'
    
    # Significant and downregulated in diestrus (upregulated in estrus)
    results_df.loc[
        (results_df['significant']) & (results_df['fold_change'] <= 1/fc_threshold),
        'category'
    ] = f'Estrus enriched (FC≤{1/fc_threshold:.2f})'
    
    # Significant but modest fold change
    results_df.loc[
        (results_df['significant']) & 
        (results_df['fold_change'] > 1/fc_threshold) & 
        (results_df['fold_change'] < fc_threshold),
        'category'
    ] = 'Significant (modest FC)'
    
    # Color mapping
    colors = {
        'Not significant': 'lightgray',
        'Significant (modest FC)': 'gold',
        f'Estrus enriched (FC≤{1/fc_threshold:.2f})': 'red',
        f'Diestrus enriched (FC≥{fc_threshold})': 'slateblue'
    }
    
    # Plot points by category
    for category in colors.keys():
        subset = results_df[results_df['category'] == category]
        ax.scatter(
            subset['log2_fold_change'],
            subset['-log10_pval'],
            c=colors[category],
            label=f"{category} (n={len(subset)})",
            alpha=0.6,
            s=50,
            edgecolors='black',
            linewidth=0.5
        )
    
    # Add threshold lines
    ax.axhline(-np.log10(fdr_threshold), color='red', linestyle='--', 
               linewidth=1.5, alpha=0.7, label=f'FDR = {fdr_threshold}')
    ax.axvline(np.log2(fc_threshold), color='blue', linestyle='--', 
               linewidth=1.5, alpha=0.7, label=f'FC = {fc_threshold}')
    ax.axvline(-np.log2(fc_threshold), color='blue', linestyle='--', 
               linewidth=1.5, alpha=0.7)
    
    # Label top significant modules
    n_labels = 15
    top_modules = results_df[results_df['significant']].nlargest(n_labels, '-log10_pval')
    
    for _, row in top_modules.iterrows():
        label = row.get('description', row['module'])
        if pd.notna(label) and len(str(label)) > 50:
            label = str(label)[:47] + '...'
        
        ax.annotate(
            label,
            xy=(row['log2_fold_change'], row['-log10_pval']),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=7,
            alpha=0.8,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3)
        )
    
    ax.set_xlabel(r'$\log_2$(Fold Change: Diestrus/Estrus)', fontsize=12, fontweight='bold')
    ax.set_ylabel(r'$-\log_{10}$(P-value)', fontsize=12, fontweight='bold')
    ax.set_title(f'Metabolic Flux Differential Analysis: {cell_type}\n(Estrus vs. Diestrus)', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.3)
    
    # Set x-axis limits if specified
    if xlim is not None:
        ax.set_xlim(xlim)
    
    # Add text box with summary stats
    n_sig_up = len(results_df[
        (results_df['significant']) & (results_df['fold_change'] >= fc_threshold)
    ])
    n_sig_down = len(results_df[
        (results_df['significant']) & (results_df['fold_change'] <= 1/fc_threshold)
    ])
    
    textstr = f'Total modules: {len(results_df)}\n'
    textstr += f'Significant: {results_df["significant"].sum()}\n'
    textstr += f'↑ Diestrus: {n_sig_up}\n'
    textstr += f'↑ Estrus: {n_sig_down}'
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.98, 0.02, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved volcano plot to {output_path}")


def plot_significant_modules_summary(all_results_df, output_path, fdr_threshold=0.05):
    """
    Create bar chart showing number of significant modules per cell type.
    
    Parameters
    ----------
    all_results_df : pd.DataFrame
        Combined results from all cell types
    output_path : Path
        Output file path
    fdr_threshold : float
        FDR threshold used
    """
    # Count significant modules per cell type
    summary = all_results_df.groupby('cell_type').agg({
        'significant': 'sum',
        'module': 'count'
    }).rename(columns={'module': 'total_modules', 'significant': 'n_significant'})
    
    summary['pct_significant'] = (summary['n_significant'] / summary['total_modules'] * 100)
    summary = summary.sort_values('n_significant', ascending=False)
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = np.arange(len(summary.index))
    bars = ax.bar(x, summary['n_significant'], color='steelblue', edgecolor='black', alpha=0.8)
    
    # Color bars by percentage
    colors = plt.cm.RdYlGn(summary['pct_significant'] / summary['pct_significant'].max())
    for bar, color in zip(bars, colors):
        bar.set_facecolor(color)
    
    # Add value labels on bars
    for i, (idx, row) in enumerate(summary.iterrows()):
        ax.text(i, row['n_significant'] + 0.5, 
                f"{int(row['n_significant'])}\n({row['pct_significant']:.1f}%)",
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax.set_xticks(x)
    ax.set_xticklabels(summary.index, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Number of Significant Modules', fontsize=12, fontweight='bold')
    ax.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax.set_title(f'Differential Metabolic Modules by Cell Type\n(FDR < {fdr_threshold})', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.grid(axis='y', alpha=0.3)
    
    # Add horizontal line for average
    mean_sig = summary['n_significant'].mean()
    ax.axhline(mean_sig, color='red', linestyle='--', linewidth=2, alpha=0.7,
               label=f'Mean = {mean_sig:.1f}')
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Saved summary bar chart to {output_path}")


def plot_clustered_heatmap_by_stage(flux_df, output_path, annotations=None, top_n=168):
    """
    Create clustered heatmap comparing estrus vs diestrus cell type metabolic profiles.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux data in long format with 'reaction_id', 'cell_type', 'stage', 'flux' columns
    output_path : Path
        Output file path
    annotations : pd.DataFrame, optional
        Module annotations for labeling
    top_n : int
        Number of most variable modules to display
    """
    print(f"\nCreating stage comparison heatmap (top {top_n} most variable modules)...")
    
    # Calculate mean flux per cell type per module per stage
    flux_summary = flux_df.groupby(['stage', 'cell_type', 'reaction_id'])['flux'].mean().reset_index()
    
    # Get stages
    stages = sorted(flux_summary['stage'].unique(), reverse=True)  # estrus first
    if len(stages) < 2:
        print(f"  Warning: Only {len(stages)} stage(s) found")
        return
    
    estrus_stage = 'estrus' if 'estrus' in stages else stages[0]
    diestrus_stage = 'diestrus' if 'diestrus' in stages else stages[1]
    
    # Create pivot tables for each stage
    estrus_data = flux_summary[flux_summary['stage'] == estrus_stage].pivot(
        index='reaction_id',
        columns='cell_type',
        values='flux'
    ).fillna(0)
    
    diestrus_data = flux_summary[flux_summary['stage'] == diestrus_stage].pivot(
        index='reaction_id',
        columns='cell_type',
        values='flux'
    ).fillna(0)
    
    # Select top N most variable modules (by variance across all cells and stages)
    all_data = pd.concat([estrus_data, diestrus_data], axis=0)
    module_variance = all_data.var(axis=1)
    top_modules = module_variance.nlargest(top_n).index
    
    estrus_data = estrus_data.loc[top_modules]
    diestrus_data = diestrus_data.loc[top_modules]
    
    # Only use cell types present in both stages
    common_cell_types = sorted(set(estrus_data.columns) & set(diestrus_data.columns))
    
    if len(common_cell_types) == 0:
        print("  Warning: No cell types found in both stages")
        return
    
    estrus_data = estrus_data[common_cell_types]
    diestrus_data = diestrus_data[common_cell_types]
    
    print(f"  Comparing {len(common_cell_types)} cell types present in both stages")
    
    # Transpose so cell types are on axes (estrus x-axis, diestrus y-axis)
    # Compute correlation between estrus and diestrus profiles for each cell type
    # This shows metabolic similarity between stages
    
    # Create correlation/comparison matrix
    # For each pair of cell types, compute correlation across modules
    from scipy.stats import pearsonr
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist, squareform
    
    cell_types = common_cell_types
    n_cells = len(cell_types)
    
    # Create a matrix showing correlation between estrus and diestrus for each cell type pair
    comparison_matrix = np.zeros((n_cells, n_cells))
    
    for i, ct_diestrus in enumerate(cell_types):
        for j, ct_estrus in enumerate(cell_types):
            # Correlation between diestrus profile of ct_i and estrus profile of ct_j
            diestrus_profile = diestrus_data[ct_diestrus].values
            estrus_profile = estrus_data[ct_estrus].values
            
            # Pearson correlation
            corr, _ = pearsonr(diestrus_profile, estrus_profile)
            comparison_matrix[i, j] = corr
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(14, 12))
    
    sns.heatmap(
        comparison_matrix,
        xticklabels=cell_types,
        yticklabels=cell_types,
        cmap='RdBu_r',
        center=0,
        vmin=-1,
        vmax=1,
        cbar_kws={'label': 'Pearson Correlation'},
        ax=ax,
        square=True,
        linewidths=0.5,
        linecolor='white'
    )
    
    ax.set_xlabel('Estrus Cell Types', fontsize=12, fontweight='bold')
    ax.set_ylabel('Diestrus Cell Types', fontsize=12, fontweight='bold')
    ax.set_title(f'Metabolic Profile Correlation: Estrus vs Diestrus\n(Based on Top {top_n} Variable Modules)',
                 fontsize=14, fontweight='bold', pad=20)
    
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved stage comparison heatmap to {output_path}")


def plot_module_correlation_network(flux_df, output_path, annotations=None, min_corr=0.7, top_n=30, stage=None, pathway_colors=None):
    """
    Create correlation network showing co-variation between metabolic modules.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux data in long format with 'reaction_id', 'cell_type', 'stage', 'flux' columns
    output_path : Path
        Output file path
    annotations : pd.DataFrame, optional
        Module annotations for labeling
    min_corr : float
        Minimum correlation threshold for displaying edges (default: 0.7)
    top_n : int
        Number of most variable modules to include (default: 30)
    stage : str, optional
        Specific stage to analyze (e.g., 'estrus' or 'diestrus'). If None, uses all stages.
    pathway_colors : dict, optional
        Predefined color mapping for pathways to ensure consistency across plots
    """
    import networkx as nx
    
    stage_label = f" ({stage})" if stage else " (all stages)"
    print(f"\nCreating module correlation network{stage_label} (top {top_n} modules, correlation ≥ {min_corr})...")
    
    # Filter by stage if specified
    if stage is not None:
        flux_df = flux_df[flux_df['stage'] == stage].copy()
    
    # Create pivot table with modules as columns and cell types as rows
    pivot = flux_df.pivot_table(
        index='cell_type',
        columns='reaction_id',
        values='flux',
        aggfunc='mean'
    ).fillna(0)
    
    # Select top N most variable modules
    module_variance = pivot.var(axis=0)
    top_modules = module_variance.nlargest(top_n).index
    pivot_top = pivot[top_modules]
    
    # Calculate correlation matrix between modules
    corr_matrix = pivot_top.corr()
    
    # Create network graph
    G = nx.Graph()
    
    # Add nodes with module names
    for module in corr_matrix.index:
        if annotations is not None:
            import re
            match = re.search(r'(M_\d+)', str(module))
            module_id = match.group(1) if match else module
            if module_id in annotations.index:
                label = annotations.loc[module_id, 'pathway']
                G.add_node(module, label=label, module_id=module_id)
            else:
                G.add_node(module, label=module, module_id=module)
        else:
            G.add_node(module, label=module, module_id=module)
    
    # Add edges for correlations above threshold
    for i, mod1 in enumerate(corr_matrix.index):
        for j, mod2 in enumerate(corr_matrix.columns):
            if i < j:  # Only upper triangle
                corr = corr_matrix.loc[mod1, mod2]
                if abs(corr) >= min_corr:
                    G.add_edge(mod1, mod2, weight=abs(corr), sign=np.sign(corr))
    
    if len(G.edges()) == 0:
        print(f"  Warning: No correlations above threshold {min_corr} found")
        print(f"  Max correlation: {corr_matrix.where(~np.eye(len(corr_matrix), dtype=bool)).max().max():.2f}")
        return
    
    print(f"  Network has {len(G.nodes())} nodes and {len(G.edges())} edges")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 14))
    
    # Use spring layout for positioning
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    # Get node colors based on pathway
    if annotations is not None:
        pathways = [G.nodes[node].get('label', '') for node in G.nodes()]
        unique_pathways = sorted(set(pathways))
        if pathway_colors is None:
            pathway_colors = dict(zip(unique_pathways, plt.cm.tab20(np.linspace(0, 1, len(unique_pathways)))))
        node_colors = [pathway_colors.get(G.nodes[node].get('label', ''), 'gray') for node in G.nodes()]
    else:
        node_colors = 'steelblue'
        unique_pathways = []
    
    # Draw edges with width proportional to correlation strength
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    signs = [G[u][v]['sign'] for u, v in edges]
    
    # Positive correlations in slateblue, negative in red
    edge_colors = ['slateblue' if s > 0 else 'red' for s in signs]
    
    nx.draw_networkx_edges(
        G, pos,
        width=[w * 3 for w in weights],
        alpha=0.6,
        edge_color=edge_colors,
        ax=ax
    )
    
    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=800,
        alpha=0.9,
        edgecolors='black',
        linewidths=1.5,
        ax=ax
    )
    
    # Draw labels
    labels = {node: G.nodes[node].get('module_id', node) for node in G.nodes()}
    nx.draw_networkx_labels(
        G, pos,
        labels=labels,
        font_size=7,
        font_weight='bold',
        ax=ax
    )
    
    # Add legend for pathways if available
    if annotations is not None and len(unique_pathways) <= 15:
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                     markerfacecolor=pathway_colors[pathway], 
                                     markersize=10, label=pathway, markeredgecolor='black')
                          for pathway in unique_pathways]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=8, 
                 title='Metabolic Pathway', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    
    title = f'Metabolic Module Correlation Network'
    if stage:
        title += f' - {stage.capitalize()}'
    title += f'\n(Top {top_n} Variable Modules, |r| ≥ {min_corr})'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved correlation network to {output_path}")
    
    # Print network statistics
    print(f"  Network statistics:")
    print(f"    - Nodes: {len(G.nodes())}")
    print(f"    - Edges: {len(G.edges())}")
    print(f"    - Average degree: {sum(dict(G.degree()).values()) / len(G.nodes()):.2f}")
    if len(G.edges()) > 0:
        avg_corr = np.mean([G[u][v]['weight'] for u, v in G.edges()])
        print(f"    - Average correlation: {avg_corr:.3f}")


def plot_overall_flux_comparison(flux_df, output_path, annotations=None, top_n=50):
    """
    Create scatter plot comparing mean flux between estrus and diestrus (pooled across cell types).
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux data in long format with 'reaction_id', 'cell_type', 'stage', 'flux' columns
    output_path : Path
        Output file path
    annotations : pd.DataFrame, optional
        Module annotations for labeling and coloring by pathway
    top_n : int
        Number of modules with largest flux differences to label
    """
    print("\nCreating overall flux comparison plot (estrus vs diestrus)...")
    
    # Calculate mean flux per module per stage (pooling cell types)
    flux_by_stage = flux_df.groupby(['stage', 'reaction_id'])['flux'].mean().reset_index()
    
    # Get stages
    stages = sorted(flux_by_stage['stage'].unique(), reverse=True)  # estrus first
    if len(stages) < 2:
        print(f"  Warning: Only {len(stages)} stage(s) found")
        return
    
    estrus_stage = 'estrus' if 'estrus' in stages else stages[0]
    diestrus_stage = 'diestrus' if 'diestrus' in stages else stages[1]
    
    # Pivot to get estrus and diestrus in separate columns
    flux_pivot = flux_by_stage.pivot(index='reaction_id', columns='stage', values='flux').fillna(0)
    
    if estrus_stage not in flux_pivot.columns or diestrus_stage not in flux_pivot.columns:
        print("  Error: Could not find both stages in data")
        return
    
    # Add pathway information if available
    if annotations is not None:
        import re
        def extract_module_id(name):
            match = re.search(r'(M_\d+)', str(name))
            return match.group(1) if match else name
        
        flux_pivot['module_id'] = [extract_module_id(m) for m in flux_pivot.index]
        flux_pivot = flux_pivot.merge(
            annotations[['pathway', 'description']],
            left_on='module_id',
            right_index=True,
            how='left'
        )
        flux_pivot['pathway'] = flux_pivot['pathway'].fillna('Unknown')
    
    # Calculate absolute difference for labeling
    flux_pivot['abs_diff'] = np.abs(flux_pivot[estrus_stage] - flux_pivot[diestrus_stage])
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Color by pathway if available
    if annotations is not None:
        pathways = flux_pivot['pathway'].unique()
        pathway_colors = dict(zip(sorted(pathways), plt.cm.tab20(np.linspace(0, 1, len(pathways)))))
        colors = [pathway_colors[p] for p in flux_pivot['pathway']]
        
        # Scatter plot with pathway colors
        for pathway in sorted(pathways):
            mask = flux_pivot['pathway'] == pathway
            ax.scatter(
                flux_pivot.loc[mask, estrus_stage],
                flux_pivot.loc[mask, diestrus_stage],
                c=[pathway_colors[pathway]],
                label=pathway,
                alpha=0.7,
                s=80,
                edgecolors='black',
                linewidth=0.5
            )
    else:
        ax.scatter(
            flux_pivot[estrus_stage],
            flux_pivot[diestrus_stage],
            alpha=0.7,
            s=80,
            c='steelblue',
            edgecolors='black',
            linewidth=0.5
        )
    
    # Add diagonal line (x=y)
    max_val = max(flux_pivot[estrus_stage].max(), flux_pivot[diestrus_stage].max())
    min_val = min(flux_pivot[estrus_stage].min(), flux_pivot[diestrus_stage].min())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, linewidth=2, label='x=y (no change)')
    
    # Label top N most different modules
    top_diff = flux_pivot.nlargest(top_n, 'abs_diff')
    for idx, row in top_diff.iterrows():
        if annotations is not None and 'description' in row:
            label = row['description']
            if pd.notna(label) and len(str(label)) > 30:
                label = str(label)[:27] + '...'
        else:
            label = idx
        
        ax.annotate(
            label,
            xy=(row[estrus_stage], row[diestrus_stage]),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=6,
            alpha=0.7,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3)
        )
    
    ax.set_xlabel('Mean Flux in Estrus (mmol/gDW/h)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean Flux in Diestrus (mmol/gDW/h)', fontsize=12, fontweight='bold')
    ax.set_title('Overall Metabolic Flux Comparison\n(Pooled Across All Cell Types)',
                 fontsize=14, fontweight='bold', pad=20)
    
    if annotations is not None:
        ax.legend(loc='upper left', fontsize=8, bbox_to_anchor=(1.02, 1), borderaxespad=0, title='Pathway')
    else:
        ax.legend()
    
    ax.grid(alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved overall flux comparison to {output_path}")
    
    # Print statistics
    corr = flux_pivot[[estrus_stage, diestrus_stage]].corr().iloc[0, 1]
    print(f"  Overall correlation: {corr:.3f}")
    n_higher_estrus = (flux_pivot[estrus_stage] > flux_pivot[diestrus_stage]).sum()
    n_higher_diestrus = (flux_pivot[diestrus_stage] > flux_pivot[estrus_stage]).sum()
    print(f"  Modules with higher flux in estrus: {n_higher_estrus}")
    print(f"  Modules with higher flux in diestrus: {n_higher_diestrus}")


def plot_overall_celltype_comparison(flux_df, output_path, top_n=10):
    """
    Create scatter plot comparing mean flux between estrus and diestrus for each cell type (pooled across modules).
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux data in long format with 'reaction_id', 'cell_type', 'stage', 'flux' columns
    output_path : Path
        Output file path
    top_n : int
        Number of cell types with largest flux differences to label
    """
    print("\nCreating overall cell type flux comparison plot (estrus vs diestrus)...")
    
    # Calculate mean flux per cell type per stage (pooling modules)
    flux_by_stage = flux_df.groupby(['stage', 'cell_type'])['flux'].mean().reset_index()
    
    # Get stages
    stages = sorted(flux_by_stage['stage'].unique(), reverse=True)  # estrus first
    if len(stages) < 2:
        print(f"  Warning: Only {len(stages)} stage(s) found")
        return
    
    estrus_stage = 'estrus' if 'estrus' in stages else stages[0]
    diestrus_stage = 'diestrus' if 'diestrus' in stages else stages[1]
    
    # Pivot to get estrus and diestrus in separate columns
    flux_pivot = flux_by_stage.pivot(index='cell_type', columns='stage', values='flux')
    
    if estrus_stage not in flux_pivot.columns or diestrus_stage not in flux_pivot.columns:
        print("  Error: Could not find both stages in data")
        return
    
    # Only keep cell types found in both stages (remove rows with NaN)
    flux_pivot = flux_pivot.dropna()
    print(f"  Cell types found in both stages: {len(flux_pivot)}")
    
    # Calculate absolute difference for labeling
    flux_pivot['abs_diff'] = np.abs(flux_pivot[estrus_stage] - flux_pivot[diestrus_stage])
    
    # Determine which cell types are above/below the diagonal
    # Below diagonal: estrus > diestrus (higher in estrus) -> red
    # Above diagonal: diestrus > estrus (higher in diestrus) -> slateblue
    flux_pivot['color'] = flux_pivot.apply(
        lambda row: 'red' if row[estrus_stage] > row[diestrus_stage] else 'slateblue',
        axis=1
    )
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Scatter plot with colors based on position relative to diagonal
    for color in ['red', 'slateblue']:
        mask = flux_pivot['color'] == color
        label = 'Higher in Estrus' if color == 'red' else 'Higher in Diestrus'
        ax.scatter(
            flux_pivot.loc[mask, estrus_stage],
            flux_pivot.loc[mask, diestrus_stage],
            alpha=0.7,
            s=120,
            c=color,
            edgecolors='black',
            linewidth=0.8,
            label=label
        )
    
    # Add diagonal line (x=y)
    max_val = max(flux_pivot[estrus_stage].max(), flux_pivot[diestrus_stage].max())
    min_val = min(flux_pivot[estrus_stage].min(), flux_pivot[diestrus_stage].min())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, linewidth=2, label='x=y (no change)')
    
    # Label top N most different cell types (or all if there are fewer than top_n)
    n_to_label = min(top_n, len(flux_pivot))
    top_diff = flux_pivot.nlargest(n_to_label, 'abs_diff')
    for idx, row in top_diff.iterrows():
        label = str(idx)
        
        ax.annotate(
            label,
            xy=(row[estrus_stage], row[diestrus_stage]),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=8,
            alpha=0.8,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.4)
        )
    
    ax.set_xlabel('Mean Flux in Estrus (mmol/gDW/h)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean Flux in Diestrus (mmol/gDW/h)', fontsize=12, fontweight='bold')
    ax.set_title('Overall Cell Type Flux Comparison\n(Pooled Across All Modules)',
                 fontsize=14, fontweight='bold', pad=20)
    
    ax.legend(loc='best', fontsize=10)
    ax.grid(alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved overall cell type flux comparison to {output_path}")
    
    # Print statistics
    corr = flux_pivot[[estrus_stage, diestrus_stage]].corr().iloc[0, 1]
    print(f"  Overall correlation: {corr:.3f}")
    n_higher_estrus = (flux_pivot[estrus_stage] > flux_pivot[diestrus_stage]).sum()
    n_higher_diestrus = (flux_pivot[diestrus_stage] > flux_pivot[estrus_stage]).sum()
    print(f"  Cell types with higher flux in estrus: {n_higher_estrus}")
    print(f"  Cell types with higher flux in diestrus: {n_higher_diestrus}")


def plot_flux_distribution_by_stage(flux_df, output_path):
    """
    Create comprehensive flux distribution plots for each stage.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux data in long format with 'reaction_id', 'cell_type', 'stage', 'flux' columns
    output_path : Path
        Output file path
    """
    print("\nCreating flux distribution plots by stage...")
    
    # Get stages (estrus first)
    stages = sorted(flux_df['stage'].unique(), reverse=True)  # estrus before diestrus alphabetically reversed
    
    # Filter out zero/negative flux values for log scale
    flux_df_pos = flux_df[flux_df['flux'] > 0].copy()
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Color scheme
    colors = {'estrus': 'red', 'diestrus': 'slateblue'}
    
    # Calculate statistics
    stage_stats = []
    for stage in stages:
        stage_data = flux_df_pos[flux_df_pos['stage'] == stage]['flux']
        mean_val = stage_data.mean()
        median_val = stage_data.median()
        std_val = stage_data.std()
        stage_stats.append((stage, mean_val, median_val, std_val))
        print(f"  {stage.capitalize()}: mean={mean_val:.4f}, median={median_val:.4f}, std={std_val:.4f}")
    
    # 1. Overlapping KDE plots (log scale)
    ax1 = axes[0, 0]
    for stage in stages:
        stage_data = flux_df_pos[flux_df_pos['stage'] == stage]['flux']
        log_data = np.log10(stage_data)
        ax1.hist(log_data, bins=50, alpha=0.5, label=stage.capitalize(), 
                color=colors[stage], edgecolor='black', linewidth=0.5, density=True)
    
    ax1.set_xlabel(r'$\log_{10}$(Flux)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Density', fontsize=11, fontweight='bold')
    ax1.set_title('Flux Distribution (Log Scale)', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(alpha=0.3)
    
    # 2. Box plots (log scale)
    ax2 = axes[0, 1]
    box_data = [np.log10(flux_df_pos[flux_df_pos['stage'] == stage]['flux'].values) for stage in stages]
    bp = ax2.boxplot(box_data, tick_labels=[s.capitalize() for s in stages], 
                     patch_artist=True, widths=0.6, showfliers=False)
    
    for patch, stage in zip(bp['boxes'], stages):
        patch.set_facecolor(colors[stage])
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
        patch.set_linewidth(1.5)
    
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp[element], color='black', linewidth=2)
    
    ax2.set_ylabel(r'$\log_{10}$(Flux)', fontsize=11, fontweight='bold')
    ax2.set_title('Flux Distribution Box Plots', fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. ECDF (Empirical Cumulative Distribution Function)
    ax3 = axes[1, 0]
    for stage in stages:
        stage_data = np.sort(flux_df_pos[flux_df_pos['stage'] == stage]['flux'].values)
        ecdf = np.arange(1, len(stage_data) + 1) / len(stage_data)
        ax3.plot(stage_data, ecdf, label=stage.capitalize(), 
                color=colors[stage], linewidth=2.5, alpha=0.8)
    
    ax3.set_xlabel('Flux (mmol/gDW/h)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Cumulative Probability', fontsize=11, fontweight='bold')
    ax3.set_xscale('log')
    ax3.set_title('Cumulative Distribution Function', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(alpha=0.3)
    
    # 4. Summary statistics table
    ax4 = axes[1, 1]
    ax4.axis('tight')
    ax4.axis('off')
    
    # Create table data
    table_data = []
    headers = ['Stage', 'Mean', 'Median', 'Std Dev', 'Min', 'Max', 'N']
    for stage in stages:
        stage_data = flux_df_pos[flux_df_pos['stage'] == stage]['flux']
        row = [
            stage.capitalize(),
            f'{stage_data.mean():.4f}',
            f'{stage_data.median():.4f}',
            f'{stage_data.std():.4f}',
            f'{stage_data.min():.2e}',
            f'{stage_data.max():.2f}',
            f'{len(stage_data):,}'
        ]
        table_data.append(row)
    
    table = ax4.table(cellText=table_data, colLabels=headers,
                     cellLoc='center', loc='center',
                     colColours=['lightgray']*len(headers))
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2.5)
    
    # Color cells by stage
    for i, stage in enumerate(stages):
        table[(i+1, 0)].set_facecolor(colors[stage])
        table[(i+1, 0)].set_alpha(0.3)
    
    ax4.set_title('Summary Statistics', fontsize=12, fontweight='bold', pad=20)
    
    plt.suptitle('Overall Flux Distribution by Stage\n(Pooled Across All Cell Types and Modules)',
                fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved comprehensive flux distribution plot to {output_path}")


def perform_pca_biplot(flux_df, output_path):
    """
    Perform PCA and create biplot comparing estrus vs diestrus profiles.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux data in long format with 'reaction_id', 'cell_type', 'stage', 'flux' columns
    output_path : Path
        Output file path
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    print("\nPerforming PCA for stage comparison...")
    
    # Create pivot table with modules as features and samples as rows
    pivot = flux_df.pivot_table(
        index=['cell_type', 'stage'],
        columns='reaction_id',
        values='flux',
        aggfunc='mean'
    ).fillna(0)
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(pivot.values)
    
    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    # Create biplot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Extract cell types and stages from multi-index
    cell_types = pivot.index.get_level_values('cell_type')
    stages = pivot.index.get_level_values('stage')
    
    # Define markers for stages
    stage_markers = {'estrus': 'o', 'diestrus': 's'}
    stage_colors = {'estrus': 'red', 'diestrus': 'slateblue'}
    
    # Get unique cell types and assign colors
    unique_cell_types = sorted(cell_types.unique())
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_cell_types)))
    cell_type_colors = dict(zip(unique_cell_types, colors))
    
    # Plot each combination of cell type and stage
    for cell_type in unique_cell_types:
        color = cell_type_colors[cell_type]
        for stage in stages.unique():
            mask = (cell_types == cell_type) & (stages == stage)
            marker = stage_markers.get(stage, 'o')
            
            # Only add label for first stage to avoid duplicate legend entries
            label = cell_type if stage == sorted(stages.unique())[0] else None
            
            ax.scatter(
                X_pca[mask, 0],
                X_pca[mask, 1],
                label=label,
                marker=marker,
                color=color,
                s=100,
                alpha=0.8,
                edgecolor='black',
                linewidth=0.5
            )
    
    # Add custom legend for markers (estrus first)
    legend_elements = [Line2D([0], [0], marker='o', color='w', 
                              markerfacecolor='red', markersize=10, 
                              label='Estrus', markeredgecolor='black'),
                      Line2D([0], [0], marker='s', color='w', 
                              markerfacecolor='slateblue', markersize=10, 
                              label='Diestrus', markeredgecolor='black')]
    
    # Get main legend and add marker legend
    leg1 = ax.legend(loc='upper left', fontsize=8, title='Cell Type', 
                     bbox_to_anchor=(1.02, 1), borderaxespad=0)
    ax.add_artist(leg1)
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9, 
             title='Stage', bbox_to_anchor=(1.02, 0.3), borderaxespad=0)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', fontsize=12, fontweight='bold')
    ax.set_title('PCA Biplot: Estrus vs Diestrus Metabolic Profiles', fontsize=14, fontweight='bold', pad=20)
    
    plt.legend(title='Stage')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved PCA biplot to {output_path}")

def perform_overall_differential_analysis(flux_df, annotations=None, fdr_threshold=0.05):
    """
    Perform differential analysis across ALL cells (ignoring cell type).
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux data in long format
    annotations : pd.DataFrame, optional
        Module annotations
    fdr_threshold : float
        FDR threshold
        
    Returns
    -------
    pd.DataFrame
        Results table with statistics for each module
    """
    print("\nPerforming overall differential analysis (all cell types combined)...")
    
    # Get stages
    stages = flux_df['stage'].unique()
    if len(stages) < 2:
        print(f"  Warning: Only {len(stages)} stage(s) found")
        return None
    
    estrus_stage = 'estrus' if 'estrus' in stages else stages[0]
    diestrus_stage = 'diestrus' if 'diestrus' in stages else stages[1]
    
    # Get all unique modules
    modules = flux_df['reaction_id'].unique()
    
    results = []
    
    for module in modules:
        module_data = flux_df[flux_df['reaction_id'] == module]
        
        estrus_flux = module_data[module_data['stage'] == estrus_stage]['flux'].values
        diestrus_flux = module_data[module_data['stage'] == diestrus_stage]['flux'].values
        
        # Skip if insufficient data
        if len(estrus_flux) < 10 or len(diestrus_flux) < 10:
            continue
        
        # Calculate statistics
        mean_estrus = np.mean(estrus_flux)
        mean_diestrus = np.mean(diestrus_flux)
        
        # Fold change
        epsilon = 1e-10
        fold_change = (mean_diestrus + epsilon) / (mean_estrus + epsilon)
        fold_change = np.clip(fold_change, epsilon, None)
        log2_fc = np.log2(fold_change)
        
        # Mann-Whitney U test
        statistic, pval = stats.mannwhitneyu(
            estrus_flux, diestrus_flux, 
            alternative='two-sided'
        )
        
        # Effect size
        effect_size = calculate_effect_size(diestrus_flux, estrus_flux)
        
        results.append({
            'module': module,
            'mean_estrus': mean_estrus,
            'mean_diestrus': mean_diestrus,
            'fold_change': fold_change,
            'log2_fold_change': log2_fc,
            'p_value': pval,
            'effect_size': effect_size,
            'n_estrus': len(estrus_flux),
            'n_diestrus': len(diestrus_flux)
        })
    
    if len(results) == 0:
        print("  No results")
        return None
    
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction
    reject, pvals_corrected, _, _ = multipletests(
        results_df['p_value'], 
        alpha=fdr_threshold, 
        method='fdr_bh'
    )
    
    results_df['p_adjusted'] = pvals_corrected
    results_df['significant'] = reject
    
    # Add annotations
    if annotations is not None:
        import re
        
        def extract_module_id(name):
            match = re.search(r'(M_\d+)', str(name))
            return match.group(1) if match else name
        
        results_df['module_id'] = results_df['module'].apply(extract_module_id)
        results_df = results_df.merge(
            annotations[['pathway', 'description']], 
            left_on='module_id', 
            right_index=True, 
            how='left'
        )
    
    results_df = results_df.sort_values('p_value')
    
    n_sig = results_df['significant'].sum()
    print(f"  Found {n_sig} significant modules (FDR < {fdr_threshold})")
    
    return results_df


def main():
    parser = argparse.ArgumentParser(
        description='Perform statistical analysis on flux results'
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input flux results CSV file'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory for results and figures'
    )
    parser.add_argument(
        '--config', '-c',
        default='config/analysis_config.yaml',
        help='Configuration file (default: config/analysis_config.yaml)'
    )
    parser.add_argument(
        '--fdr',
        type=float,
        default=0.05,
        help='FDR threshold for significance (default: 0.05)'
    )
    parser.add_argument(
        '--fc-threshold',
        type=float,
        default=1.5,
        help='Fold-change threshold for highlighting (default: 1.5)'
    )
    parser.add_argument(
        '--xlim',
        type=float,
        nargs=2,
        default=(-10, 10),
        metavar=('MIN', 'MAX'),
        help='X-axis limits for volcano plots (default: -10 10, use --xlim None None for auto)'
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
    print(f"\nLoading flux results from {args.input}")
    flux_df = pd.read_csv(args.input, index_col=0)
    
    # Load cell metadata
    metadata_path = Path(args.input).parent / 'cell_metadata.csv'
    if not metadata_path.exists():
        print(f"Error: Cell metadata file not found at {metadata_path}")
        return 1
    
    cell_metadata = pd.read_csv(metadata_path, index_col=0)
    
    # Use refined_cell_type if available
    if 'refined_cell_type' in cell_metadata.columns:
        celltype_col = 'refined_cell_type'
        print("Using refined_cell_type from subclustering")
    else:
        celltype_col = 'cell_type'
        print("Using basic cell_type (no subclustering available)")
    
    # Check for stage column
    if 'stage' not in cell_metadata.columns:
        print("Error: No 'stage' column found in metadata. Cannot perform stage comparison.")
        return 1
    
    # Merge metadata
    metadata_cols = [celltype_col, 'stage']
    flux_df = flux_df.join(cell_metadata[metadata_cols], how='left')
    flux_df.rename(columns={celltype_col: 'cell_type'}, inplace=True)
    
    # Filter out unwanted cell types
    flux_df = flux_df[~flux_df['cell_type'].isin(['Ignore', 'Inconclusive'])]
    
    # Get module columns
    module_cols = [col for col in flux_df.columns if col.startswith('M_')]
    
    # Reshape to long format
    flux_df = flux_df[module_cols + ['cell_type', 'stage']].melt(
        id_vars=['cell_type', 'stage'],
        value_vars=module_cols,
        var_name='reaction_id',
        value_name='flux'
    )
    
    print(f"\nAnalyzing {flux_df['cell_type'].nunique()} cell types")
    print(f"Stages: {flux_df['stage'].unique()}")
    
    # Create output directories
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tables_dir = output_dir / 'tables'
    tables_dir.mkdir(exist_ok=True)
    
    figures_dir = output_dir / 'figures'
    figures_dir.mkdir(exist_ok=True)
    
    fmt = config['visualisation']['figure_format']
    
    # Perform overall differential analysis (all cells combined)
    print("\n" + "="*60)
    print("OVERALL DIFFERENTIAL ANALYSIS (ALL CELLS)")
    print("="*60)
    
    overall_results = perform_overall_differential_analysis(
        flux_df,
        annotations=annotations,
        fdr_threshold=args.fdr
    )
    
    if overall_results is not None:
        # Save overall results
        overall_file = tables_dir / 'differential_flux_overall.csv'
        overall_results.to_csv(overall_file, index=False)
        print(f"✓ Saved overall results to {overall_file}")
        
        # Create overall volcano plot
        plot_volcano(
            overall_results,
            'All Cell Types Combined',
            figures_dir / f'volcano_overall.{fmt}',
            fdr_threshold=args.fdr,
            fc_threshold=args.fc_threshold,
            xlim=args.xlim
        )
    
    # Perform differential analysis for each cell type
    print("\n" + "="*60)
    print("DIFFERENTIAL FLUX ANALYSIS: ESTRUS vs. DIESTRUS")
    print("="*60)
    
    all_results = []
    
    for cell_type in sorted(flux_df['cell_type'].unique()):
        # Perform analysis
        results_df = differential_flux_analysis_by_stage(
            flux_df, 
            cell_type, 
            annotations=annotations,
            fdr_threshold=args.fdr
        )
        
        if results_df is None or len(results_df) == 0:
            continue
        
        all_results.append(results_df)
        
        # Sanitize cell type name for file paths
        safe_name = cell_type.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "").replace("+", "plus")
        
        # Save results table
        output_file = tables_dir / f'differential_flux_{safe_name}.csv'
        results_df.to_csv(output_file, index=False)
        print(f"  Saved results to {output_file}")
        
        # Create volcano plot
        volcano_file = figures_dir / f'volcano_{safe_name}.{fmt}'
        plot_volcano(
            results_df, 
            cell_type, 
            volcano_file,
            fdr_threshold=args.fdr,
            fc_threshold=args.fc_threshold,
            xlim=args.xlim
        )
    
    # Combine all results
    if len(all_results) > 0:
        combined_results = pd.concat(all_results, ignore_index=True)
        combined_file = tables_dir / 'differential_flux_all_cell_types.csv'
        combined_results.to_csv(combined_file, index=False)
        print(f"\n✓ Saved combined results to {combined_file}")
        
        # Create summary bar chart
        plot_significant_modules_summary(
            combined_results,
            figures_dir / f'significant_modules_summary.{fmt}',
            fdr_threshold=args.fdr
        )
        
        # Summary statistics
        print("\n" + "="*60)
        print("SUMMARY STATISTICS")
        print("="*60)
        
        summary = combined_results.groupby('cell_type').agg({
            'significant': 'sum',
            'module': 'count'
        }).rename(columns={'module': 'total_modules', 'significant': 'significant_modules'})
        
        summary['pct_significant'] = (summary['significant_modules'] / summary['total_modules'] * 100).round(2)
        
        print("\nSignificant modules per cell type:")
        print(summary.to_string())
        
        # Most altered modules across cell types
        print("\nTop 10 most consistently altered modules (by average -log10(p)):")
        combined_results['-log10_pval'] = -np.log10(combined_results['p_value'] + 1e-300)
        top_modules = combined_results.groupby('module').agg({
            '-log10_pval': 'mean',
            'significant': 'sum',
            'description': 'first'
        }).sort_values('-log10_pval', ascending=False).head(10)
        
        print(top_modules.to_string())
    
    # Create PCA plot
    print("\n" + "="*60)
    print("PERFORMING PCA AND CREATING BIPLOT")
    print("="*60)

    perform_pca_biplot(
        flux_df,
        figures_dir / f'pca_biplot.{fmt}'
    )

    # Create clustered heatmap
    print("\n" + "="*60)
    print("GENERATING CLUSTERED HEATMAP")
    print("="*60)
    
    plot_clustered_heatmap_by_stage(
        flux_df,
        figures_dir / f'clustered_heatmap.{fmt}',
        annotations=annotations,
        top_n=50
    )
    
    # Create consistent pathway color mapping for network plots
    if annotations is not None:
        unique_pathways = sorted(annotations['pathway'].unique())
        pathway_colors = dict(zip(unique_pathways, plt.cm.tab20(np.linspace(0, 1, len(unique_pathways)))))
    else:
        pathway_colors = None
    
    # Create correlation networks by stage
    print("\n" + "="*60)
    print("GENERATING MODULE CORRELATION NETWORKS BY STAGE")
    print("="*60)
    
    stages = flux_df['stage'].unique()
    for stage in sorted(stages, reverse=True):  # estrus first
        plot_module_correlation_network(
            flux_df,
            figures_dir / f'correlation_network_{stage}.{fmt}',
            annotations=annotations,
            min_corr=0.7,
            top_n=30,
            stage=stage,
            pathway_colors=pathway_colors
        )
    
    # Create overall flux comparison plot
    print("\n" + "="*60)
    print("GENERATING OVERALL FLUX COMPARISON")
    print("="*60)
    
    plot_overall_flux_comparison(
        flux_df,
        figures_dir / f'overall_flux_comparison.{fmt}',
        annotations=annotations,
        top_n=50
    )
    
    # Create overall cell type comparison plot
    print("\n" + "="*60)
    print("GENERATING OVERALL CELL TYPE COMPARISON")
    print("="*60)
    
    plot_overall_celltype_comparison(
        flux_df,
        figures_dir / f'overall_celltype_comparison.{fmt}',
        top_n=10
    )
    
    # Create flux distribution violin plot
    print("\n" + "="*60)
    print("GENERATING FLUX DISTRIBUTION VIOLIN PLOT")
    print("="*60)
    
    plot_flux_distribution_by_stage(
        flux_df,
        figures_dir / f'flux_distribution_by_stage.{fmt}'
    )
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nResults saved to: {output_dir}")
    print(f"  - Tables: {tables_dir}")
    print(f"  - Figures: {figures_dir}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
