#!/usr/bin/env python3
"""
Flux analysis script for metabolic modeling.

This script performs flux balance analysis on single-cell data using
gene expression constraints.
"""

import argparse
import sys
from pathlib import Path

import cobra
import scanpy as sc
import numpy as np
import pandas as pd
import yaml
from tqdm import tqdm


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def load_metabolic_model(model_name, reference_path='data/reference'):
    """
    Load genome-scale metabolic model.
    
    Parameters
    ----------
    model_name : str
        Name of the model (e.g., 'recon3d', 'iMM1415')
    reference_path : str or Path
        Path to reference data directory
        
    Returns
    -------
    model : cobra.Model
        Metabolic model
    """
    reference_path = Path(reference_path)
    
    print(f"Loading metabolic model: {model_name}")
    
    # Try to load from file
    model_file = reference_path / f"{model_name}.xml"
    if model_file.exists():
        model = cobra.io.read_sbml_model(str(model_file))
    else:
        # Try loading from BiGG database
        try:
            print(f"Model file not found locally. Attempting to load from BiGG database...")
            model = cobra.io.load_model(model_name)
        except Exception as e:
            raise FileNotFoundError(
                f"Could not load model {model_name}. "
                f"Please download it to {model_file} or use a BiGG model name. "
                f"Error: {str(e)}"
            )
    
    print(f"Model loaded: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, {len(model.genes)} genes")
    return model


def map_genes_to_reactions(model, gene_expression):
    """
    Map gene expression to reaction bounds.
    
    Parameters
    ----------
    model : cobra.Model
        Metabolic model
    gene_expression : pd.Series
        Gene expression values (gene symbols as index)
        
    Returns
    -------
    reaction_expression : pd.Series
        Reaction expression scores
    """
    reaction_expression = pd.Series(index=[r.id for r in model.reactions], dtype=float)
    
    for reaction in model.reactions:
        if len(reaction.genes) == 0:
            # No gene association - assume active
            reaction_expression[reaction.id] = 1.0
        else:
            # Evaluate GPR rule
            gene_values = []
            for gene in reaction.genes:
                if gene.name in gene_expression.index:
                    gene_values.append(gene_expression[gene.name])
                else:
                    gene_values.append(0.0)
            
            if len(gene_values) > 0:
                # NOTE: This is a simplified approach using mean expression.
                # For production use, implement proper GPR Boolean logic:
                # - AND relationships: use minimum expression
                # - OR relationships: use maximum expression
                # This would require parsing reaction.gene_reaction_rule
                reaction_expression[reaction.id] = np.mean(gene_values)
            else:
                reaction_expression[reaction.id] = 0.0
    
    return reaction_expression


def constrain_model_by_expression(model, gene_expression, threshold=0.1):
    """
    Constrain model flux bounds based on gene expression.
    
    Parameters
    ----------
    model : cobra.Model
        Metabolic model
    gene_expression : pd.Series
        Gene expression values
    threshold : float
        Expression threshold for reaction activity
        
    Returns
    -------
    model : cobra.Model
        Constrained model
    """
    # Get reaction expression scores
    reaction_expression = map_genes_to_reactions(model, gene_expression)
    
    # Constrain reactions based on expression
    for reaction_id, expr_value in reaction_expression.items():
        reaction = model.reactions.get_by_id(reaction_id)
        
        if expr_value < threshold:
            # Low expression - constrain flux to zero
            reaction.lower_bound = 0
            reaction.upper_bound = 0
        # else: keep original bounds
    
    return model


def perform_fba(model, objective=None):
    """
    Perform flux balance analysis.
    
    Parameters
    ----------
    model : cobra.Model
        Metabolic model
    objective : str, optional
        Objective function (reaction ID)
        
    Returns
    -------
    solution : cobra.Solution
        FBA solution
    """
    if objective:
        model.objective = objective
    
    solution = model.optimize()
    
    return solution


def analyze_cell_type(adata, cell_type, model_template, config):
    """
    Perform FBA for a specific cell type.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cell type annotations
    cell_type : str
        Cell type label
    model_template : cobra.Model
        Template metabolic model
    config : dict
        Configuration parameters
        
    Returns
    -------
    flux_df : pd.DataFrame
        Flux distribution for all reactions
    """
    # Get cells of this type
    cells_mask = adata.obs['cell_type'] == cell_type
    if cells_mask.sum() == 0:
        print(f"Warning: No cells found for cell type {cell_type}")
        return None
    
    # Calculate mean expression for this cell type
    if 'log1p_norm' in adata.layers:
        expression_data = adata[cells_mask].layers['log1p_norm']
    else:
        expression_data = adata[cells_mask].X
    
    mean_expression = np.array(expression_data.mean(axis=0)).flatten()
    gene_expression = pd.Series(mean_expression, index=adata.var_names)
    
    # Create model copy and constrain
    model = model_template.copy()
    threshold = config['flux_analysis']['flux_threshold']
    model = constrain_model_by_expression(model, gene_expression, threshold)
    
    # Set solver
    model.solver = config['flux_analysis']['solver']
    
    # Perform FBA
    solution = perform_fba(model, config['flux_analysis']['objective_function'])
    
    if solution.status == 'optimal':
        flux_df = pd.DataFrame({
            'reaction_id': [r.id for r in model.reactions],
            'reaction_name': [r.name for r in model.reactions],
            'flux': [solution.fluxes[r.id] for r in model.reactions],
            'cell_type': cell_type
        })
        return flux_df
    else:
        print(f"Warning: FBA did not converge for {cell_type} (status: {solution.status})")
        return None


def aggregate_fluxes_by_pathway(flux_df, pathway_mapping=None):
    """
    Aggregate fluxes by metabolic pathway.
    
    Parameters
    ----------
    flux_df : pd.DataFrame
        Flux distribution
    pathway_mapping : dict, optional
        Mapping of reactions to pathways
        
    Returns
    -------
    pathway_df : pd.DataFrame
        Pathway-level flux aggregation
    """
    # TODO: Implement pathway mapping
    # This requires a pathway annotation file
    print("Pathway aggregation not yet implemented")
    return flux_df


def main():
    parser = argparse.ArgumentParser(
        description='Perform flux balance analysis on single-cell data'
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input processed data file (.h5ad) with cell type annotations'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output CSV file for flux results'
    )
    parser.add_argument(
        '--config', '-c',
        default='config/analysis_config.yaml',
        help='Configuration file (default: config/analysis_config.yaml)'
    )
    parser.add_argument(
        '--model', '-m',
        help='Metabolic model name (overrides config)'
    )
    parser.add_argument(
        '--cell-types',
        nargs='+',
        help='Specific cell types to analyze (default: all)'
    )
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility (default, may be overridden by config)
    np.random.seed(42)
    
    # Load configuration
    print(f"Loading configuration from {args.config}")
    config = load_config(args.config)
    
    # Update random seed from config if specified
    if 'compute' in config and 'random_state' in config['compute']:
        np.random.seed(config['compute']['random_state'])
    
    # Load data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Check for cell type annotations
    if 'cell_type' not in adata.obs.columns:
        raise ValueError(
            "Cell type annotations not found. "
            "Please run cell type annotation first."
        )
    
    # Load metabolic model
    model_name = args.model if args.model else config['flux_analysis']['model']
    model = load_metabolic_model(model_name)
    
    # Get cell types to analyze
    if args.cell_types:
        cell_types = args.cell_types
    else:
        cell_types = adata.obs['cell_type'].unique()
    
    print(f"\nAnalyzing {len(cell_types)} cell types: {', '.join(cell_types)}")
    
    # Analyze each cell type
    all_results = []
    for cell_type in tqdm(cell_types, desc="Analyzing cell types"):
        flux_df = analyze_cell_type(adata, cell_type, model, config)
        if flux_df is not None:
            all_results.append(flux_df)
    
    # Combine results
    if len(all_results) > 0:
        results_df = pd.concat(all_results, ignore_index=True)
        
        # Save results
        print(f"\nSaving results to {args.output}")
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(args.output, index=False)
        
        print("\nFlux analysis complete!")
        print(f"Analyzed {len(cell_types)} cell types")
        print(f"Total flux predictions: {len(results_df)}")
    else:
        print("\nNo successful flux analyses. Please check your data and model.")
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
