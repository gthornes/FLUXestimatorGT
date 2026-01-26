#!/usr/bin/env python3
"""
Annotate scFEA metabolic modules with pathway information.

This script creates a mapping between module IDs (M_1, M_2, etc.) 
and their biological pathway functions based on the genes in each module.
"""

import pandas as pd
import sys
from pathlib import Path


def load_module_annotations_from_scfea(scfea_path):
    """
    Load module annotations from scFEA data directory.
    
    Parameters
    ----------
    scfea_path : Path
        Path to scFEA directory
        
    Returns
    -------
    pd.DataFrame
        Module annotations with compound information
    """
    info_file = scfea_path / 'data' / 'Human_M168_information.symbols.csv'
    
    df = pd.read_csv(info_file)
    
    # Create readable descriptions
    df['description'] = df.apply(
        lambda row: f"{row['Compound_IN_name']} → {row['Compound_OUT_name']}", 
        axis=1
    )
    
    # Map supermodule IDs to pathway names (based on manual inspection)
    supermodule_map = {
        1: 'Central carbon metabolism',
        2: 'Amino acid metabolism', 
        3: 'Lipid metabolism',
        4: 'TCA cycle extensions',
        5: 'Nucleotide metabolism',
        6: 'Cofactor biosynthesis',
        7: 'Secondary metabolism'
    }
    
    df['pathway'] = df['Supermodule_id'].map(supermodule_map).fillna('Other metabolism')
    
    # Set module_id as index
    df.index = df['Unnamed: 0']
    
    return df[['description', 'pathway', 'Compound_IN_name', 'Compound_OUT_name', 'Supermodule_id']]

def load_module_genes(scfea_path):
    """
    Load the gene composition of each module from scFEA data.
    
    Parameters
    ----------
    scfea_path : Path
        Path to scFEA directory
    
    Returns
    -------
    pd.DataFrame
        Module genes
    """
    module_file = scfea_path / 'data' / 'module_gene_complete_mouse_m168.csv'
    
    if not module_file.exists():
        print(f"Warning: Module file not found at {module_file}")
        return None
    
    df = pd.read_csv(module_file)
    
    # Extract genes for each module
    module_info = []
    for idx, row in df.iterrows():
        module_id = row.iloc[0]
        genes = [g for g in row.iloc[1:] if pd.notna(g)]
        module_info.append({
            'module_id': module_id,
            'n_genes': len(genes),
            'genes': ', '.join(genes[:5]) + ('...' if len(genes) > 5 else '')
        })
    
    return pd.DataFrame(module_info)


def annotate_flux_results(flux_file, output_file=None, scfea_path=None):
    """
    Add pathway annotations to flux results.
    
    Parameters
    ----------
    flux_file : str or Path
        Path to flux results CSV
    output_file : str or Path, optional
        Output path for annotated results
    scfea_path : str or Path, optional
        Path to scFEA directory
    """
    # Load annotations from scFEA if available
    annotations = load_module_annotations_from_scfea(Path(scfea_path))
    
    # Load flux results
    flux_df = pd.read_csv(flux_file, index_col=0)
    
    # Get module columns
    module_cols = [col for col in flux_df.columns if col.startswith('M_')]
    
    print(f"Flux results: {flux_df.shape[0]} cells × {len(module_cols)} modules")
    print(f"Annotations available for {len(annotations)} modules")
    
    # Calculate summary statistics per module
    module_stats = []
    for module in module_cols:
        if module in flux_df.columns:
            values = flux_df[module].dropna()
            stats = {
                'module_id': module,
                'mean_flux': values.abs().mean(),
                'median_flux': values.abs().median(),
                'max_flux': values.abs().max(),
                'std_flux': values.abs().std(),
                'n_cells': len(values)
            }
            
            # Add annotation if available
            if module in annotations.index:
                stats['pathway'] = annotations.loc[module, 'pathway']
                stats['description'] = annotations.loc[module, 'description']
                if 'Compound_IN_name' in annotations.columns:
                    stats['substrate'] = annotations.loc[module, 'Compound_IN_name']
                    stats['product'] = annotations.loc[module, 'Compound_OUT_name']
            else:
                stats['pathway'] = 'Unknown'
                stats['description'] = 'Not yet annotated'
                stats['substrate'] = 'N/A'
                stats['product'] = 'N/A'
            
            module_stats.append(stats)
    
    result_df = pd.DataFrame(module_stats)
    result_df = result_df.sort_values('mean_flux', ascending=False)
    
    # Display top modules
    print("\nTop 20 most active metabolic modules:")
    display_cols = ['module_id', 'mean_flux', 'pathway', 'description']
    print(result_df[display_cols].head(20).to_string(index=False))
    
    if output_file:
        # Handle directory vs file path
        output_path = Path(output_file)
        if output_path.is_dir() or (not output_path.suffix and not output_path.exists()):
            # It's a directory - create a filename
            output_path.mkdir(parents=True, exist_ok=True)
            output_path = output_path / 'module_annotations.csv'
        else:
            # Ensure parent directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)
        
        result_df.to_csv(output_path, index=False)
        print(f"\nFull annotated results saved to {output_path}")
    
    return result_df


def create_annotated_flux_file(flux_file, annotations, output_file=None):
    """
    Create a version of the flux file with annotated column names.
    
    Parameters
    ----------
    flux_file : str or Path
        Path to original flux results CSV
    annotations : pd.DataFrame
        Module annotations (indexed by module_id)
    output_file : str or Path, optional
        Output path for annotated file. If None, creates one based on input name.
    
    Returns
    -------
    pd.DataFrame
        Flux dataframe with annotated column names
    """
    print(f"\nCreating annotated flux file with readable column names...")
    
    # Load the original flux data
    flux_df = pd.read_csv(flux_file, index_col=0)
    
    # Get module columns
    module_cols = [col for col in flux_df.columns if col.startswith('M_')]
    
    # Create column name mapping
    col_mapping = {}
    annotated_count = 0
    for module_id in module_cols:
        if module_id in annotations.index:
            pathway = annotations.loc[module_id, 'pathway']
            desc = annotations.loc[module_id, 'description']
            # Create readable column name: Module_ID | Description | [Pathway]
            new_name = f"{module_id} | {desc} | [{pathway}]"
            col_mapping[module_id] = new_name
            annotated_count += 1
        else:
            col_mapping[module_id] = module_id
    
    # Rename columns
    flux_df_annotated = flux_df.rename(columns=col_mapping)
    
    # Determine output path
    if output_file is None:
        flux_path = Path(flux_file)
        output_file = flux_path.parent / f"{flux_path.stem}_annotated.csv"
    else:
        output_file = Path(output_file)
        # Handle directory vs file path
        if output_file.is_dir() or (not output_file.suffix and not output_file.exists()):
            # It's a directory - create a filename based on input
            flux_path = Path(flux_file)
            output_file = output_file / f"{flux_path.stem}_annotated.csv"
    
    # Ensure parent directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Save to new file
    flux_df_annotated.to_csv(output_file)
    
    print(f"Saved annotated flux file to {output_file}")
    print(f"  - Total modules: {len(module_cols)}")
    print(f"  - Annotated: {annotated_count}")
    print(f"  - Not annotated: {len(module_cols) - annotated_count}")
    
    return flux_df_annotated


def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Annotate scFEA metabolic modules with pathway information'
    )
    parser.add_argument(
        '--flux', '-f',
        required=True,
        help='Flux results CSV file'
    )
    parser.add_argument(
        '--output', '-o',
        help='Output CSV file for annotated modules'
    )
    parser.add_argument(
        '--scfea-path',
        default='../scFEA',
        help='Path to scFEA directory (default: ../scFEA)'
    )
    parser.add_argument(
        '--create-annotated-file',
        action='store_true',
        help='Create a flux file with annotated column names'
    )
    parser.add_argument(
        '--annotated-output',
        help='Output path for annotated flux file (default: adds _annotated suffix)'
    )
    
    args = parser.parse_args()
    
    # Load annotations
    annotations = load_module_annotations_from_scfea(
        Path(args.scfea_path) if args.scfea_path else None
    )
    
    # Create summary table
    annotate_flux_results(
        flux_file=args.flux,
        output_file=args.output,
        scfea_path=Path(args.scfea_path) if args.scfea_path else None
    )
    
    # Optionally create annotated flux file
    if args.create_annotated_file:
        create_annotated_flux_file(
            flux_file=args.flux,
            annotations=annotations,
            output_file=args.annotated_output
        )


if __name__ == '__main__':
    sys.exit(main())
