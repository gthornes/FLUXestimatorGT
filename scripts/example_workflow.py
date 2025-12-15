#!/usr/bin/env python3
"""
Example workflow script demonstrating the complete analysis pipeline.

This script provides a template for running the complete analysis from
raw data to flux estimation.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

import yaml


def print_step(step_num, title):
    """Print a formatted step header."""
    print("\n" + "="*60)
    print(f"Step {step_num}: {title}")
    print("="*60 + "\n")


def main():
    """Run example workflow."""
    
    print("""
╔═══════════════════════════════════════════════════════════╗
║   FLUXestimator Single-Cell Metabolic Analysis Pipeline   ║
║                    Example Workflow                       ║
╚═══════════════════════════════════════════════════════════╝
    """)
    
    print("This script demonstrates the complete analysis pipeline.")
    print("For actual analysis, you would run each step with your data.\n")
    
    # Load configuration
    config_path = Path(__file__).parent.parent / 'config' / 'analysis_config.yaml'
    
    if config_path.exists():
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        print(f"✓ Configuration loaded from {config_path}")
    else:
        print(f"✗ Configuration file not found: {config_path}")
        return 1
    
    # Pipeline steps
    print_step(1, "Data Preprocessing")
    print("Command:")
    print("  python scripts/preprocessing.py \\")
    print("    --input data/raw/your_data.h5ad \\")
    print("    --output data/processed/preprocessed_data.h5ad")
    print("\nThis step:")
    print("  - Loads raw scRNA-seq data")
    print("  - Performs quality control")
    print("  - Filters cells and genes")
    print("  - Normalizes and scales data")
    print("  - Performs PCA")
    
    print_step(2, "Cell Type Annotation")
    print("Use Jupyter notebook:")
    print("  jupyter notebook notebooks/02_cell_type_annotation.ipynb")
    print("\nThis step:")
    print("  - Computes neighborhood graph")
    print("  - Performs clustering (Leiden)")
    print("  - Creates UMAP visualization")
    print("  - Identifies marker genes")
    print("  - Annotates cell types")
    
    print_step(3, "Metabolic Flux Analysis")
    print("Command:")
    print("  python scripts/flux_analysis.py \\")
    print("    --input data/processed/annotated_data.h5ad \\")
    print("    --output results/tables/flux_results.csv")
    print("\nThis step:")
    print("  - Loads metabolic model")
    print("  - Maps gene expression to reactions")
    print("  - Performs flux balance analysis (FBA)")
    print("  - Generates cell type-specific flux predictions")
    
    print_step(4, "Visualization")
    print("Command:")
    print("  python scripts/visualization.py \\")
    print("    --input results/tables/flux_results.csv \\")
    print("    --output results/figures/")
    print("\nThis step:")
    print("  - Creates flux heatmaps")
    print("  - Generates distribution plots")
    print("  - Produces pathway comparison figures")
    print("  - Summarizes metabolic activity")
    
    print("\n" + "="*60)
    print("Pipeline Overview Complete")
    print("="*60)
    
    print("\n📚 Documentation:")
    print("  - Setup guide: docs/setup.md")
    print("  - Methodology: docs/methodology.md")
    print("  - Notebooks: notebooks/")
    
    print("\n🔧 Configuration:")
    print("  Edit parameters in: config/analysis_config.yaml")
    
    print("\n📂 Directory Structure:")
    print("  data/raw/       - Place your raw data here")
    print("  data/processed/ - Processed data files")
    print("  results/        - Output figures and tables")
    
    print("\n✨ Happy analyzing!\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
