# FLUXestimatorGT

A repository for a two-month lab project using FLUXestimator to infer metabolic information from single-cell transcriptomic data on reproductive mouse tissue.

## Project Overview

This project aims to analyze metabolic flux in reproductive mouse tissue at the single-cell level by integrating:
- Single-cell RNA sequencing (scRNA-seq) data
- FLUXestimator for metabolic flux inference
- Cell type-specific metabolic profiling

## Features

- **Single-cell data processing**: Tools for quality control, normalization, and cell type annotation
- **Metabolic flux estimation**: Integration with FLUXestimator for flux balance analysis
- **Visualization**: Comprehensive plotting utilities for metabolic pathway analysis
- **Reproducible workflows**: Jupyter notebooks documenting the complete analysis pipeline

## Project Structure

```
FLUXestimatorGT/
├── data/                    # Data directory (not tracked in git)
│   ├── raw/                # Raw scRNA-seq data
│   ├── processed/          # Processed data files
│   └── reference/          # Reference files (metabolic models, annotations)
├── notebooks/              # Jupyter notebooks for analysis
│   ├── 01_preprocessing.ipynb
│   ├── 02_cell_type_annotation.ipynb
│   ├── 03_flux_estimation.ipynb
│   └── 04_visualization.ipynb
├── scripts/                # Python scripts for processing
│   ├── preprocessing.py
│   ├── flux_analysis.py
│   └── visualization.py
├── results/                # Output directory for results
│   ├── figures/
│   └── tables/
├── docs/                   # Documentation
│   ├── methodology.md
│   └── setup.md
├── config/                 # Configuration files
│   └── analysis_config.yaml
├── requirements.txt        # Python dependencies
└── README.md              # This file
```

## Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Setup

1. Clone the repository:
```bash
git clone https://github.com/gthornes/FLUXestimatorGT.git
cd FLUXestimatorGT
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Quick Start

1. Place your scRNA-seq data in the `data/raw/` directory
2. Run the preprocessing notebook: `notebooks/01_preprocessing.ipynb`
3. Annotate cell types: `notebooks/02_cell_type_annotation.ipynb`
4. Estimate metabolic fluxes: `notebooks/03_flux_estimation.ipynb`
5. Visualize results: `notebooks/04_visualization.ipynb`

### Running Scripts

```bash
# Preprocess data
python scripts/preprocessing.py --input data/raw/sample.h5ad --output data/processed/sample_processed.h5ad

# Run flux analysis
python scripts/flux_analysis.py --input data/processed/sample_processed.h5ad --output results/flux_results.csv

# Generate visualizations
python scripts/visualization.py --input results/flux_results.csv --output results/figures/
```

## Data Requirements

This project expects single-cell RNA-seq data in one of the following formats:
- AnnData (.h5ad) - preferred format
- 10X Genomics output (matrix.mtx, genes.tsv, barcodes.tsv)
- CSV/TSV count matrices

## Methodology

The analysis pipeline includes:

1. **Quality Control**: Filtering low-quality cells and genes
2. **Normalization**: Log-normalization and scaling
3. **Cell Type Annotation**: Clustering and marker-based identification
4. **Metabolic Flux Estimation**: Using FLUXestimator with genome-scale metabolic models
5. **Differential Flux Analysis**: Comparing metabolic states across cell types
6. **Visualization**: Pathway enrichment and flux distribution plots

For detailed methodology, see [docs/methodology.md](docs/methodology.md)

## Dependencies

Key packages:
- scanpy: Single-cell analysis
- anndata: Annotated data structures
- numpy/pandas: Data manipulation
- matplotlib/seaborn: Visualization
- cobra: Constraint-based metabolic modeling
- Additional packages listed in requirements.txt

## Citation

If you use this code, please cite:
- FLUXestimator: [Original FLUXestimator publication]
- Scanpy: Wolf, F. A. et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology.

## License

This project is licensed for academic and research use.

## Contact

For questions or issues, please open an issue on GitHub or contact the repository owner.

## Acknowledgments

This project was developed as part of a two-month laboratory project focusing on reproductive biology and metabolic analysis.