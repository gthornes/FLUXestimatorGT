# Setup Guide

## System Requirements

### Hardware
- RAM: Minimum 16GB, recommended 32GB or more for large datasets
- Storage: At least 50GB free space for data and results
- CPU: Multi-core processor recommended for parallel processing
- GPU: Optional but recommended for faster scFEA training

### Software
- Operating System: Linux, macOS, or Windows
- Python: Version 3.8 or higher
- Git: For version control

## Installation Steps

### 1. Clone the Repository

```bash
git clone https://github.com/gthornes/FLUXestimatorGT.git
cd FLUXestimatorGT
```

### 2. Set Up Python Environment

#### Option A: Using venv (recommended)

```bash
# Create virtual environment
python3 -m venv venv

# Activate virtual environment
# On Linux/macOS:
source venv/bin/activate
# On Windows:
venv\Scripts\activate
```

#### Option B: Using conda

```bash
# Create conda environment
conda create -n fluxestimator python=3.9
conda activate fluxestimator
```

### 3. Install Dependencies

```bash
# Install required packages
pip install -r requirements.txt

# Verify installation
python -c "import scanpy; import pandas; import torch; print('All packages installed successfully!')"
```

### 4. Install Jupyter (if not already installed)

```bash
pip install jupyter notebook
python -m ipykernel install --user --name=fluxestimator
```

### 5. Install scFEA

scFEA (Single-Cell Flux Estimation Analysis) is the core tool for metabolic flux prediction:

```bash
# Clone scFEA repository (as sibling to FLUXestimator)
cd ..
git clone https://github.com/changwn/scFEA.git
cd scFEA

# Install scFEA dependencies
pip install torch torchvision
pip install scikit-learn pandas numpy
pip install python-magic-bin  # Windows only

# Verify scFEA data files exist
ls data/  # Should contain module_gene_m168.csv, cmMat_c70_m168.csv, etc.

# Return to FLUXestimator directory
cd ../FLUXestimatorGT
```

### 6. Download Reference Data

#### scFEA Metabolic Gene List

The scFEA metabolic gene list should already be in `data/reference/scFEA_metabolic_genes.csv`. This contains the 719 curated metabolic genes used by scFEA.

#### Mouse Metabolic Model (Optional for annotations)

```bash
# Create reference directory if it doesn't exist
mkdir -p data/reference

# Optional: Download mouse metabolic model (iMM1415) for pathway annotations
# Available from: http://bigg.ucsd.edu/models/iMM1415
```

### 7. Test Installation

Run a quick test to ensure everything is working:

```bash
# Test scanpy
python -c "import scanpy as sc; print(f'Scanpy version: {sc.__version__}')"

# Test PyTorch (required for scFEA)
python -c "import torch; print(f'PyTorch version: {torch.__version__}')"

# Test data loading
python -c "import anndata; print('AnnData is working!')"

# Check scFEA installation
python -c "import sys; sys.path.append('../scFEA/src'); from scFEA import *; print('scFEA is accessible')"
```

### 8. Launch Jupyter Notebook

```bash
# Start Jupyter notebook server
jupyter notebook

# Navigate to notebooks/ directory in the browser interface
```

## Data Setup

### Preparing Your Data

1. **Place raw data files in `data/raw/`**:
   - 10X Genomics output: matrix.mtx, features.tsv, barcodes.tsv
   - Or AnnData file: your_data.h5ad
   - Or Loom file: your_data.loom

2. **Data format expectations**:
   - Counts should be raw (not normalized)
   - Genes as rows, cells as columns (or in var/obs for AnnData)
   - Include gene symbols or Ensembl IDs

### Example Data Structure

```
data/
├── raw/
│   ├── filtered_feature_bc_matrix_estrus/
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   └── filtered_feature_bc_matrix_diestrus/
│       ├── matrix.mtx.gz
│       ├── features.tsv.gz
│       └── barcodes.tsv.gz
├── processed/
│   ├── annotated_data.h5ad
│   ├── preprocessed_data.h5ad
│   └── expr_mtx_scFEA_genes.csv
└── reference/
    ├── iMM1415.xml
    └── scFEA_metabolic_genes.csv
```

## Configuration

### Edit Analysis Parameters

Edit `config/analysis_config.yaml` to customize analysis parameters:

```yaml
preprocessing:
  min_genes_per_cell: 200      # Adjust based on your data
  max_mito_percent: 10         # Adjust based on tissue type
  n_top_genes: 3000            # Number of highly variable genes
  
clustering:
  resolution: 1.0              # Leiden clustering resolution
  n_neighbors: 15              # Number of neighbors for graph
  
flux_analysis:
  scfea_epochs: 100            # Training epochs for scFEA
  scfea_lr: 0.008              # Learning rate
```

### scFEA Configuration

scFEA uses graph neural networks to estimate metabolic flux. Key parameters:

- **Module genes**: 168 metabolic modules covering major pathways
- **Training**: Default 100 epochs with learning rate 0.008
- **Input**: Log-normalized expression of 719 metabolic genes
- **Output**: Flux predictions for 168 metabolic modules per cell

## Troubleshooting

### Common Issues

#### 1. Memory Errors
If you encounter memory errors with large datasets:
```python
# In your scripts/notebooks, use:
sc.settings.n_jobs = 4  # Reduce parallelization
# Or subset data for testing
adata_subset = adata[:5000, :]  # Use first 5000 cells
```

#### 2. scFEA Execution Issues
If scFEA fails to run:
```bash
# Ensure you're in scFEA directory when running
cd ../scFEA
python src/scFEA.py --data_dir data --input_dir ../FLUXestimator/results/scfea --test_file expression_input.csv

# Check that all data files exist in scFEA/data/
ls data/module_gene_m168.csv
ls data/cmMat_c70_m168.csv

# On Windows, set environment variable if you get library errors:
$env:KMP_DUPLICATE_LIB_OK="TRUE"
```

#### 3. PyTorch/CUDA Issues
If GPU acceleration isn't working:
```python
# Check CUDA availability
import torch
print(f"CUDA available: {torch.cuda.is_available()}")

# Force CPU mode if needed (slower but more compatible)
# scFEA will automatically use CPU if CUDA is unavailable
```

#### 4. Import Errors
If imports fail after installation:
```bash
# Ensure virtual environment is activated
which python  # Should point to venv/bin/python

# Reinstall problematic package
pip install --force-reinstall package-name
```

#### 4. Jupyter Kernel Issues
If kernel doesn't show up:
```bash
python -m ipykernel install --user --name=fluxestimator --display-name="FLUXestimator"
```

### Getting Help

1. Check the [documentation](docs/methodology.md)
2. Review example notebooks in `notebooks/`
3. Open an issue on GitHub
4. Check package documentation:
   - [Scanpy tutorials](https://scanpy-tutorials.readthedocs.io/)
   - [scFEA repository](https://github.com/changwn/scFEA)
   - [PyTorch documentation](https://pytorch.org/docs/)

## Next Steps

After setup is complete:

1. Review the methodology documentation: `docs/methodology.md`
2. Start with the preprocessing notebook: `notebooks/01_preprocessing.ipynb`
3. Follow the cell type annotation: `notebooks/02_cell_type_annotation.ipynb`
4. Run metabolic flux estimation: `notebooks/03_flux_estimation.ipynb`
5. Use annotation script: `scripts/annotate_modules.py` to add pathway descriptions
6. Generate visualizations: `scripts/visualisation.py`

## Updating

To update the repository and dependencies:

```bash
# Update repository
git pull origin main

# Update dependencies
pip install --upgrade -r requirements.txt

# Update scFEA
cd ../scFEA
git pull origin master
cd ../FLUXestimatorGT
```

## Uninstallation

To remove the environment:

```bash
# Deactivate virtual environment
deactivate

# Remove virtual environment directory
rm -rf venv/

# Or for conda:
conda env remove -n fluxestimator
```

## Additional Resources

- [Scanpy Tutorial](https://scanpy-tutorials.readthedocs.io/)
- [scFEA Paper](https://genome.cshlp.org/content/31/10/1867) - Alghamdi et al. (2021)
- [Single-cell Best Practices](https://www.sc-best-practices.org/)
- [PyTorch Tutorials](https://pytorch.org/tutorials/)
- [Graph Neural Networks](https://distill.pub/2021/gnn-intro/)

## Contact

For setup issues, please open an issue on the GitHub repository.
