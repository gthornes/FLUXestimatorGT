# Setup Guide

## System Requirements

### Hardware
- RAM: Minimum 8GB, recommended 16GB or more for large datasets
- Storage: At least 50GB free space for data and results
- CPU: Multi-core processor recommended for parallel processing

### Software
- Operating System: Linux, macOS, or Windows with WSL2
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
python -c "import scanpy; import cobra; import pandas; print('All packages installed successfully!')"
```

### 4. Install Jupyter (if not already installed)

```bash
pip install jupyter notebook
python -m ipykernel install --user --name=fluxestimator
```

### 5. Download Reference Data

#### Metabolic Models

Download genome-scale metabolic models:

```bash
# Create reference directory if it doesn't exist
mkdir -p data/reference

# Example: Download mouse metabolic model (iMM1415)
# You can manually download from:
# http://bigg.ucsd.edu/models/iMM1415
# Or use wget/curl if available:
# wget -O data/reference/iMM1415.xml http://bigg.ucsd.edu/models/iMM1415/download
```

#### Gene Annotations

Download mouse gene annotations (if needed):

```bash
# Example: Download from Ensembl
# wget -O data/reference/Mus_musculus.GRCm39.gtf.gz \
#   ftp://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz
# gunzip data/reference/Mus_musculus.GRCm39.gtf.gz
```

### 6. Test Installation

Run a quick test to ensure everything is working:

```bash
# Test scanpy
python -c "import scanpy as sc; print(f'Scanpy version: {sc.__version__}')"

# Test cobra
python -c "import cobra; print(f'COBRApy version: {cobra.__version__}')"

# Test data loading
python -c "import anndata; print('AnnData is working!')"
```

### 7. Launch Jupyter Notebook

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
│   ├── sample1/
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   └── sample2.h5ad
├── processed/
│   └── (will be created during analysis)
└── reference/
    ├── iMM1415.xml
    └── gene_annotations.gtf
```

## Configuration

### Edit Analysis Parameters

Edit `config/analysis_config.yaml` to customize analysis parameters:

```yaml
preprocessing:
  min_genes_per_cell: 200  # Adjust based on your data
  max_mito_percent: 10     # Adjust based on tissue type
  
flux_analysis:
  model: "recon3d"         # Choose your metabolic model
  solver: "glpk"           # Change if you have commercial solver
```

### Solvers for Optimization

For better performance, consider installing commercial solvers:

#### CPLEX (Free for academics)
```bash
# After obtaining CPLEX from IBM Academic Initiative
pip install cplex
```

#### Gurobi (Free for academics)
```bash
# After obtaining Gurobi license
pip install gurobipy
```

## Troubleshooting

### Common Issues

#### 1. Memory Errors
If you encounter memory errors with large datasets:
```python
# In your scripts/notebooks, use:
sc.settings.n_jobs = 4  # Reduce parallelization
# Or process data in batches
```

#### 2. Solver Issues
If optimization fails:
```python
# Switch solver in config or code:
model.solver = 'glpk'  # or 'cplex', 'gurobi'
```

#### 3. Import Errors
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
   - [COBRApy documentation](https://cobrapy.readthedocs.io/)

## Next Steps

After setup is complete:

1. Review the methodology documentation: `docs/methodology.md`
2. Start with the preprocessing notebook: `notebooks/01_preprocessing.ipynb`
3. Follow the analysis pipeline through the numbered notebooks
4. Customize scripts in `scripts/` for your specific needs

## Updating

To update the repository and dependencies:

```bash
# Update repository
git pull origin main

# Update dependencies
pip install --upgrade -r requirements.txt
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
- [COBRApy Documentation](https://cobrapy.readthedocs.io/)
- [Single-cell Best Practices](https://www.sc-best-practices.org/)
- [FBA Tutorial](https://opencobra.github.io/cobrapy/getting_started.html)

## Contact

For setup issues, please open an issue on the GitHub repository.
