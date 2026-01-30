# Methodology

## Overview

This document describes the computational methodology for inferring metabolic flux from single-cell transcriptomic data using FLUXestimator with scFEA (Single-Cell Flux Estimation Analysis). The pipeline integrates standard single-cell RNA-seq preprocessing with graph neural network-based flux prediction.

## Workflow

### 1. Data Acquisition and Preprocessing

#### 1.1 Data Loading
- Load raw single-cell RNA-seq count matrices
- Support for multiple input formats (10X, h5ad, loom)
- Read cell and gene metadata
- Concatenate multiple samples (e.g., estrus and diestrus stages)

#### 1.2 Quality Control
**Cell-level filtering:**
- Remove cells with too few genes (< 200 genes)
- Remove cells with too many genes (> 6000 genes, potential doublets)
- Filter cells with high mitochondrial content (> 10%)
- Remove potential doublets using Scrublet

**Gene-level filtering:**
- Remove genes expressed in fewer than 3 cells
- Filter mitochondrial genes (mt-) for downstream analysis
- Retain ribosomal genes (they provide biological information)

#### 1.3 Normalization
- Library size normalization (counts per 10,000)
- Log transformation: log1p(normalized_counts)
- Store in layers: 'log1p_norm' for flux analysis, 'counts' for raw data
- Identify highly variable genes (HVGs) for clustering (top 3000 genes)
- Scale HVGs to unit variance and zero mean for PCA

### 2. Dimensionality Reduction and Clustering

#### 2.1 PCA
- Perform PCA on highly variable genes
- Use top 50 principal components for downstream analysis
- Visualize variance explained

#### 2.2 Batch Correction (if multiple samples)
- Apply Harmony for batch correction when combining samples
- Integrates datasets while preserving biological variation
- Particularly important for estrus/diestrus comparisons

#### 2.3 Neighbor Graph Construction
- Construct k-nearest neighbor graph (k=15)
- Use Euclidean distance in PCA/Harmony space

#### 2.4 Clustering
- Apply Leiden algorithm for community detection
- Default resolution: 1.0-1.2 (adjust based on dataset)
- Produces clusters ranging from broad cell types (low res) to subtypes (high res)

#### 2.5 UMAP Visualization
- Compute 2D UMAP embedding for visualization
- Parameters: min_dist=0.5, n_neighbors=15

### 3. Cell Type Annotation

#### 3.1 Marker Gene Identification
- Perform differential expression analysis (Wilcoxon rank-sum test)
- Identify top marker genes per cluster
- Filter by log-fold change (> 0.5) and adjusted p-value (< 0.05)

#### 3.2 Manual Annotation
- Review marker gene expression patterns
- Assign cell type labels based on known markers
- Validate with literature for reproductive/ovarian tissue

#### 3.3 Subclustering (Optional)
- Subset specific cell type clusters
- Re-run preprocessing at higher resolution (0.4-0.5)
- Identify cell subtypes (e.g., vascular vs lymphatic endothelial cells)
- Merge refined annotations back into main dataset as 'refined_cell_type'

#### 3.4 Marker Visualization
- Create overlay plots showing marker gene expression on UMAP
- Generate dynamic grids displaying multiple markers simultaneously
- Facilitates identification of cluster identity

### 4. Metabolic Flux Estimation with scFEA

#### 4.1 Expression Matrix Preparation
- Extract expression of 719 scFEA metabolic genes from dataset
- Use log-normalized expression (log1p_norm layer)
- Export as genes (rows) × cells (columns) CSV format
- Save cell metadata including cell_type, refined_cell_type, and stage

#### 4.2 scFEA Architecture
scFEA uses a graph neural network to estimate cell-specific metabolic flux:

**Key Features:**
- **Input**: Expression of 719 metabolic genes per cell
- **Model**: Graph convolutional network trained on metabolic network structure
- **Modules**: 168 metabolic modules/reactions covering major pathways
- **Training**: 100 epochs (default), learning rate 0.008
- **Output**: Flux predictions for 168 modules per cell

**Advantages over traditional FBA:**
- No need for optimization at single-cell level
- Handles sparse scRNA-seq data effectively
- Learns from metabolic network topology
- Accounts for enzyme stoichiometry and constraints
- Faster than iterative constraint-based methods

#### 4.3 Running scFEA
```bash
# Execute from scFEA directory
python src/scFEA.py \
    --data_dir data \
    --input_dir ../FLUXestimator/results/scfea/[suffix] \
    --res_dir ../FLUXestimator/results/scfea/[suffix] \
    --test_file expression_input.csv \
    --moduleGene_file module_gene_complete_mouse_m168.csv
```

**Output files:**
- `expression_input_module168_cell[N]_*.csv` - Flux predictions per cell
- Time stamps and training logs

#### 4.4 Module Annotation
Use `annotate_modules.py` to add human-readable pathway descriptions:

```bash
python scripts/annotate_modules.py \
    --flux results/scfea/[suffix]/[flux_file].csv \
    -o results/tables/[suffix] \
    --create-annotated-file
```

This replaces module IDs (M_001, M_002) with pathway descriptions like "Acetyl-CoA → (E,E)-Farnesyl-PP" for interpretability.

### 5. Differential Flux Analysis

#### 5.1 Cell Type-Specific Metabolism
- Calculate mean/median flux per metabolic module for each cell type
- Compare subclustered cell types (e.g., different immune subtypes)
- Identify cell type-specific metabolic signatures

#### 5.2 Stage Comparison (Estrus vs Diestrus)
- Calculate mean flux per module for each stage
- Perform t-tests for each module across stages
- Apply Benjamini-Hochberg FDR correction for multiple testing
- Calculate log2 fold changes

#### 5.3 Cell Type-Stage Interactions
- For each cell type, compare flux between stages
- Requires minimum 20 cells per cell type per stage
- Identifies cell type-specific metabolic responses to hormonal stage

#### 5.4 Pathway-Level Analysis
- Aggregate related modules into metabolic pathways
- Calculate pathway activity scores
- Visualize differential pathway activity

### 6. Visualization and Interpretation

#### 6.1 Flux Distribution Plots
- **Violin plots**: Show distribution of high-activity modules across cell types
- **Heatmaps**: Display top variable modules across cell types
- **Paired stage plots**: Compare estrus vs diestrus side-by-side for each cell type
- Color schemes: Mean/median flux intensity, stage-specific coloring

#### 6.2 Differential Module Analysis
- **Volcano plots**: Log2 fold change vs -log10(p-value)
- **Bar charts**: Top differentially active modules
- **Heatmaps**: Estrus vs diestrus for top modules

#### 6.3 Cell Type Metabolic Profiles
- Bar plots of overall metabolic activity per cell type
- Pathway breadth (number of active reactions)
- Peak metabolic flux

#### 6.4 Network Visualization
- Pathway comparison plots showing relative activity
- Flux vs gene count correlations
- Module-level heatmaps with hierarchical clustering

### 7. Validation and Quality Checks

#### 7.1 Technical Validation
- Check that scFEA training converged (loss decreasing)
- Verify cell type labels are preserved through pipeline
- Ensure subclustering doesn't create artifacts
- Validate that stage labels match experimental conditions

#### 7.2 Biological Validation
- Compare results with known metabolic phenotypes
- Validate against literature for reproductive cells
- Check that estrus/diestrus differences align with hormonal changes
- Cross-reference with proteomic/metabolomic data if available

#### 7.3 Statistical Validation
- Ensure sufficient cell numbers per comparison (n≥20)
- Check for batch effects between samples
- Verify multiple testing correction is applied
- Assess clustering stability across parameter ranges

## Tools and Packages

### Primary Tools
- **scanpy**: Single-cell data processing, clustering, and visualization
- **scFEA**: Graph neural network-based metabolic flux estimation
- **PyTorch**: Deep learning framework for scFEA

### Supporting Tools
- **anndata**: Data structure for single-cell data
- **pandas/numpy**: Data manipulation and analysis
- **matplotlib/seaborn**: Visualization
- **scipy/statsmodels**: Statistical analysis
- **harmony-pytorch**: Batch correction for multi-sample datasets

## Pipeline Structure

### Notebooks
1. **01_preprocessing.ipynb**: QC, normalization, clustering
2. **02_cell_type_annotation.ipynb**: Cell type identification, subclustering
3. **03_flux_estimation.ipynb**: scFEA execution, differential analysis, visualization

### Scripts
- **annotate_modules.py**: Add pathway descriptions to flux results
- **visualisation.py**: Generate publication-quality figures
- **preprocessing.py**: Reusable preprocessing functions

## Best Practices

1. **Reproducibility**: 
   - Set random seeds for all stochastic operations
   - Document scFEA version and parameters
   - Save intermediate results (preprocessed data, annotated data)

2. **Documentation**: 
   - Keep detailed logs of parameter choices
   - Note any manual annotations or decisions
   - Track software versions in requirements.txt

3. **Version Control**: 
   - Commit changes regularly
   - Use meaningful commit messages
   - Tag releases with date suffixes

4. **Validation**: 
   - Always validate results against biological knowledge
   - Check for batch effects between samples
   - Compare with published datasets when possible

5. **Visualization**: 
   - Create comprehensive figures for all major results
   - Use consistent color schemes across plots
   - Include statistical annotations (p-values, fold changes)

6. **Subclustering**:
   - Only subcluster cell types with sufficient cells (>100)
   - Use higher resolution (0.4-0.5) for subclustering
   - Validate subtypes with marker genes
   - Preserve original annotations while adding refined ones

## Key Metabolic Pathways Analyzed

Focus on pathways relevant to reproductive biology:

### Energy Metabolism
- Glycolysis and gluconeogenesis
- TCA cycle (citric acid cycle)
- Oxidative phosphorylation (OXPHOS)
- Pentose phosphate pathway

### Lipid Metabolism
- Fatty acid synthesis and oxidation
- Steroid hormone biosynthesis
- Cholesterol metabolism
- Prostaglandin synthesis

### Amino Acid Metabolism
- Amino acid biosynthesis and degradation
- One-carbon metabolism (folate/methionine cycles)
- Glutamine and glutamate metabolism

### Nucleotide Metabolism
- Purine and pyrimidine synthesis
- DNA/RNA precursor production

### Specialized Pathways
- Reactive oxygen species (ROS) metabolism
- NAD+ biosynthesis
- Cofactor metabolism

## Common Pitfalls

1. **Over-filtering**: Too stringent QC can remove rare cell types
2. **Under-filtering**: Doublets and low-quality cells affect clustering
3. **Normalization**: Incorrect normalization affects flux predictions
4. **Resolution**: Clustering resolution affects cell type granularity
5. **Batch effects**: Failure to correct for batch can confound stage comparisons
6. **Expression thresholds**: Low gene coverage may miss metabolic activity
7. **Sample size**: Insufficient cells per comparison reduces statistical power
8. **Module interpretation**: Need pathway annotations for biological interpretation
9. **Subclustering artifacts**: Too aggressive subclustering can split true cell types
10. **Stage labels**: Ensure experimental metadata is correctly assigned

## Differences from Traditional FBA

**scFEA Advantages:**
- No need for cell-specific optimization
- Handles sparse single-cell data natively
- Learns from metabolic network topology
- Predicts flux for individual cells
- Much faster than iterative constraint-based methods

**Traditional FBA Limitations:**
- Requires extensive constraints and assumptions
- Difficult to apply at single-cell level
- Computationally expensive for thousands of cells
- May not handle sparse expression well
- Needs carefully chosen objective functions

## References

1. **scFEA**: Alghamdi, N. et al. (2021) "A graph neural network model to estimate cell-wise metabolic flux using single-cell RNA-seq data" *Genome Research* 31(10):1867-1884
2. **Scanpy**: Wolf, F.A. et al. (2018) "SCANPY: large-scale single-cell gene expression data analysis" *Genome Biology* 19:15
3. **Leiden clustering**: Traag, V.A. et al. (2019) "From Louvain to Leiden: guaranteeing well-connected communities" *Scientific Reports* 9:5233
4. **Harmony**: Korsunsky, I. et al. (2019) "Fast, sensitive and accurate integration of single-cell data with Harmony" *Nature Methods* 16:1289-1296
5. **Graph Neural Networks**: Kipf, T.N. & Welling, M. (2017) "Semi-Supervised Classification with Graph Convolutional Networks" *ICLR*

## Contact

For methodology questions, please open an issue on the GitHub repository.
