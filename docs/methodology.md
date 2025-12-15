# Methodology

## Overview

This document describes the computational methodology for inferring metabolic information from single-cell transcriptomic data on reproductive mouse tissue using FLUXestimator.

## Workflow

### 1. Data Acquisition and Preprocessing

#### 1.1 Data Loading
- Load raw single-cell RNA-seq count matrices
- Support for multiple input formats (10X, h5ad, loom)
- Read cell and gene metadata

#### 1.2 Quality Control
**Cell-level filtering:**
- Remove cells with too few genes (< 200 genes)
- Remove cells with too many genes (> 5000 genes, potential doublets)
- Filter cells with high mitochondrial content (> 10%)
- Remove potential doublets using Scrublet

**Gene-level filtering:**
- Remove genes expressed in fewer than 3 cells
- Filter mitochondrial and ribosomal genes if needed

#### 1.3 Normalization
- Library size normalization (counts per 10,000)
- Log transformation: log1p(normalized_counts)
- Scale data to unit variance and zero mean
- Identify highly variable genes (HVGs) for downstream analysis

### 2. Dimensionality Reduction and Clustering

#### 2.1 PCA
- Perform PCA on highly variable genes
- Use top 50 principal components for clustering
- Visualize variance explained

#### 2.2 Neighbor Graph Construction
- Construct k-nearest neighbor graph (k=15)
- Use Euclidean distance in PCA space

#### 2.3 Clustering
- Apply Leiden algorithm for community detection
- Test multiple resolutions (0.4, 0.6, 0.8, 1.0)
- Validate clusters using silhouette scores

#### 2.4 UMAP Visualization
- Compute 2D UMAP embedding
- Parameters: min_dist=0.5, spread=1.0, n_neighbors=15

### 3. Cell Type Annotation

#### 3.1 Marker Gene Identification
- Perform differential expression analysis (t-test/Wilcoxon)
- Identify top marker genes per cluster
- Filter by log-fold change (> 0.5) and adjusted p-value (< 0.05)

#### 3.2 Automated Annotation
- Match marker genes to known cell type signatures
- Use reference datasets (e.g., Mouse Cell Atlas)
- Score each cluster against reference cell types

#### 3.3 Manual Curation
- Review automated annotations
- Validate with known markers for reproductive tissue:
  - Sertoli cells: Sox9, Amh, Dhh
  - Leydig cells: Cyp17a1, Hsd3b1, Star
  - Spermatogonia: Utf1, Zbtb16, Gfra1
  - Spermatocytes: Sycp3, Tex14, Prdm9
  - Spermatids: Tnp1, Prm1, Prm2

### 4. Metabolic Flux Estimation

#### 4.1 Model Selection
- Use genome-scale metabolic model (e.g., Recon3D for mouse)
- Load model using COBRApy
- Define biomass objective function

#### 4.2 Gene Expression Integration
- Map gene expression to metabolic reactions
- Apply Gene-Protein-Reaction (GPR) rules
- Convert expression levels to flux bounds

**Mapping strategy:**
- For AND relationships: minimum expression
- For OR relationships: maximum expression
- Threshold low expression values

#### 4.3 Flux Balance Analysis (FBA)
- For each cell type:
  1. Extract mean gene expression profile
  2. Constrain reaction bounds based on expression
  3. Solve FBA optimization problem
  4. Extract flux distribution

**Optimization:**
- Maximize biomass production or other objective
- Subject to stoichiometric constraints
- Apply thermodynamic feasibility constraints

#### 4.4 Cell-Specific Flux Estimation
- Option 1: Aggregate by cell type (faster)
- Option 2: Single-cell level (computationally intensive)
- Store flux predictions for all reactions

### 5. Differential Flux Analysis

#### 5.1 Statistical Testing
- Compare flux distributions between cell types
- Use Wilcoxon rank-sum test or t-test
- Correct for multiple testing (Benjamini-Hochberg)

#### 5.2 Pathway-Level Analysis
- Aggregate fluxes by metabolic pathway
- Calculate pathway activity scores
- Identify differentially active pathways

#### 5.3 Key Metabolic Pathways
Focus on pathways relevant to reproductive biology:
- Energy metabolism (glycolysis, TCA, OXPHOS)
- Steroid hormone biosynthesis
- Lipid metabolism
- Amino acid metabolism
- One-carbon metabolism (folate cycle)

### 6. Visualization and Interpretation

#### 6.1 Cell Type-Specific Metabolism
- UMAP plots colored by pathway activity
- Heatmaps of flux distributions
- Violin plots for key reactions

#### 6.2 Pathway Analysis
- Bar plots of pathway enrichment
- Network visualization of metabolic connections
- Correlation analysis between pathways

#### 6.3 Metabolite-Reaction Networks
- Visualize flux through specific pathways
- Highlight differentially active reactions
- Map to biological processes

### 7. Validation and Quality Checks

#### 7.1 Technical Validation
- Check flux consistency with mass balance
- Verify no negative fluxes where inappropriate
- Ensure thermodynamic feasibility

#### 7.2 Biological Validation
- Compare with known metabolic phenotypes
- Validate against literature for reproductive cells
- Cross-reference with proteomic/metabolomic data if available

## Tools and Packages

### Primary Tools
- **scanpy**: Single-cell data processing and analysis
- **COBRApy**: Constraint-based metabolic modeling
- **FLUXestimator**: Flux estimation from expression data

### Supporting Tools
- **anndata**: Data structure for single-cell data
- **pandas/numpy**: Data manipulation
- **matplotlib/seaborn**: Visualization
- **scipy/statsmodels**: Statistical analysis

## Best Practices

1. **Reproducibility**: Set random seeds for all stochastic operations
2. **Documentation**: Keep detailed logs of parameter choices
3. **Version Control**: Track software versions and dependencies
4. **Validation**: Always validate results against biological knowledge
5. **Visualization**: Create comprehensive figures for all major results

## Common Pitfalls

1. **Over-filtering**: Too stringent QC can remove rare cell types
2. **Normalization**: Choice of method affects downstream results
3. **Resolution**: Clustering resolution affects cell type granularity
4. **Model choice**: Different metabolic models give different results
5. **Expression thresholds**: Affects which reactions are active

## References

1. Wolf, F.A. et al. (2018) SCANPY: large-scale single-cell gene expression data analysis. Genome Biology.
2. Ebrahim, A. et al. (2013) COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Systems Biology.
3. Orth, J.D. et al. (2010) What is flux balance analysis? Nature Biotechnology.
4. Traag, V.A. et al. (2019) From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports.

## Contact

For methodology questions, please open an issue on the GitHub repository.
