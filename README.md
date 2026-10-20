# CRC-SingleCell-Spatial-Analysis
Colorectal Cancer: Single-Cell and Spatial Host–Bacteria Analysis
An R workflow integrating single-cell RNA sequencing and spatial transcriptomics to characterize host cell composition, bacterial signals, and spatial host–bacteria associations in colorectal cancer (CRC).
Data
- GSE132465: single-cell RNA-seq reference; 10 matched normal–tumor patient pairs are retained for analysis.
- GSE334323: spatial transcriptomics; four sections from two matched normal–tumor pairs (CRC_03 and CRC_08).
Datasets are analyzed separately, with the single-cell reference supporting spatial cell-type deconvolution using RCTD.
Workflow
Stage	Script*	Purpose
01	01_Data_Preprocessing.R	Import, metadata validation, and quality control
02	02_Cell_Annotation.R	Normalization, clustering, marker detection, and cell-type annotation
03	03_Cell_Composition.R	Paired comparisons of cell-type and subtype proportions
04	04_Spatial_Host_Profiling.R	Spatial preprocessing, host/bacterial signal separation, RCTD, and feature mapping
05	05_Bacteria_Host_Association.R	Spatial associations between bacterial signals and selected host features
06	06_Paired_Bacterial_Changes.R	Descriptive paired bacterial changes across 19 deposited bacterial groups
07	07_Bacterial_Spatial_Cooccurrence.R	Exploratory bacterial spatial co-occurrence and descriptive adjusted analyses


Uploaded filenames may include version suffixes such as (3) or (4); use the actual filenames when running scripts.
Requirements and Usage
Use R with the packages required by each script. Dependencies include Seurat, SeuratObject, Matrix, data.table, R.utils, ggplot2, patchwork, SingleR, celldex, SingleCellExperiment, BiocParallel, spacexr, RANN, SpatialPack, cowplot, and digest. RCTD requires the dmcable/spacexr API.
1. Download the processed count matrix and cell annotations for GSE132465, and GSE334323_RAW.tar for the spatial analysis.
2. Update the hard-coded input and output paths in each script to match your environment; current defaults use E:/CRC/....
3. Run stages 01–07 in order using source("<actual_script_filename>.R"). Keep upstream output directories consistent with downstream input paths.
Outputs
Stage-specific directories contain analysis tables, figures, and intermediate R objects. Spatial stages also use validation records, checksums, and caches to verify upstream inputs and support reproducibility.
Interpretation
Spatial analyses are exploratory and include only two patient pairs. Bacterial measurements represent deposited probe-derived signals. Spatially adjusted rank-test P values are approximate; depth-adjusted coefficients and normal–tumor contrasts are descriptive. Associations do not establish causality or direct microbial interactions, and broader conclusions require independent validation.
