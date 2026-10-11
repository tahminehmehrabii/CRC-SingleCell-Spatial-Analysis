# Single-Cell and Spatial Transcriptomics of Colorectal Cancer

**Cell Communication · Predictive Modeling · Bacteria–Host Associations**

This repository contains the main R analysis scripts from my thesis project on colorectal cancer (CRC). The workflow combines single-cell RNA sequencing, spatial transcriptomics, and predictive modeling to characterize cellular composition, infer cell–cell communication, and explore associations between deposited bacterial probe signals and host features.

## Research Objectives

- Characterize cell populations in paired normal and tumor samples.
- Compare cell-type and subtype composition between tissue conditions.
- Investigate inferred cell–cell communication networks.
- Evaluate expression-based classifiers using patient-level validation.
- Map selected host features and cell-type distributions in spatial tissue sections.
- Explore bacteria–host spatial associations and paired changes in bacterial probe signals.

## Datasets

| Dataset | Role |
|---|---|
| GSE132465 | Single-cell reference and model development; composition analyses use 10 matched normal–tumor patient pairs. |
| GSE144735 | External predictive evaluation using six independent KUL3 patient pairs from the same publication as the development dataset. |
| GSE334323 | Spatial analyses of CRC_03 and CRC_08, comprising four matched normal/tumor sections. |

## Analysis Workflow

| Script | Analysis |
|---|---|
| [01_Data_Preprocessing.R](01_Data_Preprocessing.R) | Raw-count import, metadata preparation, and single-cell quality control. |
| [02_Clustering_Annotation.R](02_Clustering_Annotation.R) | Normalization, PCA, clustering, UMAP, marker detection, and annotation assessment. |
| [03_Cell_Composition.R](03_Cell_Composition.R) | Paired comparisons of cell-type and subtype proportions with multiple-testing correction. |
| [05_Cell_Communication.R](05_Cell_Communication.R) | Comparative CellChat analysis of normal and tumor tissue. |
| [06_Predictive_Modeling.R](06_Predictive_Modeling.R) | Pseudobulk feature construction, classification, patient-level cross-validation, and external evaluation. |
| [07_Spatial_Transcriptomics.R](07_Spatial_Transcriptomics.R) | RCTD deconvolution, selected-feature mapping, and spatial association analyses. |
| [08_Bacteria_Host_Associations.R](08_Bacteria_Host_Associations.R) | Spatially adjusted bacteria–host rank associations and descriptive co-enrichment. |
| [09_Paired_Bacterial_Analysis.R](09_Paired_Bacterial_Analysis.R) | Paired bacterial-signal summaries and comparison of change directions across patients. |

Original stage numbering is retained. Stage 04 is not included in this repository.

## Methods

### Single-Cell Analysis

Seurat is used for preprocessing, normalization, dimensionality reduction, clustering, and marker detection. Annotation assessment includes SingleR, while downstream composition analyses retain deposited GEO cell-type and subtype labels.

Quality control retains cells with:

- More than 200 and fewer than 6,000 detected genes.
- A mitochondrial count percentage below 20%.

### Cell Composition

Cell-type and subtype proportions are compared using paired Wilcoxon tests with Benjamini–Hochberg correction. Subtype percentages are calculated relative to all cells in each sample.

### Cell–Cell Communication

CellChat models are constructed separately for normal and tumor tissue, followed by comparative analysis of inferred signaling networks and pathways.

### Predictive Modeling

Selected-gene pseudobulk log2(CPM + 1) values and gene-signature scores are evaluated using:

- Elastic-net logistic regression.
- Random forest.
- XGBoost.

Internal evaluation uses leave-one-patient-out cross-validation. Elastic-net tuning uses inner patient-level validation. Preprocessing, tuning, and model fitting use development data; external samples are evaluated without external fitting or tuning.

### Spatial Transcriptomics

RCTD estimates cell-type contributions using the single-cell reference. Spatial host expression is represented as:

ln(1 + counts / library size × 10,000)

Spatial associations include approximate rank-based tests using SpatialPack, multiple-testing correction, and sensitivity analyses.

## Software

The analyses are implemented in R.

| Category | Main packages |
|---|---|
| Single-cell analysis | Seurat, SeuratObject, Matrix, data.table, R.utils, SingleR, celldex, SingleCellExperiment, BiocParallel |
| Cell communication | CellChat |
| Predictive modeling | glmnet, ranger, xgboost, pROC |
| Spatial analysis | spacexr, SpatialPack, RANN |
| Visualization and export | ggplot2, patchwork, cowplot, openxlsx |

Each script checks its required dependencies. Some spatial dependencies are installed automatically when missing.

The spatial workflow requires a spacexr version exposing Reference, SpatialRNA, create.RCTD, and run.RCTD.

## Running the Analyses

1. Obtain the required count matrices, annotations, and spatial archive.
2. Install the packages required by the selected analysis.
3. Update input and output paths for your local environment.
4. Run upstream stages before their dependent analyses.

The scripts retain paths from the thesis working environment. Some downstream scripts expect processed objects in Core_From_Raw, while earlier scripts save to stage-specific output folders. Align these paths with the actual upstream outputs before execution.

Stage 07 requires an annotated single-cell reference and spatial data. Stage 08 requires completed Stage 07 outputs, and Stage 09 requires completed Stage 08 outputs.

This repository provides individual analysis scripts rather than an automatically configured end-to-end pipeline.

## Outputs

Depending on the stage, the scripts generate:

- Statistical tables and sample-level summaries.
- PNG figures and spatial maps.
- Processed R objects.
- Fitted model artifacts and prediction tables.
- Session information and method records.

Large input datasets and generated analysis outputs are not included in this script collection.

## Interpretation

- Spatial analyses include two matched patients; paired direction comparisons are descriptive.
- Bacterial features represent deposited probe-group annotations, not independently verified species identities, viability, or absolute abundance.
- Spatial associations and inferred communication do not establish causal or physical interactions.
- External predictive evaluation includes six patient pairs from the same publication as the development cohort and does not establish clinical validity.

## Repository Scope

This repository presents the main analytical code for academic review and research applications. Checker scripts are excluded.

The cleaned scripts preserve the original executable analysis. Full analyses have not been rerun as part of repository preparation.
