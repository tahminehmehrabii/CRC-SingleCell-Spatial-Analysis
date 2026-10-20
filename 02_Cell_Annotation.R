# CORE WORKFLOW: 02_Cell_Annotation.R
# Existing output directories are retained for compatibility.
# ============================================================
# GSE132465
# Clustering, Marker Detection, and Cell-Type Annotation
#
# INPUT:
# E:/CRC/Results/Results/01_Data_and_Preprocessing/02_seurat_qc.rds
#
# OUTPUT DIRECTORY:
# E:/CRC/Results/Results/02_Clustering_Markers_Annotation
# ============================================================


# ============================================================
# 03 — NORMALIZATION, PCA, CLUSTERING, AND UMAP
# ============================================================

local({
  
  input_root <- "E:/CRC/Results/Results/01_Data_and_Preprocessing"
  
  output_root <- "E:/CRC/Results/Results/02_Clustering_Markers_Annotation"
  
  dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (!dir.exists(output_root)) {
    stop("Could not create output directory: ", output_root)
  }
  
  out <- function(name) file.path(output_root, name)
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Install package: Seurat")
  }
  
  library(Seurat)
  
  set.seed(42)
  
  input_file <- file.path(
    input_root,
    "02_seurat_qc.rds"
  )
  
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  obj <- readRDS(input_file)
  
  DefaultAssay(obj) <- "RNA"
  
  
  # Normalize counts
  obj <- NormalizeData(
    object = obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  
  # Select 2,000 variable genes
  obj <- FindVariableFeatures(
    object = obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  if (length(VariableFeatures(obj)) < 2000L) {
    stop("Fewer than 2,000 variable genes are available.")
  }
  
  
  # Scale variable genes
  obj <- ScaleData(
    object = obj,
    features = VariableFeatures(obj),
    verbose = FALSE
  )
  
  
  # PCA
  obj <- RunPCA(
    object = obj,
    features = VariableFeatures(obj),
    npcs = 30,
    verbose = FALSE
  )
  
  if (ncol(Embeddings(obj, reduction = "pca")) < 30L) {
    stop("Fewer than 30 principal components were calculated.")
  }
  
  
  # Elbow plot
  pdf(
    out("03_elbow.pdf"),
    width = 7,
    height = 5
  )
  
  print(
    ElbowPlot(
      obj,
      ndims = 30
    )
  )
  
  dev.off()
  
  
  # Nearest-neighbor graph
  obj <- FindNeighbors(
    object = obj,
    dims = 1:20,
    verbose = FALSE
  )
  
  
  # Clustering
  obj <- FindClusters(
    object = obj,
    resolution = 0.5,
    verbose = FALSE
  )
  
  
  # UMAP
  obj <- RunUMAP(
    object = obj,
    dims = 1:20,
    seed.use = 42,
    verbose = FALSE
  )
  
  
  # Cluster UMAP
  pdf(
    out("03_umap_clusters.pdf"),
    width = 8,
    height = 6
  )
  
  print(
    DimPlot(
      obj,
      reduction = "umap",
      label = TRUE
    )
  )
  
  dev.off()
  
  
  # Save clustered object
  saveRDS(
    obj,
    out("03_clustered.rds")
  )
  
  
  cat("\n")
  cat("========================================\n")
  cat("STEP 03 COMPLETE\n")
  cat("========================================\n")
  cat("Cells:", ncol(obj), "\n")
  cat(
    "Variable genes:",
    length(VariableFeatures(obj)),
    "\n"
  )
  cat(
    "PCA components:",
    ncol(Embeddings(obj, "pca")),
    "\n"
  )
  cat(
    "Clusters:",
    length(unique(as.character(Idents(obj)))),
    "\n"
  )
  cat(
    "Saved:",
    out("03_clustered.rds"),
    "\n"
  )
  cat("========================================\n\n")
  
})

invisible(gc())



# ============================================================
# 04 — CLUSTER MARKER DETECTION
# ============================================================

local({
  
  output_root <- "E:/CRC/Results/Results/02_Clustering_Markers_Annotation"
  
  if (!dir.exists(output_root)) {
    stop("Output directory does not exist: ", output_root)
  }
  
  out <- function(name) file.path(output_root, name)
  
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Install package: Seurat")
  }
  
  library(Seurat)
  
  
  input_file <- out("03_clustered.rds")
  
  if (!file.exists(input_file)) {
    stop("Missing input: ", input_file)
  }
  
  
  obj <- readRDS(input_file)
  
  DefaultAssay(obj) <- "RNA"
  
  
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    stop(
      "Cluster identities are missing from object metadata."
    )
  }
  
  Idents(obj) <- "seurat_clusters"
  
  
  # Positive markers for each cluster
  markers <- FindAllMarkers(
    object = obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 2,
    test.use = "wilcox"
  )
  
  
  if (!nrow(markers)) {
    stop("No cluster markers were returned.")
  }
  
  
  if (!"gene" %in% colnames(markers)) {
    stop(
      "FindAllMarkers did not return a gene column; ",
      "inspect its output."
    )
  }
  
  
  if (
    anyNA(markers$gene) ||
    !all(
      as.character(markers$gene) %in%
      rownames(obj[["RNA"]])
    )
  ) {
    
    stop(
      "Some marker gene names do not match the RNA assay. ",
      "Files were not overwritten."
    )
  }
  
  
  # Save marker table
  write.csv(
    markers,
    out("04_cluster_markers.csv"),
    row.names = FALSE
  )
  
  
  # Save Seurat object
  saveRDS(
    obj,
    out("04_markers_object.rds")
  )
  
  
  cat("\n")
  cat("========================================\n")
  cat("STEP 04 COMPLETE\n")
  cat("========================================\n")
  
  cat(
    "Cells:",
    ncol(obj),
    "\n"
  )
  
  cat(
    "Clusters:",
    length(unique(as.character(Idents(obj)))),
    "\n"
  )
  
  cat(
    "Marker rows:",
    nrow(markers),
    "\n"
  )
  
  cat(
    "Unique gene names:",
    length(unique(markers$gene)),
    "\n"
  )
  
  cat(
    "Valid gene names:",
    sum(markers$gene %in% rownames(obj[["RNA"]])),
    "\n"
  )
  
  cat(
    "Saved:",
    out("04_cluster_markers.csv"),
    "\n"
  )
  
  cat(
    "Saved:",
    out("04_markers_object.rds"),
    "\n"
  )
  
  cat("========================================\n\n")
  
})

invisible(gc())



# ============================================================
# 05 — CELL-TYPE ANNOTATION AND VALIDATION
# ============================================================

local({
  
  required_packages <- c(
    "Seurat",
    "SingleR",
    "celldex",
    "SingleCellExperiment",
    "BiocParallel"
  )
  
  
  for (p in required_packages) {
    
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(
        "Install required package: ",
        p
      )
    }
  }
  
  
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(BiocParallel)
  
  
  output_root <- "E:/CRC/Results/Results/02_Clustering_Markers_Annotation"
  
  if (!dir.exists(output_root)) {
    stop("Output directory does not exist: ", output_root)
  }
  
  out <- function(name) file.path(output_root, name)
  
  
  input_file <- out("04_markers_object.rds")
  
  if (!file.exists(input_file)) {
    stop(
      "Missing input: ",
      input_file
    )
  }
  
  
  obj <- readRDS(input_file)
  
  DefaultAssay(obj) <- "RNA"
  
  
  # GEO annotations
  obj$GEO_CellType <-
    if ("Cell_type" %in% colnames(obj@meta.data)) {
      
      as.character(obj$Cell_type)
      
    } else {
      
      "Unknown"
    }
  
  
  obj$GEO_CellSubtype <-
    if ("Cell_subtype" %in% colnames(obj@meta.data)) {
      
      as.character(obj$Cell_subtype)
      
    } else {
      
      "Unknown"
    }
  
  
  # Human Primary Cell Atlas reference
  ref <- celldex::HumanPrimaryCellAtlasData()
  
  
  # Convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(
    obj,
    assay = "RNA"
  )
  
  
  # SingleR annotation
  pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.main,
    BPPARAM = SerialParam()
  )
  
  
  # Pruned labels
  labels <- pred$pruned.labels
  
  labels[is.na(labels)] <- "Unknown"
  
  names(labels) <- rownames(pred)
  
  
  obj$SingleR_CellType <-
    unname(
      labels[colnames(obj)]
    )
  
  
  if (anyNA(obj$SingleR_CellType)) {
    stop("Unmatched SingleR cell IDs.")
  }
  
  
  # Annotation comparison table
  comparison <- data.frame(
    cell = colnames(obj),
    sample = obj$sample_id,
    GEO = obj$GEO_CellType,
    SingleR = obj$SingleR_CellType,
    stringsAsFactors = FALSE
  )
  
  
  write.csv(
    comparison,
    out("05_annotation_comparison.csv"),
    row.names = FALSE
  )
  
  
  # Use GEO annotation for downstream analyses
  obj$analysis_cell_type <- obj$GEO_CellType
  
  
  if (
    anyNA(obj$analysis_cell_type) ||
    all(obj$analysis_cell_type == "Unknown")
  ) {
    
    stop(
      "GEO cell types missing; ",
      "review annotation before analysis."
    )
  }
  
  
  # Canonical marker genes
  canonical <- c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    "MS4A1",
    "CD79A",
    "CD68",
    "LYZ",
    "NKG7",
    "GNLY",
    "PECAM1",
    "VWF"
  )
  
  
  features <- intersect(
    canonical,
    rownames(obj)
  )
  
  
  if (!length(features)) {
    stop("Canonical marker genes absent.")
  }
  
  
  # DotPlot
  pdf(
    out("05_canonical_dotplot.pdf"),
    width = 13,
    height = 6
  )
  
  print(
    DotPlot(
      obj,
      features = features,
      group.by = "analysis_cell_type"
    ) +
      RotatedAxis()
  )
  
  dev.off()
  
  
  # FeaturePlots
  pdf(
    out("05_canonical_featureplots.pdf"),
    width = 12,
    height = 8
  )
  
  for (gene in features) {
    
    print(
      FeaturePlot(
        obj,
        features = gene
      )
    )
  }
  
  dev.off()
  
  
  # Save SingleR prediction
  saveRDS(
    pred,
    out("05_singler_prediction.rds")
  )
  
  
  # Save annotated Seurat object
  saveRDS(
    obj,
    out("05_annotated.rds")
  )
  
  
  cat("\n")
  cat("========================================\n")
  cat("STEP 05 COMPLETE\n")
  cat("========================================\n")
  
  cat(
    "Cells:",
    ncol(obj),
    "\n"
  )
  
  cat(
    "GEO cell types:",
    length(unique(obj$GEO_CellType)),
    "\n"
  )
  
  cat(
    "SingleR cell types:",
    length(unique(obj$SingleR_CellType)),
    "\n"
  )
  
  cat(
    "Canonical markers found:",
    length(features),
    "of",
    length(canonical),
    "\n"
  )
  
  cat(
    "Saved:",
    out("05_annotation_comparison.csv"),
    "\n"
  )
  
  cat(
    "Saved:",
    out("05_canonical_dotplot.pdf"),
    "\n"
  )
  
  cat(
    "Saved:",
    out("05_canonical_featureplots.pdf"),
    "\n"
  )
  
  cat(
    "Saved:",
    out("05_singler_prediction.rds"),
    "\n"
  )
  
  cat(
    "Saved:",
    out("05_annotated.rds"),
    "\n"
  )
  
  cat("========================================\n")
  cat(
    "ALL OUTPUTS SAVED IN:\n",
    output_root,
    "\n"
  )
  cat("========================================\n\n")
  
})

invisible(gc())