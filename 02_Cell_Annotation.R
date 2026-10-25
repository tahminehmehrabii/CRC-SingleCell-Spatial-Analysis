# GSE132465: clustering, cluster markers and reference-based annotation comparison.
# Published GEO labels remain the downstream cell-type annotations.
local({
  input_dir <- Sys.getenv("CRC_SCRNA_OUTPUT", "E:/CRC/Results/Results/01_Data_and_Preprocessing")
  output_dir <- Sys.getenv("CRC_ANNOTATION_OUTPUT", "E:/CRC/Results/Results/02_Clustering_Markers_Annotation")
  packages <- c("Seurat", "SingleR", "celldex", "SingleCellExperiment", "BiocParallel")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(Seurat))
  set.seed(42)

  input_file <- file.path(input_dir, "02_seurat_qc.rds")
  if (!file.exists(input_file)) stop("Missing input: ", input_file)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
  out <- function(name) file.path(output_dir, name)
  save_pdf <- function(name, draw, width, height) {
    grDevices::pdf(out(name), width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw()
  }

  # Normalize, reduce dimensions and cluster using uncorrected RNA-derived PCs.
  obj <- readRDS(input_file)
  if (!inherits(obj, "Seurat") || !"RNA" %in% Assays(obj)) stop("Expected a Seurat object with an RNA assay")
  DefaultAssay(obj) <- "RNA"
  needed <- c("sample_id", "Cell_type")
  if (!all(needed %in% names(obj[[]]))) stop("Missing metadata: ", paste(setdiff(needed, names(obj[[]])), collapse = ", "))
  geo_labels <- as.character(obj$Cell_type)
  if (anyNA(geo_labels) || any(!nzchar(trimws(geo_labels))) || all(geo_labels == "Unknown")) stop("Published cell-type labels are unavailable")

  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  if (length(VariableFeatures(obj)) < 2000L) stop("Fewer than 2,000 variable genes")
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
  if (ncol(Embeddings(obj, "pca")) < 30L) stop("Fewer than 30 principal components")
  save_pdf("03_elbow.pdf", function() print(ElbowPlot(obj, ndims = 30)), 7, 5)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, seed.use = 42, verbose = FALSE)
  save_pdf("03_umap_clusters.pdf", function() print(DimPlot(obj, reduction = "umap", label = TRUE)), 8, 6)
  saveRDS(obj, out("03_clustered.rds"))

  # Positive cluster markers; cell-level tests are descriptive, not patient-level DE.
  if (!"seurat_clusters" %in% names(obj[[]])) stop("Missing cluster identities")
  Idents(obj) <- "seurat_clusters"
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                            logfc.threshold = 2, test.use = "wilcox")
  if (!nrow(markers)) stop("No cluster markers returned")
  if (!"gene" %in% names(markers) || anyNA(markers$gene) ||
      !all(markers$gene %in% rownames(obj[["RNA"]]))) stop("Invalid marker gene identifiers")
  utils::write.csv(markers, out("04_cluster_markers.csv"), row.names = FALSE)
  saveRDS(obj, out("04_markers_object.rds"))

  # celldex may download/cache the HPCA reference on first use.
  obj$GEO_CellType <- geo_labels
  obj$GEO_CellSubtype <- if ("Cell_subtype" %in% names(obj[[]])) as.character(obj$Cell_subtype) else "Unknown"
  ref <- celldex::HumanPrimaryCellAtlasData()
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  pred <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main,
                          BPPARAM = BiocParallel::SerialParam())
  if (is.null(rownames(pred)) || anyDuplicated(rownames(pred))) stop("Invalid SingleR cell identifiers")
  matched <- match(colnames(obj), rownames(pred))
  if (anyNA(matched)) stop("Unmatched SingleR cell identifiers")
  labels <- as.character(pred$pruned.labels[matched])
  labels[is.na(labels)] <- "Unknown"
  obj$SingleR_CellType <- labels
  obj$analysis_cell_type <- obj$GEO_CellType
  comparison <- data.frame(cell = colnames(obj), sample = obj$sample_id,
                           GEO = obj$GEO_CellType, SingleR = obj$SingleR_CellType)
  utils::write.csv(comparison, out("05_annotation_comparison.csv"), row.names = FALSE)
  rm(sce, ref)
  invisible(gc())

  canonical <- c("EPCAM", "KRT8", "KRT18", "CD3D", "CD3E", "CD4", "CD8A",
                 "MS4A1", "CD79A", "CD68", "LYZ", "NKG7", "GNLY", "PECAM1", "VWF")
  features <- intersect(canonical, rownames(obj))
  if (!length(features)) stop("No canonical markers available")
  save_pdf("05_canonical_dotplot.pdf", function() {
    print(DotPlot(obj, features = features, group.by = "analysis_cell_type") + RotatedAxis())
  }, 13, 6)
  save_pdf("05_canonical_featureplots.pdf", function() {
    for (gene in features) print(FeaturePlot(obj, features = gene))
  }, 12, 8)
  saveRDS(pred, out("05_singler_prediction.rds"))
  saveRDS(obj, out("05_annotated.rds"))
  writeLines(capture.output(utils::sessionInfo()), out("02_sessionInfo.txt"))
  cat("\nClustering and annotation complete\nCells:", ncol(obj),
      "\nClusters:", length(unique(as.character(Idents(obj)))),
      "\nMarker rows:", nrow(markers),
      "\nPublished cell types:", length(unique(obj$GEO_CellType)),
      "\nSingleR labels (including Unknown):", length(unique(obj$SingleR_CellType)),
      "\nCanonical markers:", length(features), "of", length(canonical),
      "\nOutput:", output_dir, "\n")
})
invisible(gc())
