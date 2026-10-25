# GSE132465: import 10 matched tumor-normal pairs and apply cell-level QC.
# Requires locally downloaded GEO files.
local({
  input_root <- Sys.getenv("CRC_SCRNA_INPUT", "E:/CRC/Results/GSE132465_scRNA")
  output_dir <- Sys.getenv("CRC_SCRNA_OUTPUT", "E:/CRC/Results/Results/01_Data_and_Preprocessing")
  patients <- sprintf("SMC%02d", 1:10)
  block_rows <- 100L

  packages <- c("Matrix", "data.table", "R.utils", "Seurat", "ggplot2")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required packages: ", paste(missing, collapse = ", "))

  count_file <- file.path(input_root, "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
  annotation_file <- file.path(input_root, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")
  missing <- c(count_file, annotation_file)[!file.exists(c(count_file, annotation_file))]
  if (length(missing)) stop("Missing input: ", paste(missing, collapse = "; "))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
  out <- function(name) file.path(output_dir, name)

  # Select the prespecified matched samples using published annotations.
  metadata <- as.data.frame(data.table::fread(annotation_file, sep = "\t", check.names = FALSE))
  required <- c("Index", "Sample", "Patient", "Class")
  if (!all(required %in% names(metadata))) stop("Missing metadata: ", paste(setdiff(required, names(metadata)), collapse = ", "))
  if (anyNA(metadata$Index) || any(!nzchar(metadata$Index)) || anyDuplicated(metadata$Index)) stop("Invalid cell IDs")
  rownames(metadata) <- as.character(metadata$Index)
  wanted <- as.vector(outer(patients, c("-N", "-T"), paste0))
  metadata <- metadata[metadata$Sample %in% wanted, , drop = FALSE]
  if (anyNA(metadata[, required])) stop("Missing sample, patient or tissue labels")
  if (!all(wanted %in% metadata$Sample)) stop("One or more matched samples are missing")
  expected_patient <- sub("-[NT]$", "", as.character(metadata$Sample))
  expected_class <- ifelse(grepl("-N$", metadata$Sample), "Normal", "Tumor")
  if (any(metadata$Patient != expected_patient) || any(metadata$Class != expected_class)) stop("Inconsistent patient/tissue labels")

  # Stream gene blocks rather than importing the full dense matrix.
  read_sparse_counts <- function(path, selected, block_rows) {
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    header <- readLines(con, n = 1L, warn = FALSE)
    if (!length(header)) stop("Empty count file")
    header <- gsub('^"|"$', "", strsplit(header, "\t", fixed = TRUE)[[1L]])
    cells <- header[-1L]
    if (anyDuplicated(cells)) stop("Duplicate matrix cell IDs")
    indices <- match(selected, cells)
    if (anyNA(indices)) stop("Selected cells absent from count matrix")
    blocks <- list()
    genes <- character()
    repeat {
      lines <- readLines(con, n = block_rows, warn = FALSE)
      if (!length(lines)) break
      fields <- strsplit(lines, "\t", fixed = TRUE)
      if (any(lengths(fields) != length(header))) stop("Malformed count row")
      gene <- gsub('^"|"$', "", vapply(fields, `[`, character(1), 1L))
      values <- lapply(fields, function(z) as.numeric(z[indices + 1L]))
      invalid <- vapply(values, function(v) any(!is.finite(v)) || any(v < 0) || any(v != floor(v)), logical(1))
      if (any(invalid)) stop("Invalid UMI counts")
      nz <- lapply(values, function(v) which(v != 0))
      blocks[[length(blocks) + 1L]] <- Matrix::sparseMatrix(
        i = rep.int(seq_along(nz), lengths(nz)),
        j = unlist(nz, use.names = FALSE),
        x = unlist(Map(function(v, j) v[j], values, nz), use.names = FALSE),
        dims = c(length(lines), length(selected))
      )
      genes <- c(genes, gene)
      if (length(blocks) %% 25L == 0L) message("Imported genes: ", length(genes))
    }
    if (!length(blocks)) stop("No count rows")
    if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) stop("Invalid gene identifiers")
    mat <- do.call(rbind, blocks)
    dimnames(mat) <- list(genes, selected)
    mat
  }

  counts <- read_sparse_counts(count_file, rownames(metadata), block_rows)
  stopifnot(identical(colnames(counts), rownames(metadata)))
  utils::write.csv(as.data.frame(table(metadata$Patient, metadata$Class)), out("01_pair_counts.csv"), row.names = FALSE)
  utils::write.csv(data.frame(genes = nrow(counts), cells = ncol(counts)), out("01_import_summary.csv"), row.names = FALSE)
  saveRDS(counts, out("01_counts_paired.rds"))
  saveRDS(metadata, out("01_metadata_paired.rds"))

  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata,
                                  project = "GSE132465_CRC", min.cells = 3, min.features = 200)
  rm(counts, metadata)
  invisible(gc())
  obj$sample_id <- as.character(obj$Sample)
  obj$patient_id <- sub("-[NT]$", "", obj$sample_id)
  obj$group <- as.character(obj$Class)
  if (!any(grepl("^MT-", rownames(obj)))) stop("No MT- genes for mitochondrial QC")
  obj$percent.mt <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")

  save_qc <- function(object, phase) {
    utils::write.csv(object[[]], out(paste0("02_qc_", phase, ".csv")), row.names = TRUE)
    grDevices::pdf(out(paste0("02_qc_", phase, ".pdf")), width = 12, height = 5)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(Seurat::VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
  }
  # Before-QC records follow initial Seurat inclusion, not raw input.
  before <- ncol(obj)
  save_qc(obj, "before")
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
  if (!ncol(obj)) stop("No cells remain after QC")
  save_qc(obj, "after")
  utils::write.csv(data.frame(Stage = c("Before_QC", "After_QC"), Cells = c(before, ncol(obj))),
                   out("02_qc_summary.csv"), row.names = FALSE)
  saveRDS(obj, out("02_seurat_qc.rds"))
  writeLines(capture.output(utils::sessionInfo()), out("01_sessionInfo.txt"))
  cat("\nPreprocessing complete\nCells before additional QC:", before,
      "\nCells after QC:", ncol(obj), "\nCells removed:", before - ncol(obj),
      "\nRetained genes:", nrow(obj), "\nMaximum mitochondrial percentage:", max(obj$percent.mt),
      "\nOutput:", output_dir, "\n")
})
invisible(gc())
