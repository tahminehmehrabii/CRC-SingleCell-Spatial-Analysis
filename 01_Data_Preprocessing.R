# Data Preprocessing
# GSE132465 colorectal cancer scRNA-seq; raw-count import and quality control.

local({

  # Input and output paths

  input_root <- "E:/CRC/Results/GSE132465_scRNA"

  root <- "E:/CRC/Results/Results/01_Data_and_Preprocessing"

  dir.create(
    root,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(root)) {
    stop("Cannot create directory: ", root)
  }

  out <- function(name) file.path(root, name)

  # Required packages

  for (p in c(
    "Matrix",
    "data.table",
    "R.utils",
    "Seurat",
    "ggplot2"
  )) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Install package: ", p)
    }
  }

  library(Matrix)
  library(data.table)

  # Input files

  count_file <- file.path(
    input_root,
    "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
  )

  annotation_file <- file.path(
    input_root,
    "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
  )

  input_files <- c(count_file, annotation_file)

  if (!all(file.exists(input_files))) {
    stop(
      "Missing GEO input files: ",
      paste(
        input_files[!file.exists(input_files)],
        collapse = "; "
      )
    )
  }

  # Metadata and paired sample selection

  metadata <- as.data.frame(
    fread(
      annotation_file,
      sep = "\t",
      check.names = FALSE
    )
  )

  required <- c("Index", "Sample", "Patient", "Class")

  if (!all(required %in% names(metadata))) {
    stop(
      "Missing metadata columns: ",
      paste(
        setdiff(required, names(metadata)),
        collapse = ", "
      )
    )
  }

  if (
    anyNA(metadata$Index) ||
    anyDuplicated(metadata$Index)
  ) {
    stop("Invalid or duplicated cell IDs")
  }

  rownames(metadata) <- as.character(metadata$Index)

  patients <- sprintf("SMC%02d", 1:10)

  wanted <- as.vector(
    outer(patients, c("-N", "-T"), paste0)
  )

  metadata <- metadata[
    as.character(metadata$Sample) %in% wanted,
    ,
    drop = FALSE
  ]

  for (patient in patients) {

    observed <- unique(as.character(metadata$Sample))

    if (!all(
      paste0(patient, c("-N", "-T")) %in% observed
    )) {
      stop("Missing pair for ", patient)
    }

    patient_ids <- unique(
      as.character(
        metadata$Patient[
          as.character(metadata$Sample) %in%
            paste0(patient, c("-N", "-T"))
        ]
      )
    )

    if (length(patient_ids) != 1L) {
      stop("Inconsistent patient IDs for ", patient)
    }
  }

  # Import selected cells into a sparse matrix

  read_sparse_counts <- function(
    path,
    selected,
    block_rows = 100L
  ) {

    con <- gzfile(path, "rt")
    on.exit(close(con))

    header <- strsplit(
      readLines(con, n = 1L),
      "\t",
      fixed = TRUE
    )[[1L]]

    all_cells <- header[-1L]

    if (anyDuplicated(all_cells)) {
      stop("Duplicate matrix cell IDs")
    }

    indices <- match(selected, all_cells)

    if (anyNA(indices)) {
      stop(
        "Selected annotation cells absent from count matrix: ",
        sum(is.na(indices))
      )
    }

    blocks <- list()
    genes <- character()
    b <- 0L

    repeat {

      lines <- readLines(con, n = block_rows)

      if (!length(lines)) break

      fields <- strsplit(lines, "\t", fixed = TRUE)

      if (any(lengths(fields) != length(header))) {
        stop("Malformed count matrix row")
      }

      gene <- vapply(
        fields,
        `[`,
        character(1),
        1L
      )

      values <- lapply(
        fields,
        function(z) as.numeric(z[indices + 1L])
      )

      invalid_counts <- vapply(
        values,
        function(v) {
          anyNA(v) ||
            any(!is.finite(v)) ||
            any(v < 0) ||
            any(v != floor(v))
        },
        logical(1)
      )

      if (any(invalid_counts)) {
        stop("Invalid UMI counts")
      }

      nz <- lapply(
        values,
        function(v) which(v != 0)
      )

      ii <- rep.int(seq_along(nz), lengths(nz))
      jj <- unlist(nz, use.names = FALSE)

      xx <- unlist(
        Map(
          function(v, j) v[j],
          values,
          nz
        ),
        use.names = FALSE
      )

      b <- b + 1L

      blocks[[b]] <- sparseMatrix(
        i = ii,
        j = jj,
        x = xx,
        dims = c(length(lines), length(selected))
      )

      genes <- c(genes, gene)

      if (b %% 25L == 0L) {
        message("Imported genes: ", length(genes))
      }
    }

    if (!length(blocks)) {
      stop("Empty matrix")
    }

    mat <- do.call(rbind, blocks)

    if (
      anyNA(genes) ||
      any(!nzchar(genes)) ||
      anyDuplicated(genes)
    ) {
      stop(
        "Blank/duplicate gene identifiers; inspect raw matrix"
      )
    }

    rownames(mat) <- genes
    colnames(mat) <- selected

    mat
  }

  con <- gzfile(count_file, "rt")
  header <- readLines(con, n = 1L)
  close(con)

  matrix_cells <- strsplit(
    header,
    "\t",
    fixed = TRUE
  )[[1L]][-1L]

  selected <- intersect(
    rownames(metadata),
    matrix_cells
  )

  if (!length(selected)) {
    stop("No matching selected cells")
  }

  if (length(selected) != nrow(metadata)) {
    stop(
      "Annotation cells are absent from count matrix; inspect inputs"
    )
  }

  metadata <- metadata[selected, , drop = FALSE]

  counts <- read_sparse_counts(
    count_file,
    selected
  )

  stopifnot(
    identical(colnames(counts), rownames(metadata))
  )

  # Save import outputs

  write.csv(
    as.data.frame(
      table(metadata$Patient, metadata$Class)
    ),
    out("01_pair_counts.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      genes = nrow(counts),
      cells = ncol(counts)
    ),
    out("01_import_summary.csv"),
    row.names = FALSE
  )

  saveRDS(
    counts,
    out("01_counts_paired.rds")
  )

  saveRDS(
    metadata,
    out("01_metadata_paired.rds")
  )

})

invisible(gc())

local({

  # Load imported data

  library(Seurat)

  root <- "E:/CRC/Results/Results/01_Data_and_Preprocessing"

  if (!dir.exists(root)) {
    stop("Directory does not exist: ", root)
  }

  out <- function(name) file.path(root, name)

  counts <- readRDS(
    out("01_counts_paired.rds")
  )

  metadata <- readRDS(
    out("01_metadata_paired.rds")
  )

  stopifnot(
    identical(colnames(counts), rownames(metadata))
  )

  # Create seurat object

  obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    project = "GSE132465_CRC",
    min.cells = 3,
    min.features = 200
  )

  obj$sample_id <- as.character(obj$Sample)

  obj$patient_id <- sub(
    "-[NT]$",
    "",
    obj$sample_id
  )

  if (any(
    as.character(obj$Patient) != obj$patient_id
  )) {
    stop("Patient metadata mismatch")
  }

  obj$group <- as.character(obj$Class)

  if (
    !all(obj$group %in% c("Normal", "Tumor")) ||
    any(
      obj$sample_id != paste0(
        obj$patient_id,
        ifelse(obj$group == "Normal", "-N", "-T")
      )
    )
  ) {
    stop("Unexpected tissue/patient labels")
  }

  if (!any(grepl("^MT-", rownames(obj)))) {
    stop(
      "No MT- genes: mitochondrial QC cannot be computed"
    )
  }

  obj$percent.mt <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )

  # Save qc metadata and plots before filtering

  write.csv(
    obj@meta.data,
    out("02_qc_before.csv"),
    row.names = TRUE
  )

  pdf(
    out("02_qc_before.pdf"),
    width = 12,
    height = 5
  )

  print(
    VlnPlot(
      obj,
      features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt"
      ),
      ncol = 3,
      pt.size = 0
    )
  )

  dev.off()

  before <- ncol(obj)

  # Quality-control filtering

  obj <- subset(
    obj,
    subset =
      nFeature_RNA > 200 &
      nFeature_RNA < 6000 &
      percent.mt < 20
  )

  if (ncol(obj) == 0) {
    stop("No cells remain after QC")
  }

  # Save qc metadata and plots after filtering

  write.csv(
    obj@meta.data,
    out("02_qc_after.csv"),
    row.names = TRUE
  )

  qc_summary <- data.frame(
    Stage = c("Before_QC", "After_QC"),
    Cells = c(before, ncol(obj))
  )

  write.csv(
    qc_summary,
    out("02_qc_summary.csv"),
    row.names = FALSE
  )

  pdf(
    out("02_qc_after.pdf"),
    width = 12,
    height = 5
  )

  print(
    VlnPlot(
      obj,
      features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt"
      ),
      ncol = 3,
      pt.size = 0
    )
  )

  dev.off()

  # Save final seurat object and session information

  saveRDS(
    obj,
    out("02_seurat_qc.rds")
  )

  writeLines(
    capture.output(sessionInfo()),
    out("01_sessionInfo.txt")
  )

  cat("\nData and Preprocessing completed.\n")
  cat("Cells before QC:", before, "\n")
  cat("Cells after QC:", ncol(obj), "\n")
  cat("Cells removed:", before - ncol(obj), "\n")
  cat("Genes in final Seurat object:", nrow(obj), "\n")
  cat(
    "Maximum mitochondrial percentage:",
    max(obj$percent.mt),
    "\n"
  )
  cat("Output folder:", root, "\n")

})

invisible(gc())
