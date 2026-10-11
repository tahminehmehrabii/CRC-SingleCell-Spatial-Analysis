# Spatial Transcriptomics
# GSE334323: two paired patients; RCTD deconvolution and selected-feature mapping.
# Spatial expression is ln(1 + counts/library size * 10000), not ML log2-CPM.
# Spatial rank tests are approximate; deposited bacterial probe groups are not verified taxonomy.

# Spatial analysis of selected features

# Configuration: edit paths here only

output_root <- "E:/CRC/Results/Results/07-Spatial_Analysis_of_Selected_Features"

spatial_config <- list(
  data_root = "E:/CRC/Results/GSE334323_ST",
  archive = "E:/CRC/Results/GSE334323_ST/GSE334323_RAW.tar",
  reference_file = "E:/CRC/Results/GSE132465_scRNA/Core_From_Raw/05_annotated.rds",
  raw_dir = file.path(output_root, "Matched_RAW_Extracted"),
  analysis_dir = output_root,
  figure_dir = file.path(output_root, "Figures"),
  seed = 42L,
  max_cores = 8L
)

# Install missing cran dependencies for spatial tests and excel export

local({
  cran_packages <- c("SpatialPack", "openxlsx")
  missing <- cran_packages[!vapply(cran_packages, requireNamespace,
                                   logical(1), quietly = TRUE)]
  if (length(missing)) {
    lib <- .libPaths()[1L]
    if (!dir.exists(lib)) dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    if (file.access(lib, 2L) != 0L) {
      lib <- file.path(path.expand("~"), "R", "CRC_packages",
                       paste(R.version$major, strsplit(R.version$minor, ".", fixed=TRUE)[[1L]][1L], sep="."))
      dir.create(lib, recursive=TRUE, showWarnings=FALSE)
      if (!dir.exists(lib) || file.access(lib, 2L) != 0L)
        stop("Cannot create a writable R package library: ", lib)
    }
    .libPaths(unique(c(lib, .libPaths())))
    old_timeout <- getOption("timeout")
    options(timeout=max(600, old_timeout))
    tryCatch({
      cat("Installing missing package(s): ", paste(missing, collapse=", "), "\n", sep="")
      install.packages(missing, lib=lib, repos="https://cloud.r-project.org")
    }, error=function(e) stop("Package installation failed: ", conditionMessage(e),
                              "\nCheck the internet connection and R installation permissions.", call.=FALSE),
    finally=options(timeout=old_timeout))
  }
  unavailable <- cran_packages[!vapply(cran_packages, requireNamespace,
                                       logical(1), quietly=TRUE)]
  if (length(unavailable)) stop("Package installation did not complete: ",
                                paste(unavailable, collapse=", "),
                                "\nReview the installation messages above; analysis has not started.", call.=FALSE)
  cat("Spatial-test dependencies: OK\n")
})

run_spatial_workflow <- function(config) {
  previous_options <- options(width = 220, stringsAsFactors = FALSE)
  on.exit(options(previous_options), add = TRUE)
  required_packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "spacexr", "SpatialPack", "openxlsx")
  missing_packages <- required_packages[!vapply(required_packages,
                                                requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop("Install missing package(s) before running: ",
         paste(missing_packages, collapse = ", "),
         "\nFor SpatialPack/openxlsx: install.packages(c('SpatialPack','openxlsx')).\nFor spacexr: remotes::install_github('dmcable/spacexr').",
         call. = FALSE)
  }
  rctd_api <- c("Reference", "SpatialRNA", "create.RCTD", "run.RCTD")
  if (!all(rctd_api %in% getNamespaceExports("spacexr"))) {
    stop("This script requires the dmcable/spacexr GitHub API; installed spacexr has a different API.",
         call. = FALSE)
  }
  for (f in c(config$archive, config$reference_file)) {
    if (!file.exists(f)) stop("Missing input: ", f, call. = FALSE)
  }
  for (p in c(config$raw_dir, config$analysis_dir, config$figure_dir)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(p)) stop("Cannot create directory: ", p, call. = FALSE)
  }
  set.seed(config$seed)

  # Raw spatial data, qc, rctd and feature analyses
  spatial_stage_01 <- function(config) {

    # Inspect raw tar
    # Extract matched samples only
    # Load spatial matrices
    # Host / bacterial feature separation
    # Spatial qc
    # Host normalization
    # Rctd deconvolution
    # Ml-selected genes
    # Ml signature scores
    # Matched normal vs tumor
    # Feature ↔ bacterial load
    # Feature ↔ rctd
    # Bacteria ↔ rctd
    # Spatial maps
    # Final matched-pair consensus

    invisible(gc())

    set.seed(config$seed)

    # Paths

    root <- config$data_root

    archive <- config$archive

    reference_file <- config$reference_file

    raw_dir <- config$raw_dir

    outdir <- config$analysis_dir

    dir.create(
      raw_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    dir.create(
      outdir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    out <- function(x) {
      file.path(outdir, x)
    }

    if (!file.exists(archive)) {
      stop(
        "RAW archive not found:\n",
        archive
      )
    }

    if (!file.exists(reference_file)) {
      stop(
        "scRNA reference not found:\n",
        reference_file
      )
    }

    # Packages

    cran_packages <- c(
      "Seurat",
      "SeuratObject",
      "Matrix",
      "ggplot2"
    )

    for (p in cran_packages) {

      if (!requireNamespace(
        p,
        quietly = TRUE
      )) {

        stop(
          "Install package first: ",
          p
        )
      }
    }

    suppressPackageStartupMessages(
      library(Seurat)
    )

    suppressPackageStartupMessages(
      library(SeuratObject)
    )

    suppressPackageStartupMessages(
      library(Matrix)
    )

    suppressPackageStartupMessages(
      library(ggplot2)
    )

    suppressPackageStartupMessages(
      library(spacexr)
    )

    # Matched sample metadata

    meta <- data.frame(

      GSM = c(
        "GSM9785428",
        "GSM9785429",
        "GSM9785432",
        "GSM9785433"
      ),

      section = c(
        "CRC_03_Normal",
        "CRC_03_Tumor",
        "CRC_08_Normal",
        "CRC_08_Tumor"
      ),

      patient = c(
        "CRC_03",
        "CRC_03",
        "CRC_08",
        "CRC_08"
      ),

      tissue = c(
        "Normal",
        "Tumor",
        "Normal",
        "Tumor"
      ),

      pair_status = "Matched",

      stringsAsFactors = FALSE
    )

    matched_patients <- c(
      "CRC_03",
      "CRC_08"
    )

    matched_sections <- meta$section

    write.csv(
      meta,
      out("00_Matched_sample_metadata.csv"),
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("MATCHED SAMPLES ONLY\n")
    cat("============================================================\n")

    print(
      meta,
      row.names = FALSE
    )

    # First: inspect raw tar

    cat("\n")
    cat("============================================================\n")
    cat("INSPECTING GSE334323_RAW.tar\n")
    cat("============================================================\n")

    tar_entries <- utils::untar(
      archive,
      list = TRUE
    )

    if (!length(tar_entries)) {
      stop(
        "The TAR archive appears empty."
      )
    }

    cat(
      "\nTotal entries in TAR:",
      length(tar_entries),
      "\n"
    )

    cat(
      "\nFirst 100 entries:\n"
    )

    print(
      head(
        tar_entries,
        100
      )
    )

    matched_entry_mask <- Reduce(
      `|`,
      lapply(
        c(
          meta$GSM,
          meta$section
        ),
        function(x) {
          grepl(
            x,
            tar_entries,
            fixed = TRUE
          )
        }
      )
    )

    matched_tar_entries <- tar_entries[
      matched_entry_mask
    ]

    cat(
      "\nEntries directly matching CRC_03 / CRC_08:\n"
    )

    print(
      matched_tar_entries
    )

    write.csv(
      data.frame(
        archive_entry = tar_entries,
        matched =
          matched_entry_mask,
        stringsAsFactors = FALSE
      ),
      out("01_RAW_TAR_inventory.csv"),
      row.names = FALSE
    )

    # Security check before extraction

    unsafe <- grepl(
      "(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",
      tar_entries
    )

    if (any(unsafe)) {

      stop(
        "Unsafe path detected inside TAR archive."
      )
    }

    # Extract raw archive

    cat(
      "\nExtracting RAW archive...\n"
    )

    utils::untar(
      archive,
      exdir = raw_dir
    )

    all_files <- list.files(
      raw_dir,
      recursive = TRUE,
      full.names = TRUE
    )

    if (!length(all_files)) {
      stop(
        "No files were extracted from TAR."
      )
    }

    cat(
      "Files extracted:",
      length(all_files),
      "\n"
    )

    write.csv(
      data.frame(
        path = all_files,
        size_bytes =
          file.info(
            all_files
          )$size
      ),
      out("02_Extracted_file_inventory.csv"),
      row.names = FALSE
    )

    # Helper functions

    read_tsv_auto <- function(f) {

      con <- if (
        grepl(
          "\\.gz$",
          f,
          ignore.case = TRUE
        )
      ) {

        gzfile(
          f,
          "rt"
        )

      } else {

        file(
          f,
          "rt"
        )
      }

      on.exit(
        close(con)
      )

      read.delim(
        con,
        header = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    read_matrix_auto <- function(f) {

      con <- if (
        grepl(
          "\\.gz$",
          f,
          ignore.case = TRUE
        )
      ) {

        gzfile(
          f,
          "rt"
        )

      } else {

        file(
          f,
          "rt"
        )
      }

      on.exit(
        close(con)
      )

      as(
        Matrix::readMM(con),
        "dgCMatrix"
      )
    }

    choose_one <- function(
    files,
    pattern,
    description,
    section
    ) {

      z <- files[
        grepl(
          pattern,
          basename(files),
          ignore.case = TRUE
        )
      ]

      if (length(z) != 1) {

        stop(
          "\n",
          section,
          ": expected exactly one ",
          description,
          " but found ",
          length(z),
          "\nCandidates:\n",
          paste(
            z,
            collapse = "\n"
          )
        )
      }

      z
    }

    sanitize_taxon <- function(x) {

      x <- trimws(
        as.character(x)
      )

      x[
        is.na(x) |
          !nzchar(x)
      ] <- "Unlabeled_bacterial_feature"

      x
    }

    bacterial_probe_names <- function(features) {
      if (ncol(features) < 2L) stop("Expected feature IDs and names in columns 1 and 2.")
      trimws(as.character(features[[2]]))
    }
    probe_group <- function(x) {
      trimws(sub("_probe[0-9]+$", "", x, ignore.case=TRUE))
    }

    bacterial_groups_allowed <- c("Alistipes", "Bacteroides_Phocaeicola", "Butyricimonas", "Campylobacter", "Coprobacter", "Dialister", "Enterobacteriaceae", "Fusobacterium", "Helicobacter", "Intestinimonas", "Lachnospiraceae", "Leptotrichia", "Neisseria", "Parabacteroides", "Parvimonas", "Peptostreptococcus", "Porphyromonas", "Prevotellaceae", "Streptococcaceae")
    valid_probe <- function(x) {
      !is.na(x) & grepl("_probe[0-9]+$", x, ignore.case=TRUE) &
        tolower(probe_group(x)) %in% tolower(bacterial_groups_allowed)
    }

    classify_bacterial_features <- function(features, host_mask) {

      valid_probe(bacterial_probe_names(features))
    }
    extract_taxon_label <- function(features_subset) {
      nm <- bacterial_probe_names(features_subset)
      if (!all(valid_probe(nm))) stop("Invalid bacterial probe names; inspect the deposited feature annotation.")
      probe_group(nm)
    }

    safe_cor <- function(
    x,
    y
    ) {

      ok <-
        is.finite(x) &
        is.finite(y)

      if (sum(ok) < 5) {

        return(
          c(
            rho = NA_real_,
            p = NA_real_,
            n = sum(ok)
          )
        )
      }

      x <- x[ok]
      y <- y[ok]

      if (
        length(unique(x)) < 2 ||
        length(unique(y)) < 2
      ) {

        return(
          c(
            rho = NA_real_,
            p = NA_real_,
            n = length(x)
          )
        )
      }

      ct <- suppressWarnings(
        cor.test(
          x,
          y,
          method = "spearman",
          exact = FALSE
        )
      )

      c(
        rho =
          unname(
            ct$estimate
          ),
        p =
          ct$p.value,
        n =
          length(x)
      )
    }

    # Load scrna-seq reference

    cat("\n")
    cat("============================================================\n")
    cat("LOADING scRNA-seq REFERENCE\n")
    cat("============================================================\n")

    sc_reference <- readRDS(
      reference_file
    )

    if (!inherits(
      sc_reference,
      "Seurat"
    )) {

      stop(
        "Reference file is not a Seurat object."
      )
    }

    DefaultAssay(
      sc_reference
    ) <- "RNA"

    # Get reference counts

    get_seurat_counts <- function(
    obj,
    assay = "RNA"
    ) {

      assay_obj <- obj[[assay]]

      if (inherits(
        assay_obj,
        "Assay5"
      )) {

        layers <- SeuratObject::Layers(
          assay_obj,
          search = "^counts"
        )

        if (!length(layers)) {
          stop(
            "No counts layer in reference."
          )
        }

        if (length(layers) > 1) {

          obj <- SeuratObject::JoinLayers(
            obj,
            assay = assay,
            layers = layers,
            new = "counts"
          )
        }

        x <- SeuratObject::GetAssayData(
          obj,
          assay = assay,
          layer = "counts"
        )

      } else {

        x <- Seurat::GetAssayData(
          obj,
          assay = assay,
          slot = "counts"
        )
      }

      as(
        x,
        "dgCMatrix"
      )
    }

    reference_counts_all <-
      get_seurat_counts(
        sc_reference
      )

    reference_genes <-
      rownames(
        reference_counts_all
      )

    host_symbols <-
      reference_genes

    # Reference cell-type labels

    candidate_celltype_columns <- c(
      "analysis_cell_type",
      "Cell_type",
      "cell_type",
      "GEO_CellType",
      "SingleR_CellType"
    )

    available_celltype_columns <- intersect(
      candidate_celltype_columns,
      colnames(
        sc_reference@meta.data
      )
    )

    if (length(
      available_celltype_columns
    )) {

      celltype_column <-
        available_celltype_columns[1]

      reference_celltype <-
        as.character(
          sc_reference@meta.data[
            ,
            celltype_column
          ]
        )

    } else {

      celltype_column <-
        "active.ident"

      reference_celltype <-
        as.character(
          Idents(
            sc_reference
          )
        )
    }

    names(
      reference_celltype
    ) <- colnames(
      sc_reference
    )

    cat(
      "RCTD cell-type column:",
      celltype_column,
      "\n"
    )

    # Clean / downsample reference

    valid_ref <- names(
      reference_celltype
    )[
      !is.na(
        reference_celltype
      ) &
        nzchar(
          reference_celltype
        ) &
        !tolower(
          reference_celltype
        ) %in%
        c(
          "unknown",
          "unassigned",
          "na"
        )
    ]

    ref_table <- table(
      reference_celltype[
        valid_ref
      ]
    )

    valid_types <- names(
      ref_table
    )[
      ref_table >= 25
    ]

    valid_ref <- valid_ref[
      reference_celltype[
        valid_ref
      ] %in%
        valid_types
    ]

    set.seed(config$seed)

    ref_split <- split(
      valid_ref,
      reference_celltype[
        valid_ref
      ]
    )

    reference_keep <- unlist(
      lapply(
        ref_split,
        function(z) {

          if (length(z) > 500) {

            sample(
              z,
              500
            )

          } else {

            z
          }
        }
      ),
      use.names = FALSE
    )

    reference_counts <-
      reference_counts_all[
        ,
        reference_keep,
        drop = FALSE
      ]

    reference_labels <- factor(
      reference_celltype[
        reference_keep
      ]
    )

    reference_nUMI <-
      Matrix::colSums(
        reference_counts
      )

    good_reference <-
      reference_nUMI > 0

    reference_counts <-
      reference_counts[
        ,
        good_reference,
        drop = FALSE
      ]

    reference_labels <-
      droplevels(
        reference_labels[
          good_reference
        ]
      )

    reference_nUMI <-
      reference_nUMI[
        good_reference
      ]

    rctd_reference <- spacexr::Reference(
      reference_counts,
      reference_labels,
      reference_nUMI
    )

    cat(
      "Reference cells:",
      ncol(reference_counts),
      "\n"
    )

    cat(
      "Reference cell types:",
      nlevels(reference_labels),
      "\n"
    )

    print(
      sort(
        table(reference_labels),
        decreasing = TRUE
      )
    )

    # Ml genes and signatures

    ml_genes <- c(
      "KRT18",
      "CEACAM5",
      "CD68",
      "EPCAM",
      "KRT8",
      "GNLY",
      "LYZ",
      "NKG7",
      "CD79A",
      "MS4A1",
      "CD4",
      "CD3E",
      "PECAM1",
      "VWF",
      "CD8A",
      "CD3D"
    )

    ml_direction <- c(
      KRT18   = "Higher_in_Tumor",
      CEACAM5 = "Higher_in_Tumor",
      CD68    = "Higher_in_Tumor",
      EPCAM   = "Higher_in_Tumor",
      KRT8    = "Higher_in_Tumor",
      GNLY    = "Higher_in_Tumor",
      LYZ     = "Higher_in_Tumor",
      NKG7    = "Higher_in_Tumor",
      CD79A   = "Lower_in_Tumor",
      MS4A1   = "Lower_in_Tumor",
      CD4     = "Higher_in_Tumor",
      CD3E    = "Lower_in_Tumor",
      PECAM1  = "Lower_in_Tumor",
      VWF     = "Lower_in_Tumor",
      CD8A    = "Lower_in_Tumor",
      CD3D    = "Lower_in_Tumor"
    )

    signatures <- list(

      Epithelial = c(
        "EPCAM",
        "KRT8",
        "KRT18",
        "CEACAM5"
      ),

      Myeloid = c(
        "CD68",
        "LYZ"
      ),

      Cytotoxic = c(
        "NKG7",
        "GNLY"
      ),

      B_cell = c(
        "MS4A1",
        "CD79A"
      ),

      T_cell = c(
        "CD3D",
        "CD3E",
        "CD4",
        "CD8A"
      ),

      Endothelial = c(
        "PECAM1",
        "VWF"
      )
    )

    signature_direction <- c(
      Signature__Epithelial =
        "Higher_in_Tumor",
      Signature__Myeloid =
        "Higher_in_Tumor",
      Signature__Cytotoxic =
        "Higher_in_Tumor",
      Signature__B_cell =
        "Lower_in_Tumor",
      Signature__T_cell =
        "Lower_in_Tumor",
      Signature__Endothelial =
        "Lower_in_Tumor"
    )

    # Storage

    host_norm_list <- list()

    spot_expression_list <- list()

    bacterial_spot_list <- list()

    bacterial_taxa_list <- list()

    rctd_long_list <- list()

    rctd_status_list <- list()

    qc_list <- list()

    availability_list <- list()

    manifest_list <- list()

    # Process four matched sections

    for (
      i in seq_len(
        nrow(meta)
      )
    ) {

      gsm <- meta$GSM[i]
      ss <- meta$section[i]
      patient_i <- meta$patient[i]
      tissue_i <- meta$tissue[i]

      cat("\n")
      cat("============================================================\n")
      cat(
        "PROCESSING: ",
        ss,
        " | ",
        gsm,
        "\n",
        sep = ""
      )
      cat("============================================================\n")

      # Find files associated with this sample

      sample_files <- all_files[
        grepl(
          gsm,
          all_files,
          fixed = TRUE
        ) |
          grepl(
            ss,
            all_files,
            fixed = TRUE
          )
      ]

      if (!length(sample_files)) {

        stop(
          "No files found for ",
          ss,
          " / ",
          gsm
        )
      }

      # Extract nested archives

      nested <- sample_files[
        grepl(
          "\\.(tar|tar\\.gz|tgz|zip)$",
          sample_files,
          ignore.case = TRUE
        )
      ]

      if (length(nested)) {

        sample_extract_dir <- file.path(
          raw_dir,
          paste0(
            "nested_",
            ss
          )
        )

        dir.create(
          sample_extract_dir,
          recursive = TRUE,
          showWarnings = FALSE
        )

        for (nf in nested) {

          cat(
            "Extracting nested archive: ",
            basename(nf),
            "\n",
            sep = ""
          )

          if (grepl(
            "\\.zip$",
            nf,
            ignore.case = TRUE
          )) {

            entries_nested <- unzip(
              nf,
              list = TRUE
            )$Name

            unsafe_nested <- grepl(
              "(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",
              entries_nested
            )

            if (any(unsafe_nested)) {
              stop(
                "Unsafe path in nested ZIP."
              )
            }

            unzip(
              nf,
              exdir =
                sample_extract_dir
            )

          } else {

            entries_nested <- untar(
              nf,
              list = TRUE
            )

            unsafe_nested <- grepl(
              "(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",
              entries_nested
            )

            if (any(unsafe_nested)) {
              stop(
                "Unsafe path in nested TAR."
              )
            }

            untar(
              nf,
              exdir =
                sample_extract_dir
            )
          }
        }

        nested_files <- list.files(
          sample_extract_dir,
          recursive = TRUE,
          full.names = TRUE
        )

        sample_files <- unique(
          c(
            sample_files,
            nested_files
          )
        )
      }

      # Find 10x matrix files

      filtered_files <- sample_files[
        grepl(
          "filtered",
          sample_files,
          ignore.case = TRUE
        )
      ]

      if (
        any(
          grepl(
            "matrix.*\\.mtx(\\.gz)?$",
            filtered_files,
            ignore.case = TRUE
          )
        )
      ) {

        matrix_candidates <-
          filtered_files

      } else {

        matrix_candidates <-
          sample_files
      }

      matrix_file <- choose_one(
        matrix_candidates,
        "matrix.*\\.mtx(\\.gz)?$",
        "matrix.mtx",
        ss
      )

      feature_file <- choose_one(
        matrix_candidates,
        "(features|genes).*\\.tsv(\\.gz)?$",
        "features.tsv",
        ss
      )

      barcode_file <- choose_one(
        matrix_candidates,
        "barcodes.*\\.tsv(\\.gz)?$",
        "barcodes.tsv",
        ss
      )

      position_file <- choose_one(
        sample_files,
        "tissue_positions(_list)?\\.csv(\\.gz)?$",
        "tissue_positions.csv",
        ss
      )

      manifest_list[[ss]] <- data.frame(
        section = ss,
        GSM = gsm,
        matrix = matrix_file,
        features = feature_file,
        barcodes = barcode_file,
        positions = position_file,
        stringsAsFactors = FALSE
      )

      cat(
        "Matrix: ",
        basename(matrix_file),
        "\n",
        sep = ""
      )

      cat(
        "Features: ",
        basename(feature_file),
        "\n",
        sep = ""
      )

      cat(
        "Barcodes: ",
        basename(barcode_file),
        "\n",
        sep = ""
      )

      cat(
        "Positions: ",
        basename(position_file),
        "\n",
        sep = ""
      )

      # Read matrix

      features <- read_tsv_auto(
        feature_file
      )

      barcodes <- as.character(
        read_tsv_auto(
          barcode_file
        )[[1]]
      )

      mat <- read_matrix_auto(
        matrix_file
      )

      if (
        nrow(mat) !=
        nrow(features)
      ) {

        stop(
          ss,
          ": feature/matrix row mismatch."
        )
      }

      if (
        ncol(mat) !=
        length(barcodes)
      ) {

        stop(
          ss,
          ": barcode/matrix mismatch."
        )
      }

      if (anyDuplicated(barcodes)) {

        stop(
          ss,
          ": duplicated barcodes."
        )
      }

      if (any(!is.finite(mat@x)) || any(mat@x < 0) || any(mat@x != floor(mat@x))) {
        stop(ss, ": spatial counts must be finite nonnegative integers.")
      }
      colnames(mat) <-
        barcodes

      # Feature ids / symbols

      fid <- as.character(
        features[[1]]
      )

      symbol <- if (
        ncol(features) >= 2
      ) {

        as.character(
          features[[2]]
        )

      } else {

        fid
      }

      missing_symbol <-
        is.na(symbol) |
        !nzchar(symbol)

      symbol[
        missing_symbol
      ] <- fid[
        missing_symbol
      ]

      # Identify host features

      host_mask <- grepl(
        "^ENSG[0-9]+(\\.[0-9]+)?$",
        fid
      ) |
        symbol %in%
        host_symbols

      probe_rows <- grepl("_probe[0-9]+$", symbol, ignore.case=TRUE)
      host_mask <- host_mask & !probe_rows

      host_valid <-
        host_mask &
        !is.na(symbol) &
        nzchar(symbol)

      cat(
        "Host features identified:",
        sum(host_valid),
        "\n"
      )

      if (sum(host_valid) < 1000) {

        stop(
          ss,
          ": too few human host features detected."
        )
      }

      # Identify bacterial features

      bacterial_mask <-
        classify_bacterial_features(
          features,
          host_mask
        )

      if (!any(bacterial_mask)) stop(ss, ": no valid deposited bacterial probes found; do not reuse earlier bacterial results.")

      probe_indices <- which(grepl("_probe[0-9]+$", symbol, ignore.case=TRUE))
      probe_names <- trimws(symbol[probe_indices])
      probe_groups <- probe_group(probe_names)
      included <- valid_probe(probe_names)
      viral <- grepl("^(HAdV|HBV|HCV|HPV|HSV|HTLV|EBV|VZV)([-_]|$)", probe_groups, ignore.case=TRUE)
      classification <- ifelse(included, "Bacterial", ifelse(viral, "Viral_excluded",
                                                             ifelse(grepl("gapdh|control|blank|negative|positive",probe_groups,ignore.case=TRUE),
                                                                    "Control_excluded", "Unrecognized_excluded")))
      probe_audit <- data.frame(section=ss, feature_id=fid[probe_indices],
                                feature_name=probe_names, deposited_group=probe_groups,
                                classification=classification, included_in_bacterial_total=included)
      write.csv(probe_audit, file.path(outdir,paste0("Bacterial_Probe_Audit_",ss,".csv")),row.names=FALSE)
      print(as.data.frame(table(probe_audit$classification)),row.names=FALSE)
      if(any(classification=="Unrecognized_excluded")) {
        cat("Unrecognized probes excluded; inspect the audit before reporting completeness.\n")
        print(probe_audit[classification=="Unrecognized_excluded",],row.names=FALSE)
      }

      bacterial_taxon <- rep(
        NA_character_,
        nrow(features)
      )

      if (any(
        bacterial_mask
      )) {

        bacterial_taxon[
          bacterial_mask
        ] <- extract_taxon_label(
          features[
            bacterial_mask,
            ,
            drop = FALSE
          ]
        )
      }

      cat(
        "Bacterial candidate features:",
        sum(bacterial_mask),
        "\n"
      )

      # Collapse host counts by gene symbol

      host_symbols_this <-
        symbol[
          host_valid
        ]

      unique_host_symbols <-
        unique(
          host_symbols_this
        )

      collapse_host <-
        Matrix::sparseMatrix(

          i = match(
            host_symbols_this,
            unique_host_symbols
          ),

          j = seq_along(
            host_symbols_this
          ),

          x = 1,

          dims = c(
            length(
              unique_host_symbols
            ),
            length(
              host_symbols_this
            )
          )
        )

      host_counts <- as(
        collapse_host %*%
          mat[
            host_valid,
            ,
            drop = FALSE
          ],
        "dgCMatrix"
      )

      dimnames(
        host_counts
      ) <- list(
        unique_host_symbols,
        barcodes
      )

      # Collapse bacterial counts

      if (any(
        bacterial_mask
      )) {

        bt <- sanitize_taxon(
          bacterial_taxon[
            bacterial_mask
          ]
        )

        unique_bt <- unique(
          bt
        )

        collapse_bacteria <-
          Matrix::sparseMatrix(

            i = match(
              bt,
              unique_bt
            ),

            j = seq_along(bt),

            x = 1,

            dims = c(
              length(unique_bt),
              length(bt)
            )
          )

        bacterial_counts <- as(
          collapse_bacteria %*%
            mat[
              bacterial_mask,
              ,
              drop = FALSE
            ],
          "dgCMatrix"
        )

        dimnames(
          bacterial_counts
        ) <- list(
          unique_bt,
          barcodes
        )

      } else {

        bacterial_counts <- Matrix::Matrix(
          0,
          nrow = 0,
          ncol = length(barcodes),
          sparse = TRUE
        )

        colnames(
          bacterial_counts
        ) <- barcodes
      }

      # Read spatial positions

      pc <- if (
        grepl(
          "\\.gz$",
          position_file,
          ignore.case = TRUE
        )
      ) {

        gzfile(
          position_file,
          "rt"
        )

      } else {

        file(
          position_file,
          "rt"
        )
      }

      first_line <- readLines(
        pc,
        n = 1
      )

      close(pc)

      has_header <- grepl(
        "barcode",
        first_line,
        ignore.case = TRUE
      )

      pc <- if (
        grepl(
          "\\.gz$",
          position_file,
          ignore.case = TRUE
        )
      ) {

        gzfile(
          position_file,
          "rt"
        )

      } else {

        file(
          position_file,
          "rt"
        )
      }

      pos <- read.csv(
        pc,
        header = has_header,
        stringsAsFactors = FALSE
      )

      close(pc)

      if (!has_header) {

        if (ncol(pos) != 6) {

          stop(
            ss,
            ": unexpected tissue positions format."
          )
        }

        names(pos) <- c(
          "barcode",
          "in_tissue",
          "array_row",
          "array_col",
          "pxl_row_in_fullres",
          "pxl_col_in_fullres"
        )
      }

      required_positions <- c(
        "barcode",
        "in_tissue",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres"
      )

      if (
        !all(
          required_positions %in%
          names(pos)
        )
      ) {

        stop(
          ss,
          ": required spatial coordinate columns missing."
        )
      }

      if (anyDuplicated(pos$barcode)) stop(ss, ": duplicated coordinate barcodes.")
      if (anyNA(pos[, required_positions]) ||
          any(!is.finite(pos$pxl_col_in_fullres)) ||
          any(!is.finite(pos$pxl_row_in_fullres))) {
        stop(ss, ": invalid spatial positions.")
      }
      idx <- match(
        colnames(host_counts),
        pos$barcode
      )

      if (anyNA(idx)) {

        stop(
          ss,
          ": some barcodes lack spatial coordinates."
        )
      }

      pos <- pos[
        idx,
        ,
        drop = FALSE
      ]

      # Qc

      host_UMI <-
        Matrix::colSums(
          host_counts
        )

      host_genes <-
        Matrix::colSums(
          host_counts > 0
        )

      keep <-
        pos$in_tissue == 1 &
        host_UMI >= 100 &
        host_genes >= 50

      if (!any(
        keep
      )) {

        stop(
          ss,
          ": no spots passed QC."
        )
      }

      cat(
        "Imported spots:",
        ncol(host_counts),
        "\n"
      )

      cat(
        "Retained spots:",
        sum(keep),
        "\n"
      )

      qc_list[[ss]] <- data.frame(

        section = ss,
        patient = patient_i,
        tissue = tissue_i,

        imported_spots =
          ncol(host_counts),

        retained_spots =
          sum(keep),

        host_features =
          nrow(host_counts),

        bacterial_features =
          nrow(bacterial_counts),

        median_host_UMI =
          median(
            host_UMI[keep]
          ),

        median_host_genes =
          median(
            host_genes[keep]
          ),

        stringsAsFactors = FALSE
      )

      # Subset retained spots

      host_keep <-
        host_counts[
          ,
          keep,
          drop = FALSE
        ]

      bacterial_keep <-
        bacterial_counts[
          ,
          keep,
          drop = FALSE
        ]

      X <- pos$pxl_col_in_fullres[
        keep
      ]

      Y <- pos$pxl_row_in_fullres[
        keep
      ]

      # Host normalization

      lib_size <-
        Matrix::colSums(
          host_keep
        )

      host_norm <- host_keep %*%
        Matrix::Diagonal(
          x =
            10000 /
            pmax(
              lib_size,
              1
            )
        )

      host_norm@x <-
        log1p(
          host_norm@x
        )

      dimnames(
        host_norm
      ) <- dimnames(
        host_keep
      )

      host_norm_list[[ss]] <-
        host_norm

      # Save host seurat object

      spatial_obj <-
        CreateSeuratObject(
          counts = host_keep,
          assay = "Spatial",
          min.cells = 0,
          min.features = 0
        )

      spatial_obj$section <-
        ss

      spatial_obj$patient_id <-
        patient_i

      spatial_obj$group <-
        tissue_i

      spatial_obj$pair_status <-
        "Matched"

      spatial_obj$pixel_x <-
        X

      spatial_obj$pixel_y <-
        Y

      spatial_obj <- NormalizeData(
        spatial_obj,
        assay = "Spatial",
        normalization.method =
          "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
      )

      saveRDS(
        spatial_obj,
        out(
          paste0(
            "Host_",
            ss,
            ".rds"
          )
        )
      )

      # Ml gene availability

      genes_present <- intersect(
        ml_genes,
        rownames(host_norm)
      )

      genes_missing <- setdiff(
        ml_genes,
        rownames(host_norm)
      )

      availability_list[[ss]] <-
        data.frame(

          section = ss,
          patient = patient_i,
          tissue = tissue_i,

          genes_requested =
            length(ml_genes),

          genes_present =
            length(genes_present),

          missing =
            paste(
              genes_missing,
              collapse = ";"
            ),

          stringsAsFactors = FALSE
        )

      cat(
        "ML genes available:",
        length(genes_present),
        "/",
        length(ml_genes),
        "\n"
      )

      # Extract ml genes

      for (
        gene in genes_present
      ) {

        spot_expression_list[[
          length(
            spot_expression_list
          ) + 1L
        ]] <- data.frame(

          section = ss,
          patient = patient_i,
          tissue = tissue_i,

          barcode =
            colnames(
              host_norm
            ),

          X = X,
          Y = Y,

          feature = gene,

          feature_type =
            "Gene",

          expression =
            as.numeric(
              host_norm[
                gene,
              ]
            ),

          stringsAsFactors = FALSE
        )
      }

      # Ml signature scores

      for (
        signature_name in names(
          signatures
        )
      ) {

        signature_genes <- intersect(
          signatures[[
            signature_name
          ]],
          rownames(
            host_norm
          )
        )

        if (!length(
          signature_genes
        )) {

          next
        }

        signature_score <-
          Matrix::colMeans(
            host_norm[
              signature_genes,
              ,
              drop = FALSE
            ]
          )

        spot_expression_list[[
          length(
            spot_expression_list
          ) + 1L
        ]] <- data.frame(

          section = ss,
          patient = patient_i,
          tissue = tissue_i,

          barcode =
            colnames(host_norm),

          X = X,
          Y = Y,

          feature =
            paste0(
              "Signature__",
              signature_name
            ),

          feature_type =
            "Signature",

          expression =
            as.numeric(
              signature_score
            ),

          stringsAsFactors = FALSE
        )
      }

      # Bacterial signal

      if (
        nrow(
          bacterial_keep
        ) > 0
      ) {

        total_bacteria <-
          Matrix::colSums(
            bacterial_keep
          )

        detected_taxa <-
          Matrix::colSums(
            bacterial_keep > 0
          )

        bacterial_spot_list[[ss]] <-
          data.frame(

            section = ss,
            patient = patient_i,
            tissue = tissue_i,
            pair_status = "Matched",

            barcode =
              colnames(
                bacterial_keep
              ),

            X = X,
            Y = Y,

            bacterial_total_raw =
              total_bacteria,

            bacterial_detected_taxa =
              detected_taxa,

            bacterial_positive =
              total_bacteria > 0,

            stringsAsFactors = FALSE
          )

        bacterial_taxa_list[[ss]] <-
          data.frame(

            section = ss,
            patient = patient_i,
            tissue = tissue_i,

            taxon =
              rownames(
                bacterial_keep
              ),

            total_counts =
              as.numeric(
                Matrix::rowSums(
                  bacterial_keep
                )
              ),

            detected_spots =
              as.numeric(
                Matrix::rowSums(
                  bacterial_keep > 0
                )
              ),

            stringsAsFactors = FALSE
          )

      } else {

        bacterial_spot_list[[ss]] <-
          data.frame(

            section = ss,
            patient = patient_i,
            tissue = tissue_i,
            pair_status = "Matched",

            barcode =
              colnames(
                host_keep
              ),

            X = X,
            Y = Y,

            bacterial_total_raw = 0,
            bacterial_detected_taxa = 0,
            bacterial_positive = FALSE,

            stringsAsFactors = FALSE
          )
      }

      # Rctd

      common_genes <- intersect(
        rownames(
          host_keep
        ),
        rownames(
          reference_counts
        )
      )

      cat(
        "RCTD shared genes:",
        length(common_genes),
        "\n"
      )

      if (length(
        common_genes
      ) < 500) {

        warning(
          ss,
          ": low RCTD gene overlap."
        )
      }

      query_counts <-
        host_keep[
          common_genes,
          ,
          drop = FALSE
        ]

      query_nUMI <-
        Matrix::colSums(
          query_counts
        )

      valid_query <-
        query_nUMI > 0

      query_counts <-
        query_counts[
          ,
          valid_query,
          drop = FALSE
        ]

      query_nUMI <-
        query_nUMI[
          valid_query
        ]

      query_barcodes <-
        colnames(
          query_counts
        )

      coord_table <- data.frame(

        x =
          X[
            match(
              query_barcodes,
              colnames(
                host_keep
              )
            )
          ],

        y =
          Y[
            match(
              query_barcodes,
              colnames(
                host_keep
              )
            )
          ],

        row.names =
          query_barcodes
      )

      spatial_rna <-
        spacexr::SpatialRNA(

          coord_table,

          query_counts,

          query_nUMI
        )

      cat(
        "Running RCTD for ",
        ss,
        "...\n",
        sep = ""
      )

      rctd_result <- try(
        {

          fit <- spacexr::create.RCTD(

            spatial_rna,

            rctd_reference,

            max_cores = config$max_cores
          )

          fit <- spacexr::run.RCTD(

            fit,

            doublet_mode =
              "full"
          )

          fit
        },
        silent = TRUE
      )

      if (inherits(
        rctd_result,
        "try-error"
      )) {

        warning(
          "RCTD failed for ",
          ss,
          ": ",
          as.character(
            rctd_result
          )
        )

        rctd_status_list[[ss]] <-
          data.frame(

            section = ss,

            status = "FAILED",

            message =
              as.character(
                rctd_result
              ),

            stringsAsFactors = FALSE
          )

      } else {

        weights <-
          rctd_result@results$weights

        if (is.null(
          weights
        )) {

          warning(
            ss,
            ": RCTD returned no weights."
          )

          rctd_status_list[[ss]] <-
            data.frame(

              section = ss,

              status =
                "NO_WEIGHTS",

              message =
                "RCTD returned no weights.",

              stringsAsFactors = FALSE
            )

        } else {

          weights <-
            as.matrix(
              weights
            )

          row_overlap <- if (
            is.null(
              rownames(weights)
            )
          ) {

            0

          } else {

            sum(
              query_barcodes %in%
                rownames(weights)
            )
          }

          col_overlap <- if (
            is.null(
              colnames(weights)
            )
          ) {

            0

          } else {

            sum(
              query_barcodes %in%
                colnames(weights)
            )
          }

          if (
            col_overlap >
            row_overlap
          ) {

            weights <-
              t(weights)
          }

          if (is.null(
            rownames(weights)
          )) {

            if (
              nrow(weights) ==
              length(query_barcodes)
            ) {

              rownames(
                weights
              ) <- query_barcodes

            } else {

              stop(
                ss,
                ": cannot align RCTD weights."
              )
            }
          }

          retained_barcodes <-
            query_barcodes[
              query_barcodes %in%
                rownames(
                  weights
                )
            ]

          if (!length(
            retained_barcodes
          )) {

            stop(
              ss,
              ": zero overlap between RCTD weights and query barcodes."
            )
          }

          weights <-
            weights[
              retained_barcodes,
              ,
              drop = FALSE
            ]

          weight_sum <-
            rowSums(
              weights
            )

          weights_norm <-
            weights

          good_weight <-
            is.finite(weight_sum) &
            weight_sum > 0

          weights_norm[
            good_weight,
          ] <-
            weights[
              good_weight,
              ,
              drop = FALSE
            ] /
            weight_sum[
              good_weight
            ]

          weights_norm[
            !good_weight,
          ] <- 0

          coord_used <-
            coord_table[
              retained_barcodes,
              ,
              drop = FALSE
            ]

          long_parts <-
            lapply(

              colnames(
                weights_norm
              ),

              function(cell_type) {

                data.frame(

                  section = ss,
                  patient = patient_i,
                  tissue = tissue_i,
                  pair_status =
                    "Matched",

                  barcode =
                    retained_barcodes,

                  X =
                    coord_used[
                      retained_barcodes,
                      "x"
                    ],

                  Y =
                    coord_used[
                      retained_barcodes,
                      "y"
                    ],

                  cell_type =
                    cell_type,

                  proportion =
                    as.numeric(
                      weights_norm[
                        ,
                        cell_type
                      ]
                    ),

                  stringsAsFactors =
                    FALSE
                )
              }
            )

          rctd_long_list[[ss]] <-
            do.call(
              rbind,
              long_parts
            )

          rctd_status_list[[ss]] <-
            data.frame(

              section = ss,

              status = "OK",

              query_spots =
                length(
                  query_barcodes
                ),

              retained_spots =
                length(
                  retained_barcodes
                ),

              retention_percent =
                100 *
                length(
                  retained_barcodes
                ) /
                length(
                  query_barcodes
                ),

              message = "",

              stringsAsFactors = FALSE
            )

          write.csv(
            data.frame(
              barcode =
                rownames(
                  weights_norm
                ),
              weights_norm,
              check.names = FALSE
            ),
            out(
              paste0(
                "RCTD_weights_",
                ss,
                ".csv"
              )
            ),
            row.names = FALSE
          )
        }
      }

      rm(
        mat,
        host_counts,
        bacterial_counts,
        host_keep,
        bacterial_keep,
        spatial_obj
      )

      invisible(
        gc()
      )
    }

    # Save primary matched outputs

    manifest <- do.call(
      rbind,
      manifest_list
    )

    write.csv(
      manifest,
      out("03_Input_manifest_MATCHED.csv"),
      row.names = FALSE
    )

    qc_summary <- do.call(
      rbind,
      qc_list
    )

    write.csv(
      qc_summary,
      out("04_QC_summary_MATCHED.csv"),
      row.names = FALSE
    )

    availability <- do.call(
      rbind,
      availability_list
    )

    write.csv(
      availability,
      out("05_ML_gene_availability_MATCHED.csv"),
      row.names = FALSE
    )

    bacterial <- do.call(
      rbind,
      bacterial_spot_list
    )

    rownames(
      bacterial
    ) <- NULL

    bacterial$log_bacterial <-
      log1p(
        bacterial$bacterial_total_raw
      )

    write.csv(
      bacterial,
      out("06_Bacterial_signal_All_spots_MATCHED.csv"),
      row.names = FALSE
    )

    valid_taxa <- bacterial_taxa_list[
      vapply(
        bacterial_taxa_list,
        function(x) {
          !is.null(x) &&
            nrow(x) > 0
        },
        logical(1)
      )
    ]

    if (length(valid_taxa)) {

      bacterial_taxa <-
        do.call(
          rbind,
          valid_taxa
        )

    } else {

      bacterial_taxa <-
        data.frame()
    }

    write.csv(
      bacterial_taxa,
      out("07_Bacterial_taxa_MATCHED.csv"),
      row.names = FALSE
    )

    # Save rctd output

    status_columns <- unique(unlist(lapply(rctd_status_list, names), use.names = FALSE))
    rctd_status <- do.call(rbind, lapply(rctd_status_list, function(z) {
      for (nm in setdiff(status_columns, names(z))) z[[nm]] <- NA
      z[, status_columns, drop = FALSE]
    }))

    write.csv(
      rctd_status,
      out("08_RCTD_status_MATCHED.csv"),
      row.names = FALSE
    )

    if (nrow(rctd_status) != 4L || any(rctd_status$status != "OK")) {
      stop("All four sections must complete RCTD before proceeding; inspect 08_RCTD_status_MATCHED.csv.")
    }

    successful_rctd <- rctd_long_list[
      vapply(
        rctd_long_list,
        function(x) {
          !is.null(x) &&
            nrow(x) > 0
        },
        logical(1)
      )
    ]

    if (!length(
      successful_rctd
    )) {

      stop(
        paste0(
          "RCTD failed for all four sections. ",
          "Inspect 08_RCTD_status_MATCHED.csv."
        )
      )
    }

    rctd <- do.call(
      rbind,
      successful_rctd
    )

    rownames(
      rctd
    ) <- NULL

    write.csv(
      rctd,
      out("09_RCTD_All_spot_celltype_proportions_MATCHED.csv"),
      row.names = FALSE
    )

    # Combine ml features

    spot_expr <- do.call(
      rbind,
      spot_expression_list
    )

    rownames(
      spot_expr
    ) <- NULL

    if (
      !all(
        spot_expr$patient %in%
        matched_patients
      )
    ) {

      stop(
        "Non-matched patient found in feature table."
      )
    }

    write.csv(
      spot_expr,
      out("10_ML_features_All_spots_MATCHED.csv"),
      row.names = FALSE
    )

    # Add bacterial signal to ml features

    spot_expr_bact <- merge(

      spot_expr,

      bacterial[
        ,
        c(
          "section",
          "patient",
          "tissue",
          "barcode",
          "bacterial_total_raw",
          "bacterial_detected_taxa",
          "bacterial_positive",
          "log_bacterial"
        ),
        drop = FALSE
      ],

      by = c(
        "section",
        "patient",
        "tissue",
        "barcode"
      ),

      all.x = TRUE,

      sort = FALSE
    )

    spot_expr_bact$bacterial_total_raw[
      is.na(
        spot_expr_bact$bacterial_total_raw
      )
    ] <- 0

    spot_expr_bact$log_bacterial[
      is.na(
        spot_expr_bact$log_bacterial
      )
    ] <- 0

    write.csv(
      spot_expr_bact,
      out("11_ML_features_with_bacteria_MATCHED.csv"),
      row.names = FALSE
    )

    section_feature_groups <- split(

      spot_expr,

      interaction(
        spot_expr$section,
        spot_expr$feature,
        drop = TRUE
      )
    )

    section_summary <- do.call(

      rbind,

      lapply(
        section_feature_groups,
        function(z) {

          data.frame(

            section =
              unique(
                z$section
              ),

            patient =
              unique(
                z$patient
              ),

            tissue =
              unique(
                z$tissue
              ),

            feature =
              unique(
                z$feature
              ),

            feature_type =
              unique(
                z$feature_type
              ),

            mean_expression =
              mean(
                z$expression,
                na.rm = TRUE
              ),

            median_expression =
              median(
                z$expression,
                na.rm = TRUE
              ),

            detected_percent =
              100 *
              mean(
                z$expression > 0,
                na.rm = TRUE
              ),

            n_spots =
              nrow(z),

            stringsAsFactors = FALSE
          )
        }
      )
    )

    rownames(
      section_summary
    ) <- NULL

    write.csv(
      section_summary,
      out("12_Feature_section_summary_MATCHED.csv"),
      row.names = FALSE
    )

    # True matched normal vs tumor

    normal <- section_summary[
      section_summary$tissue ==
        "Normal",
      ,
      drop = FALSE
    ]

    tumor <- section_summary[
      section_summary$tissue ==
        "Tumor",
      ,
      drop = FALSE
    ]

    normal <- normal[
      ,
      c(
        "patient",
        "feature",
        "feature_type",
        "mean_expression",
        "median_expression",
        "detected_percent"
      )
    ]

    tumor <- tumor[
      ,
      c(
        "patient",
        "feature",
        "feature_type",
        "mean_expression",
        "median_expression",
        "detected_percent"
      )
    ]

    names(normal)[4:6] <- c(
      "Normal_mean",
      "Normal_median",
      "Normal_detected_percent"
    )

    names(tumor)[4:6] <- c(
      "Tumor_mean",
      "Tumor_median",
      "Tumor_detected_percent"
    )

    paired <- merge(

      normal,

      tumor,

      by = c(
        "patient",
        "feature",
        "feature_type"
      ),

      all = FALSE
    )

    paired$Delta_mean <-
      paired$Tumor_mean -
      paired$Normal_mean

    paired$Delta_median <-
      paired$Tumor_median -
      paired$Normal_median

    paired$Observed_ST_direction <-
      ifelse(

        paired$Delta_mean > 0,

        "Higher_in_Tumor",

        ifelse(
          paired$Delta_mean < 0,
          "Lower_in_Tumor",
          "No_change"
        )
      )

    paired$Expected_ML_direction <-
      NA_character_

    gene_rows <-
      paired$feature_type ==
      "Gene"

    signature_rows <-
      paired$feature_type ==
      "Signature"

    paired$Expected_ML_direction[
      gene_rows
    ] <- unname(
      ml_direction[
        paired$feature[
          gene_rows
        ]
      ]
    )

    paired$Expected_ML_direction[
      signature_rows
    ] <- unname(
      signature_direction[
        paired$feature[
          signature_rows
        ]
      ]
    )

    paired$Direction_agreement <-
      paired$Expected_ML_direction ==
      paired$Observed_ST_direction

    write.csv(
      paired,
      out("13_Matched_pair_ML_ST_validation.csv"),
      row.names = FALSE
    )

    # Consistency across both patients

    pair_consistency <- do.call(

      rbind,

      lapply(
        split(
          paired,
          paired$feature
        ),
        function(z) {

          crc03 <-
            z[
              z$patient ==
                "CRC_03",
              ,
              drop = FALSE
            ]

          crc08 <-
            z[
              z$patient ==
                "CRC_08",
              ,
              drop = FALSE
            ]

          data.frame(

            feature =
              unique(
                z$feature
              ),

            feature_type =
              unique(
                z$feature_type
              ),

            Expected_ML_direction =
              unique(
                z$Expected_ML_direction
              )[1],

            CRC03_direction =
              if (nrow(crc03)) {
                crc03$Observed_ST_direction[1]
              } else {
                NA_character_
              },

            CRC08_direction =
              if (nrow(crc08)) {
                crc08$Observed_ST_direction[1]
              } else {
                NA_character_
              },

            Same_direction_in_both_pairs =
              length(
                unique(
                  z$Observed_ST_direction
                )
              ) == 1,

            ML_agreement_pairs =
              sum(
                z$Direction_agreement,
                na.rm = TRUE
              ),

            Mean_delta =
              mean(
                z$Delta_mean,
                na.rm = TRUE
              ),

            stringsAsFactors = FALSE
          )
        }
      )
    )

    rownames(
      pair_consistency
    ) <- NULL

    write.csv(
      pair_consistency,
      out("14_Matched_pair_direction_consistency.csv"),
      row.names = FALSE
    )

    # Feature ↔ total bacterial load

    feature_bacteria_groups <- split(

      spot_expr_bact,

      interaction(
        spot_expr_bact$section,
        spot_expr_bact$feature,
        drop = TRUE
      )
    )

    feature_bacteria_cor <- do.call(

      rbind,

      lapply(
        feature_bacteria_groups,
        function(z) {

          cc <- safe_cor(
            z$expression,
            z$log_bacterial
          )

          data.frame(

            section =
              unique(
                z$section
              ),

            patient =
              unique(
                z$patient
              ),

            tissue =
              unique(
                z$tissue
              ),

            feature =
              unique(
                z$feature
              ),

            feature_type =
              unique(
                z$feature_type
              ),

            spearman_rho =
              unname(
                cc["rho"]
              ),

            p_value =
              unname(
                cc["p"]
              ),

            n_spots =
              unname(
                cc["n"]
              ),

            stringsAsFactors = FALSE
          )
        }
      )
    )

    feature_bacteria_cor$FDR <-
      p.adjust(
        feature_bacteria_cor$p_value,
        method = "BH"
      )

    feature_bacteria_cor <-
      feature_bacteria_cor[
        order(
          abs(
            feature_bacteria_cor$spearman_rho
          ),
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ]

    write.csv(
      feature_bacteria_cor,
      out("15_Feature_TotalBacteria_correlations_MATCHED.csv"),
      row.names = FALSE
    )

    # Feature ↔ rctd cell-type proportions

    rctd_feature_data <- merge(

      spot_expr,

      rctd[
        ,
        c(
          "section",
          "barcode",
          "cell_type",
          "proportion"
        ),
        drop = FALSE
      ],

      by = c(
        "section",
        "barcode"
      ),

      all = FALSE,

      sort = FALSE
    )

    if (!nrow(
      rctd_feature_data
    )) {

      stop(
        "No overlap between ML-feature spots and RCTD."
      )
    }

    feature_rctd_groups <- split(

      rctd_feature_data,

      interaction(
        rctd_feature_data$section,
        rctd_feature_data$feature,
        rctd_feature_data$cell_type,
        drop = TRUE
      )
    )

    feature_rctd_cor <- do.call(

      rbind,

      lapply(
        feature_rctd_groups,
        function(z) {

          cc <- safe_cor(
            z$expression,
            z$proportion
          )

          data.frame(

            section =
              unique(
                z$section
              ),

            patient =
              unique(
                z$patient
              ),

            tissue =
              unique(
                z$tissue
              ),

            feature =
              unique(
                z$feature
              ),

            feature_type =
              unique(
                z$feature_type
              ),

            cell_type =
              unique(
                z$cell_type
              ),

            spearman_rho =
              unname(
                cc["rho"]
              ),

            p_value =
              unname(
                cc["p"]
              ),

            n_spots =
              unname(
                cc["n"]
              ),

            stringsAsFactors = FALSE
          )
        }
      )
    )

    feature_rctd_cor$FDR <-
      p.adjust(
        feature_rctd_cor$p_value,
        method = "BH"
      )

    feature_rctd_cor <-
      feature_rctd_cor[
        order(
          abs(
            feature_rctd_cor$spearman_rho
          ),
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ]

    write.csv(
      feature_rctd_cor,
      out("16_Feature_RCTD_correlations_MATCHED.csv"),
      row.names = FALSE
    )

    # Bacteria ↔ rctd

    bacteria_rctd_data <- merge(

      bacterial[
        ,
        c(
          "section",
          "patient",
          "tissue",
          "barcode",
          "log_bacterial"
        ),
        drop = FALSE
      ],

      rctd[
        ,
        c(
          "section",
          "barcode",
          "cell_type",
          "proportion"
        ),
        drop = FALSE
      ],

      by = c(
        "section",
        "barcode"
      ),

      all = FALSE,

      sort = FALSE
    )

    bacteria_rctd_groups <- split(

      bacteria_rctd_data,

      interaction(
        bacteria_rctd_data$section,
        bacteria_rctd_data$cell_type,
        drop = TRUE
      )
    )

    bacteria_rctd_cor <- do.call(

      rbind,

      lapply(
        bacteria_rctd_groups,
        function(z) {

          cc <- safe_cor(
            z$log_bacterial,
            z$proportion
          )

          data.frame(

            section =
              unique(
                z$section
              ),

            patient =
              unique(
                z$patient
              ),

            tissue =
              unique(
                z$tissue
              ),

            cell_type =
              unique(
                z$cell_type
              ),

            spearman_rho =
              unname(
                cc["rho"]
              ),

            p_value =
              unname(
                cc["p"]
              ),

            n_spots =
              unname(
                cc["n"]
              ),

            stringsAsFactors = FALSE
          )
        }
      )
    )

    bacteria_rctd_cor$FDR <-
      p.adjust(
        bacteria_rctd_cor$p_value,
        method = "BH"
      )

    bacteria_rctd_cor <-
      bacteria_rctd_cor[
        order(
          abs(
            bacteria_rctd_cor$spearman_rho
          ),
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ]

    write.csv(
      bacteria_rctd_cor,
      out("17_Bacteria_RCTD_correlations_MATCHED.csv"),
      row.names = FALSE
    )

    for (
      feature_i in unique(
        spot_expr$feature
      )
    ) {

      d <- spot_expr[
        spot_expr$feature == feature_i,
        ,
        drop = FALSE
      ]

      d$section <- factor(
        d$section,
        levels = matched_sections
      )

      p <- ggplot(
        d,
        aes(
          x = X,
          y = Y,
          color = expression
        )
      ) +

        geom_point(
          size = 0.45,
          alpha = 0.9,
          stroke = 0
        ) +

        scale_y_reverse() +

        coord_equal() +

        facet_wrap(
          ~ section,
          nrow = 1,
          scales = "fixed"
        ) +

        labs(
          title = feature_i,
          x = NULL,
          y = NULL,
          color = "Expression"
        ) +

        theme_void(
          base_family = "Arial"
        ) +

        theme(
          plot.title =
            element_text(
              face = "bold",
              hjust = 0.5
            ),

          strip.text =
            element_text(
              face = "bold"
            )
        )

      safe_name <- gsub(
        "[^A-Za-z0-9_]+",
        "_",
        feature_i
      )

      ggsave(
        filename =
          out(
            paste0(
              "Spatial_",
              safe_name,
              "_MATCHED.png"
            )
          ),
        plot = p,
        width = 12,
        height = 4,
        units = "in",
        dpi = 600,
        bg = "white"
      )
    }

    bacterial$section <- factor(
      bacterial$section,
      levels = matched_sections
    )

    p_bacteria <- ggplot(
      bacterial,
      aes(
        x = X,
        y = Y,
        color = log_bacterial
      )
    ) +

      geom_point(
        size = 0.45,
        alpha = 0.9,
        stroke = 0
      ) +

      scale_y_reverse() +

      coord_equal() +

      facet_wrap(
        ~ section,
        nrow = 1,
        scales = "fixed"
      ) +

      labs(
        title = "Total bacterial signal",
        x = NULL,
        y = NULL,
        color = "log1p bacterial counts"
      ) +

      theme_void(
        base_family = "Arial"
      ) +

      theme(
        plot.title =
          element_text(
            face = "bold",
            hjust = 0.5
          ),

        strip.text =
          element_text(
            face = "bold"
          )
      )

    ggsave(
      filename =
        out(
          "18_Bacterial_signal_spatial_MATCHED.png"
        ),
      plot = p_bacteria,
      width = 12,
      height = 4,
      units = "in",
      dpi = 600,
      bg = "white"
    )

    # Final consensus table

    bacteria_feature_summary <- do.call(

      rbind,

      lapply(
        split(
          feature_bacteria_cor,
          feature_bacteria_cor$feature
        ),

        function(z) {

          valid_rho <- z$spearman_rho[
            is.finite(
              z$spearman_rho
            )
          ]

          data.frame(

            feature =
              unique(
                z$feature
              ),

            median_bacteria_rho =
              if (length(valid_rho)) {
                median(valid_rho)
              } else {
                NA_real_
              },

            max_abs_bacteria_rho =
              if (length(valid_rho)) {
                max(
                  abs(
                    valid_rho
                  )
                )
              } else {
                NA_real_
              },

            significant_bacteria_sections =
              sum(
                z$FDR < 0.05,
                na.rm = TRUE
              ),

            stringsAsFactors = FALSE
          )
        }
      )
    )

    rownames(
      bacteria_feature_summary
    ) <- NULL

    # Strongest rctd association per feature

    rctd_feature_summary_list <- lapply(

      split(
        feature_rctd_cor,
        feature_rctd_cor$feature
      ),

      function(z) {

        feature_name <- unique(
          z$feature
        )[1]

        z <- z[
          is.finite(
            z$spearman_rho
          ),
          ,
          drop = FALSE
        ]

        if (!nrow(z)) {

          return(
            data.frame(
              feature = feature_name,
              strongest_cell_type = NA_character_,
              strongest_RCTD_rho = NA_real_,
              stringsAsFactors = FALSE
            )
          )
        }

        med <- aggregate(
          spearman_rho ~ cell_type,
          data = z,
          FUN = function(x) {
            median(
              x,
              na.rm = TRUE
            )
          }
        )

        med <- med[
          order(
            abs(
              med$spearman_rho
            ),
            decreasing = TRUE
          ),
          ,
          drop = FALSE
        ]

        data.frame(

          feature =
            feature_name,

          strongest_cell_type =
            med$cell_type[1],

          strongest_RCTD_rho =
            med$spearman_rho[1],

          stringsAsFactors = FALSE
        )
      }
    )

    rctd_feature_summary <- do.call(
      rbind,
      rctd_feature_summary_list
    )

    rownames(
      rctd_feature_summary
    ) <- NULL

    # Final matched-pair consensus

    consensus <- merge(
      pair_consistency,
      bacteria_feature_summary,
      by = "feature",
      all.x = TRUE
    )

    consensus <- merge(
      consensus,
      rctd_feature_summary,
      by = "feature",
      all.x = TRUE
    )

    consensus$Spatial_ML_support <- ifelse(

      consensus$ML_agreement_pairs == 2,

      "Supported_in_both_pairs",

      ifelse(

        consensus$ML_agreement_pairs == 1,

        "Supported_in_one_pair",

        "Not_supported"
      )
    )

    consensus <- consensus[
      order(
        consensus$ML_agreement_pairs,
        consensus$Same_direction_in_both_pairs,
        abs(
          consensus$Mean_delta
        ),
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]

    write.csv(
      consensus,
      out(
        "19_FINAL_ML_ST_Bacteria_RCTD_Consensus_MATCHED.csv"
      ),
      row.names = FALSE
    )

    cat("\n\n")
    cat("============================================================\n")
    cat("MATCHED ST ANALYSIS COMPLETED\n")
    cat("============================================================\n")

    cat("\nPatients analyzed:\n")

    print(
      unique(
        spot_expr$patient
      )
    )

    cat("\nSections analyzed:\n")

    print(
      unique(
        spot_expr$section
      )
    )

    cat("\nQC SUMMARY:\n")

    print(
      qc_summary,
      row.names = FALSE
    )

    cat("\nML GENE AVAILABILITY:\n")

    print(
      availability,
      row.names = FALSE
    )

    cat("\nRCTD STATUS:\n")

    print(
      rctd_status,
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("MATCHED NORMAL vs TUMOR VALIDATION\n")
    cat("============================================================\n")

    print(
      paired[
        order(
          paired$feature,
          paired$patient
        ),
      ],
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("CONSISTENCY BETWEEN CRC_03 AND CRC_08\n")
    cat("============================================================\n")

    print(
      pair_consistency[
        order(
          pair_consistency$ML_agreement_pairs,
          abs(
            pair_consistency$Mean_delta
          ),
          decreasing = TRUE
        ),
      ],
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("TOP FEATURE ↔ BACTERIA ASSOCIATIONS\n")
    cat("============================================================\n")

    print(
      head(
        feature_bacteria_cor,
        30
      ),
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("TOP FEATURE ↔ RCTD ASSOCIATIONS\n")
    cat("============================================================\n")

    print(
      head(
        feature_rctd_cor,
        30
      ),
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("BACTERIA ↔ RCTD ASSOCIATIONS\n")
    cat("============================================================\n")

    print(
      head(
        bacteria_rctd_cor,
        30
      ),
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("FINAL CONSENSUS\n")
    cat("============================================================\n")

    print(
      consensus,
      row.names = FALSE
    )

    cat("\nOutput directory:\n")

    cat(
      outdir,
      "\n"
    )

    cat("\n")
    cat("============================================================\n")
    cat("IMPORTANT: ONLY CRC_03 AND CRC_08 MATCHED PAIRS WERE USED.\n")
    cat("CRC_06 AND CRC_07 WERE NOT USED.\n")
    cat("============================================================\n")

    # Session info

    writeLines(
      capture.output(
        sessionInfo()
      ),
      out(
        "20_sessionInfo.txt"
      )
    )

    cat(
      "\nSession info saved to:\n",
      out("20_sessionInfo.txt"),
      "\n",
      sep = ""
    )

    cat("\n")
    cat("============================================================\n")
    cat("ALL ANALYSES FINISHED SUCCESSFULLY\n")
    cat("============================================================\n")

    invisible(
      gc()
    )
  }
  cat("\n[01/06] Raw Spatial Data, QC, RCTD and Feature Analyses\n")
  spatial_stage_01(config)
  invisible(gc())

  # Final gene and signature panel
  spatial_stage_02 <- function(config) {

    invisible(gc())

    options(
      width = 220,
      stringsAsFactors = FALSE
    )

    # Paths

    root <- config$data_root

    outdir <- config$analysis_dir

    out <- function(x) {
      file.path(outdir, x)
    }

    if (!dir.exists(outdir)) {
      stop(
        "Output directory not found:\n",
        outdir
      )
    }

    # Input files

    spot_file <- out(
      "10_ML_features_All_spots_MATCHED.csv"
    )

    paired_file <- out(
      "13_Matched_pair_ML_ST_validation.csv"
    )

    consistency_file <- out(
      "14_Matched_pair_direction_consistency.csv"
    )

    required_files <- c(
      spot_file,
      paired_file,
      consistency_file
    )

    missing_files <- required_files[
      !file.exists(required_files)
    ]

    if (length(missing_files)) {
      stop(
        "Missing required files:\n",
        paste(
          missing_files,
          collapse = "\n"
        )
      )
    }

    # Read saved results

    cat("\nReading saved ST results...\n")

    spot_expr <- read.csv(
      spot_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    paired <- read.csv(
      paired_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    consistency <- read.csv(
      consistency_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    cat(
      "spot_expr rows:",
      nrow(spot_expr),
      "\n"
    )

    cat(
      "paired rows:",
      nrow(paired),
      "\n"
    )

    cat(
      "consistency rows:",
      nrow(consistency),
      "\n"
    )

    # Final 15 genes

    final_genes <- c(
      "CEACAM5",
      "CD68",
      "EPCAM",
      "KRT8",
      "GNLY",
      "LYZ",
      "NKG7",
      "CD79A",
      "MS4A1",
      "CD4",
      "CD3E",
      "PECAM1",
      "VWF",
      "CD8A",
      "CD3D"
    )

    # Final 6 signatures

    final_signatures <- c(
      "Signature__Epithelial",
      "Signature__Myeloid",
      "Signature__Cytotoxic",
      "Signature__B_cell",
      "Signature__T_cell",
      "Signature__Endothelial"
    )

    # Final 21 features

    final_features <- c(
      final_genes,
      final_signatures
    )

    cat("\n========================================\n")
    cat("FINAL PANEL\n")
    cat("========================================\n")

    cat(
      "Final genes:",
      length(final_genes),
      "\n"
    )

    cat(
      "Final signatures:",
      length(final_signatures),
      "\n"
    )

    cat(
      "Final features:",
      length(final_features),
      "\n"
    )

    # Check krt18

    cat("\nKRT18 check:\n")

    if ("KRT18" %in% spot_expr$feature) {

      cat(
        "WARNING: KRT18 exists in spot_expr\n"
      )

    } else {

      cat(
        "KRT18 is absent from ST data — removed from final panel\n"
      )
    }

    # Filter spot-level data

    spot_expr_final <- spot_expr[
      spot_expr$feature %in%
        final_features,
      ,
      drop = FALSE
    ]

    rownames(
      spot_expr_final
    ) <- NULL

    observed_spot_features <- sort(
      unique(
        spot_expr_final$feature
      )
    )

    cat(
      "\nUnique features in final spot data:",
      length(observed_spot_features),
      "\n"
    )

    print(
      observed_spot_features
    )

    # Check final feature set

    missing_spot_features <- setdiff(
      final_features,
      observed_spot_features
    )

    extra_spot_features <- setdiff(
      observed_spot_features,
      final_features
    )

    cat(
      "\nMissing final features:\n"
    )

    if (length(missing_spot_features)) {

      print(
        missing_spot_features
      )

    } else {

      cat(
        "None\n"
      )
    }

    cat(
      "\nUnexpected extra features:\n"
    )

    if (length(extra_spot_features)) {

      print(
        extra_spot_features
      )

    } else {

      cat(
        "None\n"
      )
    }

    if (
      length(
        observed_spot_features
      ) != 21
    ) {

      warning(
        "Final spot-level panel does not contain exactly 21 features."
      )

    } else {

      cat(
        "\nPASS: Final spot-level panel = 21 features\n"
      )
    }

    # Filter matched validation results

    paired_final <- paired[
      paired$feature %in%
        final_features,
      ,
      drop = FALSE
    ]

    rownames(
      paired_final
    ) <- NULL

    cat(
      "\nMatched validation rows:",
      nrow(paired_final),
      "\n"
    )

    cat(
      "Matched validation unique features:",
      length(
        unique(
          paired_final$feature
        )
      ),
      "\n"
    )

    # Filter consistency results

    consistency_final <- consistency[
      consistency$feature %in%
        final_features,
      ,
      drop = FALSE
    ]

    rownames(
      consistency_final
    ) <- NULL

    cat(
      "\nConsistency rows:",
      nrow(consistency_final),
      "\n"
    )

    cat(
      "Consistency unique features:",
      length(
        unique(
          consistency_final$feature
        )
      ),
      "\n"
    )

    # Define 8 priority features

    priority_features <- c(
      "CEACAM5",
      "EPCAM",
      "KRT8",
      "MS4A1",
      "VWF",
      "PECAM1",
      "Signature__Epithelial",
      "Signature__Endothelial"
    )

    # Extract priority results

    priority_results <- consistency_final[
      consistency_final$feature %in%
        priority_features,
      ,
      drop = FALSE
    ]

    priority_results <- priority_results[
      match(
        priority_features,
        priority_results$feature
      ),
      ,
      drop = FALSE
    ]

    priority_results <- priority_results[
      !is.na(
        priority_results$feature
      ),
      ,
      drop = FALSE
    ]

    rownames(
      priority_results
    ) <- NULL

    cat("\n")
    cat("========================================\n")
    cat("8 PRIORITY FEATURES\n")
    cat("========================================\n")

    print(
      priority_results,
      row.names = FALSE
    )

    # Check priority features

    missing_priority <- setdiff(
      priority_features,
      priority_results$feature
    )

    cat(
      "\nPriority features found:",
      nrow(priority_results),
      "/",
      length(priority_features),
      "\n"
    )

    if (length(missing_priority)) {

      cat(
        "\nMissing priority features:\n"
      )

      print(
        missing_priority
      )

    } else {

      cat(
        "\nPASS: All 8 priority features found\n"
      )
    }

    # Verify both-pair support

    if (
      "ML_agreement_pairs" %in%
      names(priority_results)
    ) {

      cat(
        "\nML agreement pairs:\n"
      )

      print(
        priority_results[
          ,
          c(
            "feature",
            "ML_agreement_pairs"
          ),
          drop = FALSE
        ],
        row.names = FALSE
      )

      if (
        all(
          priority_results$ML_agreement_pairs ==
          2
        )
      ) {

        cat(
          "\nPASS: All 8 features supported in both matched pairs\n"
        )

      } else {

        warning(
          "One or more priority features are not supported in both pairs."
        )
      }
    }

    # Show directions

    direction_columns <- intersect(
      c(
        "feature",
        "Expected_ML_direction",
        "CRC03_direction",
        "CRC08_direction",
        "Same_direction_in_both_pairs",
        "ML_agreement_pairs",
        "Mean_delta"
      ),
      names(
        priority_results
      )
    )

    cat(
      "\n========================================\n"
    )

    cat(
      "PRIORITY FEATURE DIRECTIONS\n"
    )

    cat(
      "========================================\n"
    )

    print(
      priority_results[
        ,
        direction_columns,
        drop = FALSE
      ],
      row.names = FALSE
    )

    cat(
      "\n========================================\n"
    )

    cat(
      "FINAL 21-FEATURE PANEL\n"
    )

    cat(
      "========================================\n"
    )

    final_panel_table <- data.frame(

      feature =
        final_features,

      type =
        c(
          rep(
            "Gene",
            length(
              final_genes
            )
          ),

          rep(
            "Signature",
            length(
              final_signatures
            )
          )
        ),

      priority =
        final_features %in%
        priority_features,

      stringsAsFactors = FALSE
    )

    print(
      final_panel_table,
      row.names = FALSE
    )

    # Save final panel

    write.csv(
      final_panel_table,
      out(
        "21_Final_ST_21_Feature_Panel.csv"
      ),
      row.names = FALSE
    )

    # Save spot-level data

    write.csv(
      spot_expr_final,
      out(
        "22_Final_21_Features_All_spots.csv"
      ),
      row.names = FALSE
    )

    # Save matched validation

    write.csv(
      paired_final,
      out(
        "23_Final_21_Features_Matched_validation.csv"
      ),
      row.names = FALSE
    )

    # Save consistency table

    write.csv(
      consistency_final,
      out(
        "24_Final_21_Features_Consistency.csv"
      ),
      row.names = FALSE
    )

    # Save 8 priority features

    write.csv(
      priority_results,
      out(
        "25_Priority_8_Features_Supported_Both_Pairs.csv"
      ),
      row.names = FALSE
    )

    # Final check

    cat("\n")
    cat("========================================\n")
    cat("FINAL CHECK\n")
    cat("========================================\n")

    cat(
      "15 genes:",
      length(final_genes),
      "\n"
    )

    cat(
      "6 signatures:",
      length(final_signatures),
      "\n"
    )

    cat(
      "21 final features:",
      length(final_features),
      "\n"
    )

    cat(
      "Observed final ST features:",
      length(
        unique(
          spot_expr_final$feature
        )
      ),
      "\n"
    )

    cat(
      "Priority features:",
      nrow(
        priority_results
      ),
      "\n"
    )

    if (
      length(
        unique(
          spot_expr_final$feature
        )
      ) == 21 &&
      nrow(
        priority_results
      ) == 8
    ) {

      cat(
        "\n",
        "========================================\n",
        "FINAL STATUS: PASS\n",
        "21-feature panel created successfully.\n",
        "8 priority features extracted successfully.\n",
        "KRT18 excluded.\n",
        "========================================\n",
        sep = ""
      )

    } else {

      cat(
        "\n",
        "========================================\n",
        "FINAL STATUS: REVIEW REQUIRED\n",
        "Check missing features above.\n",
        "========================================\n",
        sep = ""
      )
    }

    # Saved files

    cat(
      "\nSaved files:\n"
    )

    cat(
      "21_Final_ST_21_Feature_Panel.csv\n"
    )

    cat(
      "22_Final_21_Features_All_spots.csv\n"
    )

    cat(
      "23_Final_21_Features_Matched_validation.csv\n"
    )

    cat(
      "24_Final_21_Features_Consistency.csv\n"
    )

    cat(
      "25_Priority_8_Features_Supported_Both_Pairs.csv\n"
    )

    cat(
      "\nOutput directory:\n",
      outdir,
      "\n",
      sep = ""
    )

    invisible(
      gc()
    )
  }
  cat("\n[02/06] Final Gene and Signature Panel\n")
  spatial_stage_02(config)
  invisible(gc())

  # Result review and feature summaries
  spatial_stage_03 <- function(config) {

    invisible(gc())

    options(
      width = 220,
      stringsAsFactors = FALSE
    )

    # Paths

    root <- config$data_root

    outdir <- config$analysis_dir

    out <- function(x) {
      file.path(outdir, x)
    }

    if (!dir.exists(outdir)) {
      stop(
        "Output directory not found:\n",
        outdir
      )
    }

    # Input files

    panel_file <- out(
      "21_Final_ST_21_Feature_Panel.csv"
    )

    spot_file <- out(
      "22_Final_21_Features_All_spots.csv"
    )

    paired_file <- out(
      "23_Final_21_Features_Matched_validation.csv"
    )

    consistency_file <- out(
      "24_Final_21_Features_Consistency.csv"
    )

    priority_file <- out(
      "25_Priority_8_Features_Supported_Both_Pairs.csv"
    )

    required_files <- c(
      panel_file,
      spot_file,
      paired_file,
      consistency_file,
      priority_file
    )

    missing_files <- required_files[
      !file.exists(required_files)
    ]

    if (length(missing_files)) {
      stop(
        "Missing final files:\n",
        paste(
          missing_files,
          collapse = "\n"
        )
      )
    }

    # Read files

    cat("\nReading final files...\n")

    panel <- read.csv(
      panel_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    spot <- read.csv(
      spot_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    paired <- read.csv(
      paired_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    consistency <- read.csv(
      consistency_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    priority <- read.csv(
      priority_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    cat("\nLoaded files:\n")
    cat("Panel rows:", nrow(panel), "\n")
    cat("Spot rows:", nrow(spot), "\n")
    cat("Paired rows:", nrow(paired), "\n")
    cat("Consistency rows:", nrow(consistency), "\n")
    cat("Priority rows:", nrow(priority), "\n")

    # Expected features

    final_genes <- c(
      "CEACAM5",
      "CD68",
      "EPCAM",
      "KRT8",
      "GNLY",
      "LYZ",
      "NKG7",
      "CD79A",
      "MS4A1",
      "CD4",
      "CD3E",
      "PECAM1",
      "VWF",
      "CD8A",
      "CD3D"
    )

    final_signatures <- c(
      "Signature__Epithelial",
      "Signature__Myeloid",
      "Signature__Cytotoxic",
      "Signature__B_cell",
      "Signature__T_cell",
      "Signature__Endothelial"
    )

    final_features <- c(
      final_genes,
      final_signatures
    )

    priority_features <- c(
      "CEACAM5",
      "EPCAM",
      "KRT8",
      "MS4A1",
      "VWF",
      "PECAM1",
      "Signature__Epithelial",
      "Signature__Endothelial"
    )

    # Basic checks

    cat("\n")
    cat("============================================================\n")
    cat("FINAL FILE CHECKS\n")
    cat("============================================================\n")

    check_table <- data.frame(
      Check = character(0),
      Result = character(0),
      stringsAsFactors = FALSE
    )

    add_check <- function(name, condition) {

      check_table <<- rbind(
        check_table,
        data.frame(
          Check = name,
          Result = ifelse(
            isTRUE(condition),
            "PASS",
            "FAIL"
          ),
          stringsAsFactors = FALSE
        )
      )
    }

    add_check(
      "Panel has 21 rows",
      nrow(panel) == 21
    )

    add_check(
      "Panel has 15 genes",
      sum(panel$type == "Gene") == 15
    )

    add_check(
      "Panel has 6 signatures",
      sum(panel$type == "Signature") == 6
    )

    add_check(
      "Spot table has 21 unique features",
      length(
        unique(
          spot$feature
        )
      ) == 21
    )

    add_check(
      "KRT18 absent",
      !"KRT18" %in%
        unique(
          spot$feature
        )
    )

    add_check(
      "Consistency has 21 features",
      length(
        unique(
          consistency$feature
        )
      ) == 21
    )

    add_check(
      "Priority table has 8 features",
      length(
        unique(
          priority$feature
        )
      ) == 8
    )

    add_check(
      "All priority features present",
      setequal(
        unique(
          priority$feature
        ),
        priority_features
      )
    )

    if (
      "ML_agreement_pairs" %in%
      names(priority)
    ) {

      add_check(
        "All 8 priority features supported in both pairs",
        all(
          priority$ML_agreement_pairs ==
            2
        )
      )
    }

    print(
      check_table,
      row.names = FALSE
    )

    # Panel overview

    cat("\n")
    cat("============================================================\n")
    cat("FINAL 21-FEATURE PANEL\n")
    cat("============================================================\n")

    print(
      panel,
      row.names = FALSE
    )

    # Spot counts

    cat("\n")
    cat("============================================================\n")
    cat("SPOT COUNTS BY SECTION\n")
    cat("============================================================\n")

    spot_counts <- aggregate(
      barcode ~
        patient +
        tissue +
        section,
      data = spot,
      FUN = function(x) {
        length(
          unique(x)
        )
      }
    )

    names(
      spot_counts
    )[4] <- "n_spots"

    print(
      spot_counts,
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("SECTION-LEVEL FEATURE SUMMARY\n")
    cat("============================================================\n")

    section_summary <- aggregate(
      expression ~
        patient +
        tissue +
        section +
        feature +
        feature_type,
      data = spot,
      FUN = mean
    )

    names(
      section_summary
    )[
      names(section_summary) ==
        "expression"
    ] <- "mean_expression"

    detected_summary <- aggregate(
      expression ~
        patient +
        tissue +
        section +
        feature +
        feature_type,
      data = spot,
      FUN = function(x) {
        100 *
          mean(
            x > 0
          )
      }
    )

    names(
      detected_summary
    )[
      names(detected_summary) ==
        "expression"
    ] <- "detected_percent"

    section_summary <- merge(
      section_summary,
      detected_summary,
      by = c(
        "patient",
        "tissue",
        "section",
        "feature",
        "feature_type"
      ),
      all = TRUE,
      sort = FALSE
    )

    print(
      section_summary[
        order(
          section_summary$feature,
          section_summary$patient,
          section_summary$tissue
        ),
        ,
        drop = FALSE
      ],
      row.names = FALSE
    )

    # Priority 8 features

    cat("\n")
    cat("============================================================\n")
    cat("8 PRIORITY FEATURES\n")
    cat("============================================================\n")

    priority_order <- match(
      priority_features,
      priority$feature
    )

    priority_order <- priority_order[
      !is.na(
        priority_order
      )
    ]

    priority_sorted <- priority[
      priority_order,
      ,
      drop = FALSE
    ]

    print(
      priority_sorted,
      row.names = FALSE
    )

    priority_spot <- spot[
      spot$feature %in%
        priority_features,
      ,
      drop = FALSE
    ]

    priority_section_summary <- aggregate(
      expression ~
        patient +
        tissue +
        section +
        feature,
      data = priority_spot,
      FUN = mean
    )

    names(
      priority_section_summary
    )[
      names(priority_section_summary) ==
        "expression"
    ] <- "mean_expression"

    cat("\nPriority mean expression by section:\n")

    print(
      priority_section_summary[
        order(
          priority_section_summary$feature,
          priority_section_summary$patient,
          priority_section_summary$tissue
        ),
        ,
        drop = FALSE
      ],
      row.names = FALSE
    )

    # Matched normal vs tumor

    cat("\n")
    cat("============================================================\n")
    cat("MATCHED NORMAL VS TUMOR RESULTS\n")
    cat("============================================================\n")

    print(
      paired[
        order(
          paired$feature,
          paired$patient
        ),
        ,
        drop = FALSE
      ],
      row.names = FALSE
    )

    # Priority matched results

    priority_paired <- paired[
      paired$feature %in%
        priority_features,
      ,
      drop = FALSE
    ]

    priority_paired <- priority_paired[
      order(
        match(
          priority_paired$feature,
          priority_features
        ),
        priority_paired$patient
      ),
      ,
      drop = FALSE
    ]

    cat("\n")
    cat("============================================================\n")
    cat("PRIORITY 8 — MATCHED NORMAL VS TUMOR\n")
    cat("============================================================\n")

    print(
      priority_paired,
      row.names = FALSE
    )

    # Feature ranking for article

    cat("\n")
    cat("============================================================\n")
    cat("FEATURE RANKING FOR ARTICLE\n")
    cat("============================================================\n")

    article_rank <- consistency

    article_rank$abs_mean_delta <- abs(
      article_rank$Mean_delta
    )

    article_rank$priority_feature <-
      article_rank$feature %in%
      priority_features

    article_rank <- article_rank[
      order(
        article_rank$ML_agreement_pairs,
        article_rank$Same_direction_in_both_pairs,
        article_rank$abs_mean_delta,
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]

    cat(
      "\nTop features based on:\n",
      "1. agreement in both matched patients\n",
      "2. same direction in both patients\n",
      "3. magnitude of Tumor-Normal difference\n",
      sep = ""
    )

    print(
      article_rank,
      row.names = FALSE
    )

    # Top priority features by effect size

    priority_rank <- article_rank[
      article_rank$feature %in%
        priority_features,
      ,
      drop = FALSE
    ]

    priority_rank <- priority_rank[
      order(
        priority_rank$abs_mean_delta,
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]

    cat("\n")
    cat("============================================================\n")
    cat("PRIORITY FEATURES RANKED BY |MEAN DELTA|\n")
    cat("============================================================\n")

    print(
      priority_rank,
      row.names = FALSE
    )

    # Simple article result groups

    tumor_high_priority <- priority_rank[
      priority_rank$CRC03_direction ==
        "Higher_in_Tumor" &
        priority_rank$CRC08_direction ==
        "Higher_in_Tumor",
      ,
      drop = FALSE
    ]

    tumor_low_priority <- priority_rank[
      priority_rank$CRC03_direction ==
        "Lower_in_Tumor" &
        priority_rank$CRC08_direction ==
        "Lower_in_Tumor",
      ,
      drop = FALSE
    ]

    cat("\n")
    cat("============================================================\n")
    cat("PRIORITY FEATURES — HIGHER IN TUMOR\n")
    cat("============================================================\n")

    print(
      tumor_high_priority,
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("PRIORITY FEATURES — LOWER IN TUMOR\n")
    cat("============================================================\n")

    print(
      tumor_low_priority,
      row.names = FALSE
    )

    # Figure candidates

    cat("\n")
    cat("============================================================\n")
    cat("ARTICLE FIGURE CANDIDATES\n")
    cat("============================================================\n")

    figure_candidates <- data.frame(

      Rank = seq_len(
        nrow(
          priority_rank
        )
      ),

      Feature =
        priority_rank$feature,

      Direction_CRC03 =
        priority_rank$CRC03_direction,

      Direction_CRC08 =
        priority_rank$CRC08_direction,

      Mean_delta =
        priority_rank$Mean_delta,

      Suggested_use =
        ifelse(
          priority_rank$feature %in%
            c(
              "CEACAM5",
              "EPCAM",
              "KRT8",
              "Signature__Epithelial"
            ),

          "Tumor-high epithelial result",

          "Tumor-low vascular/immune result"
        ),

      stringsAsFactors = FALSE
    )

    print(
      figure_candidates,
      row.names = FALSE
    )

    # Top 4 simple figure features

    top4 <- head(
      priority_rank,
      4
    )

    cat("\n")
    cat("============================================================\n")
    cat("TOP 4 FEATURES FOR MAIN FIGURES\n")
    cat("============================================================\n")

    print(
      top4,
      row.names = FALSE
    )

    # Check spatial pngs

    cat("\n")
    cat("============================================================\n")
    cat("SPATIAL PNG CHECK — PRIORITY FEATURES\n")
    cat("============================================================\n")

    priority_png <- paste0(
      "Spatial_",
      gsub(
        "[^A-Za-z0-9_]+",
        "_",
        priority_features
      ),
      "_MATCHED.png"
    )

    png_table <- data.frame(

      feature =
        priority_features,

      file =
        priority_png,

      exists =
        file.exists(
          file.path(
            outdir,
            priority_png
          )
        ),

      stringsAsFactors = FALSE
    )

    print(
      png_table,
      row.names = FALSE
    )

    # Save review outputs

    write.csv(
      section_summary,
      out(
        "26_Final_21_Features_Section_Summary.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      priority_section_summary,
      out(
        "27_Priority_8_Section_Summary.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      priority_paired,
      out(
        "28_Priority_8_Matched_Normal_Tumor.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      article_rank,
      out(
        "29_Final_21_Features_Article_Ranking.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      priority_rank,
      out(
        "30_Priority_8_Article_Ranking.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      figure_candidates,
      out(
        "31_Article_Figure_Candidates.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      top4,
      out(
        "32_Top4_Main_Figure_Features.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      png_table,
      out(
        "33_Priority_8_Spatial_PNG_Check.csv"
      ),
      row.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("FINAL ARTICLE REVIEW SUMMARY\n")
    cat("============================================================\n")

    cat(
      "\nFinal feature panel:",
      length(
        unique(
          spot$feature
        )
      ),
      "features\n"
    )

    cat(
      "Priority features:",
      length(
        priority_features
      ),
      "\n"
    )

    cat(
      "Priority features higher in Tumor in both pairs:",
      nrow(
        tumor_high_priority
      ),
      "\n"
    )

    cat(
      "Priority features lower in Tumor in both pairs:",
      nrow(
        tumor_low_priority
      ),
      "\n"
    )

    cat(
      "\nTumor-high priority features:\n"
    )

    print(
      tumor_high_priority$feature
    )

    cat(
      "\nTumor-low priority features:\n"
    )

    print(
      tumor_low_priority$feature
    )

    cat(
      "\nTop 4 main-figure candidates:\n"
    )

    print(
      top4$feature
    )

    cat(
      "\nPriority spatial figures available:",
      sum(
        png_table$exists
      ),
      "/",
      nrow(
        png_table
      ),
      "\n"
    )

    cat(
      "\nSaved review files:\n"
    )

    cat(
      "26_Final_21_Features_Section_Summary.csv\n"
    )

    cat(
      "27_Priority_8_Section_Summary.csv\n"
    )

    cat(
      "28_Priority_8_Matched_Normal_Tumor.csv\n"
    )

    cat(
      "29_Final_21_Features_Article_Ranking.csv\n"
    )

    cat(
      "30_Priority_8_Article_Ranking.csv\n"
    )

    cat(
      "31_Article_Figure_Candidates.csv\n"
    )

    cat(
      "32_Top4_Main_Figure_Features.csv\n"
    )

    cat(
      "33_Priority_8_Spatial_PNG_Check.csv\n"
    )

    cat("\n")
    cat("============================================================\n")
    cat("REVIEW COMPLETED\n")
    cat("============================================================\n")

    invisible(
      gc()
    )
  }
  cat("\n[03/06] Result Review and Feature Summaries\n")
  spatial_stage_03(config)
  invisible(gc())

  # Combined and individual spatial png maps
  spatial_stage_04 <- function(config) {

    invisible(gc())

    options(
      width = 220,
      stringsAsFactors = FALSE
    )

    # Paths

    root <- config$data_root

    outdir <- config$analysis_dir

    figdir <- config$figure_dir

    dir.create(
      figdir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    out <- function(x) {
      file.path(outdir, x)
    }

    figout <- function(x) {
      file.path(figdir, x)
    }

    cat(
      "\nFigure output directory:\n",
      figdir,
      "\n",
      sep = ""
    )

    # Packages

    packages <- c(
      "ggplot2"
    )

    for (p in packages) {

      if (!requireNamespace(
        p,
        quietly = TRUE
      )) {

        install.packages(
          p,
          repos = "https://cloud.r-project.org"
        )
      }
    }

    suppressPackageStartupMessages(
      library(ggplot2)
    )

    # Input files

    rctd_file <- out(
      "09_RCTD_All_spot_celltype_proportions_MATCHED.csv"
    )

    feature_file <- out(
      "22_Final_21_Features_All_spots.csv"
    )

    bacteria_file <- out(
      "06_Bacterial_signal_All_spots_MATCHED.csv"
    )

    required_files <- c(
      rctd_file,
      feature_file,
      bacteria_file
    )

    missing_files <- required_files[
      !file.exists(
        required_files
      )
    ]

    if (length(
      missing_files
    )) {

      stop(
        "Missing required input files:\n",
        paste(
          missing_files,
          collapse = "\n"
        )
      )
    }

    # Read data

    cat(
      "\nReading spatial data...\n"
    )

    rctd <- read.csv(
      rctd_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    features <- read.csv(
      feature_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    bacteria <- read.csv(
      bacteria_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    cat(
      "RCTD rows:",
      nrow(rctd),
      "\n"
    )

    cat(
      "Feature rows:",
      nrow(features),
      "\n"
    )

    cat(
      "Bacterial rows:",
      nrow(bacteria),
      "\n"
    )

    # Required columns check

    required_rctd <- c(
      "section",
      "X",
      "Y",
      "cell_type",
      "proportion"
    )

    required_features <- c(
      "section",
      "X",
      "Y",
      "feature",
      "expression"
    )

    required_bacteria <- c(
      "section",
      "X",
      "Y",
      "bacterial_total_raw"
    )

    if (!all(
      required_rctd %in%
      names(rctd)
    )) {

      stop(
        "Missing RCTD columns: ",
        paste(
          setdiff(
            required_rctd,
            names(rctd)
          ),
          collapse = ", "
        )
      )
    }

    if (!all(
      required_features %in%
      names(features)
    )) {

      stop(
        "Missing feature columns: ",
        paste(
          setdiff(
            required_features,
            names(features)
          ),
          collapse = ", "
        )
      )
    }

    if (!all(
      required_bacteria %in%
      names(bacteria)
    )) {

      stop(
        "Missing bacterial columns: ",
        paste(
          setdiff(
            required_bacteria,
            names(bacteria)
          ),
          collapse = ", "
        )
      )
    }

    # Section order

    section_order <- c(
      "CRC_03_Normal",
      "CRC_03_Tumor",
      "CRC_08_Normal",
      "CRC_08_Tumor"
    )

    rctd$section <- factor(
      rctd$section,
      levels = section_order
    )

    features$section <- factor(
      features$section,
      levels = section_order
    )

    bacteria$section <- factor(
      bacteria$section,
      levels = section_order
    )

    # Normalize coordinates within each section

    normalize_xy <- function(df) {

      parts <- split(
        df,
        df$section
      )

      parts <- lapply(
        parts,
        function(z) {

          xr <- range(
            z$X,
            na.rm = TRUE
          )

          yr <- range(
            z$Y,
            na.rm = TRUE
          )

          if (
            diff(xr) == 0
          ) {

            z$X_plot <- 0.5

          } else {

            z$X_plot <- (
              z$X - xr[1]
            ) /
              diff(xr)
          }

          if (
            diff(yr) == 0
          ) {

            z$Y_plot <- 0.5

          } else {

            z$Y_plot <- (
              z$Y - yr[1]
            ) /
              diff(yr)
          }

          z
        }
      )

      result <- do.call(
        rbind,
        parts
      )

      rownames(result) <- NULL

      result
    }

    rctd <- normalize_xy(
      rctd
    )

    features <- normalize_xy(
      features
    )

    bacteria <- normalize_xy(
      bacteria
    )

    # Cell types

    cell_types <- sort(
      unique(
        as.character(
          rctd$cell_type
        )
      )
    )

    cat(
      "\nRCTD cell types:\n"
    )

    print(
      cell_types
    )

    cat(
      "\nNumber of cell types:",
      length(cell_types),
      "\n"
    )

    # Priority features

    priority_features <- c(
      "CEACAM5",
      "EPCAM",
      "KRT8",
      "MS4A1",
      "VWF",
      "PECAM1",
      "Signature__Epithelial",
      "Signature__Endothelial"
    )

    missing_priority <- setdiff(
      priority_features,
      unique(
        features$feature
      )
    )

    if (length(
      missing_priority
    )) {

      stop(
        "Missing priority features: ",
        paste(
          missing_priority,
          collapse = ", "
        )
      )
    }

    priority_data <- features[
      features$feature %in%
        priority_features,
      ,
      drop = FALSE
    ]

    priority_data$feature <- factor(
      priority_data$feature,
      levels = priority_features
    )

    # Top 4 features

    main_features <- c(
      "CEACAM5",
      "Signature__Epithelial",
      "EPCAM",
      "KRT8"
    )

    main_data <- features[
      features$feature %in%
        main_features,
      ,
      drop = FALSE
    ]

    main_data$feature <- factor(
      main_data$feature,
      levels = main_features
    )

    # Bacterial log signal

    if (
      !"log_bacterial" %in%
      names(bacteria)
    ) {

      bacteria$log_bacterial <- log1p(
        bacteria$bacterial_total_raw
      )
    }

    # Common publication theme

    publication_theme <- theme_void(
      base_size = 8
    ) +

      theme(

        plot.title = element_text(
          face = "bold",
          size = 13,
          hjust = 0.5,
          margin = margin(
            b = 8
          )
        ),

        strip.text.x = element_text(
          face = "bold",
          size = 7,
          margin = margin(
            b = 3
          )
        ),

        strip.text.y = element_text(
          face = "bold",
          size = 8,
          angle = 0
        ),

        legend.title = element_text(
          size = 8
        ),

        legend.text = element_text(
          size = 7
        ),

        legend.key.height = grid::unit(
          0.35,
          "cm"
        ),

        panel.spacing = grid::unit(
          0.10,
          "cm"
        ),

        plot.margin = margin(
          5,
          5,
          5,
          5
        )
      )

    # Figure a

    cat(
      "\nCreating Figure A...\n"
    )

    p_rctd <- ggplot(
      rctd,
      aes(
        x = X_plot,
        y = Y_plot,
        color = proportion
      )
    ) +

      geom_point(
        size = 0.32,
        alpha = 0.95,
        stroke = 0
      ) +

      scale_y_reverse() +

      coord_fixed(
        ratio = 1,
        expand = FALSE
      ) +

      facet_grid(
        rows = vars(section),
        cols = vars(cell_type)
      ) +

      scale_color_viridis_c(
        option = "magma",
        direction = 1
      ) +

      labs(
        title =
          "Spatial distribution of RCTD-inferred cell-type proportions",
        color =
          "Cell-type\nproportion"
      ) +

      publication_theme

    ggsave(
      filename = figout(
        "Figure_A_RCTD_CellTypes.png"
      ),
      plot = p_rctd,
      width = max(
        12,
        length(cell_types) * 2
      ),
      height = 9,
      units = "in",
      dpi = 600,
      bg = "white"
    )

    # Figure b

    cat(
      "Creating Figure B...\n"
    )

    p_features <- ggplot(
      priority_data,
      aes(
        x = X_plot,
        y = Y_plot,
        color = expression
      )
    ) +

      geom_point(
        size = 0.32,
        alpha = 0.95,
        stroke = 0
      ) +

      scale_y_reverse() +

      coord_fixed(
        ratio = 1,
        expand = FALSE
      ) +

      facet_grid(
        rows = vars(section),
        cols = vars(feature)
      ) +

      scale_color_viridis_c(
        option = "magma",
        direction = 1
      ) +

      labs(
        title =
          "Spatial expression of validated ML-associated features",
        color =
          "Expression"
      ) +

      publication_theme

    ggsave(
      filename = figout(
        "Figure_B_Priority8_Features.png"
      ),
      plot = p_features,
      width = 18,
      height = 9,
      units = "in",
      dpi = 600,
      bg = "white"
    )

    # Figure b2

    cat(
      "Creating Figure B2...\n"
    )

    p_main <- ggplot(
      main_data,
      aes(
        x = X_plot,
        y = Y_plot,
        color = expression
      )
    ) +

      geom_point(
        size = 0.38,
        alpha = 0.95,
        stroke = 0
      ) +

      scale_y_reverse() +

      coord_fixed(
        ratio = 1,
        expand = FALSE
      ) +

      facet_grid(
        rows = vars(section),
        cols = vars(feature)
      ) +

      scale_color_viridis_c(
        option = "magma"
      ) +

      labs(
        title =
          "Spatial distribution of strongest replicated tumor-associated features",
        color =
          "Expression"
      ) +

      publication_theme

    ggsave(
      filename = figout(
        "Figure_B2_Top4_Features.png"
      ),
      plot = p_main,
      width = 11,
      height = 9,
      units = "in",
      dpi = 600,
      bg = "white"
    )

    # Figure c

    cat(
      "Creating Figure C...\n"
    )

    p_bacteria <- ggplot(
      bacteria,
      aes(
        x = X_plot,
        y = Y_plot,
        color = log_bacterial
      )
    ) +

      geom_point(
        size = 0.55,
        alpha = 0.95,
        stroke = 0
      ) +

      scale_y_reverse() +

      coord_fixed(
        ratio = 1,
        expand = FALSE
      ) +

      facet_wrap(
        ~ section,
        nrow = 1
      ) +

      scale_color_viridis_c(
        option = "magma",
        direction = 1
      ) +

      labs(
        title =
          "Spatial distribution of bacterial signal",
        color =
          "log1p\nbacterial counts"
      ) +

      publication_theme

    ggsave(
      filename = figout(
        "Figure_C_Bacterial_Load.png"
      ),
      plot = p_bacteria,
      width = 12,
      height = 3.6,
      units = "in",
      dpi = 600,
      bg = "white"
    )

    # Individual section figures

    cat(
      "\nCreating individual section figures...\n"
    )

    for (
      ss in section_order
    ) {

      d1 <- rctd[
        rctd$section == ss,
        ,
        drop = FALSE
      ]

      p1 <- ggplot(
        d1,
        aes(
          x = X_plot,
          y = Y_plot,
          color = proportion
        )
      ) +

        geom_point(
          size = 0.45,
          alpha = 0.95,
          stroke = 0
        ) +

        scale_y_reverse() +

        coord_fixed(
          ratio = 1,
          expand = FALSE
        ) +

        facet_wrap(
          ~ cell_type,
          nrow = 1
        ) +

        scale_color_viridis_c(
          option = "magma"
        ) +

        labs(
          title =
            paste0(
              ss,
              " — RCTD cell-type proportions"
            ),
          color =
            "Proportion"
        ) +

        publication_theme

      ggsave(
        filename = figout(
          paste0(
            ss,
            "_A_CellTypes.png"
          )
        ),
        plot = p1,
        width = max(
          12,
          length(cell_types) * 2
        ),
        height = 3.5,
        units = "in",
        dpi = 600,
        bg = "white"
      )

      d2 <- priority_data[
        priority_data$section == ss,
        ,
        drop = FALSE
      ]

      p2 <- ggplot(
        d2,
        aes(
          x = X_plot,
          y = Y_plot,
          color = expression
        )
      ) +

        geom_point(
          size = 0.45,
          alpha = 0.95,
          stroke = 0
        ) +

        scale_y_reverse() +

        coord_fixed(
          ratio = 1,
          expand = FALSE
        ) +

        facet_wrap(
          ~ feature,
          nrow = 1
        ) +

        scale_color_viridis_c(
          option = "magma"
        ) +

        labs(
          title =
            paste0(
              ss,
              " — Priority spatial features"
            ),
          color =
            "Expression"
        ) +

        publication_theme

      ggsave(
        filename = figout(
          paste0(
            ss,
            "_B_Priority8_Features.png"
          )
        ),
        plot = p2,
        width = 18,
        height = 3.6,
        units = "in",
        dpi = 600,
        bg = "white"
      )

      d3 <- bacteria[
        bacteria$section == ss,
        ,
        drop = FALSE
      ]

      p3 <- ggplot(
        d3,
        aes(
          x = X_plot,
          y = Y_plot,
          color = log_bacterial
        )
      ) +

        geom_point(
          size = 0.65,
          alpha = 0.95,
          stroke = 0
        ) +

        scale_y_reverse() +

        coord_fixed(
          ratio = 1,
          expand = FALSE
        ) +

        scale_color_viridis_c(
          option = "magma"
        ) +

        labs(
          title =
            paste0(
              ss,
              " — Bacterial signal"
            ),
          color =
            "log1p bacterial counts"
        ) +

        publication_theme

      ggsave(
        filename = figout(
          paste0(
            ss,
            "_C_Bacterial_Load.png"
          )
        ),
        plot = p3,
        width = 4.5,
        height = 4,
        units = "in",
        dpi = 600,
        bg = "white"
      )
    }

    # Final check

    png_files <- list.files(
      figdir,
      pattern = "\\.png$",
      full.names = FALSE
    )

    cat("\n")
    cat("============================================================\n")
    cat("PUBLICATION SPATIAL FIGURES COMPLETED\n")
    cat("============================================================\n")

    cat(
      "\nOutput directory:\n",
      figdir,
      "\n",
      sep = ""
    )

    cat(
      "\nPNG files created:",
      length(png_files),
      "\n"
    )

    cat(
      "\nCreated PNG files:\n"
    )

    print(
      png_files
    )

    cat(
      "\nMain figures:\n"
    )

    cat(
      "1. Figure_A_RCTD_CellTypes.png\n"
    )

    cat(
      "2. Figure_B_Priority8_Features.png\n"
    )

    cat(
      "3. Figure_B2_Top4_Features.png\n"
    )

    cat(
      "4. Figure_C_Bacterial_Load.png\n"
    )

    cat(
      "\nDONE\n"
    )

    invisible(
      gc()
    )
  }
  cat("\n[04/06] Combined and Individual Spatial PNG Maps\n")
  spatial_stage_04(config)
  invisible(gc())

  # Organize png maps by section
  spatial_stage_05 <- function(config) {

    invisible(gc())

    options(
      stringsAsFactors = FALSE,
      width = 220
    )

    # Root directory

    figdir <- config$figure_dir

    if (!dir.exists(figdir)) {
      stop(
        "Figure directory does not exist:\n",
        figdir
      )
    }

    cat(
      "\nRoot figure directory:\n",
      figdir,
      "\n",
      sep = ""
    )

    # Samples

    samples <- c(
      "CRC_03_Normal",
      "CRC_03_Tumor",
      "CRC_08_Normal",
      "CRC_08_Tumor"
    )

    # Expected original figures

    figure_info <- data.frame(

      sample = c(
        rep("CRC_03_Normal", 3),
        rep("CRC_03_Tumor", 3),
        rep("CRC_08_Normal", 3),
        rep("CRC_08_Tumor", 3)
      ),

      number = rep(
        1:3,
        times = 4
      ),

      source_name = c(

        "CRC_03_Normal_A_CellTypes.png",
        "CRC_03_Normal_B_Priority8_Features.png",
        "CRC_03_Normal_C_Bacterial_Load.png",

        "CRC_03_Tumor_A_CellTypes.png",
        "CRC_03_Tumor_B_Priority8_Features.png",
        "CRC_03_Tumor_C_Bacterial_Load.png",

        "CRC_08_Normal_A_CellTypes.png",
        "CRC_08_Normal_B_Priority8_Features.png",
        "CRC_08_Normal_C_Bacterial_Load.png",

        "CRC_08_Tumor_A_CellTypes.png",
        "CRC_08_Tumor_B_Priority8_Features.png",
        "CRC_08_Tumor_C_Bacterial_Load.png"
      ),

      final_title = c(

        "RCTD_CellTypes",
        "Priority_Features",
        "Bacterial_Signal",

        "RCTD_CellTypes",
        "Priority_Features",
        "Bacterial_Signal",

        "RCTD_CellTypes",
        "Priority_Features",
        "Bacterial_Signal",

        "RCTD_CellTypes",
        "Priority_Features",
        "Bacterial_Signal"
      ),

      stringsAsFactors = FALSE
    )

    # Find all png files recursively

    all_png <- list.files(
      figdir,
      pattern = "\\.png$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    cat(
      "\nTotal PNG files found recursively:",
      length(all_png),
      "\n"
    )

    if (!length(all_png)) {
      stop(
        "No PNG files were found anywhere under:\n",
        figdir
      )
    }

    # Show all found png files

    cat("\n")
    cat("============================================================\n")
    cat("ALL PNG FILES FOUND\n")
    cat("============================================================\n")

    print(
      all_png
    )

    # Find each expected source file

    find_source_file <- function(target_name) {

      hits <- all_png[
        tolower(
          basename(all_png)
        ) ==
          tolower(target_name)
      ]

      if (length(hits) == 0) {

        return(
          NA_character_
        )
      }

      if (length(hits) == 1) {

        return(
          hits[1]
        )
      }

      direct_hit <- hits[
        normalizePath(dirname(hits), winslash = "/", mustWork = FALSE) ==
          normalizePath(
            figdir,
            winslash = "/",
            mustWork = FALSE
          )
      ]

      if (length(direct_hit)) {

        return(
          direct_hit[1]
        )
      }

      hits[1]
    }

    figure_info$source_path <- vapply(
      figure_info$source_name,
      find_source_file,
      character(1)
    )

    figure_info$source_exists <- !is.na(
      figure_info$source_path
    ) &
      file.exists(
        figure_info$source_path
      )

    cat("\n")
    cat("============================================================\n")
    cat("SOURCE SEARCH RESULT\n")
    cat("============================================================\n")

    print(
      figure_info[
        ,
        c(
          "sample",
          "source_name",
          "source_exists",
          "source_path"
        )
      ],
      row.names = FALSE
    )

    # If exact names are not found:

    for (i in seq_len(
      nrow(figure_info)
    )) {

      if (!figure_info$source_exists[i]) {

        ss <- figure_info$sample[i]
        tt <- figure_info$final_title[i]

        pattern_terms <- switch(

          tt,

          "RCTD_CellTypes" =
            c(
              ss,
              "CellTypes"
            ),

          "Priority_Features" =
            c(
              ss,
              "Priority"
            ),

          "Bacterial_Signal" =
            c(
              ss,
              "Bacterial"
            )
        )

        hits <- all_png

        for (term in pattern_terms) {

          hits <- hits[
            grepl(
              term,
              basename(hits),
              fixed = TRUE,
              ignore.case = TRUE
            )
          ]
        }

        if (length(hits) >= 1) {

          figure_info$source_path[i] <- hits[1]

          figure_info$source_exists[i] <- TRUE
        }
      }
    }

    # Final source check

    cat("\n")
    cat("============================================================\n")
    cat("FINAL SOURCE CHECK\n")
    cat("============================================================\n")

    print(
      figure_info[
        ,
        c(
          "sample",
          "final_title",
          "source_exists",
          "source_path"
        )
      ],
      row.names = FALSE
    )

    missing_rows <- which(
      !figure_info$source_exists
    )

    if (length(
      missing_rows
    )) {

      cat(
        "\nFiles still missing:\n"
      )

      print(
        figure_info[
          missing_rows,
          c(
            "sample",
            "source_name"
          ),
          drop = FALSE
        ],
        row.names = FALSE
      )

      stop(
        "\nSome source images could not be located.\n",
        "Check the printed ALL PNG FILES FOUND list above."
      )
    }

    # Create sample folders

    for (ss in samples) {

      sample_dir <- file.path(
        figdir,
        ss
      )

      dir.create(
        sample_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }

    # Final file names

    figure_info$final_name <- sprintf(
      "Figure_%02d_%s.png",
      figure_info$number,
      figure_info$final_title
    )

    figure_info$output_path <- file.path(
      figdir,
      figure_info$sample,
      figure_info$final_name
    )

    # Copy + rename

    figure_info$copied <- FALSE

    for (i in seq_len(
      nrow(
        figure_info
      )
    )) {

      src <- figure_info$source_path[i]

      dst <- figure_info$output_path[i]

      src_norm <- normalizePath(
        src,
        winslash = "/",
        mustWork = FALSE
      )

      dst_norm <- normalizePath(
        dst,
        winslash = "/",
        mustWork = FALSE
      )

      if (
        identical(
          src_norm,
          dst_norm
        )
      ) {

        figure_info$copied[i] <- TRUE

      } else {

        figure_info$copied[i] <- file.copy(
          from = src,
          to = dst,
          overwrite = TRUE
        )
      }
    }

    # Check outputs

    figure_info$output_exists <- file.exists(
      figure_info$output_path
    )

    cat("\n")
    cat("============================================================\n")
    cat("FINAL ORGANIZATION CHECK\n")
    cat("============================================================\n")

    print(
      figure_info[
        ,
        c(
          "sample",
          "final_name",
          "copied",
          "output_exists"
        )
      ],
      row.names = FALSE
    )

    # Verify 3 figures per sample

    cat("\n")
    cat("============================================================\n")
    cat("FILES PER SAMPLE\n")
    cat("============================================================\n")

    for (ss in samples) {

      sample_dir <- file.path(
        figdir,
        ss
      )

      pngs <- list.files(
        sample_dir,
        pattern = "\\.png$",
        full.names = FALSE,
        ignore.case = TRUE
      )

      cat(
        "\n",
        ss,
        " — ",
        length(pngs),
        " PNG files\n",
        sep = ""
      )

      print(
        pngs
      )
    }

    cat("\n")
    cat("============================================================\n")

    if (
      all(
        figure_info$output_exists
      )
    ) {

      cat(
        "FINAL STATUS: PASS\n"
      )

      cat(
        "All 12 images organized successfully.\n"
      )

    } else {

      cat(
        "FINAL STATUS: REVIEW REQUIRED\n"
      )

      cat(
        "Some output images are missing.\n"
      )
    }

    cat(
      "============================================================\n"
    )

    # Final folder structure

    cat(
      "\nFinal folders:\n"
    )

    for (ss in samples) {

      cat(
        "\n",
        file.path(
          figdir,
          ss
        ),
        "\n",
        sep = ""
      )
    }

    cat(
      "\nDONE\n"
    )

    invisible(
      gc()
    )
  }
  cat("\n[05/06] Organize PNG Maps by Section\n")
  spatial_stage_05(config)
  invisible(gc())

  # Final highlighted spatial png maps
  spatial_stage_06 <- function(config) {

    invisible(gc())

    options(
      width = 220,
      stringsAsFactors = FALSE
    )

    # Paths

    st_root <- config$data_root

    st_outdir <- config$analysis_dir

    figdir <- config$figure_dir

    dir.create(
      figdir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    out <- function(x) {
      file.path(
        st_outdir,
        x
      )
    }

    # Packages

    if (!requireNamespace(
      "ggplot2",
      quietly = TRUE
    )) {
      install.packages(
        "ggplot2",
        repos = "https://cloud.r-project.org"
      )
    }

    suppressPackageStartupMessages(
      library(ggplot2)
    )

    # Input files

    rctd_file <- out(
      "09_RCTD_All_spot_celltype_proportions_MATCHED.csv"
    )

    feature_file <- out(
      "22_Final_21_Features_All_spots.csv"
    )

    bacteria_file <- out(
      "06_Bacterial_signal_All_spots_MATCHED.csv"
    )

    required_files <- c(
      rctd_file,
      feature_file,
      bacteria_file
    )

    missing_files <- required_files[
      !file.exists(
        required_files
      )
    ]

    if (length(
      missing_files
    )) {
      stop(
        "Missing required files:\n",
        paste(
          missing_files,
          collapse = "\n"
        )
      )
    }

    # Read data

    cat("\nReading data...\n")

    rctd <- read.csv(
      rctd_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    features <- read.csv(
      feature_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    bacteria <- read.csv(
      bacteria_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # Required column check

    rctd_required <- c(
      "section",
      "X",
      "Y",
      "cell_type",
      "proportion"
    )

    feature_required <- c(
      "section",
      "X",
      "Y",
      "feature",
      "expression"
    )

    bacteria_required <- c(
      "section",
      "X",
      "Y",
      "bacterial_total_raw"
    )

    if (!all(
      rctd_required %in%
      names(rctd)
    )) {
      stop(
        "Missing RCTD columns: ",
        paste(
          setdiff(
            rctd_required,
            names(rctd)
          ),
          collapse = ", "
        )
      )
    }

    if (!all(
      feature_required %in%
      names(features)
    )) {
      stop(
        "Missing feature columns: ",
        paste(
          setdiff(
            feature_required,
            names(features)
          ),
          collapse = ", "
        )
      )
    }

    if (!all(
      bacteria_required %in%
      names(bacteria)
    )) {
      stop(
        "Missing bacterial columns: ",
        paste(
          setdiff(
            bacteria_required,
            names(bacteria)
          ),
          collapse = ", "
        )
      )
    }

    # Sample order

    samples <- c(
      "CRC_03_Normal",
      "CRC_03_Tumor",
      "CRC_08_Normal",
      "CRC_08_Tumor"
    )

    rctd$section <- factor(
      rctd$section,
      levels = samples
    )

    features$section <- factor(
      features$section,
      levels = samples
    )

    bacteria$section <- factor(
      bacteria$section,
      levels = samples
    )

    # Create sample folders

    for (ss in samples) {
      dir.create(
        file.path(
          figdir,
          ss
        ),
        recursive = TRUE,
        showWarnings = FALSE
      )
    }

    # Normalize coordinates

    normalize_xy <- function(df) {

      parts <- split(
        df,
        df$section
      )

      parts <- lapply(
        parts,
        function(z) {

          xr <- range(
            z$X,
            na.rm = TRUE
          )

          yr <- range(
            z$Y,
            na.rm = TRUE
          )

          if (
            is.finite(diff(xr)) &&
            diff(xr) > 0
          ) {

            z$X_plot <- (
              z$X - xr[1]
            ) /
              diff(xr)

          } else {

            z$X_plot <- 0.5
          }

          if (
            is.finite(diff(yr)) &&
            diff(yr) > 0
          ) {

            z$Y_plot <- (
              z$Y - yr[1]
            ) /
              diff(yr)

          } else {

            z$Y_plot <- 0.5
          }

          z
        }
      )

      ans <- do.call(
        rbind,
        parts
      )

      rownames(ans) <- NULL

      ans
    }

    rctd <- normalize_xy(
      rctd
    )

    features <- normalize_xy(
      features
    )

    bacteria <- normalize_xy(
      bacteria
    )

    # Bacterial log signal

    if (!"log_bacterial" %in%
        names(bacteria)) {

      bacteria$log_bacterial <- log1p(
        bacteria$bacterial_total_raw
      )
    }

    # Priority features

    priority_features <- c(
      "CEACAM5",
      "EPCAM",
      "KRT8",
      "MS4A1",
      "VWF",
      "PECAM1",
      "Signature__Epithelial",
      "Signature__Endothelial"
    )

    missing_priority <- setdiff(
      priority_features,
      unique(
        features$feature
      )
    )

    if (length(
      missing_priority
    )) {
      stop(
        "Missing priority features: ",
        paste(
          missing_priority,
          collapse = ", "
        )
      )
    }

    features <- features[
      features$feature %in%
        priority_features,
      ,
      drop = FALSE
    ]

    features$feature <- factor(
      features$feature,
      levels = priority_features
    )

    # Visualization settings

    RCTD_THRESHOLD <- 0.10

    EXPRESSION_QUANTILE <- 0.75

    BACTERIA_QUANTILE <- 0.75

    LOW_SPOT_COLOR <- "#1F1F1F"

    LOW_SPOT_ALPHA <- 0.60

    # Feature threshold

    add_expression_threshold <- function(df) {

      groups <- split(
        df,
        interaction(
          df$section,
          df$feature,
          drop = TRUE
        )
      )

      groups <- lapply(
        groups,
        function(z) {

          positive <- z$expression[
            is.finite(
              z$expression
            ) &
              z$expression > 0
          ]

          if (!length(
            positive
          )) {

            cutoff <- Inf

          } else {

            cutoff <- as.numeric(
              quantile(
                positive,
                probs =
                  EXPRESSION_QUANTILE,
                na.rm = TRUE,
                names = FALSE
              )
            )
          }

          z$highlight <-
            is.finite(
              z$expression
            ) &
            z$expression >= cutoff &
            z$expression > 0

          z
        }
      )

      ans <- do.call(
        rbind,
        groups
      )

      rownames(ans) <- NULL

      ans
    }

    features <- add_expression_threshold(
      features
    )

    # Bacterial threshold

    add_bacteria_threshold <- function(df) {

      groups <- split(
        df,
        df$section
      )

      groups <- lapply(
        groups,
        function(z) {

          positive <- z$log_bacterial[
            is.finite(
              z$log_bacterial
            ) &
              z$log_bacterial > 0
          ]

          if (!length(
            positive
          )) {

            cutoff <- Inf

          } else {

            cutoff <- as.numeric(
              quantile(
                positive,
                probs =
                  BACTERIA_QUANTILE,
                na.rm = TRUE,
                names = FALSE
              )
            )
          }

          z$highlight <-
            is.finite(
              z$log_bacterial
            ) &
            z$log_bacterial >= cutoff &
            z$log_bacterial > 0

          z
        }
      )

      ans <- do.call(
        rbind,
        groups
      )

      rownames(ans) <- NULL

      ans
    }

    bacteria <- add_bacteria_threshold(
      bacteria
    )

    # Rctd highlight

    rctd$highlight <-
      is.finite(
        rctd$proportion
      ) &
      rctd$proportion >=
      RCTD_THRESHOLD

    # Publication theme

    publication_theme <- theme_void(
      base_size = 10
    ) +

      theme(

        plot.title = element_text(
          face = "bold",
          size = 15,
          hjust = 0.5,
          margin = margin(
            b = 12
          )
        ),

        strip.text = element_text(
          face = "bold",
          size = 10
        ),

        legend.title = element_text(
          size = 9
        ),

        legend.text = element_text(
          size = 8
        ),

        legend.key.height = grid::unit(
          0.5,
          "cm"
        ),

        panel.spacing = grid::unit(
          0.14,
          "cm"
        ),

        plot.margin = margin(
          8,
          8,
          8,
          8
        )
      )

    # High-spot colour scale

    high_colour_scale <- scale_color_gradientn(
      colours = c(
        "#FFD54F",
        "#FF9800",
        "#F4511E",
        "#B71C1C"
      )
    )

    # Generate figures

    for (ss in samples) {

      cat(
        "\nProcessing: ",
        ss,
        "\n",
        sep = ""
      )

      sample_dir <- file.path(
        figdir,
        ss
      )

      d1 <- rctd[
        rctd$section == ss,
        ,
        drop = FALSE
      ]

      low1 <- d1[
        !d1$highlight,
        ,
        drop = FALSE
      ]

      high1 <- d1[
        d1$highlight,
        ,
        drop = FALSE
      ]

      p1 <- ggplot() +

        geom_point(
          data = low1,
          aes(
            x = X_plot,
            y = Y_plot
          ),
          color = LOW_SPOT_COLOR,
          size = 0.35,
          alpha = LOW_SPOT_ALPHA,
          stroke = 0
        ) +

        geom_point(
          data = high1,
          aes(
            x = X_plot,
            y = Y_plot,
            color = proportion,
            size = proportion
          ),
          alpha = 0.95,
          stroke = 0
        ) +

        scale_size_continuous(
          range = c(
            0.55,
            1.55
          ),
          guide = "none"
        ) +

        high_colour_scale +

        scale_y_reverse() +

        coord_fixed(
          ratio = 1,
          expand = FALSE
        ) +

        facet_wrap(
          ~ cell_type,
          nrow = 1
        ) +

        labs(
          title = paste0(
            ss,
            " — RCTD cell-type proportions"
          ),
          color = "Proportion"
        ) +

        publication_theme

      ggsave(
        filename = file.path(
          sample_dir,
          "Figure_01_RCTD_CellTypes.png"
        ),
        plot = p1,
        width = max(
          13,
          length(
            unique(
              d1$cell_type
            )
          ) * 2.2
        ),
        height = 4.5,
        units = "in",
        dpi = 600,
        bg = "white"
      )

      d2 <- features[
        features$section == ss,
        ,
        drop = FALSE
      ]

      low2 <- d2[
        !d2$highlight,
        ,
        drop = FALSE
      ]

      high2 <- d2[
        d2$highlight,
        ,
        drop = FALSE
      ]

      p2 <- ggplot() +

        geom_point(
          data = low2,
          aes(
            x = X_plot,
            y = Y_plot
          ),
          color = LOW_SPOT_COLOR,
          size = 0.30,
          alpha = LOW_SPOT_ALPHA,
          stroke = 0
        ) +

        geom_point(
          data = high2,
          aes(
            x = X_plot,
            y = Y_plot,
            color = expression,
            size = expression
          ),
          alpha = 0.95,
          stroke = 0
        ) +

        scale_size_continuous(
          range = c(
            0.55,
            1.55
          ),
          guide = "none"
        ) +

        high_colour_scale +

        scale_y_reverse() +

        coord_fixed(
          ratio = 1,
          expand = FALSE
        ) +

        facet_wrap(
          ~ feature,
          nrow = 1
        ) +

        labs(
          title = paste0(
            ss,
            " — Priority spatial features"
          ),
          color = "Expression"
        ) +

        publication_theme

      ggsave(
        filename = file.path(
          sample_dir,
          "Figure_02_Priority_Features.png"
        ),
        plot = p2,
        width = 18,
        height = 4.5,
        units = "in",
        dpi = 600,
        bg = "white"
      )

      d3 <- bacteria[
        bacteria$section == ss,
        ,
        drop = FALSE
      ]

      low3 <- d3[
        !d3$highlight,
        ,
        drop = FALSE
      ]

      high3 <- d3[
        d3$highlight,
        ,
        drop = FALSE
      ]

      p3 <- ggplot() +

        geom_point(
          data = low3,
          aes(
            x = X_plot,
            y = Y_plot
          ),
          color = LOW_SPOT_COLOR,
          size = 0.45,
          alpha = LOW_SPOT_ALPHA,
          stroke = 0
        ) +

        geom_point(
          data = high3,
          aes(
            x = X_plot,
            y = Y_plot,
            color = log_bacterial,
            size = log_bacterial
          ),
          alpha = 0.95,
          stroke = 0
        ) +

        scale_size_continuous(
          range = c(
            0.8,
            2.2
          ),
          guide = "none"
        ) +

        high_colour_scale +

        scale_y_reverse() +

        coord_fixed(
          ratio = 1,
          expand = FALSE
        ) +

        labs(
          title = paste0(
            ss,
            " — Bacterial signal"
          ),
          color =
            "log1p bacterial counts"
        ) +

        publication_theme

      ggsave(
        filename = file.path(
          sample_dir,
          "Figure_03_Bacterial_Signal.png"
        ),
        plot = p3,
        width = 5.5,
        height = 5.2,
        units = "in",
        dpi = 600,
        bg = "white"
      )
    }

    # Final check

    cat("\n")
    cat("============================================================\n")
    cat("FINAL OUTPUT CHECK\n")
    cat("============================================================\n")

    all_ok <- TRUE

    for (ss in samples) {

      sample_dir <- file.path(
        figdir,
        ss
      )

      expected <- c(
        "Figure_01_RCTD_CellTypes.png",
        "Figure_02_Priority_Features.png",
        "Figure_03_Bacterial_Signal.png"
      )

      exists <- file.exists(
        file.path(
          sample_dir,
          expected
        )
      )

      cat(
        "\n",
        ss,
        "\n",
        sep = ""
      )

      print(
        data.frame(
          Figure = expected,
          Exists = exists,
          stringsAsFactors = FALSE
        ),
        row.names = FALSE
      )

      if (!all(
        exists
      )) {
        all_ok <- FALSE
      }
    }

    cat("\n")
    cat("============================================================\n")

    if (all_ok) {

      cat(
        "FINAL STATUS: PASS\n"
      )

      cat(
        "All 12 PNG figures created successfully.\n"
      )

    } else {

      cat(
        "FINAL STATUS: REVIEW REQUIRED\n"
      )
    }

    cat(
      "============================================================\n"
    )

    cat(
      "\nLow/zero spot color: ",
      LOW_SPOT_COLOR,
      "\n",
      sep = ""
    )

    cat(
      "High spots: yellow -> orange -> red\n"
    )

    cat(
      "High spots are larger.\n"
    )

    cat(
      "\nDONE\n"
    )

    invisible(
      gc()
    )
  }
  cat("\n[06/06] Final Highlighted Spatial PNG Maps\n")
  spatial_stage_06(config)
  invisible(gc())

  cat("\n07 - Spatial Analysis of Selected Features: COMPLETE\n")
  cat("Analysis outputs: ", config$analysis_dir, "\n", sep = "")
  cat("PNG maps: ", config$figure_dir, "\n", sep = "")
  cat("Cohort: two matched patients, four tissue sections.\n")
  cat("Final sample folders: CRC_03_Normal, CRC_03_Tumor, CRC_08_Normal, CRC_08_Tumor.\n")
  cat("Twelve organized PNG maps: three per section (RCTD, selected features, bacterial signal).\n")
  invisible(TRUE)
}

# Spatially adjusted rank-correlation inference and supplementary workbook

run_spatial_adjustment <- function(config) {
  root <- config$analysis_dir
  if(!requireNamespace("SpatialPack",quietly=TRUE) || !requireNamespace("openxlsx",quietly=TRUE))
    stop("Install packages: install.packages(c('SpatialPack','openxlsx'))")
  cert_path <- file.path(root,"Bacterial_extraction_VALIDATED.csv")
  if(!file.exists(cert_path)) stop("Completed bacteria-only extraction certificate missing.")
  cert <- read.csv(cert_path,stringsAsFactors=FALSE)
  if(nrow(cert)!=1L || !"Extraction_version" %in% names(cert) ||
     cert$Extraction_version != "feature_name_column2_bacteria_only_v2") stop("Unsupported bacterial extraction certificate.")
  output_dir <- "E:/CRC/Results/Results/Supplementary"
  workbook <- file.path(output_dir,"Supplementary Table S15.xlsx")
  read_input <- function(n,req) {
    f <- file.path(root,n)
    if(!file.exists(f)) stop("Missing spatial inference input: ",f)
    d <- read.csv(f,check.names=FALSE,stringsAsFactors=FALSE)
    if(!all(req %in% names(d))) stop("Unexpected columns in ",n,": need ",paste(req,collapse=", "))
    d
  }
  b <- read_input("06_Bacterial_signal_All_spots_MATCHED.csv",c("section","patient","tissue","barcode","X","Y","bacterial_total_raw"))
  f <- read_input("22_Final_21_Features_All_spots.csv",c("section","barcode","feature","expression","X","Y"))
  r <- read_input("09_RCTD_All_spot_celltype_proportions_MATCHED.csv",c("section","barcode","cell_type","proportion","X","Y"))
  key <- function(d) paste(d$section,d$barcode,sep="::")
  if(anyDuplicated(key(b)) || anyDuplicated(paste(key(f),f$feature)) || anyDuplicated(paste(key(r),r$cell_type)))
    stop("Duplicate specimen/spot/feature identifiers.")
  sections <- c("CRC_03_Normal","CRC_03_Tumor","CRC_08_Normal","CRC_08_Tumor")
  if(!setequal(unique(b$section),sections)) stop("Expected four matched sections for spatial testing.")
  feature_names <- sort(unique(f$feature)); type_names <- sort(unique(r$cell_type))
  if(length(feature_names)!=21L || length(type_names)!=6L) stop("Expected 21 features and six RCTD types.")
  results <- list(); diagnostics <- list(); counter <- 0L

  bins <- c(13L,8L,20L)
  test_one <- function(z,meta,nclass) {
    ok <- is.finite(z$x) & is.finite(z$y) & is.finite(z$X) & is.finite(z$Y)
    z <- z[ok,,drop=FALSE]; n <- nrow(z)
    rho <- if(n>=3L && sd(z$x)>0 && sd(z$y)>0) cor(z$x,z$y,method="spearman") else NA_real_
    base <- data.frame(meta,N_spots=n,Spearman_rho=rho,Nclass=nclass,
                       Effective_sample_size=NA_real_,Degrees_freedom=NA_real_,P_spatial=NA_real_,
                       Status="NOT_ESTIMABLE",Warnings="",stringsAsFactors=FALSE)
    if(n<30L || !is.finite(rho)) {base$Warnings <- "Fewer than 30 valid spots or a constant variable";return(base)}
    coords <- as.matrix(z[,c("X","Y")]);storage.mode(coords) <- "double"
    if(anyDuplicated(data.frame(coords))) stop("Duplicate coordinates within a section; inspect spot mapping.")
    xr <- rank(z$x,ties.method="average"); yr <- rank(z$y,ties.method="average")
    warn <- character()
    ans <- tryCatch(withCallingHandlers(SpatialPack::modified.ttest(xr,yr,coords,nclass=nclass),
                                        warning=function(w){warn <<- c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
    if(inherits(ans,"error")) {base$Warnings <- conditionMessage(ans);return(base)}
    if(!is.finite(ans$corr) || abs(ans$corr-rho)>1e-8) stop("Rank correlation differs from Spearman; inspect alignment.")
    base$Effective_sample_size <- ans$ESS; base$Degrees_freedom <- ans$dof
    base$Warnings <- paste(unique(warn),collapse=" | ")
    if(!is.finite(ans$ESS) || !is.finite(ans$dof) || ans$dof<=0 ||
       !is.finite(ans$p.value) || ans$p.value<0 || ans$p.value>1) {
      base$Warnings <- paste(base$Warnings,"Invalid effective sample size, degrees of freedom or P value")
      return(base)
    }
    base$P_spatial <- ans$p.value
    base$Status <- if(length(warn)) "REVIEW_WARNING" else "OK"

    if(length(warn)) base$P_spatial <- NA_real_
    diagnostics[[length(diagnostics)+1L]] <<- data.frame(meta,Nclass=nclass,
                                                         Lag_class=seq_along(ans$upper.bounds),Upper_distance=ans$upper.bounds,
                                                         Pair_count=ans$card,Moran_rank_x=ans$imoran[,1],Moran_rank_y=ans$imoran[,2])
    base
  }

  align <- function(base,d,label,value) {
    rows <- which(!is.na(d[[label]]) & d[[label]] == value)
    z <- d[rows,,drop=FALSE]
    if (!nrow(z)) stop("Missing variable in section: ", value)
    if (anyNA(z$barcode) || anyDuplicated(z$barcode))
      stop("Missing or duplicated variable barcodes: ", value)
    idx <- match(base$barcode,z$barcode)
    matched <- which(!is.na(idx))
    if(length(matched)) {
      delta <- abs(base$X[matched]-z$X[idx[matched]])+abs(base$Y[matched]-z$Y[idx[matched]])
      if(any(!is.finite(delta)) || any(delta>1e-6)) stop("Inconsistent spot coordinates: ",value)
    }
    original_idx <- rows[idx]
    if (length(matched) &&
        (!identical(as.character(d$barcode[original_idx[matched]]),
                    as.character(base$barcode[matched])) ||
         any(d[[label]][original_idx[matched]] != value)))
      stop("Variable/barcode alignment failed: ", value)
    original_idx
  }
  add_pair <- function(z,meta) {
    for(k in bins) {
      counter <<- counter+1L
      if(counter %% 30L==1L) cat("Spatial tests: ",counter," / 1836; ",meta$Section,"; ",meta$Variable_x," / ",meta$Variable_y,"\n",sep="")
      results[[length(results)+1L]] <<- test_one(z,meta,k)
    }
  }
  for(ss in sections) {
    base <- b[b$section==ss,,drop=FALSE];base <- base[order(base$barcode),]
    fs <- f[f$section==ss,,drop=FALSE];rs <- r[r$section==ss,,drop=FALSE]
    meta <- function(family,x,y) data.frame(Family=family,Section=ss,Patient=base$patient[1],Tissue=base$tissue[1],Variable_x=x,Variable_y=y)
    for(fn in feature_names) {
      ix <- align(base,fs,"feature",fn); xv <- fs$expression[ix]
      add_pair(data.frame(X=base$X,Y=base$Y,x=xv,y=log1p(base$bacterial_total_raw)),meta("Feature_Bacteria",fn,"Total_bacterial_signal"))
      for(ct in type_names) {
        iy <- align(base,rs,"cell_type",ct)
        add_pair(data.frame(X=base$X,Y=base$Y,x=xv,y=rs$proportion[iy]),meta("Feature_RCTD",fn,ct))
      }
    }
    for(ct in type_names) {
      iy <- align(base,rs,"cell_type",ct)
      add_pair(data.frame(X=base$X,Y=base$Y,x=log1p(base$bacterial_total_raw),y=rs$proportion[iy]),meta("Bacteria_RCTD","Total_bacterial_signal",ct))
    }
  }
  all_results <- do.call(rbind,results)
  all_results$Alignment_version <- "parent_row_barcode_match_v3"
  all_results$FDR_spatial <- NA_real_
  all_results$FDR_family_planned_tests <- NA_integer_
  for(k in bins) for(fam in unique(all_results$Family)) {
    ids <- which(all_results$Family==fam & all_results$Nclass==k)
    good <- ids[is.finite(all_results$P_spatial[ids]) & all_results$Status[ids]=="OK"]
    all_results$FDR_family_planned_tests[ids] <- length(ids)

    if(length(good)) all_results$FDR_spatial[good] <- p.adjust(all_results$P_spatial[good],"BH",n=length(ids))
  }
  primary <- all_results[all_results$Nclass==13L,,drop=FALSE]
  sensitivity <- all_results[all_results$Nclass!=13L,,drop=FALSE]
  stopifnot(nrow(primary)==612L,nrow(sensitivity)==1224L)
  summary <- do.call(rbind,lapply(split(primary,primary$Family),function(z) data.frame(
    Family=z$Family[1],Planned_tests=nrow(z),Valid_tests=sum(z$Status=="OK" & is.finite(z$P_spatial)),
    Not_estimable_or_review=sum(z$Status!="OK"),FDR_below_0_05=sum(z$FDR_spatial<.05,na.rm=TRUE))))
  notes <- data.frame(Note=c(
    "Alignment revision: parent_row_barcode_match_v3; each variable is matched to its own original rows. Earlier spatial correlations must be replaced.",
    "Primary analysis: SpatialPack modified.ttest on within-section average ranks; nclass=13. Correlation equals Spearman rho; spatial P is approximate.",
    "Sensitivity analyses: nclass=8 and 20; same data and hypotheses; do not select the most favorable P/FDR.",
    "612 primary comparisons: Feature_Bacteria 84, Feature_RCTD 504, Bacteria_RCTD 24, spanning four sections.",
    "BH adjustment separately across each entire comparison family (all four sections), and separately for each nclass. Planned family size retained for missing tests.",
    "Coordinates use original X/Y pixel coordinates in common units; no separate axis scaling. Matching uses section and barcode; coordinate equality checked.",
    "Bacteria-only total comes from 19 deposited bacterial groups; viral and control probes excluded. Labels are annotation-based, not proof of viable bacteria.",
    "Unestimable or warning-bearing tests have no reportable spatial P/FDR. No fallback to naive independent-spot P values.",
    "Inference assumes an adequately stationary/isotropic spatial covariance model within each section; distance-bin Moran diagnostics are provided. Tissue heterogeneity can violate these assumptions.",
    "Spatial rank testing is an approximate adaptation of the modified Pearson test, especially with many ties/zero counts; interpret as exploratory sensitivity evidence.",
    "Only two independent patients. Section-level P/FDR does not establish patient-level tumor-normal differences, causality, or physical microbial interaction.",
    "Original spot-independent P/FDR tables remain for provenance and should not support claims of spatially adjusted significance.",
    "References: Clifford P, Richardson S, Hemon D. Biometrics. 1989;45:123-134. Dutilleul P. Biometrics. 1993;49:305-314.",
    "Documentation: https://search.r-project.org/CRAN/refmans/SpatialPack/html/modified.ttest.html"))
  session <- data.frame(Information=c(paste("R",R.version.string),paste("SpatialPack",packageVersion("SpatialPack")),
                                      paste("openxlsx",packageVersion("openxlsx")),capture.output(sessionInfo())))
  write.csv(primary,file.path(root,"33_Spatial_Adjusted_Correlations.csv"),row.names=FALSE,na="")
  write.csv(sensitivity,file.path(root,"34_Spatial_Adjusted_Sensitivity.csv"),row.names=FALSE,na="")
  write.csv(summary,file.path(root,"35_Spatial_Adjusted_Summary.csv"),row.names=FALSE,na="")
  diag <- if(length(diagnostics)) do.call(rbind,diagnostics) else data.frame(Note="No estimable spatial tests")
  write.csv(diag,file.path(root,"36_Spatial_Lag_Diagnostics.csv"),row.names=FALSE,na="")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  wb <- openxlsx::createWorkbook()
  sheets <- list(Primary_Spatial_Tests=primary,Sensitivity_8_20=sensitivity,Summary=summary,
                 Spatial_Lag_Diagnostics=diag,Methods_and_Limits=notes,Session=session)
  for(nm in names(sheets)) {
    openxlsx::addWorksheet(wb,nm);openxlsx::writeDataTable(wb,nm,sheets[[nm]])
    openxlsx::freezePane(wb,nm,firstRow=TRUE)
    openxlsx::setColWidths(wb,nm,cols=seq_len(ncol(sheets[[nm]])),widths=20)
  }
  openxlsx::saveWorkbook(wb,workbook,overwrite=TRUE)
  cat("\nSPATIALLY ADJUSTED RANK-CORRELATION SUMMARY\n")
  print(summary,row.names=FALSE)
  print(primary,row.names=FALSE)
  cat("\nWorkbook: ",workbook,"\nTwo patients only; tests are approximate and section-specific.\n",sep="")
}

audit_file <- file.path(spatial_config$analysis_dir, "Bacterial_extraction_VALIDATED.csv")
run_spatial_workflow(spatial_config)
write.csv(data.frame(Extraction_version="feature_name_column2_bacteria_only_v2",
                     Feature_name_column=2L, Feature_type_column=3L, Bacterial_group_allowlist_size=19L,
                     Interpretation="Bacterial-only deposited probe groups; viruses and controls excluded; not independently verified taxonomy",
                     Completed=as.character(Sys.time())), audit_file, row.names=FALSE)

run_spatial_adjustment(spatial_config)
