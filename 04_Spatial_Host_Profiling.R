# Corrected stage 04: host marker/QC variable collision fixed; validated RCTD caches retained.
# Stage 04 — necessary spatial preprocessing, RCTD and descriptive maps.
# Uses supplied audited raw extraction; no classifier, trajectory or CellChat is run.
# Output directory is explicitly stage 04. Downstream readers must use this same path.
output_root <- "E:/CRC/Results/Results/04-Spatial_Analysis_of_Selected_Features"

spatial_config <- list(
  data_root = "E:/CRC/Results/GSE334323_ST",
  archive = "E:/CRC/Results/GSE334323_ST/GSE334323_RAW.tar",
  reference_file = "E:/CRC/Results/Results/02_Clustering_Markers_Annotation/05_annotated.rds",
  raw_dir = file.path(output_root, "Matched_RAW_Extracted"),
  analysis_dir = output_root,
  figure_dir = file.path(output_root, "Figures"),
  seed = 42L,
  max_cores = 2L
)

run_spatial_workflow <- function(config) {
  previous_options <- options(width = 220, stringsAsFactors = FALSE)
  on.exit(options(previous_options), add = TRUE)
  required_packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "spacexr", "digest")
  missing_packages <- required_packages[!vapply(required_packages,
                                                requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop("Install missing package(s) before running: ",
         paste(missing_packages, collapse = ", "),
         "\nFor spacexr: remotes::install_github('dmcable/spacexr').",
         call. = FALSE)
  }
  if(requireNamespace("future",quietly=TRUE)) future::plan("sequential")
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
  complete <- file.path(config$analysis_dir,"Core_Spatial_Complete.rds")
  if(file.exists(complete)) unlink(complete)
  set.seed(config$seed)
  
  spatial_stage_01 <- function(config) {
    # ============================================================
    # GSE334323 — MATCHED-PAIR SPATIAL TRANSCRIPTOMICS
    # END-TO-END FROM RAW TAR
    #
    # STARTING INPUT:
    #   E:/CRC/Results/GSE334323_ST/GSE334323_RAW.tar
    #
    # scRNA-seq reference:
    #   E:/CRC/Results/Results/02_Clustering_Markers_Annotation/05_annotated.rds
    #
    # ONLY MATCHED PAIRS:
    #   CRC_03_Normal
    #   CRC_03_Tumor
    #   CRC_08_Normal
    #   CRC_08_Tumor
    #
    # Raw extraction, host/bacterial separation, QC, normalization and RCTD only.
    invisible(gc())
    
    set.seed(config$seed)
    
    # ============================================================
    # 01.01. PATHS
    # ============================================================
    
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
    
    # ============================================================
    # 01.02. PACKAGES
    # ============================================================
    
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
    
    # ============================================================
    # 01.03. MATCHED SAMPLE METADATA
    # ============================================================
    
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
    
    # ============================================================
    # 01.04. FIRST: INSPECT RAW TAR
    # ============================================================
    
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
    
    # Check which entries belong to matched GSMs / sections
    
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
    
    # ============================================================
    # 01.05. SECURITY CHECK BEFORE EXTRACTION
    # ============================================================
    
    unsafe <- grepl(
      "(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",
      tar_entries
    )
    
    if (any(unsafe)) {
      
      stop(
        "Unsafe path detected inside TAR archive."
      )
    }
    
    # ============================================================
    # 01.06. EXTRACT RAW ARCHIVE
    # ============================================================
    
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
    
    # ============================================================
    # 01.07. HELPER FUNCTIONS
    # ============================================================
    
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
    
    # Feature column 2 is the deposited name; column 3 is assay type.
    bacterial_probe_names <- function(features) {
      if (ncol(features) < 2L) stop("Expected feature IDs and names in columns 1 and 2.")
      trimws(as.character(features[[2]]))
    }
    probe_group <- function(x) {
      trimws(sub("_probe[0-9]+$", "", x, ignore.case=TRUE))
    }
    # Explicit deposited bacterial probe-group allowlist (19 groups).
    # Other probe groups, including viruses and controls, NEVER enter bacterial totals.
    bacterial_groups_allowed <- c("Alistipes", "Bacteroides_Phocaeicola", "Butyricimonas", "Campylobacter", "Coprobacter", "Dialister", "Enterobacteriaceae", "Fusobacterium", "Helicobacter", "Intestinimonas", "Lachnospiraceae", "Leptotrichia", "Neisseria", "Parabacteroides", "Parvimonas", "Peptostreptococcus", "Porphyromonas", "Prevotellaceae", "Streptococcaceae")
    valid_probe <- function(x) {
      !is.na(x) & grepl("_probe[0-9]+$", x, ignore.case=TRUE) &
        tolower(probe_group(x)) %in% tolower(bacterial_groups_allowed)
    }
    
    classify_bacterial_features <- function(features, host_mask) {
      # Probe naming identifies deposited candidates, not independently verified taxa.
      valid_probe(bacterial_probe_names(features))
    }
    extract_taxon_label <- function(features_subset) {
      nm <- bacterial_probe_names(features_subset)
      if (!all(valid_probe(nm))) stop("Invalid bacterial probe names; inspect the deposited feature annotation.")
      probe_group(nm)
    }
    
    # 01.08. LOAD scRNA-seq REFERENCE
    # ============================================================
    
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
    
    # ============================================================
    # 01.09. GET REFERENCE COUNTS
    # ============================================================
    
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
    
    # ============================================================
    # 01.10. REFERENCE CELL-TYPE LABELS
    # ============================================================
    
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
    
    # ============================================================
    # 01.11. CLEAN / DOWNSAMPLE REFERENCE
    # ============================================================
    
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
    
    # ============================================================
    # 01.12. HOST MARKER GENES AND SIGNATURES
    # ============================================================
    
    host_marker_genes <- c(
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
    
    # 01.13. STORAGE
    # ============================================================
    
    
    spot_expression_list <- list()
    
    bacterial_spot_list <- list()
    
    bacterial_taxa_list <- list()
    
    rctd_long_list <- list()
    
    rctd_status_list <- list()
    
    qc_list <- list()
    
    availability_list <- list()
    
    manifest_list <- list()
    
    # ============================================================
    # 01.14. PROCESS FOUR MATCHED SECTIONS
    # ============================================================
    
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
      
      # ----------------------------------------------------------
      # 14A. FIND FILES ASSOCIATED WITH THIS SAMPLE
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14B. EXTRACT NESTED ARCHIVES
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14C. FIND 10X MATRIX FILES
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14D. READ MATRIX
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14E. FEATURE IDS / SYMBOLS
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14F. IDENTIFY HOST FEATURES
      # ----------------------------------------------------------
      
      host_mask <- grepl(
        "^ENSG[0-9]+(\\.[0-9]+)?$",
        fid
      ) |
        symbol %in%
        host_symbols
      
      # Exclude all probe features (including GAPDH controls) from host RNA.
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
      
      # ----------------------------------------------------------
      # 14G. IDENTIFY BACTERIAL FEATURES
      # ----------------------------------------------------------
      
      bacterial_mask <-
        classify_bacterial_features(
          features,
          host_mask
        )
      
      if (!any(bacterial_mask)) stop(ss, ": no valid deposited bacterial probes found; do not reuse earlier bacterial results.")
      
      # Audit all probe rows, including excluded viral/control/unrecognized names.
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
      
      # ----------------------------------------------------------
      # 14H. COLLAPSE HOST COUNTS BY GENE SYMBOL
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14I. COLLAPSE BACTERIAL COUNTS
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14J. READ SPATIAL POSITIONS
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14K. QC
      # ----------------------------------------------------------
      
      host_UMI <-
        Matrix::colSums(
          host_counts
        )
      
      host_detected_genes <-
        Matrix::colSums(
          host_counts > 0
        )
      
      keep <-
        pos$in_tissue == 1 &
        host_UMI >= 100 &
        host_detected_genes >= 50
      
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
            host_detected_genes[keep]
          ),
        
        stringsAsFactors = FALSE
      )
      
      # ----------------------------------------------------------
      # 14L. SUBSET RETAINED SPOTS
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14M. HOST NORMALIZATION
      # ----------------------------------------------------------
      
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
      # ----------------------------------------------------------
      # 14O. HOST MARKER GENE AVAILABILITY
      # ----------------------------------------------------------
      
      # Names of requested markers are separate from per-spot QC counts.
      if (!is.character(host_marker_genes) || length(host_marker_genes)!=16L ||
          anyNA(host_marker_genes) || anyDuplicated(host_marker_genes))
        stop("Internal error: requested host marker names were overwritten.")
      genes_present <- intersect(
        host_marker_genes,
        rownames(host_norm)
      )
      
      genes_missing <- setdiff(
        host_marker_genes,
        rownames(host_norm)
      )
      
      availability_list[[ss]] <-
        data.frame(
          
          section = ss,
          patient = patient_i,
          tissue = tissue_i,
          
          genes_requested =
            length(host_marker_genes),
          
          genes_present =
            length(genes_present),
          
          missing =
            paste(
              genes_missing,
              collapse = ";"
            ),
          
          stringsAsFactors = FALSE
        )
      
      # Fail before expensive RCTD if a required standalone gene is genuinely absent.
      # KRT18 is optional as a signature member, but is not a standalone final-panel gene.
      required_panel_genes <- setdiff(host_marker_genes,"KRT18")
      missing_required <- setdiff(required_panel_genes,rownames(host_norm))
      if(length(missing_required)) {
        write.csv(do.call(rbind,availability_list),out("05_Host_gene_availability_MATCHED.csv"),row.names=FALSE)
        stop(ss,": required host markers missing: ",paste(missing_required,collapse=", "),
             ". Inspect deposited host symbols; no zero-value imputation was performed.")
      }
      cat(
        "Host-marker genes available:",
        length(genes_present),
        "/",
        length(host_marker_genes),
        "\n"
      )
      
      # ----------------------------------------------------------
      # 14P. EXTRACT HOST MARKER GENES
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14Q. HOST SIGNATURE SCORES
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14R. BACTERIAL SIGNAL
      # ----------------------------------------------------------
      
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
      
      # ----------------------------------------------------------
      # 14S. RCTD
      # ----------------------------------------------------------
      
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
      
      rctd_cache <- out(paste0("RCTD_Core_",ss,".rds"))
      rctd_signature <- digest::digest(list(query=query_counts,coordinates=coord_table,
                                            reference=reference_counts,labels=reference_labels,seed=config$seed,
                                            cores=config$max_cores,mode="full",spacexr=as.character(packageVersion("spacexr"))),algo="xxhash64")
      cached_fit <- if(file.exists(rctd_cache)) tryCatch(readRDS(rctd_cache),error=function(e) NULL) else NULL
      if(!is.null(cached_fit) && identical(cached_fit$signature,rctd_signature) &&
         is.matrix(cached_fit$weights) && all(is.finite(cached_fit$weights))) {
        message("Reusing checked RCTD fit: ",ss)
        rctd_result <- list(cached_weights=cached_fit$weights)
      } else {
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
        
        if(!inherits(rctd_result,"try-error") && !is.null(rctd_result@results$weights)) {
          tmp <- paste0(rctd_cache,".writing")
          saveRDS(list(signature=rctd_signature,weights=as.matrix(rctd_result@results$weights)),tmp)
          if(file.exists(rctd_cache)) {
            if(!file.copy(rctd_cache,paste0(rctd_cache,".previous"),overwrite=TRUE)) stop("Cannot back up RCTD cache")
            if(!file.remove(rctd_cache)) stop("Cannot replace RCTD cache")
          }
          if(!file.rename(tmp,rctd_cache)) stop("Could not finalize RCTD checkpoint")
        }
      }
      rm(cached_fit); invisible(gc())
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
        
        weights <- if(is.list(rctd_result) && !is.null(rctd_result$cached_weights))
          rctd_result$cached_weights else rctd_result@results$weights
        
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
          
          # Determine whether barcodes are rows or columns
          
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
        bacterial_keep
      )
      
      rm(list=intersect(c("fit","rctd_result","weights","weights_norm","host_norm","query_counts","spatial_rna"),ls()))
      invisible(gc())
    }
    
    # ============================================================
    # 01.15. SAVE PRIMARY MATCHED OUTPUTS
    # ============================================================
    
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
      out("05_Host_gene_availability_MATCHED.csv"),
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
    
    # ============================================================
    # 01.16. SAVE RCTD OUTPUT
    # ============================================================
    
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
    
    # ============================================================
    # 01.17. COMBINE HOST FEATURES
    # ============================================================
    
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
      out("10_Host_features_All_spots_MATCHED.csv"),
      row.names = FALSE
    )
    
    # ============================================================
    
    writeLines(capture.output(sessionInfo()),out("20_sessionInfo.txt"))
    invisible(TRUE)
  }
  spatial_stage_01(config)
  invisible(gc())
  # Final fixed panel: 15 genes plus six available-gene-average signatures.
  # This is a reused exploratory marker panel, not a newly preregistered panel.
  final_genes <- c("CEACAM5","CD68","EPCAM","KRT8","GNLY","LYZ","NKG7","CD79A",
                   "MS4A1","CD4","CD3E","PECAM1","VWF","CD8A","CD3D")
  signatures <- paste0("Signature__",c("Epithelial","Myeloid","Cytotoxic","B_cell","T_cell","Endothelial"))
  final_features <- c(final_genes,signatures)
  f <- read.csv(file.path(config$analysis_dir,"10_Host_features_All_spots_MATCHED.csv"),stringsAsFactors=FALSE)
  f <- f[f$feature %in% final_features,,drop=FALSE]
  if(!setequal(unique(f$feature),final_features)) stop("Missing fixed-panel features; inspect host gene availability")
  if(anyDuplicated(paste(f$section,f$barcode,f$feature))) stop("Duplicate host feature/spot rows")
  expected_sections <- c("CRC_03_Normal","CRC_03_Tumor","CRC_08_Normal","CRC_08_Tumor")
  if(!setequal(f$section,expected_sections) || any(!is.finite(f$expression) | f$expression<0))
    stop("Invalid final host feature values/sections")
  for(ss in expected_sections) {
    z <- f[f$section==ss,,drop=FALSE]
    if(!setequal(z$feature,final_features)) stop(ss,": incomplete final 21-feature panel")
    bars <- unique(z$barcode)
    for(nm in final_features) if(!setequal(z$barcode[z$feature==nm],bars))
      stop(ss,": incomplete host feature spot grid: ",nm)
  }
  write.csv(data.frame(feature=final_features,feature_type=c(rep("Gene",15),rep("Signature",6)),
                       Selection="Reused exploratory host marker panel; no classifier or predicted direction"),
            file.path(config$analysis_dir,"21_Final_ST_21_Feature_Panel.csv"),row.names=FALSE)
  write.csv(f,file.path(config$analysis_dir,"22_Final_21_Features_All_spots.csv"),row.names=FALSE)
  summary <- do.call(rbind,lapply(split(f,interaction(f$section,f$feature,drop=TRUE)),function(z)
    data.frame(section=z$section[1],patient=z$patient[1],tissue=z$tissue[1],feature=z$feature[1],
               N_spots=nrow(z),Mean=mean(z$expression),Median=median(z$expression),Detected_percent=100*mean(z$expression>0))))
  write.csv(summary,file.path(config$analysis_dir,"26_Final_21_Features_Section_Summary.csv"),row.names=FALSE)
  b <- read.csv(file.path(config$analysis_dir,"06_Bacterial_signal_All_spots_MATCHED.csv"))
  r <- read.csv(file.path(config$analysis_dir,"09_RCTD_All_spot_celltype_proportions_MATCHED.csv"))
  key <- function(x) paste(x$section,x$barcode,sep="::")
  if(anyDuplicated(key(b)) || anyDuplicated(paste(key(r),r$cell_type))) stop("Duplicate core spot IDs")
  for(z in list(f,r)) {
    ix <- match(key(z),key(b))
    if(anyNA(ix) || any(abs(z$X-b$X[ix])>1e-6 | abs(z$Y-b$Y[ix])>1e-6)) stop("Core barcode/coordinate mismatch")
  }
  # Only the essential maps are produced: cell proportions, fixed host features and total bacterial signal.
  library(ggplot2)
  map_plot <- function(z,val,facet=NULL,title) {
    z$Value <- z[[val]]
    p <- ggplot(z,aes(X,-Y,colour=Value))+geom_point(size=.45)+coord_equal()+
      scale_colour_viridis_c(option="C",na.value="grey80")+theme_void(base_family="Arial")+
      labs(title=title,colour=NULL)+theme(plot.title=element_text(size=11),strip.text=element_text(size=8))
    if(!is.null(facet)) p <- p+facet_wrap(stats::as.formula(paste("~",facet)),ncol=3,scales="fixed")
    p
  }
  for(ss in unique(b$section)) {
    bb <- b[b$section==ss,];ff <- f[f$section==ss,];rr <- r[r$section==ss,]
    dir.create(file.path(config$figure_dir,ss),recursive=TRUE,showWarnings=FALSE)
    p <- map_plot(rr,"proportion","cell_type",paste(ss,"RCTD proportions"))
    ggsave(file.path(config$figure_dir,ss,"Core_RCTD.png"),p,width=9,height=6,dpi=300,bg="white")
    p <- map_plot(ff,"expression","feature",paste(ss,"host marker expression; shared display scale"))
    ggsave(file.path(config$figure_dir,ss,"Core_Host_Features.png"),p,width=9,height=16,dpi=300,bg="white",limitsize=FALSE)
    bb$Log_count <- log1p(bb$bacterial_total_raw)
    p <- map_plot(bb,"Log_count",title=paste(ss,"ln(1 + total bacterial probe counts)"))
    ggsave(file.path(config$figure_dir,ss,"Core_Bacterial_Signal.png"),p,width=5,height=5,dpi=300,bg="white")
  }
  write.csv(data.frame(Extraction_version="feature_name_column2_bacteria_only_v2",
                       Feature_name_column=2L,Feature_type_column=3L,Bacterial_group_allowlist_size=19L,
                       Interpretation="Deposited bacteria-only probe labels, not confirmed taxa/viability",
                       Completed=as.character(Sys.time())),file.path(config$analysis_dir,"Bacterial_extraction_VALIDATED.csv"),row.names=FALSE)
  writeLines(c("Fixed exploratory marker panel; no classifier and no assumed tumor-change direction.",
               "Stage 04 fix: character host_marker_genes is distinct from numeric per-spot host_detected_genes.",
               "Standalone genes exclude KRT18; the original epithelial signature still uses KRT18 when available.",
               "Signature means use available genes. Shared display scales are for plotting, not evidence of comparable biology.",
               "Only two patients; descriptive spatial maps. Statistical spatial tests are run in stage 05.",
               "Old independent-spot P/FDR and classifier-direction validation outputs are not used."),
             file.path(config$analysis_dir,"Core_Methods_Notes.txt"))
  files <- c("04_QC_summary_MATCHED.csv","05_Host_gene_availability_MATCHED.csv",
             "06_Bacterial_signal_All_spots_MATCHED.csv","08_RCTD_status_MATCHED.csv",
             "09_RCTD_All_spot_celltype_proportions_MATCHED.csv","21_Final_ST_21_Feature_Panel.csv",
             "22_Final_21_Features_All_spots.csv","26_Final_21_Features_Section_Summary.csv",
             "Bacterial_extraction_VALIDATED.csv","20_sessionInfo.txt","Core_Methods_Notes.txt")
  saveRDS(list(revision="core-v1",host_marker_fix="marker-qc-separation-v1",files=files,md5=unname(tools::md5sum(file.path(config$analysis_dir,files))),
               reference_file=config$reference_file,reference_md5=unname(tools::md5sum(config$reference_file)),
               completed=Sys.time()),complete)
  cat("Stage 04 complete. Run 04_Checker.R, then 05_Bacteria_Host_Association.R.\n")
}
run_spatial_workflow(spatial_config)
