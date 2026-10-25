# GSE334323: matched spatial preprocessing, RCTD and descriptive host/bacterial maps.
# Requires local GEO files and the dmcable/spacexr API. No installation or downloads.
local({
  output_root <- Sys.getenv("CRC_SPATIAL_OUTPUT", "E:/CRC/Results/Results/04-Spatial_Analysis_of_Selected_Features")
  config <- list(
    archive = Sys.getenv("CRC_SPATIAL_ARCHIVE", "E:/CRC/Results/GSE334323_ST/GSE334323_RAW.tar"),
    reference_file = file.path(Sys.getenv("CRC_ANNOTATION_OUTPUT", "E:/CRC/Results/Results/02_Clustering_Markers_Annotation"), "05_annotated.rds"),
    raw_dir = file.path(output_root, "Matched_RAW_Extracted"),
    figure_dir = file.path(output_root, "Figures"), seed = 42L, max_cores = 2L)
  packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "spacexr", "digest")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required packages: ", paste(missing, collapse = ", "))
  if (!all(c("Reference", "SpatialRNA", "create.RCTD", "run.RCTD") %in% getNamespaceExports("spacexr")))
    stop("Requires the dmcable/spacexr API")
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan("sequential")
    on.exit(future::plan(old_plan), add = TRUE)
  }
  for (f in c(config$archive, config$reference_file)) if (!file.exists(f)) stop("Missing input: ", f)
  for (d in c(output_root, config$raw_dir, config$figure_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(d)) stop("Cannot create directory: ", d)
  }
  out <- function(name) file.path(output_root, name)
  write_tab <- function(z, name) utils::write.csv(z, out(name), row.names = FALSE)
  complete <- out("Core_Spatial_Complete.rds")
  if (file.exists(complete) && !file.remove(complete)) stop("Cannot invalidate completion record")
  set.seed(config$seed)
  meta <- data.frame(GSM = c("GSM9785428", "GSM9785429", "GSM9785432", "GSM9785433"),
    section = c("CRC_03_Normal", "CRC_03_Tumor", "CRC_08_Normal", "CRC_08_Tumor"),
    patient = rep(c("CRC_03", "CRC_08"), each = 2), tissue = rep(c("Normal", "Tumor"), 2), pair_status = "Matched")
  write_tab(meta, "00_Matched_sample_metadata.csv")

  extract_archive <- function(path, destination) {
    zip <- grepl("\\.zip$", path, ignore.case = TRUE)
    entries <- if (zip) utils::unzip(path, list = TRUE)$Name else utils::untar(path, list = TRUE)
    safe_names <- gsub("\\", "/", entries, fixed = TRUE)
    if (!length(entries) || any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", safe_names))) stop("Empty archive or unsafe archive paths")
    if (zip) utils::unzip(path, exdir = destination) else utils::untar(path, exdir = destination)
    entries
  }
  entries <- extract_archive(config$archive, config$raw_dir)
  matched <- Reduce(`|`, lapply(c(meta$GSM, meta$section), function(x) grepl(x, entries, fixed = TRUE)))
  write_tab(data.frame(archive_entry = entries, matched = matched), "01_RAW_TAR_inventory.csv")
  all_files <- list.files(config$raw_dir, recursive = TRUE, full.names = TRUE)
  if (!length(all_files)) stop("No extracted files")
  write_tab(data.frame(path = all_files, size_bytes = file.info(all_files)$size), "02_Extracted_file_inventory.csv")
  read_text <- function(path, reader, ...) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con), add = TRUE)
    reader(con, ...)
  }
  read_tsv <- function(path) read_text(path, utils::read.delim, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  choose <- function(files, pattern, section) {
    hits <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
    if (length(hits) != 1L) stop(section, ": expected one file matching ", pattern, "; found ", length(hits))
    hits
  }
  collapse_counts <- function(mat, labels) {
    groups <- unique(labels)
    a <- Matrix::sparseMatrix(i = match(labels, groups), j = seq_along(labels), x = 1,
                              dims = c(length(groups), length(labels)))
    result <- as(a %*% mat, "dgCMatrix")
    dimnames(result) <- list(groups, colnames(mat))
    result
  }
  allowed <- c("Alistipes", "Bacteroides_Phocaeicola", "Butyricimonas", "Campylobacter", "Coprobacter",
    "Dialister", "Enterobacteriaceae", "Fusobacterium", "Helicobacter", "Intestinimonas", "Lachnospiraceae",
    "Leptotrichia", "Neisseria", "Parabacteroides", "Parvimonas", "Peptostreptococcus", "Porphyromonas",
    "Prevotellaceae", "Streptococcaceae")
  probe_group <- function(x) trimws(sub("_probe[0-9]+$", "", x, ignore.case = TRUE))
  valid_probe <- function(x) !is.na(x) & grepl("_probe[0-9]+$", x, ignore.case = TRUE) &
    tolower(probe_group(x)) %in% tolower(allowed)
  markers <- c("KRT18", "CEACAM5", "CD68", "EPCAM", "KRT8", "GNLY", "LYZ", "NKG7", "CD79A", "MS4A1",
               "CD4", "CD3E", "PECAM1", "VWF", "CD8A", "CD3D")
  signatures <- list(Epithelial = c("EPCAM", "KRT8", "KRT18", "CEACAM5"), Myeloid = c("CD68", "LYZ"),
    Cytotoxic = c("NKG7", "GNLY"), B_cell = c("MS4A1", "CD79A"), T_cell = c("CD3D", "CD3E", "CD4", "CD8A"),
    Endothelial = c("PECAM1", "VWF"))

  # Reference counts and published broad labels; at most 500 cells per eligible type.
  ref_obj <- readRDS(config$reference_file)
  if (!inherits(ref_obj, "Seurat")) stop("Invalid reference object")
  DefaultAssay(ref_obj) <- "RNA"
  if (inherits(ref_obj[["RNA"]], "Assay5")) {
    layers <- SeuratObject::Layers(ref_obj[["RNA"]], search = "^counts")
    if (!length(layers)) stop("No reference counts")
    if (length(layers) > 1L) ref_obj <- SeuratObject::JoinLayers(ref_obj, assay = "RNA", layers = layers, new = "counts")
    counts <- SeuratObject::GetAssayData(ref_obj, assay = "RNA", layer = "counts")
  } else counts <- SeuratObject::GetAssayData(ref_obj, assay = "RNA", slot = "counts")
  counts <- as(counts, "dgCMatrix")
  host_symbols <- rownames(counts)
  choices <- intersect(c("analysis_cell_type", "Cell_type", "cell_type", "GEO_CellType", "SingleR_CellType"), names(ref_obj[[]]))
  label_column <- if (length(choices)) choices[1] else "active.ident"
  labels <- if (length(choices)) as.character(ref_obj[[]][[label_column]]) else as.character(Idents(ref_obj))
  names(labels) <- colnames(ref_obj)
  valid <- names(labels)[!is.na(labels) & nzchar(labels) & !tolower(labels) %in% c("unknown", "unassigned", "na")]
  eligible <- names(which(table(labels[valid]) >= 25L))
  valid <- valid[labels[valid] %in% eligible]
  set.seed(config$seed)
  selected <- unlist(lapply(split(valid, labels[valid]), function(x) if (length(x) > 500L) sample(x, 500L) else x), use.names = FALSE)
  if (!length(selected)) stop("No eligible reference cells")
  reference_counts <- counts[, selected, drop = FALSE]
  reference_labels <- factor(labels[selected])
  reference_nUMI <- Matrix::colSums(reference_counts)
  good <- reference_nUMI > 0
  reference_counts <- reference_counts[, good, drop = FALSE]
  reference_labels <- droplevels(reference_labels[good])
  reference_nUMI <- reference_nUMI[good]
  reference <- spacexr::Reference(reference_counts, reference_labels, reference_nUMI)
  cat("Reference label source:", label_column, "\n")
  print(table(reference_labels))
  rm(ref_obj, counts)
  invisible(gc())
  host_list <- bacteria_list <- taxa_list <- rctd_list <- qc_list <- availability_list <- manifest_list <- status_list <- list()

  for (i in seq_len(nrow(meta))) {
    ss <- meta$section[i]; gsm <- meta$GSM[i]
    message("Processing ", ss)
    files <- all_files[grepl(gsm, all_files, fixed = TRUE) | grepl(ss, all_files, fixed = TRUE)]
    nested <- files[grepl("\\.(tar|tar\\.gz|tgz|zip)$", files, ignore.case = TRUE)]
    if (length(nested)) {
      destination <- file.path(config$raw_dir, paste0("nested_", ss))
      dir.create(destination, recursive = TRUE, showWarnings = FALSE)
      for (nf in nested) extract_archive(nf, destination)
      files <- unique(c(files, list.files(destination, recursive = TRUE, full.names = TRUE)))
    }
    filtered <- files[grepl("filtered", files, ignore.case = TRUE)]
    candidates <- if (any(grepl("matrix.*\\.mtx(\\.gz)?$", filtered, ignore.case = TRUE))) filtered else files
    mf <- choose(candidates, "matrix.*\\.mtx(\\.gz)?$", ss)
    ff <- choose(candidates, "(features|genes).*\\.tsv(\\.gz)?$", ss)
    bf <- choose(candidates, "barcodes.*\\.tsv(\\.gz)?$", ss)
    pf <- choose(files, "tissue_positions(_list)?\\.csv(\\.gz)?$", ss)
    manifest_list[[ss]] <- data.frame(section = ss, GSM = gsm, matrix = mf, features = ff, barcodes = bf, positions = pf)
    features <- read_tsv(ff)
    barcodes <- as.character(read_tsv(bf)[[1]])
    mat <- as(read_text(mf, Matrix::readMM), "dgCMatrix")
    if (ncol(features) < 2L || nrow(mat) != nrow(features) || ncol(mat) != length(barcodes) ||
        anyNA(barcodes) || anyDuplicated(barcodes) || any(!nzchar(barcodes))) stop(ss, ": invalid matrix annotations")
    if (any(!is.finite(mat@x)) || any(mat@x < 0) || any(mat@x != floor(mat@x))) stop(ss, ": invalid UMI counts")
    colnames(mat) <- barcodes
    fid <- as.character(features[[1]]); symbol <- as.character(features[[2]])
    missing_symbol <- is.na(symbol) | !nzchar(symbol)
    symbol[missing_symbol] <- fid[missing_symbol]
    probe <- !is.na(symbol) & grepl("_probe[0-9]+$", symbol, ignore.case = TRUE)
    host <- (grepl("^ENSG[0-9]+(\\.[0-9]+)?$", fid) | symbol %in% host_symbols) & !probe & !is.na(symbol) & nzchar(symbol)
    if (sum(host) < 1000L) stop(ss, ": insufficient host features")
    bacterial <- valid_probe(trimws(as.character(features[[2]])))
    if (!any(bacterial)) stop(ss, ": no allowed bacterial probes")
    nm <- trimws(symbol[probe]); group <- probe_group(nm); included <- valid_probe(nm)
    viral <- grepl("^(HAdV|HBV|HCV|HPV|HSV|HTLV|EBV|VZV)([-_]|$)", group, ignore.case = TRUE)
    classification <- ifelse(included, "Bacterial", ifelse(viral, "Viral_excluded",
      ifelse(grepl("gapdh|control|blank|negative|positive", group, ignore.case = TRUE), "Control_excluded", "Unrecognized_excluded")))
    write_tab(data.frame(section = ss, feature_id = fid[probe], feature_name = nm, deposited_group = group,
                         classification = classification, included_in_bacterial_total = included), paste0("Bacterial_Probe_Audit_", ss, ".csv"))
    host_counts <- collapse_counts(mat[host, , drop = FALSE], symbol[host])
    bacterial_counts <- collapse_counts(mat[bacterial, , drop = FALSE], probe_group(trimws(as.character(features[[2]][bacterial]))))
    first <- read_text(pf, readLines, n = 1L, warn = FALSE)
    has_header <- grepl("barcode", first, ignore.case = TRUE)
    pos <- read_text(pf, utils::read.csv, header = has_header, stringsAsFactors = FALSE)
    if (!has_header) {
      if (ncol(pos) != 6L) stop(ss, ": invalid positions format")
      names(pos) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    }
    needed <- c("barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres")
    if (!all(needed %in% names(pos))) stop(ss, ": missing coordinates")
    if (anyNA(pos[needed]) || anyDuplicated(pos$barcode) || any(!is.finite(pos$pxl_row_in_fullres)) ||
        any(!is.finite(pos$pxl_col_in_fullres))) stop(ss, ": invalid coordinates")
    ix <- match(barcodes, pos$barcode)
    if (anyNA(ix)) stop(ss, ": missing spot coordinates")
    pos <- pos[ix, ]
    umi <- Matrix::colSums(host_counts); detected <- Matrix::colSums(host_counts > 0)
    keep <- pos$in_tissue == 1 & umi >= 100 & detected >= 50
    if (!any(keep)) stop(ss, ": no spots pass QC")
    qc_list[[ss]] <- data.frame(section = ss, patient = meta$patient[i], tissue = meta$tissue[i],
      imported_spots = ncol(host_counts), retained_spots = sum(keep), host_features = nrow(host_counts),
      bacterial_features = nrow(bacterial_counts), median_host_UMI = median(umi[keep]), median_host_genes = median(detected[keep]))
    host_keep <- host_counts[, keep, drop = FALSE]; bacterial_keep <- bacterial_counts[, keep, drop = FALSE]
    X <- pos$pxl_col_in_fullres[keep]; Y <- pos$pxl_row_in_fullres[keep]
    norm <- host_keep %*% Matrix::Diagonal(x = 10000 / Matrix::colSums(host_keep))
    norm@x <- log1p(norm@x); dimnames(norm) <- dimnames(host_keep)
    present <- intersect(markers, rownames(norm)); absent <- setdiff(markers, rownames(norm))
    availability_list[[ss]] <- data.frame(section = ss, patient = meta$patient[i], tissue = meta$tissue[i],
      genes_requested = length(markers), genes_present = length(present), missing = paste(absent, collapse = ";"))
    if (length(setdiff(setdiff(markers, "KRT18"), present))) {
      write_tab(do.call(rbind, availability_list), "05_Host_gene_availability_MATCHED.csv")
      stop(ss, ": required host markers absent; no zero imputation")
    }
    spot_rows <- function() data.frame(section = ss, patient = meta$patient[i], tissue = meta$tissue[i],
                                       barcode = colnames(host_keep), X = X, Y = Y)
    host_list[[ss]] <- do.call(rbind, c(lapply(present, function(g)
      data.frame(spot_rows(), feature = g, feature_type = "Gene", expression = as.numeric(norm[g, ]))),
      lapply(names(signatures), function(n) {
        g <- intersect(signatures[[n]], rownames(norm))
        if (!length(g)) stop(ss, ": empty signature ", n)
        data.frame(spot_rows(), feature = paste0("Signature__", n), feature_type = "Signature", expression = as.numeric(Matrix::colMeans(norm[g, , drop = FALSE])))
      })))
    total <- Matrix::colSums(bacterial_keep)
    bacteria_list[[ss]] <- data.frame(spot_rows(), pair_status = "Matched", bacterial_total_raw = total,
      bacterial_detected_taxa = Matrix::colSums(bacterial_keep > 0), bacterial_positive = total > 0)
    taxa_list[[ss]] <- data.frame(section = ss, patient = meta$patient[i], tissue = meta$tissue[i],
      taxon = rownames(bacterial_keep), total_counts = as.numeric(Matrix::rowSums(bacterial_keep)),
      detected_spots = as.numeric(Matrix::rowSums(bacterial_keep > 0)))

    # Cache signatures match the original workflow, including reference selection.
    common <- intersect(rownames(host_keep), rownames(reference_counts))
    if (length(common) < 500L) warning(ss, ": low reference gene overlap")
    query <- host_keep[common, , drop = FALSE]
    query_umi <- Matrix::colSums(query); good <- query_umi > 0
    query <- query[, good, drop = FALSE]; query_umi <- query_umi[good]
    coords <- data.frame(x = X[good], y = Y[good], row.names = colnames(query))
    cache <- out(paste0("RCTD_Core_", ss, ".rds"))
    signature <- digest::digest(list(query = query, coordinates = coords, reference = reference_counts,
      labels = reference_labels, seed = config$seed, cores = config$max_cores, mode = "full",
      spacexr = as.character(utils::packageVersion("spacexr"))), algo = "xxhash64")
    cached <- if (file.exists(cache)) tryCatch(readRDS(cache), error = function(e) NULL) else NULL
    attempt <- tryCatch({
      if (!is.null(cached) && identical(cached$signature, signature) && is.matrix(cached$weights) &&
          all(is.finite(cached$weights))) {
        message("Reusing RCTD cache: ", ss)
        weights <- cached$weights
      } else {
        query_obj <- spacexr::SpatialRNA(coords, query, query_umi)
        fit <- spacexr::create.RCTD(query_obj, reference, max_cores = config$max_cores)
        fit <- spacexr::run.RCTD(fit, doublet_mode = "full")
        weights <- as.matrix(fit@results$weights)
        if (!length(weights) || any(!is.finite(weights))) stop("Invalid RCTD weights")
        tmp <- paste0(cache, ".writing")
        saveRDS(list(signature = signature, weights = weights), tmp)
        if (file.exists(cache)) {
          if (!file.copy(cache, paste0(cache, ".previous"), overwrite = TRUE)) stop("Cannot back up cache")
          if (!file.remove(cache)) stop("Cannot replace cache")
        }
        if (!file.rename(tmp, cache)) stop("Cannot finalize cache")
      }
      if (sum(colnames(query) %in% colnames(weights)) > sum(colnames(query) %in% rownames(weights))) weights <- t(weights)
      if (is.null(rownames(weights)) || is.null(colnames(weights)) || anyDuplicated(rownames(weights)) ||
          anyDuplicated(colnames(weights)) || any(!is.finite(weights)) || any(weights < 0)) stop("Invalid RCTD weight identifiers/values")
      bars <- colnames(query)[colnames(query) %in% rownames(weights)]
      if (!length(bars)) stop("No RCTD spot matches")
      weights <- weights[bars, , drop = FALSE]
      sums <- rowSums(weights); valid <- is.finite(sums) & sums > 0
      weights[valid, ] <- weights[valid, , drop = FALSE] / sums[valid]
      weights[!valid, ] <- 0
      write_tab(data.frame(barcode = rownames(weights), weights, check.names = FALSE), paste0("RCTD_weights_", ss, ".csv"))
      rctd_list[[ss]] <- do.call(rbind, lapply(colnames(weights), function(ct) data.frame(
        section = ss, patient = meta$patient[i], tissue = meta$tissue[i], pair_status = "Matched", barcode = bars,
        X = coords[bars, "x"], Y = coords[bars, "y"], cell_type = ct, proportion = as.numeric(weights[, ct]))))
      data.frame(section = ss, status = "OK", query_spots = ncol(query), retained_spots = length(bars),
                 retention_percent = 100 * length(bars) / ncol(query), message = "")
    }, error = function(e) data.frame(section = ss, status = "FAILED", query_spots = ncol(query),
      retained_spots = NA_integer_, retention_percent = NA_real_, message = conditionMessage(e)))
    status_list[[ss]] <- attempt
    rm(list = intersect(c("mat", "host_counts", "bacterial_counts", "host_keep", "bacterial_keep", "norm",
      "query", "query_obj", "fit", "cached", "weights"), ls()))
    invisible(gc())
  }
  write_tab(do.call(rbind, manifest_list), "03_Input_manifest_MATCHED.csv")
  write_tab(do.call(rbind, qc_list), "04_QC_summary_MATCHED.csv")
  write_tab(do.call(rbind, availability_list), "05_Host_gene_availability_MATCHED.csv")
  b <- do.call(rbind, bacteria_list); b$log_bacterial <- log1p(b$bacterial_total_raw)
  write_tab(b, "06_Bacterial_signal_All_spots_MATCHED.csv")
  write_tab(do.call(rbind, taxa_list), "07_Bacterial_taxa_MATCHED.csv")
  status <- do.call(rbind, status_list)
  write_tab(status, "08_RCTD_status_MATCHED.csv")
  writeLines(capture.output(utils::sessionInfo()), out("20_sessionInfo.txt"))
  if (nrow(status) != 4L || any(status$status != "OK")) stop("Incomplete RCTD; inspect 08_RCTD_status_MATCHED.csv")
  r <- do.call(rbind, rctd_list)
  write_tab(r, "09_RCTD_All_spot_celltype_proportions_MATCHED.csv")
  f <- do.call(rbind, host_list)
  write_tab(f, "10_Host_features_All_spots_MATCHED.csv")
  final_genes <- c("CEACAM5", "CD68", "EPCAM", "KRT8", "GNLY", "LYZ", "NKG7", "CD79A", "MS4A1", "CD4", "CD3E", "PECAM1", "VWF", "CD8A", "CD3D")
  final_features <- c(final_genes, paste0("Signature__", names(signatures)))
  f <- f[f$feature %in% final_features, , drop = FALSE]
  key <- function(z) paste(z$section, z$barcode, sep = "::")
  if (!setequal(f$section, meta$section) || any(!is.finite(f$expression)) || any(f$expression < 0) ||
      anyDuplicated(paste(key(f), f$feature)) || anyDuplicated(key(b)) || anyDuplicated(paste(key(r), r$cell_type))) stop("Invalid core spot/feature rows")
  for (ss in meta$section) {
    z <- f[f$section == ss, ]
    for (nm in final_features) if (!setequal(z$barcode[z$feature == nm], unique(z$barcode))) stop(ss, ": incomplete feature grid")
  }
  for (z in list(f, r)) {
    ix <- match(key(z), key(b))
    if (anyNA(ix) || any(!is.finite(z$X)) || any(!is.finite(z$Y)) ||
        any(abs(z$X - b$X[ix]) > 1e-6 | abs(z$Y - b$Y[ix]) > 1e-6)) stop("Core coordinate mismatch")
  }
  write_tab(data.frame(feature = final_features, feature_type = c(rep("Gene", 15), rep("Signature", 6)),
    Selection = "Reused exploratory host marker panel; no classifier or predicted direction"), "21_Final_ST_21_Feature_Panel.csv")
  write_tab(f, "22_Final_21_Features_All_spots.csv")
  summary <- do.call(rbind, lapply(split(f, interaction(f$section, f$feature, drop = TRUE)), function(z)
    data.frame(section = z$section[1], patient = z$patient[1], tissue = z$tissue[1], feature = z$feature[1],
      N_spots = nrow(z), Mean = mean(z$expression), Median = median(z$expression), Detected_percent = 100 * mean(z$expression > 0))))
  write_tab(summary, "26_Final_21_Features_Section_Summary.csv")
  map_plot <- function(z, value, facet = NULL, title) {
    z$Value <- z[[value]]
    p <- ggplot(z, aes(X, -Y, colour = Value)) + geom_point(size = 0.45) + coord_equal() +
      scale_colour_viridis_c(option = "C", na.value = "grey80") + theme_void(base_family = "Arial") +
      labs(title = title, colour = NULL) + theme(plot.title = element_text(size = 11), strip.text = element_text(size = 8))
    if (!is.null(facet)) p <- p + facet_wrap(stats::reformulate(facet), ncol = 3, scales = "fixed")
    p
  }
  for (ss in meta$section) {
    destination <- file.path(config$figure_dir, ss)
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    save_map <- function(p, name, w, h) ggsave(file.path(destination, name), p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
    save_map(map_plot(r[r$section == ss, ], "proportion", "cell_type", paste(ss, "RCTD proportions")), "Core_RCTD.png", 9, 6)
    save_map(map_plot(f[f$section == ss, ], "expression", "feature", paste(ss, "host markers; shared display scale")), "Core_Host_Features.png", 9, 16)
    bb <- b[b$section == ss, ]; bb$Log_count <- log1p(bb$bacterial_total_raw)
    save_map(map_plot(bb, "Log_count", title = paste(ss, "ln(1 + bacterial probe counts)")), "Core_Bacterial_Signal.png", 5, 5)
  }
  write_tab(data.frame(Extraction_version = "feature_name_column2_bacteria_only_v2", Feature_name_column = 2L,
    Feature_type_column = 3L, Bacterial_group_allowlist_size = 19L,
    Interpretation = "Deposited bacteria-only probe labels, not confirmed taxa/viability", Completed = as.character(Sys.time())), "Bacterial_extraction_VALIDATED.csv")
  writeLines(c(
    "Exploratory 15-gene/six-signature panel; no classifier or assumed direction.",
    "Host marker names and per-spot detected-gene counts are separate variables.",
    "KRT18 is excluded from standalone panels but used in the epithelial signature when present.",
    "Signature scores are available-gene means of log-normalized host expression.",
    "Map colour scales are shared within each section, not across sections.",
    "Bacterial values are deposited probe-group counts, not verified organisms, viability or normalized bacterial abundance.",
    "RCTD weights are reference-dependent mixture estimates; zero-weight rows do not provide resolved composition.",
    "Only two patients: maps are descriptive, not independent-spot significance tests.",
    "No classifier, trajectory or CellChat analysis is performed."), out("Core_Methods_Notes.txt"))
  files <- c("04_QC_summary_MATCHED.csv", "05_Host_gene_availability_MATCHED.csv", "06_Bacterial_signal_All_spots_MATCHED.csv",
    "08_RCTD_status_MATCHED.csv", "09_RCTD_All_spot_celltype_proportions_MATCHED.csv", "21_Final_ST_21_Feature_Panel.csv",
    "22_Final_21_Features_All_spots.csv", "26_Final_21_Features_Section_Summary.csv", "Bacterial_extraction_VALIDATED.csv",
    "20_sessionInfo.txt", "Core_Methods_Notes.txt")
  saveRDS(list(revision = "core-v1", host_marker_fix = "marker-qc-separation-v1", files = files,
    md5 = unname(tools::md5sum(out(files))), reference_file = config$reference_file,
    reference_md5 = unname(tools::md5sum(config$reference_file)), completed = Sys.time()), complete)
  cat("\nSpatial stage complete\nOutput:", output_root, "\n")
})
invisible(gc())
