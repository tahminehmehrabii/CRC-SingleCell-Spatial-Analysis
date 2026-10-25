# GSE334323: exploratory bacteria-host associations in four sections from two patients.
# SpatialPack tests on ranks are approximate; enrichment and pooled rho are descriptive.
local({
  INPUT_DIR <- Sys.getenv("CRC_SPATIAL_OUTPUT", "E:/CRC/Results/Results/04-Spatial_Analysis_of_Selected_Features")
  OUTPUT_DIR <- Sys.getenv("CRC_ASSOCIATION_OUTPUT", "E:/CRC/Results/Results/05-Bacteria_Host_Spatial_Association")
  RAW_ROOT <- Sys.getenv("CRC_SPATIAL_RAW", "E:/CRC/Results/GSE334323_ST")
  ARCHIVE <- Sys.getenv("CRC_SPATIAL_ARCHIVE", file.path(RAW_ROOT, "GSE334323_RAW.tar"))
  SECTIONS <- c("CRC_03_Normal", "CRC_03_Tumor", "CRC_08_Normal", "CRC_08_Tumor")
  GSM <- setNames(c("GSM9785428", "GSM9785429", "GSM9785432", "GSM9785433"), SECTIONS)
  TARGETS <- c("Myeloids", "Signature__Epithelial", "EPCAM")
  NCLASS <- c(13L, 8L, 20L); K <- 6L; MIN_SPOTS <- 30L
  MIN_POSITIVE <- 5L; MIN_COUNTS <- 10; HIGH_Q <- 0.75; SEED <- 42L
  TOP_N <- 10L; MAIN_EXCLUDE <- "Bacteroides_Phocaeicola"
  set.seed(SEED)
  packages <- c("Matrix", "RANN", "SpatialPack", "ggplot2", "cowplot", "digest")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required packages: ", paste(missing, collapse = ", "))
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OUTPUT_DIR)) stop("Cannot create output directory")
  out <- function(n) file.path(OUTPUT_DIR, n)
  tab <- function(d, n) utils::write.csv(d, out(n), row.names = FALSE, na = "")
  key <- function(d) paste(d$section, d$barcode, sep = "::")
  read_input <- function(n, required) {
    path <- file.path(INPUT_DIR, n)
    if (!file.exists(path)) stop("Missing Stage 04 input: ", path)
    d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(required %in% names(d))) stop("Missing columns in ", n, ": ", paste(setdiff(required, names(d)), collapse = ", "))
    d
  }
  verify_completion <- function(root, name) {
    path <- file.path(root, name)
    if (!file.exists(path)) stop("Upstream stage incomplete: ", path)
    z <- readRDS(path)
    actual <- unname(tools::md5sum(file.path(root, z$files)))
    if (!identical(z$revision, "core-v1") || !length(z$files) || anyNA(actual) || !identical(z$md5, actual)) stop("Upstream checksum mismatch")
    invisible(z)
  }
  verify_completion(INPUT_DIR, "Core_Spatial_Complete.rds")
  own_complete <- out("Core_Association_Complete.rds")
  if (file.exists(own_complete) && !file.remove(own_complete)) stop("Cannot invalidate completion record")
  cert <- read_input("Bacterial_extraction_VALIDATED.csv", "Extraction_version")
  if (nrow(cert) != 1L || cert$Extraction_version != "feature_name_column2_bacteria_only_v2") stop("Requires corrected Stage 04 bacterial extraction")
  b <- read_input("06_Bacterial_signal_All_spots_MATCHED.csv", c("section", "barcode", "X", "Y", "bacterial_total_raw"))
  f <- read_input("22_Final_21_Features_All_spots.csv", c("section", "barcode", "X", "Y", "feature", "expression"))
  r <- read_input("09_RCTD_All_spot_celltype_proportions_MATCHED.csv", c("section", "barcode", "X", "Y", "cell_type", "proportion"))
  for (n in c("b", "f", "r")) {
    z <- get(n)
    for (v in c("section", "barcode")) z[[v]] <- trimws(as.character(z[[v]]))
    if (anyNA(z[c("section", "barcode")]) || any(!nzchar(z$section)) || any(!nzchar(z$barcode)) ||
        !setequal(z$section, SECTIONS)) stop("Invalid IDs/sections: ", n)
    if (!is.numeric(z$X) || !is.numeric(z$Y) || any(!is.finite(z$X)) || any(!is.finite(z$Y))) stop("Invalid coordinates: ", n)
    assign(n, z)
  }
  if (anyDuplicated(key(b)) || anyDuplicated(paste(key(f), f$feature)) || anyDuplicated(paste(key(r), r$cell_type))) stop("Duplicate spot/variable IDs")
  if (!is.numeric(b$bacterial_total_raw) || any(!is.finite(b$bacterial_total_raw)) ||
      any(b$bacterial_total_raw < 0 | b$bacterial_total_raw != floor(b$bacterial_total_raw))) stop("Invalid bacterial counts")
  if (!is.numeric(f$expression) || !is.numeric(r$proportion)) stop("Host values must be numeric")
  if (any(r$proportion < 0 | r$proportion > 1, na.rm = TRUE)) stop("RCTD proportions outside [0,1]")
  if (!all(TARGETS[-1] %in% f$feature) || !"Myeloids" %in% r$cell_type) stop("Missing host targets")
  align <- function(z) {
    ix <- match(key(z), key(b))
    if (anyNA(ix) || any(abs(z$X - b$X[ix]) > 1e-6 | abs(z$Y - b$Y[ix]) > 1e-6)) stop("Spot/coordinate alignment mismatch")
    invisible(ix)
  }
  align(f); align(r)
  feature_names <- sort(unique(f$feature)); cell_types <- sort(unique(r$cell_type))
  tab(data.frame(Input = c("Bacterial", "Features", "RCTD"), Rows = c(nrow(b), nrow(f), nrow(r)),
    Unique_spots = c(length(unique(key(b))), length(unique(key(f))), length(unique(key(r))))), "00_Input_Alignment.csv")

  # Reconstruct bacterial probe-group counts and verify against Stage 04 totals.
  cache <- out("Input_Cache"); dir.create(cache, showWarnings = FALSE)
  raw_files <- unique(c(list.files(RAW_ROOT, recursive = TRUE, full.names = TRUE), list.files(cache, recursive = TRUE, full.names = TRUE)))
  locate <- function(ss, kind, extract = TRUE) {
    suffix <- switch(kind, features = "features\\.tsv(\\.gz)?$", barcodes = "barcodes\\.tsv(\\.gz)?$", matrix = "matrix\\.mtx(\\.gz)?$")
    hits <- raw_files[(grepl(GSM[[ss]], basename(raw_files), fixed = TRUE) | grepl(ss, basename(raw_files), fixed = TRUE)) &
                        grepl(suffix, basename(raw_files), ignore.case = TRUE)]
    if (!length(hits) && extract) {
      if (!file.exists(ARCHIVE)) stop("Missing raw archive: ", ss)
      entries <- utils::untar(ARCHIVE, list = TRUE)
      entries <- entries[grepl(GSM[[ss]], entries, fixed = TRUE)]
      safe_names <- gsub("\\", "/", entries, fixed = TRUE)
      if (!length(entries) || any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", safe_names))) stop("Missing/unsafe archive entries")
      utils::untar(ARCHIVE, files = entries, exdir = cache)
      raw_files <<- unique(c(raw_files, list.files(cache, recursive = TRUE, full.names = TRUE)))
      return(locate(ss, kind, FALSE))
    }
    if (!length(hits)) stop("Cannot locate ", kind, " for ", ss)
    hashes <- unname(tools::md5sum(hits))
    if (anyNA(hashes) || length(unique(hashes)) > 1L) stop("Conflicting/unreadable raw files: ", ss, " / ", kind)
    hits[1]
  }
  tsv <- function(path) utils::read.delim(path, header = FALSE, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  bacterial_list <- audit_list <- extraction_list <- list()
  for (ss in SECTIONS) {
    a <- read_input(paste0("Bacterial_Probe_Audit_", ss, ".csv"), c("feature_id", "feature_name", "deposited_group", "classification", "included_in_bacterial_total"))
    a <- a[tolower(as.character(a$included_in_bacterial_total)) %in% c("true", "1"), , drop = FALSE]
    if (!nrow(a) || anyNA(a[c("feature_id", "feature_name", "classification", "deposited_group")]) || any(a$classification != "Bacterial")) stop("Invalid bacterial probe audit")
    if (any(grepl("gene.expression|gapdh|control|blank|unknown|^(HAdV|HBV|HCV|HPV|HSV|HTLV|EBV|VZV)([-_]|$)", a$deposited_group, ignore.case = TRUE))) stop("Non-bacterial audited group")
    ft <- tsv(locate(ss, "features")); bar <- as.character(tsv(locate(ss, "barcodes"))[[1]])
    mm <- locate(ss, "matrix")
    con <- if (grepl("\\.gz$", mm)) gzfile(mm, "rt") else file(mm, "rt")
    mat <- tryCatch(Matrix::readMM(con), finally = close(con))
    if (nrow(mat) == length(bar) && ncol(mat) == nrow(ft)) mat <- t(mat)
    if (nrow(mat) != nrow(ft) || ncol(mat) != length(bar) || anyNA(bar) || anyDuplicated(bar)) stop("Matrix/barcode mismatch")
    mat <- as(mat, "CsparseMatrix")
    if (any(!is.finite(mat@x) | mat@x < 0 | mat@x != floor(mat@x))) stop("Invalid raw counts")
    if (ncol(ft) < 2L) stop("Missing feature-name column")
    raw_id <- paste(as.character(ft[[1]]), trimws(as.character(ft[[2]])), sep = "::")
    ai <- match(paste(a$feature_id, a$feature_name, sep = "::"), raw_id)
    if (anyNA(ai) || anyDuplicated(ai)) stop("Probe audit/raw identifiers disagree")
    bs <- b[b$section == ss, , drop = FALSE]; bix <- match(bs$barcode, bar)
    if (anyNA(bix)) stop("QC barcodes absent from raw matrix")
    groups <- sort(unique(a$deposited_group))
    gm <- vapply(groups, function(g) as.numeric(Matrix::colSums(mat[ai[a$deposited_group == g], bix, drop = FALSE])), numeric(nrow(bs)))
    if (any(rowSums(gm) != bs$bacterial_total_raw)) stop("Reconstructed totals disagree: ", ss)
    bacterial_list[[ss]] <- do.call(rbind, lapply(seq_along(groups), function(j)
      data.frame(section = ss, barcode = bs$barcode, X = bs$X, Y = bs$Y, bacterial_group = groups[j], bacterial_count = gm[, j])))
    audit_list[[ss]] <- a
    extraction_list[[ss]] <- data.frame(Section = ss, Probes = nrow(a), Groups = length(groups), QC_spots = nrow(bs),
      Reconstructed_total = sum(gm), Stage04_total = sum(bs$bacterial_total_raw), Exact_total_agreement = TRUE)
    rm(mat, gm); invisible(gc())
  }
  taxa <- do.call(rbind, bacterial_list); groups <- sort(unique(taxa$bacterial_group))
  tab(do.call(rbind, audit_list), "01_Bacterial_Probe_Audit.csv")
  tab(do.call(rbind, extraction_list), "02_Bacterial_Extraction_Check.csv")
  tab(taxa, "03_Bacterial_Groups_All_QC_Spots.csv")

  spatial_test <- function(x, y, coords, nclass) {
    ok <- is.finite(x) & is.finite(y) & apply(is.finite(coords), 1, all)
    x <- x[ok]; y <- y[ok]; coords <- coords[ok, , drop = FALSE]
    z <- data.frame(N_spots = length(x), Rho = NA_real_, ESS = NA_real_, DF = NA_real_, P_spatial = NA_real_, Status = "OK", Message = "")
    diagnostic <- NULL
    if (length(x) < MIN_SPOTS) z$Status <- "Too_few_spots"
    else if (length(unique(x)) < 2L || length(unique(y)) < 2L) z$Status <- "Constant_variable"
    else if (anyDuplicated(data.frame(coords))) z$Status <- "Duplicate_coordinates"
    else {
      z$Rho <- stats::cor(x, y, method = "spearman"); warns <- character()
      fit <- tryCatch(withCallingHandlers(SpatialPack::modified.ttest(rank(x, ties.method = "average"),
        rank(y, ties.method = "average"), coords, nclass = nclass), warning = function(w) {
          warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
        }), error = function(e) e)
      if (inherits(fit, "error")) { z$Status <- "Error"; z$Message <- conditionMessage(fit) }
      else {
        z$ESS <- as.numeric(fit$ESS)[1]; z$DF <- as.numeric(fit$dof)[1]
        if (length(warns)) { z$Status <- "Warning"; z$Message <- paste(unique(warns), collapse = " | ") }
        else if (!is.finite(z$ESS) || z$ESS <= 2 || !is.finite(z$DF) || z$DF <= 0 ||
                 length(fit$p.value) != 1L || !is.finite(fit$p.value) || fit$p.value < 0 || fit$p.value > 1) z$Status <- "Invalid_spatial_estimate"
        else z$P_spatial <- fit$p.value
        if (!is.null(fit$imoran)) {
          mo <- as.matrix(fit$imoran)
          if (ncol(mo) >= 2L && nrow(mo) == length(fit$upper.bounds)) diagnostic <- data.frame(
            Distance_upper = as.numeric(fit$upper.bounds), Moran_x = mo[, 1], Moran_y = mo[, 2], Pair_count = as.numeric(fit$card))
        }
      }
    }
    list(result = z, diagnostic = diagnostic)
  }
  overlap <- function(x, y) {
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    if (!length(x) || !any(x > 0)) return(data.frame(Bacterial_high_threshold = NA_real_, Host_high_threshold = NA_real_,
      High_bacteria = 0L, High_host = 0L, High_both = 0L, Expected_high_both = NA_real_, Enrichment = NA_real_))
    bt <- unname(stats::quantile(x[x > 0], HIGH_Q)); ht <- unname(stats::quantile(y, HIGH_Q))
    hb <- x > 0 & x >= bt; hh <- y >= ht; expected <- sum(hb) * sum(hh) / length(x)
    data.frame(Bacterial_high_threshold = bt, Host_high_threshold = ht, High_bacteria = sum(hb), High_host = sum(hh),
      High_both = sum(hb & hh), Expected_high_both = expected, Enrichment = if (expected > 0) sum(hb & hh) / expected else NA_real_)
  }
  jobs <- abundance <- nn_audit <- list()
  add_job <- function(family, ss, bac, target, z, x, y, low = FALSE, enrich = FALSE) {
    jobs[[length(jobs) + 1L]] <<- list(Family = family, Section = ss, Bacterial_group = bac, Host_target = target,
      coords = as.matrix(z[c("X", "Y")]), x = x, y = y, low = low, overlap = if (enrich) overlap(x, y) else NULL)
  }
  for (ss in SECTIONS) {
    bs <- b[b$section == ss, , drop = FALSE]
    for (ff in feature_names) {
      z <- f[f$section == ss & f$feature == ff, , drop = FALSE]
      add_job("Total_Feature", ss, "Total_bacterial_signal", ff, z, bs$bacterial_total_raw[match(z$barcode, bs$barcode)], z$expression)
    }
    for (ct in cell_types) {
      z <- r[r$section == ss & r$cell_type == ct, , drop = FALSE]
      add_job("Total_RCTD", ss, "Total_bacterial_signal", ct, z, bs$bacterial_total_raw[match(z$barcode, bs$barcode)], z$proportion)
    }
    hs <- r[r$section == ss & r$cell_type == "Myeloids", c("section", "barcode", "X", "Y", "proportion"), drop = FALSE]
    names(hs)[5] <- "Myeloids"
    for (ff in TARGETS[-1]) {
      z <- f[f$section == ss & f$feature == ff, , drop = FALSE]
      hs[[ff]] <- z$expression[match(hs$barcode, z$barcode)]
    }
    hs <- hs[order(hs$barcode), , drop = FALSE]
    if (nrow(hs) <= K || anyDuplicated(hs[c("X", "Y")])) stop("Invalid neighborhood coordinates: ", ss)
    nearest <- RANN::nn2(as.matrix(hs[c("X", "Y")]), k = K + 1L)
    ni <- t(vapply(seq_len(nrow(hs)), function(i) {
      ids <- nearest$nn.idx[i, ]; ids[ids != i][seq_len(K)]
    }, integer(K)))
    nn_audit[[ss]] <- data.frame(Section = ss, Spots = nrow(hs), K = K,
      Median_max_neighbor_distance = median(nearest$nn.dists[, K + 1L]), Max_neighbor_distance = max(nearest$nn.dists[, K + 1L]))
    for (g in groups) {
      db <- taxa[taxa$section == ss & taxa$bacterial_group == g, , drop = FALSE]
      if (!nrow(db)) stop("Bacterial group absent: ", g, " / ", ss)
      vals <- db$bacterial_count[match(hs$barcode, db$barcode)]
      if (anyNA(vals)) stop("Bacterial/RCTD barcode mismatch")
      low <- sum(vals) < MIN_COUNTS || sum(vals > 0) < MIN_POSITIVE
      abundance[[paste(ss, g)]] <- data.frame(Section = ss, Bacterial_group = g, N_spots = length(vals), Total_counts = sum(vals),
        Positive_spots = sum(vals > 0), Positive_percent = 100 * mean(vals > 0), Mean = mean(vals), Median = median(vals), Pass_information_filter = !low)
      for (target in TARGETS) {
        hv <- hs[[target]]
        local_host <- vapply(seq_len(nrow(hs)), function(i) {
          v <- hv[ni[i, ]]; if (all(is.finite(v))) mean(v) else NA_real_
        }, numeric(1))
        add_job("Taxon_SameSpot", ss, g, target, hs, vals, hv, low, TRUE)
        add_job("Taxon_Neighborhood", ss, g, target, hs, vals, local_host, low, TRUE)
      }
    }
  }
  tab(do.call(rbind, abundance), "04_Bacterial_Group_Abundance.csv")
  tab(do.call(rbind, nn_audit), "05_Neighborhood_Distance_Audit.csv")

  # Preserve planned family sizes in BH adjustment, including untestable comparisons.
  test_cache <- out("Spatial_Test_Cache"); dir.create(test_cache, showWarnings = FALSE)
  tables <- diagnostics <- list()
  for (nclass in NCLASS) {
    ans <- vector("list", length(jobs))
    for (i in seq_along(jobs)) {
      job <- jobs[[i]]
      signature <- digest::digest(list(job = job, nclass = nclass, min_spots = MIN_SPOTS,
        revision = "core-rank-v1", SpatialPack = as.character(utils::packageVersion("SpatialPack"))), algo = "xxhash64")
      cache_file <- file.path(test_cache, paste0(signature, ".rds"))
      saved <- if (file.exists(cache_file)) tryCatch(readRDS(cache_file), error = function(e) NULL) else NULL
      if (!is.null(saved) && identical(saved$signature, signature)) test <- saved$test
      else {
        if (job$low) {
          ok <- is.finite(job$x) & is.finite(job$y)
          rho <- if (sum(ok) >= MIN_SPOTS && length(unique(job$x[ok])) > 1L && length(unique(job$y[ok])) > 1L)
            stats::cor(job$x[ok], job$y[ok], method = "spearman") else NA_real_
          test <- list(result = data.frame(N_spots = sum(ok), Rho = rho, ESS = NA_real_, DF = NA_real_, P_spatial = NA_real_,
            Status = "Low_bacterial_information", Message = "Below fixed count/positive-spot criteria"), diagnostic = NULL)
        } else test <- spatial_test(job$x, job$y, job$coords, nclass)
        tmp <- paste0(cache_file, ".writing"); saveRDS(list(signature = signature, test = test), tmp)
        if (file.exists(cache_file) && !file.remove(cache_file)) stop("Cannot replace spatial checkpoint")
        if (!file.rename(tmp, cache_file)) stop("Cannot finalize spatial checkpoint")
      }
      job_meta <- data.frame(Family = job$Family, Section = job$Section, Bacterial_group = job$Bacterial_group, Host_target = job$Host_target, Nclass = nclass)
      ans[[i]] <- cbind(job_meta, test$result)
      if (!is.null(test$diagnostic)) diagnostics[[length(diagnostics) + 1L]] <- cbind(job_meta, test$diagnostic)
    }
    z <- do.call(rbind, ans); z$FDR_spatial <- NA_real_
    for (family in unique(z$Family)) {
      ii <- which(z$Family == family); ok <- ii[z$Status[ii] == "OK" & is.finite(z$P_spatial[ii])]
      if (length(ok)) z$FDR_spatial[ok] <- p.adjust(z$P_spatial[ok], "BH", n = length(ii))
    }
    tables[[as.character(nclass)]] <- z
    tab(z, paste0("Checkpoint_Nclass_", nclass, ".csv"))
  }
  all_tests <- do.call(rbind, tables); primary <- tables[["13"]]
  id <- function(z) paste(z$Family, z$Section, z$Bacterial_group, z$Host_target, sep = "::")
  enrichment <- do.call(rbind, lapply(jobs, function(j) {
    if (is.null(j$overlap)) return(NULL)
    cbind(data.frame(Family = j$Family, Section = j$Section, Bacterial_group = j$Bacterial_group, Host_target = j$Host_target), j$overlap)
  }))
  tab(primary, "06_Primary_Spatial_Adjusted_Correlations.csv")
  tab(all_tests[all_tests$Nclass != 13, , drop = FALSE], "07_Spatial_Sensitivity.csv")
  tab(enrichment, "08_Descriptive_HighSignal_Enrichment.csv")
  tab(if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame(Note = "No estimable spatial diagnostics"), "09_Spatial_Moran_Diagnostics.csv")
  summary <- do.call(rbind, lapply(split(all_tests, list(all_tests$Nclass, all_tests$Family), drop = TRUE), function(z)
    data.frame(Nclass = z$Nclass[1], Family = z$Family[1], Planned = nrow(z), Valid = sum(z$Status == "OK"),
      Significant = sum(z$FDR_spatial < 0.05, na.rm = TRUE), Invalid = nrow(z) - sum(z$Status == "OK"))))
  tab(summary, "10_Test_Summary.csv")
  stability <- primary[c("Family", "Section", "Bacterial_group", "Host_target", "Rho", "FDR_spatial")]
  for (k in c(8, 20)) {
    z <- tables[[as.character(k)]]
    stability[[paste0("FDR_nclass_", k)]] <- z$FDR_spatial[match(id(primary), id(z))]
  }
  stability$Significant_all_three <- is.finite(stability$FDR_spatial) & stability$FDR_spatial < 0.05 &
    is.finite(stability$FDR_nclass_8) & stability$FDR_nclass_8 < 0.05 & is.finite(stability$FDR_nclass_20) & stability$FDR_nclass_20 < 0.05
  tab(stability, "11_Comparison_Sensitivity_Stability.csv")
  consistency <- do.call(rbind, lapply(split(primary, paste(primary$Family, primary$Bacterial_group, primary$Host_target)), function(z)
    data.frame(Family = z$Family[1], Bacterial_group = z$Bacterial_group[1], Host_target = z$Host_target[1],
      Sections_with_rho = sum(is.finite(z$Rho)), Positive_sections = sum(z$Rho > 0, na.rm = TRUE),
      Negative_sections = sum(z$Rho < 0, na.rm = TRUE), Significant_sections = sum(z$FDR_spatial < 0.05, na.rm = TRUE),
      Median_rho = if (any(is.finite(z$Rho))) median(z$Rho, na.rm = TRUE) else NA_real_)))
  tab(consistency, "12_Direction_Consistency_Descriptive.csv")
  pooled <- do.call(rbind, lapply(split(seq_along(jobs), vapply(jobs, function(j)
    paste(j$Family, j$Bacterial_group, j$Host_target, sep = "::"), character(1))), function(ix) {
      js <- jobs[ix]; x <- unlist(lapply(js, `[[`, "x")); y <- unlist(lapply(js, `[[`, "y")); ok <- is.finite(x) & is.finite(y)
      data.frame(Family = js[[1]]$Family, Bacterial_group = js[[1]]$Bacterial_group, Host_target = js[[1]]$Host_target, N_spots = sum(ok),
        Rho = if (sum(ok) >= MIN_SPOTS && length(unique(x[ok])) > 1L && length(unique(y[ok])) > 1L)
          stats::cor(x[ok], y[ok], method = "spearman") else NA_real_)
    }))
  tab(pooled, "13_Pooled_Correlations_Descriptive_Only.csv")
  tz <- primary[grepl("^Taxon_", primary$Family), , drop = FALSE]
  ranking <- do.call(rbind, lapply(split(tz, tz$Bacterial_group), function(z)
    data.frame(Bacterial_group = z$Bacterial_group[1], Valid_tests = sum(z$Status == "OK"), Significant_tests = sum(z$FDR_spatial < 0.05, na.rm = TRUE),
      Maximum_absolute_rho = if (any(is.finite(z$Rho))) max(abs(z$Rho), na.rm = TRUE) else NA_real_)))
  ranking <- ranking[order(-ranking$Significant_tests, -ranking$Maximum_absolute_rho, na.last = TRUE), ]
  tab(ranking, "14_Exploratory_Bacterial_Ranking.csv")
  top <- head(ranking$Bacterial_group[ranking$Valid_tests > 0 & !ranking$Bacterial_group %in% MAIN_EXCLUDE], TOP_N)
  tab(data.frame(Bacterial_group = top), "15_Main_Figure_Selected_Groups.csv")

  # Main figures are selected exploratorily; stars use primary spatial FDR only.
  suppressPackageStartupMessages(library(ggplot2))
  theme_p <- theme_classic(base_family = "Arial", base_size = 9) + theme(
    plot.background = element_rect(fill = "white", colour = NA), axis.text = element_text(size = 8, colour = "black"),
    strip.text = element_text(size = 9, face = "bold"), plot.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 8))
  save_png <- function(p, name, height = 6.5) ggsave(out(name), p, width = 7.09, height = height, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
  fmt <- function(z) {
    z$Section <- factor(z$Section, levels = SECTIONS, labels = c("03 Normal", "03 Tumor", "08 Normal", "08 Tumor"))
    z$Star <- ifelse(!is.na(z$FDR_spatial) & z$FDR_spatial < 0.05, "*", "")
    z
  }
  total_plot <- function(family, title) {
    z <- fmt(primary[primary$Family == family, , drop = FALSE])
    ggplot(z, aes(Section, Host_target, fill = Rho)) + geom_tile(colour = "white") + geom_text(aes(label = Star), size = 3) +
      scale_fill_gradient2(low = "#7570B3", mid = "white", high = "#D81B60", limits = c(-1, 1), na.value = "#BDBDBD", name = "Spearman rho") +
      labs(x = NULL, y = NULL, title = title) + theme_p
  }
  totals <- cowplot::plot_grid(total_plot("Total_Feature", "Total bacterial signal vs host features"),
    total_plot("Total_RCTD", "Total bacterial signal vs RCTD"), ncol = 1, rel_heights = c(3, 1), labels = c("A", "B"), label_size = 15)
  save_png(totals, "05_Total_Bacteria_Host_Associations.png", 8)
  taxon_plot <- function(family, use, bubble = FALSE) {
    z <- tz[tz$Family == family & tz$Bacterial_group %in% use, , drop = FALSE]
    if (!nrow(z)) return(ggplot() + theme_void() + labs(title = "No eligible bacterial groups"))
    z$Enrichment <- enrichment$Enrichment[match(id(z), id(enrichment))]
    z <- fmt(z); z$Bacterial_group <- factor(z$Bacterial_group, levels = rev(use))
    if (bubble) {
      z$Display_enrichment <- pmin(3, pmax(0, z$Enrichment))
      p <- ggplot(z, aes(Section, Bacterial_group)) + facet_wrap(~Host_target, nrow = 1) +
        geom_point(aes(fill = Rho, size = Display_enrichment), shape = 21, colour = "black", stroke = 0.25) +
        scale_size_continuous(range = c(1, 5), limits = c(0, 3), name = "Enrichment\n(capped at 3)")
    } else p <- ggplot(z, aes(Section, Bacterial_group)) + facet_wrap(~Host_target, nrow = 1) + geom_tile(aes(fill = Rho), colour = "white", linewidth = 0.2)
    p + geom_text(aes(label = Star), size = 3, vjust = -0.3) +
      scale_fill_gradient2(low = "#7570B3", mid = "white", high = "#D81B60", limits = c(-1, 1), na.value = "#BDBDBD", name = "Spearman rho") +
      labs(x = NULL, y = NULL) + theme_p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  if (length(top)) save_png(cowplot::plot_grid(taxon_plot("Taxon_SameSpot", top), taxon_plot("Taxon_Neighborhood", top, TRUE),
    ncol = 1, labels = c("A", "B"), label_size = 15), "05_Main_Taxon_SameSpot_and_Neighborhood.png", 9)
  save_png(cowplot::plot_grid(taxon_plot("Taxon_SameSpot", ranking$Bacterial_group), taxon_plot("Taxon_Neighborhood", ranking$Bacterial_group),
    ncol = 1, labels = c("A", "B"), label_size = 15), "05_Supplementary_All_Bacterial_Groups.png", max(9, 0.38 * length(groups) + 3))
  for (family in c("Total_Feature", "Total_RCTD")) {
    q <- primary[primary$Family == family & primary$Status == "OK", , drop = FALSE]
    q <- head(q[order(q$FDR_spatial, -abs(q$Rho)), ], 4)
    plots <- lapply(seq_len(nrow(q)), function(i) {
      ix <- which(vapply(jobs, function(j) j$Family == family && j$Section == q$Section[i] && j$Host_target == q$Host_target[i], logical(1)))
      j <- jobs[[ix]]; z <- data.frame(Bacterial_signal = log1p(j$x), Host_value = j$y)
      ggplot(z, aes(Bacterial_signal, Host_value)) + geom_point(size = 0.4, alpha = 0.65, colour = "#0072B2") +
        labs(title = paste(q$Section[i], q$Host_target[i]), x = "log(1 + bacterial probe counts)", y = "Host value",
          subtitle = sprintf("rho = %.3f; spatial FDR = %.3g", q$Rho[i], q$FDR_spatial[i])) + theme_p
    })
    if (length(plots)) save_png(cowplot::plot_grid(plotlist = plots, ncol = 2), paste0("05_", family, "_Selected_Scatterplots.png"), 6)
  }
  writeLines(capture.output(utils::sessionInfo()), out("17_sessionInfo.txt"))
  generated <- list.files(OUTPUT_DIR, pattern = "\\.(csv|png|txt)$", full.names = FALSE)
  tab(data.frame(File = generated), "18_Output_Manifest.csv")
  completed_files <- c("02_Bacterial_Extraction_Check.csv", "03_Bacterial_Groups_All_QC_Spots.csv",
    "06_Primary_Spatial_Adjusted_Correlations.csv", "07_Spatial_Sensitivity.csv", "10_Test_Summary.csv", "17_sessionInfo.txt")
  saveRDS(list(revision = "core-v1", files = completed_files,
    md5 = unname(tools::md5sum(out(completed_files))), completed = Sys.time()), own_complete)
})
invisible(gc())
