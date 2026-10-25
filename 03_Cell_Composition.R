# GSE132465: paired cell-type/subtype composition in 10 matched patients.
# Percentages use all retained cells per sample, including for subtypes.
# GEO subtypes are published annotations, not newly inferred subclusters.
local({
  input_dir <- Sys.getenv("CRC_ANNOTATION_OUTPUT", "E:/CRC/Results/Results/02_Clustering_Markers_Annotation")
  output_dir <- Sys.getenv("CRC_COMPOSITION_OUTPUT", "E:/CRC/Results/Results/03_Cell-Type and Subtype Composition Analysis")
  packages <- c("SeuratObject", "ggplot2", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(ggplot2))
  patients <- sprintf("SMC%02d", 1:10)
  samples <- as.vector(t(outer(patients, c("-N", "-T"), paste0)))
  broad_labels <- c("B cells", "Epithelial cells", "Mast cells", "Myeloids", "Stromal cells", "T cells")
  broad_order <- c("B", "Epithelial", "Mast", "Myeloid", "Stromal", "T")
  input_file <- file.path(input_dir, "05_annotated.rds")
  if (!file.exists(input_file)) stop("Missing input: ", input_file)
  obj <- readRDS(input_file)
  if (!inherits(obj, "Seurat")) stop("Expected a Seurat object")
  meta <- obj[[]]
  if (!identical(rownames(meta), colnames(obj))) stop("Metadata order mismatch")
  rm(obj)
  invisible(gc())
  required <- c("sample_id", "patient_id", "group", "analysis_cell_type", "GEO_CellSubtype")
  if (!all(required %in% names(meta))) stop("Missing metadata: ", paste(setdiff(required, names(meta)), collapse = ", "))
  meta[required] <- lapply(meta[required], as.character)
  core <- setdiff(required, "GEO_CellSubtype")
  if (anyNA(meta[core]) || anyDuplicated(rownames(meta)) ||
      any(vapply(meta[core], function(v) any(!nzchar(trimws(v))), logical(1)))) stop("Invalid metadata")
  if (!setequal(meta$sample_id, samples) || !setequal(meta$patient_id, patients) ||
      !all(meta$group %in% c("Normal", "Tumor")) ||
      any(meta$sample_id != paste0(meta$patient_id, ifelse(meta$group == "Normal", "-N", "-T")))) stop("Invalid patient pairing")
  if (!setequal(meta$analysis_cell_type, broad_labels)) stop("Unexpected or missing broad cell types")
  meta$broad_type <- broad_order[match(meta$analysis_cell_type, broad_labels)]
  meta$subtype <- trimws(meta$GEO_CellSubtype)
  invalid <- is.na(meta$subtype) | tolower(meta$subtype) %in% c("", "unknown", "unassigned", "na", "n/a", "none")
  meta$subtype[invalid] <- NA_character_
  totals <- table(factor(meta$sample_id, levels = samples))
  if (any(totals == 0)) stop("Empty sample")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory")
  out <- function(name) file.path(output_dir, name)
  write_tab <- function(z, name) utils::write.csv(z, out(name), row.names = FALSE, na = "")
  save_plot <- function(p, name, width, height) {
    ggsave(out(paste0(name, ".png")), p, width = width, height = height,
           units = "in", device = "png", dpi = 600, bg = "white", limitsize = FALSE)
  }

  make_rows <- function(cell_samples, feature, parent, level) {
    n <- as.integer(table(factor(cell_samples, levels = samples)))
    data.frame(sample = samples, feature = feature, n_cells = n,
               percentage = 100 * n / as.numeric(totals[samples]),
               patient = sub("-[NT]$", "", samples),
               group = ifelse(grepl("-N$", samples), "Normal", "Tumor"),
               parent = parent, level = level)
  }
  broad_df <- do.call(rbind, lapply(broad_order, function(tp)
    make_rows(meta$sample_id[meta$broad_type == tp], tp, "Broad", "Broad")))
  hierarchy <- unique(meta[!is.na(meta$subtype), c("broad_type", "subtype")])
  hierarchy <- hierarchy[order(match(hierarchy$broad_type, broad_order), hierarchy$subtype), , drop = FALSE]
  rownames(hierarchy) <- NULL
  if (!nrow(hierarchy)) stop("No valid published subtypes")
  subtype_df <- do.call(rbind, lapply(seq_len(nrow(hierarchy)), function(i) {
    keep <- meta$broad_type == hierarchy$broad_type[i] & !is.na(meta$subtype) &
      meta$subtype == hierarchy$subtype[i]
    make_rows(meta$sample_id[keep], hierarchy$subtype[i], hierarchy$broad_type[i], "Subtype")
  }))
  composition <- rbind(broad_df, subtype_df)
  composition$group <- factor(composition$group, levels = c("Normal", "Tumor"))
  composition$patient <- factor(composition$patient, levels = patients)
  keys <- unique(composition[c("level", "parent", "feature")])
  # Preserve the original planned family: six broad types and 35 subtypes.
  if (nrow(keys) != 41L) stop("Expected 41 populations; review annotations before changing the test family")
  population <- function(i) composition[composition$level == keys$level[i] &
    composition$parent == keys$parent[i] & composition$feature == keys$feature[i], , drop = FALSE]
  pair_values <- function(z) {
    nn <- z[z$group == "Normal", ]; tt <- z[z$group == "Tumor", ]
    if (nrow(nn) != 10L || nrow(tt) != 10L || anyDuplicated(nn$patient) ||
        anyDuplicated(tt$patient) || !setequal(as.character(nn$patient), patients) ||
        !setequal(as.character(tt$patient), patients)) stop("Incomplete population pairing")
    nn <- nn$percentage[match(patients, as.character(nn$patient))]
    tt <- tt$percentage[match(patients, as.character(tt$patient))]
    list(normal = nn, tumor = tt, delta = tt - nn)
  }
  paired_values <- lapply(seq_len(nrow(keys)), function(i) pair_values(population(i)))
  tests <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    v <- paired_values[[i]]
    pv <- if (all(v$delta == 0)) 1 else suppressWarnings(stats::wilcox.test(
      v$tumor, v$normal, paired = TRUE, exact = FALSE, correct = TRUE)$p.value)
    data.frame(keys[i, ], n_pairs = length(v$delta), median_normal = median(v$normal),
               median_tumor = median(v$tumor), median_paired_change = median(v$delta), p_value = pv)
  }))
  tests$FDR <- p.adjust(tests$p_value, "BH")

  # Conditional exact signed-rank sensitivity, allowing tied ranks and zeros.
  exact_signed_rank <- function(delta) {
    if (any(!is.finite(delta))) stop("Nonfinite paired differences")
    delta <- delta[delta != 0]
    if (!length(delta)) return(1)
    if (length(delta) > 20L) stop("Too many pairs for sign enumeration")
    ranks <- rank(abs(delta), ties.method = "average")
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(delta)), KEEP.OUT.ATTRS = FALSE))
    mean(abs(drop(signs %*% ranks)) >= abs(sum(sign(delta) * ranks)) - 1e-10)
  }
  sensitivity <- tests
  sensitivity$n_nonzero <- vapply(paired_values, function(v) sum(v$delta != 0), integer(1))
  sensitivity$n_zero <- vapply(paired_values, function(v) sum(v$delta == 0), integer(1))
  sensitivity$tied_nonzero_magnitudes <- vapply(paired_values, function(v)
    anyDuplicated(abs(v$delta[v$delta != 0])) > 0L, logical(1))
  sensitivity$p_exact <- vapply(paired_values, function(v) exact_signed_rank(v$delta), numeric(1))
  sensitivity$FDR_exact <- p.adjust(sensitivity$p_exact, "BH")
  sensitivity$significant_primary <- sensitivity$FDR < 0.05
  sensitivity$significant_exact <- sensitivity$FDR_exact < 0.05
  sensitivity$status_changed <- sensitivity$significant_primary != sensitivity$significant_exact
  fmt <- function(x) ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x)))
  tests$label <- paste0("p = ", fmt(tests$p_value), "\nFDR = ", fmt(tests$FDR))
  write_tab(composition, "04_CellType_Subcluster_Sample_Percentages.csv")
  write_tab(tests, "04_CellType_Subcluster_Paired_Wilcoxon_FDR.csv")
  write_tab(hierarchy, "04_CellType_Subcluster_Hierarchy.csv")
  write_tab(sensitivity, "05_Composition_Exact_Sensitivity.csv")

  # Compatibility outputs: broad-only BH is distinct from joint BH over 41 tests.
  legacy_data <- broad_df
  legacy_data$cell_type <- broad_labels[match(legacy_data$feature, broad_order)]
  legacy_data <- legacy_data[c("sample", "cell_type", "n_cells", "percentage", "patient", "group")]
  legacy_tests <- tests[tests$level == "Broad", c("feature", "n_pairs", "median_paired_change", "p_value")]
  names(legacy_tests)[1] <- "cell_type"
  legacy_tests$cell_type <- broad_labels[match(legacy_tests$cell_type, broad_order)]
  legacy_tests$p_adj_BH <- p.adjust(legacy_tests$p_value, "BH")
  write_tab(legacy_data, "03_Broad_Sample_percentages.csv")
  write_tab(legacy_tests, "03_Broad_Paired_tests.csv")
  lp <- ggplot(legacy_data, aes(factor(group, levels = c("Normal", "Tumor")), percentage)) +
    geom_line(aes(group = patient), colour = "#BDBDBD", linewidth = 0.3) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.45, alpha = 0.7) + geom_point(size = 1.5) +
    facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c(Normal = "#BDBDBD", Tumor = "#208C91")) +
    labs(x = NULL, y = "Cells in sample (%)") + theme_classic(base_family = "Arial") + theme(legend.position = "none")
  save_plot(lp, "03_Broad_Composition", 9, 6)

  make_plot <- function(d, t, order, title, columns) {
    d$feature <- factor(d$feature, levels = order)
    t$feature <- factor(t$feature, levels = order)
    max_y <- aggregate(percentage ~ feature, d, max)
    names(max_y)[2] <- "max_y"
    t <- merge(t, max_y, by = "feature", all.x = TRUE, sort = FALSE)
    t$label_y <- ifelse(t$max_y > 0, 1.18 * t$max_y, 0.25)
    ggplot(d, aes(group, percentage)) +
      geom_line(aes(group = patient), colour = "#BDBDBD", linewidth = 0.35, alpha = 0.8) +
      geom_boxplot(aes(fill = group), width = 0.48, linewidth = 0.45, outlier.shape = NA, alpha = 0.75) +
      geom_point(aes(fill = group), shape = 21, size = 1.8, stroke = 0.3, colour = "black") +
      geom_text(data = t, aes(x = 1.5, y = label_y, label = label), inherit.aes = FALSE,
                size = 2.5, lineheight = 0.9, vjust = 0) +
      facet_wrap(~feature, scales = "free_y", ncol = columns) +
      scale_fill_manual(values = c(Normal = "#BDBDBD", Tumor = "#208C91")) +
      scale_y_continuous(expand = expansion(mult = c(0.03, 0.35))) +
      labs(x = NULL, y = "Cells in sample (%)", title = title) +
      theme_classic(base_family = "Arial", base_size = 9) +
      theme(legend.position = "none", axis.text = element_text(colour = "black", size = 8),
            strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9),
            plot.title = element_text(face = "bold", size = 10),
            panel.spacing.x = grid::unit(0.6, "lines"), panel.spacing.y = grid::unit(0.8, "lines"))
  }
  broad_plot <- make_plot(composition[composition$level == "Broad", ], tests[tests$level == "Broad", ],
                          broad_order, "Broad cell types", 6)
  save_plot(broad_plot, "04_Broad_CellTypes_Composition", 14, 3.4)
  plots <- list(broad_plot)
  facet_rows <- 0L
  for (parent in broad_order) {
    order <- unique(hierarchy$subtype[hierarchy$broad_type == parent])
    if (!length(order)) next
    p <- make_plot(composition[composition$level == "Subtype" & composition$parent == parent, ],
                   tests[tests$level == "Subtype" & tests$parent == parent, ],
                   order, paste(parent, "published subtypes"), min(6L, length(order)))
    nr <- ceiling(length(order) / 6)
    facet_rows <- facet_rows + nr
    plots[[length(plots) + 1L]] <- p
    save_plot(p, paste0("04_", gsub("[^A-Za-z0-9]", "_", parent), "_Subclusters_Composition"), 14, max(3.4, 3 * nr))
  }
  combined <- patchwork::wrap_plots(plots, ncol = 1) + patchwork::plot_annotation(
    title = "Cell-type and published subtype composition in matched CRC samples",
    subtitle = "Paired Wilcoxon signed-rank tests; joint BH correction across 41 populations",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 13),
                  plot.subtitle = element_text(family = "Arial", size = 9)))
  save_plot(combined, "04_CellType_Subcluster_Composition", 14, 3.5 + 3 * facet_rows)
  counts <- as.data.frame(table(meta$broad_type))
  names(counts) <- c("Broad_Cell_Type", "Cells")
  write_tab(counts, "04_Broad_CellType_Total_Cell_Counts.csv")
  counts <- as.data.frame(table(meta$broad_type, meta$subtype, useNA = "no"))
  names(counts) <- c("Broad_Cell_Type", "Subtype", "Cells")
  write_tab(counts[counts$Cells > 0, ], "04_Subcluster_Total_Cell_Counts.csv")
  writeLines(capture.output(utils::sessionInfo()), out("03_sessionInfo.txt"))
  cat("\nComposition analysis complete\nCells:", nrow(meta), "\nPatients:", length(patients),
      "\nSamples:", length(samples), "\nPublished subtypes:", nrow(hierarchy),
      "\nPrimary significant:", sum(sensitivity$significant_primary),
      "\nExact significant:", sum(sensitivity$significant_exact),
      "\nSensitivity status changes:", sum(sensitivity$status_changed), "\nOutput:", output_dir, "\n")
  print(tests, row.names = FALSE)
  print(sensitivity, row.names = FALSE)
})
invisible(gc())
