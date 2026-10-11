# Cell Communication
# CellChat models normal and tumor separately using unprojected signaling data.
# Condition-level communication estimates are exploratory.

rm(list = ls())
invisible(gc())

set.seed(42)

# Paths

input_file <- "E:/CRC/Results/GSE132465_scRNA/Core_From_Raw/05_annotated.rds"

output_root <- "E:/CRC/Results/Results/05- CellChat"

dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_root)) {
  stop("Cannot create output directory: ", output_root)
}

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

out <- function(x) {
  file.path(output_root, x)
}

# Packages

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "CellChat",
  "patchwork",
  "ggplot2",
  "Matrix"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Matrix))

PRESTO_AVAILABLE <- requireNamespace(
  "presto",
  quietly = TRUE
)

if (PRESTO_AVAILABLE) {
  message("presto detected -> fast Wilcoxon mode will be used.")
} else {
  message(
    "presto is NOT installed -> standard Wilcoxon mode will be used ",
    "(do.fast = FALSE)."
  )
}

# Font / style

FONT <- "Arial"

if (.Platform$OS.type == "windows") {
  try(
    windowsFonts(
      Arial = windowsFont("Arial")
    ),
    silent = TRUE
  )
}

WHITE <- "#FFFFFF"
BLACK <- "#000000"

CELLTYPE_COLORS <- c(
  "B cells"          = "#E69F00",
  "Epithelial cells" = "#D81B60",
  "Mast cells"       = "#7570B3",
  "Myeloids"         = "#009E73",
  "Stromal cells"    = "#56B4E9",
  "T cells"          = "#0072B2"
)

theme_pub <- function() {
  ggplot2::theme_classic(
    base_family = FONT,
    base_size = 9
  ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = WHITE,
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = WHITE,
        colour = NA
      ),
      text = ggplot2::element_text(
        family = FONT,
        colour = BLACK
      ),
      plot.title = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        size = 9
      ),
      axis.text = ggplot2::element_text(
        size = 8
      ),
      legend.title = ggplot2::element_text(
        size = 8.5,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 8
      )
    )
}

# Load annotated seurat object

obj <- readRDS(input_file)

if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.")
}

required_meta <- c(
  "sample_id",
  "patient_id",
  "group",
  "analysis_cell_type"
)

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata column(s): ",
    paste(missing_meta, collapse = ", ")
  )
}

meta_all <- obj@meta.data

if (!identical(
  rownames(meta_all),
  colnames(obj)
)) {
  stop("Metadata row names do not match Seurat cell names.")
}

meta_all$sample_id <- as.character(meta_all$sample_id)
meta_all$patient_id <- as.character(meta_all$patient_id)
meta_all$group <- as.character(meta_all$group)
meta_all$analysis_cell_type <- as.character(meta_all$analysis_cell_type)

if (!all(meta_all$group %in% c("Normal", "Tumor"))) {
  stop("group must contain only Normal and Tumor.")
}

if (anyNA(meta_all[, required_meta, drop = FALSE])) {
  stop("Missing values detected in required metadata.")
}

# Validate expected cell types

expected_cell_types <- c(
  "B cells",
  "Epithelial cells",
  "Mast cells",
  "Myeloids",
  "Stromal cells",
  "T cells"
)

observed_cell_types <- sort(
  unique(meta_all$analysis_cell_type)
)

unexpected_types <- setdiff(
  observed_cell_types,
  expected_cell_types
)

if (length(unexpected_types) > 0) {
  message(
    "Additional cell type(s) detected: ",
    paste(unexpected_types, collapse = ", ")
  )
}

if (length(intersect(
  expected_cell_types,
  observed_cell_types
)) < 2) {
  stop("Too few expected broad cell types are present.")
}

# Get normalized rna expression

DefaultAssay(obj) <- "RNA"

get_normalized_expression <- function(seurat_obj) {

  assay_obj <- seurat_obj[["RNA"]]

  if (inherits(assay_obj, "Assay5")) {

    data_layers <- SeuratObject::Layers(
      assay_obj,
      search = "^data$|^data\\."
    )

    if (length(data_layers) == 0) {

      message(
        "Normalized RNA data layer not found. Running NormalizeData..."
      )

      seurat_obj <- NormalizeData(
        seurat_obj,
        assay = "RNA",
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
      )

      assay_obj <- seurat_obj[["RNA"]]

      data_layers <- SeuratObject::Layers(
        assay_obj,
        search = "^data$|^data\\."
      )
    }

    if (length(data_layers) > 1) {

      message(
        "Multiple RNA data layers detected -> joining them into one 'data' layer."
      )

      seurat_obj <- SeuratObject::JoinLayers(
        seurat_obj,
        assay = "RNA",
        layers = data_layers,
        new = "data"
      )
    }

    expr <- SeuratObject::GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    )

  } else {

    expr <- SeuratObject::GetAssayData(
      seurat_obj,
      assay = "RNA",
      slot = "data"
    )

    if (nrow(expr) == 0 || ncol(expr) == 0) {

      message(
        "Normalized RNA data slot not found. Running NormalizeData..."
      )

      seurat_obj <- NormalizeData(
        seurat_obj,
        assay = "RNA",
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
      )

      expr <- SeuratObject::GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  }

  if (!inherits(expr, "dgCMatrix")) {
    expr <- as(expr, "dgCMatrix")
  }

  if (!identical(
    colnames(expr),
    colnames(seurat_obj)
  )) {
    stop("Expression matrix columns do not match Seurat cells.")
  }

  list(
    obj = seurat_obj,
    expr = expr
  )
}

norm_result <- get_normalized_expression(obj)

obj <- norm_result$obj
expr_all <- norm_result$expr

rm(norm_result)
invisible(gc())

# Cellchat database

CellChatDB.use <- CellChat::CellChatDB.human

# Safe overexpressed-gene wrapper

identify_overexpressed_genes_safe <- function(cellchat_obj) {

  fn <- CellChat::identifyOverExpressedGenes
  fn_formals <- names(formals(fn))

  if ("do.fast" %in% fn_formals) {

    message(
      "identifyOverExpressedGenes: do.fast = ",
      PRESTO_AVAILABLE
    )

    return(
      fn(
        cellchat_obj,
        do.fast = PRESTO_AVAILABLE
      )
    )
  }

  message(
    "CellChat version has no do.fast argument; using package defaults."
  )

  fn(cellchat_obj)
}

# Helper — run cellchat for one condition

run_cellchat_condition <- function(
    condition_name,
    expression_matrix,
    meta_data
) {

  message(
    "\n============================================================"
  )

  message(
    "\nRunning CellChat for: ",
    condition_name
  )

  message(
    "\n============================================================"
  )

  cells_keep <- rownames(meta_data)[
    meta_data$group == condition_name
  ]

  if (length(cells_keep) == 0) {
    stop(
      "No cells found for condition: ",
      condition_name
    )
  }

  expr <- expression_matrix[
    ,
    cells_keep,
    drop = FALSE
  ]

  meta <- meta_data[
    cells_keep,
    ,
    drop = FALSE
  ]

  if (!identical(
    colnames(expr),
    rownames(meta)
  )) {
    stop(
      "Expression / metadata mismatch in ",
      condition_name
    )
  }

  meta$samples <- factor(
    meta$sample_id
  )

  meta$analysis_cell_type <- factor(
    meta$analysis_cell_type,
    levels = expected_cell_types
  )

  meta$analysis_cell_type <- droplevels(
    meta$analysis_cell_type
  )

  if (
    anyNA(meta$samples) ||
    anyNA(meta$analysis_cell_type)
  ) {
    stop(
      "Invalid sample or cell-type factor in condition: ",
      condition_name
    )
  }

  cell_counts <- table(
    meta$analysis_cell_type
  )

  message(
    "\nCell counts:\n",
    paste(
      names(cell_counts),
      as.integer(cell_counts),
      sep = ": ",
      collapse = "\n"
    )
  )

  sample_counts <- table(
    meta$samples
  )

  message(
    "\nSamples:\n",
    paste(
      names(sample_counts),
      as.integer(sample_counts),
      sep = ": ",
      collapse = "\n"
    )
  )

  cellchat <- CellChat::createCellChat(
    object = expr,
    meta = meta,
    group.by = "analysis_cell_type"
  )

  cellchat@DB <- CellChatDB.use

  cellchat <- CellChat::subsetData(
    cellchat
  )

  if (
    nrow(cellchat@data.signaling) == 0 ||
    ncol(cellchat@data.signaling) == 0
  ) {
    stop(
      "CellChat data.signaling is empty for ",
      condition_name
    )
  }

  future_available <- requireNamespace(
    "future",
    quietly = TRUE
  )

  if (future_available) {

    old_plan <- future::plan()

    on.exit(
      future::plan(old_plan),
      add = TRUE
    )

    future::plan(
      "sequential"
    )
  }

  cellchat <- identify_overexpressed_genes_safe(
    cellchat
  )

  cellchat <- CellChat::identifyOverExpressedInteractions(
    cellchat
  )

  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = "triMean",
    raw.use = TRUE,
    population.size = TRUE,
    seed.use = 42L
  )

  cellchat <- CellChat::filterCommunication(
    cellchat,
    min.cells = 10
  )

  if (
    is.null(cellchat@net$prob) ||
    length(cellchat@net$prob) == 0
  ) {
    stop(
      "No communication probabilities were generated for ",
      condition_name
    )
  }

  cellchat <- CellChat::computeCommunProbPathway(
    cellchat
  )

  cellchat <- CellChat::aggregateNet(
    cellchat
  )

  cellchat <- CellChat::netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )

  cellchat
}

# Run normal

cellchat_normal <- run_cellchat_condition(
  condition_name = "Normal",
  expression_matrix = expr_all,
  meta_data = meta_all
)

saveRDS(
  cellchat_normal,
  out("04_CellChat_Normal.rds")
)

# Run tumor

cellchat_tumor <- run_cellchat_condition(
  condition_name = "Tumor",
  expression_matrix = expr_all,
  meta_data = meta_all
)

saveRDS(
  cellchat_tumor,
  out("04_CellChat_Tumor.rds")
)

# Export significant interactions

normal_interactions <- CellChat::subsetCommunication(
  cellchat_normal
)

tumor_interactions <- CellChat::subsetCommunication(
  cellchat_tumor
)

write.csv(
  normal_interactions,
  out("04_Normal_Interactions.csv"),
  row.names = FALSE
)

write.csv(
  tumor_interactions,
  out("04_Tumor_Interactions.csv"),
  row.names = FALSE
)

# Merge normal + tumor

cellchat_list <- list(
  Normal = cellchat_normal,
  Tumor = cellchat_tumor
)

cellchat_merged <- CellChat::mergeCellChat(
  cellchat_list,
  add.names = names(cellchat_list)
)

saveRDS(
  cellchat_merged,
  out("04_CellChat_Merged_Normal_vs_Tumor.rds")
)

# Comparison — number of interactions

p_compare_count <- CellChat::compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
) +
  ggplot2::labs(
    title = "Number of interactions"
  ) +
  theme_pub()

ggplot2::ggsave(
  filename = out("04_Compare_Number_of_Interactions.png"),
  plot = p_compare_count,
  width = 4.8,
  height = 4.0,
  units = "in",
  dpi = 600,
  bg = WHITE
)

# Comparison — interaction strength

p_compare_weight <- CellChat::compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
) +
  ggplot2::labs(
    title = "Interaction strength"
  ) +
  theme_pub()

ggplot2::ggsave(
  filename = out("04_Compare_Interaction_Strength.png"),
  plot = p_compare_weight,
  width = 4.8,
  height = 4.0,
  units = "in",
  dpi = 600,
  bg = WHITE
)

# Circle network helper

save_circle_network <- function(
    cellchat_obj,
    matrix_obj,
    title_text,
    file_name,
    cell_colors = CELLTYPE_COLORS
) {

  group_names <- levels(
    cellchat_obj@idents
  )

  group_size <- as.numeric(
    table(
      factor(
        cellchat_obj@idents,
        levels = group_names
      )
    )
  )

  colors_use <- cell_colors[
    group_names
  ]

  if (anyNA(colors_use)) {

    missing_names <- group_names[
      is.na(colors_use)
    ]

    fallback <- grDevices::hcl.colors(
      length(missing_names),
      "Dark 3"
    )

    colors_use[
      is.na(colors_use)
    ] <- fallback
  }

  grDevices::png(
    filename = out(file_name),
    width = 6.5,
    height = 6.5,
    units = "in",
    res = 600,
    bg = WHITE
  )

  CellChat::netVisual_circle(
    matrix_obj,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = colors_use,
    title.name = title_text
  )

  grDevices::dev.off()
}

# Normal network — number of interactions

save_circle_network(
  cellchat_obj = cellchat_normal,
  matrix_obj = cellchat_normal@net$count,
  title_text = "Normal — Number of interactions",
  file_name = "04_Normal_Network_Count.png"
)

# Tumor network — number of interactions

save_circle_network(
  cellchat_obj = cellchat_tumor,
  matrix_obj = cellchat_tumor@net$count,
  title_text = "Tumor — Number of interactions",
  file_name = "04_Tumor_Network_Count.png"
)

# Normal network — interaction strength

save_circle_network(
  cellchat_obj = cellchat_normal,
  matrix_obj = cellchat_normal@net$weight,
  title_text = "Normal — Interaction weights/strength",
  file_name = "04_Normal_Network_Strength.png"
)

# Tumor network — interaction strength

save_circle_network(
  cellchat_obj = cellchat_tumor,
  matrix_obj = cellchat_tumor@net$weight,
  title_text = "Tumor — Interaction weights/strength",
  file_name = "04_Tumor_Network_Strength.png"
)

# Signaling role scatter — normal

p_role_normal <- CellChat::netAnalysis_signalingRole_scatter(
  cellchat_normal
) +
  ggplot2::labs(
    title = "Normal"
  ) +
  theme_pub()

ggplot2::ggsave(
  filename = out("04_Normal_Signaling_Role.png"),
  plot = p_role_normal,
  width = 6.2,
  height = 5.0,
  units = "in",
  dpi = 600,
  bg = WHITE
)

# Signaling role scatter — tumor

p_role_tumor <- CellChat::netAnalysis_signalingRole_scatter(
  cellchat_tumor
) +
  ggplot2::labs(
    title = "Tumor"
  ) +
  theme_pub()

ggplot2::ggsave(
  filename = out("04_Tumor_Signaling_Role.png"),
  plot = p_role_tumor,
  width = 6.2,
  height = 5.0,
  units = "in",
  dpi = 600,
  bg = WHITE
)

# Differential network — count

grDevices::png(
  filename = out("04_Differential_Network_Count_Tumor_vs_Normal.png"),
  width = 6.5,
  height = 6.5,
  units = "in",
  res = 600,
  bg = WHITE
)

CellChat::netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "count"
)

grDevices::dev.off()

# Differential network — weight

grDevices::png(
  filename = out("04_Differential_Network_Strength_Tumor_vs_Normal.png"),
  width = 6.5,
  height = 6.5,
  units = "in",
  res = 600,
  bg = WHITE
)

CellChat::netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "weight"
)

grDevices::dev.off()

# Pathway information

normal_pathways <- cellchat_normal@netP$pathways
tumor_pathways  <- cellchat_tumor@netP$pathways

pathway_summary <- data.frame(
  pathway = sort(
    unique(
      c(
        normal_pathways,
        tumor_pathways
      )
    )
  ),
  stringsAsFactors = FALSE
)

pathway_summary$Normal <- pathway_summary$pathway %in% normal_pathways
pathway_summary$Tumor  <- pathway_summary$pathway %in% tumor_pathways

write.csv(
  pathway_summary,
  out("04_Pathway_Presence_Normal_vs_Tumor.csv"),
  row.names = FALSE
)

# Save session info

session_lines <- capture.output(
  sessionInfo()
)

writeLines(
  session_lines,
  out("04_CellChat_sessionInfo.txt")
)

cat(
  "\n============================================================\n"
)

cat(
  "04 - CellChat COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nInput:\n",
  input_file,
  "\n",
  sep = ""
)

cat(
  "\nOutput folder:\n",
  output_root,
  "\n",
  sep = ""
)

cat(
  "\npresto available: ",
  PRESTO_AVAILABLE,
  "\n",
  sep = ""
)

cat(
  "Wilcoxon mode: ",
  ifelse(
    PRESTO_AVAILABLE,
    "FAST (presto)",
    "STANDARD (do.fast = FALSE)"
  ),
  "\n",
  sep = ""
)

cat(
  "\ncomputeCommunProb raw.use: TRUE\n"
)

cat(
  "\nNormal cells: ",
  sum(meta_all$group == "Normal"),
  "\n",
  sep = ""
)

cat(
  "Tumor cells: ",
  sum(meta_all$group == "Tumor"),
  "\n",
  sep = ""
)

cat(
  "\nNormal significant LR rows: ",
  nrow(normal_interactions),
  "\n",
  sep = ""
)

cat(
  "Tumor significant LR rows: ",
  nrow(tumor_interactions),
  "\n",
  sep = ""
)

cat(
  "\nSaved merged comparison object:\n",
  out("04_CellChat_Merged_Normal_vs_Tumor.rds"),
  "\n",
  sep = ""
)

cat(
  "\nNext step:\n"
)

cat(
  "Build the final multi-panel reference-style CellChat figure.\n"
)

cat(
  "============================================================\n"
)

extra_packages <- c(
  "cowplot",
  "grid",
  "png",
  "ComplexHeatmap"
)

missing_extra <- extra_packages[
  !vapply(
    extra_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_extra) > 0) {
  stop(
    "Install missing package(s) for final figure: ",
    paste(missing_extra, collapse = ", ")
  )
}

suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(png))

REF_CELL_COLORS <- c(
  "B cells"          = "#4DAF4A",
  "Epithelial cells" = "#E41A1C",
  "Mast cells"       = "#F781BF",
  "Myeloids"         = "#984EA3",
  "Stromal cells"    = "#56B4E9",
  "T cells"          = "#377EB8"
)

get_cols <- function(cc) {

  lev <- levels(cc@idents)

  cols <- REF_CELL_COLORS[lev]

  if (anyNA(cols)) {

    fallback_names <- lev[is.na(cols)]

    fallback_cols <- grDevices::hcl.colors(
      length(fallback_names),
      "Dark 3"
    )

    names(fallback_cols) <- fallback_names

    cols[is.na(cols)] <- fallback_cols
  }

  cols <- cols[lev]
  names(cols) <- lev

  cols
}

save_network_pair <- function(
    cc,
    file
) {

  lev <- levels(cc@idents)

  colors_use <- get_cols(cc)

  if (
    is.null(names(colors_use)) ||
    !identical(
      names(colors_use),
      lev
    )
  ) {
    stop(
      "Color vector names do not match CellChat identity levels."
    )
  }

  group_size <- as.numeric(
    table(
      factor(
        cc@idents,
        levels = lev
      )
    )
  )

  grDevices::png(
    filename = file,
    width = 10.6,
    height = 4.8,
    units = "in",
    res = 600,
    bg = "white"
  )

  graphics::par(
    mfrow = c(1, 2),
    mar = c(1.0, 1.0, 2.5, 1.0),
    oma = c(0, 0, 0, 0),
    xpd = NA,
    family = FONT
  )

  CellChat::netVisual_circle(
    cc@net$count,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = colors_use,
    vertex.label.cex = 1.05,
    edge.width.max = 8,
    title.name = "Number of interactions"
  )

  CellChat::netVisual_circle(
    cc@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = colors_use,
    vertex.label.cex = 1.05,
    edge.width.max = 8,
    title.name = "Interaction weights/strength"
  )

  grDevices::dev.off()

  if (
    !file.exists(file) ||
    !is.finite(file.info(file)$size) ||
    file.info(file)$size <= 0
  ) {
    stop(
      "Failed to create network PNG: ",
      file
    )
  }

  invisible(file)
}

save_role_scatter <- function(
    cc,
    title_text,
    file
) {

  role_colors <- get_cols(cc)

  p <- CellChat::netAnalysis_signalingRole_scatter(
    cc,
    color.use = role_colors
  ) +
    ggplot2::labs(
      title = title_text,
      x = "Outgoing interaction strength",
      y = "Incoming interaction strength"
    ) +
    ggplot2::theme_classic(
      base_family = FONT,
      base_size = 9
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 10
      ),
      axis.title = ggplot2::element_text(
        size = 9
      ),
      axis.text = ggplot2::element_text(
        size = 8
      ),
      legend.title = ggplot2::element_text(
        size = 8,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 7.5
      )
    )

  ggplot2::ggsave(
    filename = file,
    plot = p,
    width = 5.0,
    height = 4.8,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  invisible(file)
}

pathway_strength <- function(cc) {

  pathways <- as.character(
    cc@netP$pathways
  )

  prob <- cc@netP$prob

  if (!length(pathways)) {
    stop(
      "No signaling pathways found."
    )
  }

  if (
    is.null(prob) ||
    length(dim(prob)) != 3
  ) {
    stop(
      "Invalid CellChat netP$prob."
    )
  }

  if (
    dim(prob)[3] != length(pathways)
  ) {
    stop(
      "Pathway names do not match netP$prob."
    )
  }

  strength <- vapply(
    seq_along(pathways),
    function(i) {
      sum(
        prob[, , i],
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  names(strength) <- pathways

  strength <- strength[
    is.finite(strength) &
      strength > 0
  ]

  sort(
    strength,
    decreasing = TRUE
  )
}

save_one_heatmap <- function(
    cc,
    pathway,
    file
) {

  grDevices::png(
    filename = file,
    width = 2.2,
    height = 3.4,
    units = "in",
    res = 450,
    bg = "white"
  )

  ht <- CellChat::netVisual_heatmap(
    cc,
    signaling = pathway,
    color.heatmap = "Reds",
    title.name = paste0(
      pathway,
      " signaling network"
    ),
    font.size = 6,
    font.size.title = 7
  )

  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right"
  )

  grDevices::dev.off()

  invisible(file)
}

save_top8_pathway_strip <- function(
    cc,
    condition_name,
    file
) {

  strengths <- pathway_strength(cc)

  top8 <- head(
    names(strengths),
    8
  )

  cat(
    "\nTop 8 pathways for ",
    condition_name,
    ":\n",
    paste0(
      seq_along(top8),
      ". ",
      top8,
      "  (",
      signif(
        strengths[top8],
        4
      ),
      ")",
      collapse = "\n"
    ),
    "\n",
    sep = ""
  )

  write.csv(
    data.frame(
      Rank = seq_along(strengths),
      Pathway = names(strengths),
      Total_communication_strength = as.numeric(strengths)
    ),
    out(
      paste0(
        "04_",
        condition_name,
        "_Pathway_Strength_Ranking.csv"
      )
    ),
    row.names = FALSE
  )

  tmp_dir <- file.path(
    output_root,
    paste0(
      "cellchat_",
      condition_name
    )
  )

  dir.create(
    tmp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  tmp_files <- character(
    length(top8)
  )

  for (i in seq_along(top8)) {

    tmp_files[i] <- file.path(
      tmp_dir,
      sprintf(
        "%02d_%s.png",
        i,
        gsub(
          "[^A-Za-z0-9_-]",
          "_",
          top8[i]
        )
      )
    )

    save_one_heatmap(
      cc = cc,
      pathway = top8[i],
      file = tmp_files[i]
    )
  }

  grobs <- lapply(
    tmp_files,
    function(f) {

      img <- png::readPNG(f)

      grid::rasterGrob(
        img,
        interpolate = TRUE
      )
    }
  )

  strip <- cowplot::plot_grid(
    plotlist = grobs,
    nrow = 1,
    align = "h"
  )

  ggplot2::ggsave(
    filename = file,
    plot = strip,
    width = 15.6,
    height = 3.6,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  invisible(
    list(
      file = file,
      top8 = top8,
      strengths = strengths
    )
  )
}

Afile <- out(
  "04_Reference_A_Normal_Networks.png"
)

Bfile <- out(
  "04_Reference_B_Normal_Centrality.png"
)

Cfile <- out(
  "04_Reference_C_Normal_Pathways.png"
)

Dfile <- out(
  "04_Reference_D_Tumor_Networks.png"
)

Efile <- out(
  "04_Reference_E_Tumor_Centrality.png"
)

Ffile <- out(
  "04_Reference_F_Tumor_Pathways.png"
)

normal_combined <- out(
  "04_CellChat_Normal_COMBINED_REFERENCE_STYLE.png"
)

tumor_combined <- out(
  "04_CellChat_Tumor_COMBINED_REFERENCE_STYLE.png"
)

final_both <- out(
  "04_CellChat_Normal_vs_Tumor_REFERENCE_STYLE_FINAL.png"
)

cat(
  "\n============================================================\n"
)

cat(
  "BUILDING REFERENCE-STYLE CELCHAT FIGURES\n"
)

cat(
  "============================================================\n"
)

save_network_pair(
  cellchat_normal,
  Afile
)

save_role_scatter(
  cellchat_normal,
  "Normal",
  Bfile
)

normal_top8 <- save_top8_pathway_strip(
  cellchat_normal,
  "Normal",
  Cfile
)

save_network_pair(
  cellchat_tumor,
  Dfile
)

save_role_scatter(
  cellchat_tumor,
  "Tumor",
  Efile
)

tumor_top8 <- save_top8_pathway_strip(
  cellchat_tumor,
  "Tumor",
  Ffile
)

max_n <- max(
  length(normal_top8$top8),
  length(tumor_top8$top8)
)

write.csv(
  data.frame(
    Normal = c(
      normal_top8$top8,
      rep(
        NA_character_,
        max_n -
          length(normal_top8$top8)
      )
    ),
    Tumor = c(
      tumor_top8$top8,
      rep(
        NA_character_,
        max_n -
          length(tumor_top8$top8)
      )
    )
  ),
  out(
    "04_Selected_Pathways_For_Reference_Figure.csv"
  ),
  row.names = FALSE
)

png_panel <- function(
    file,
    label
) {

  if (
    !file.exists(file) ||
    !is.finite(file.info(file)$size) ||
    file.info(file)$size <= 0
  ) {
    stop(
      "Invalid panel file: ",
      file
    )
  }

  img <- png::readPNG(
    file
  )

  grob <- grid::rasterGrob(
    img,
    interpolate = TRUE
  )

  cowplot::ggdraw() +
    cowplot::draw_grob(
      grob,
      0,
      0,
      1,
      1
    ) +
    cowplot::draw_label(
      label,
      x = 0.006,
      y = 0.994,
      hjust = 0,
      vjust = 1,
      fontfamily = FONT,
      fontface = "bold",
      size = 15
    )
}

A <- png_panel(
  Afile,
  "(A)"
)

B <- png_panel(
  Bfile,
  "(B)"
)

C <- png_panel(
  Cfile,
  "(C)"
)

normal_top <- cowplot::plot_grid(
  A,
  B,
  nrow = 1,
  rel_widths = c(
    2.05,
    1
  )
)

normal_final <- cowplot::plot_grid(
  normal_top,
  C,
  ncol = 1,
  rel_heights = c(
    1,
    0.58
  )
)

ggplot2::ggsave(
  filename = normal_combined,
  plot = normal_final,
  width = 15,
  height = 7.9,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

D <- png_panel(
  Dfile,
  "(D)"
)

E <- png_panel(
  Efile,
  "(E)"
)

F <- png_panel(
  Ffile,
  "(F)"
)

tumor_top <- cowplot::plot_grid(
  D,
  E,
  nrow = 1,
  rel_widths = c(
    2.05,
    1
  )
)

tumor_final <- cowplot::plot_grid(
  tumor_top,
  F,
  ncol = 1,
  rel_heights = c(
    1,
    0.58
  )
)

ggplot2::ggsave(
  filename = tumor_combined,
  plot = tumor_final,
  width = 15,
  height = 7.9,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

final_reference <- cowplot::plot_grid(
  normal_top,
  C,
  tumor_top,
  F,
  ncol = 1,
  rel_heights = c(
    1.05,
    0.72,
    1.05,
    0.72
  )
)

ggplot2::ggsave(
  filename = final_both,
  plot = final_reference,
  width = 15,
  height = 13,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

files_to_check <- c(
  Afile,
  Bfile,
  Cfile,
  Dfile,
  Efile,
  Ffile,
  normal_combined,
  tumor_combined,
  final_both
)

cat(
  "\n============================================================\n"
)

cat(
  "FINAL FIGURE FILE CHECK\n"
)

cat(
  "============================================================\n"
)

for (f in files_to_check) {

  status <- file.exists(f) &&
    is.finite(file.info(f)$size) &&
    file.info(f)$size > 0

  cat(
    if (status) {
      "[OK]   "
    } else {
      "[FAIL] "
    },
    basename(f),
    if (status) {
      paste0(
        "  ",
        round(
          file.info(f)$size /
            1024^2,
          2
        ),
        " MB"
      )
    } else {
      ""
    },
    "\n",
    sep = ""
  )
}

cat(
  "\nFinal figure:\n",
  final_both,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)
