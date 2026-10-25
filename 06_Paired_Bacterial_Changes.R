# GSE334323: descriptive bacterial changes in two matched tumor-normal pairs.
# Reuses Stage 05 spatial tests and original FDR values; no tests are refitted.
local({
  input_dir <- Sys.getenv("CRC_ASSOCIATION_OUTPUT", "E:/CRC/Results/Results/05-Bacteria_Host_Spatial_Association")
  output_dir <- Sys.getenv("CRC_PAIRED_BACTERIAL_OUTPUT", "E:/CRC/Results/Results/06-Paired_Bacterial_CRC_Association")
  sections <- c("CRC_03_Normal", "CRC_03_Tumor", "CRC_08_Normal", "CRC_08_Tumor")
  patients <- c("CRC_03", "CRC_08")
  targets <- c("EPCAM", "Signature__Epithelial", "Myeloids")
  read_tab <- function(name, required) {
    path <- file.path(input_dir, name)
    if (!file.exists(path)) stop("Missing Stage 05 input: ", path)
    d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(required %in% names(d))) stop("Missing columns in ", name, ": ", paste(setdiff(required, names(d)), collapse = ", "))
    d
  }
  key <- function(d) paste(d$section, d$barcode, sep = "::")
  completion_file <- file.path(input_dir, "Core_Association_Complete.rds")
  if (!file.exists(completion_file)) stop("Stage 05 completion record is missing")
  completion <- readRDS(completion_file)
  hashes <- unname(tools::md5sum(file.path(input_dir, completion$files)))
  if (!identical(completion$revision, "core-v1") || !length(completion$files) ||
      anyNA(hashes) || !identical(completion$md5, hashes)) stop("Stage 05 checksum mismatch")

  audit <- read_tab("02_Bacterial_Extraction_Check.csv", c("Section", "Exact_total_agreement"))
  if (nrow(audit) != 4L || anyDuplicated(audit$Section) || !setequal(audit$Section, sections) ||
      !all(tolower(as.character(audit$Exact_total_agreement)) %in% c("true", "1"))) stop("Bacterial extraction audit failed")
  taxa <- read_tab("03_Bacterial_Groups_All_QC_Spots.csv", c("section", "barcode", "bacterial_group", "bacterial_count"))
  ids <- c("section", "barcode", "bacterial_group")
  if (anyNA(taxa[ids]) || any(vapply(taxa[ids], function(x) any(!nzchar(trimws(x))), logical(1)))) stop("Invalid spot/group identifiers")
  groups <- sort(unique(taxa$bacterial_group))
  if (length(groups) != 19L || !setequal(taxa$section, sections)) stop("Expected 19 groups and four sections")
  if (anyDuplicated(paste(key(taxa), taxa$bacterial_group, sep = "::"))) stop("Duplicate group/spot rows")
  v <- taxa$bacterial_count
  if (!is.numeric(v) || any(!is.finite(v) | v < 0 | v != floor(v))) stop("Invalid bacterial counts")
  totals <- tapply(v, key(taxa), sum)
  for (ss in sections) {
    z <- taxa[taxa$section == ss, , drop = FALSE]
    bars <- unique(z$barcode)
    for (g in groups) if (!setequal(z$barcode[z$bacterial_group == g], bars)) stop("Incomplete spot grid: ", ss, " / ", g)
  }

  # Zero-inclusive raw probe summaries; fractions are compositional, not host-normalized.
  summaries <- do.call(rbind, lapply(sections, function(ss) do.call(rbind, lapply(groups, function(g) {
    z <- taxa[taxa$section == ss & taxa$bacterial_group == g, , drop = FALSE]
    x <- z$bacterial_count; total <- as.numeric(totals[key(z)])
    data.frame(Patient = sub("_(Normal|Tumor)$", "", ss), Tissue = sub("^CRC_[0-9]+_", "", ss),
      Section = ss, Bacterial_group = g, QC_spots = length(x), Total_probe_counts = sum(x),
      Mean_counts_per_QC_spot = mean(x), Median_counts_per_QC_spot = median(x),
      Positive_spots = sum(x > 0), Detection_percent = 100 * mean(x > 0),
      Fraction_of_all_bacterial_counts = if (sum(total) > 0) sum(x) / sum(total) else NA_real_)
  }))))
  paired <- do.call(rbind, lapply(patients, function(pt) do.call(rbind, lapply(groups, function(g) {
    z <- summaries[summaries$Patient == pt & summaries$Bacterial_group == g, , drop = FALSE]
    normal <- z[z$Tissue == "Normal", , drop = FALSE]
    tumor <- z[z$Tissue == "Tumor", , drop = FALSE]
    if (nrow(normal) != 1L || nrow(tumor) != 1L) stop("Incomplete section pairing")
    delta <- tumor$Mean_counts_per_QC_spot - normal$Mean_counts_per_QC_spot
    data.frame(Patient = pt, Bacterial_group = g,
      Normal_mean_counts = normal$Mean_counts_per_QC_spot, Tumor_mean_counts = tumor$Mean_counts_per_QC_spot,
      Mean_count_change = delta,
      Mean_count_ratio = if (normal$Mean_counts_per_QC_spot > 0) tumor$Mean_counts_per_QC_spot / normal$Mean_counts_per_QC_spot else NA_real_,
      Direction = if (delta > 0) "Higher_in_tumor" else if (delta < 0) "Lower_in_tumor" else "Unchanged",
      Normal_detection_percent = normal$Detection_percent, Tumor_detection_percent = tumor$Detection_percent,
      Detection_change_percentage_points = tumor$Detection_percent - normal$Detection_percent,
      Bacterial_fraction_change = tumor$Fraction_of_all_bacterial_counts - normal$Fraction_of_all_bacterial_counts)
  }))))
  consistency <- do.call(rbind, lapply(groups, function(g) {
    z <- paired[paired$Bacterial_group == g, , drop = FALSE]
    data.frame(Bacterial_group = g, Patients = nrow(z), Higher_in_tumor = sum(z$Mean_count_change > 0),
      Lower_in_tumor = sum(z$Mean_count_change < 0),
      Same_direction_both_patients = all(z$Mean_count_change > 0) || all(z$Mean_count_change < 0),
      Median_paired_mean_count_change = median(z$Mean_count_change),
      Detection_increased_both = all(z$Detection_change_percentage_points > 0),
      Detection_decreased_both = all(z$Detection_change_percentage_points < 0))
  }))

  required <- c("Family", "Section", "Bacterial_group", "Host_target", "Nclass", "Rho", "Status", "P_spatial", "FDR_spatial")
  check_tests <- function(d, classes) {
    ids <- c("Family", "Section", "Bacterial_group", "Host_target", "Status")
    if (anyNA(d[ids]) || any(vapply(d[ids], function(x) any(!nzchar(trimws(x))), logical(1)))) stop("Invalid spatial identifiers")
    test_id <- paste(d$Family, d$Section, d$Bacterial_group, d$Host_target, d$Nclass, sep = "::")
    if (anyDuplicated(test_id) || anyNA(d$Nclass) || !all(d$Nclass %in% classes)) stop("Invalid spatial test grid")
    for (nm in c("P_spatial", "FDR_spatial")) if (!is.numeric(d[[nm]]) ||
      any(!is.na(d[[nm]]) & (!is.finite(d[[nm]]) | d[[nm]] < 0 | d[[nm]] > 1))) stop("Invalid P/FDR values")
    if (!is.numeric(d$Rho) || any(!is.na(d$Rho) & (!is.finite(d$Rho) | abs(d$Rho) > 1 + 1e-8))) stop("Invalid rho")
    z <- d[d$Family == "Taxon_SameSpot" & d$Host_target %in% targets, , drop = FALSE]
    expected <- expand.grid(Section = sections, Bacterial_group = groups, Host_target = targets, Nclass = classes, stringsAsFactors = FALSE)
    test_key <- function(x) paste(x$Section, x$Bacterial_group, x$Host_target, x$Nclass, sep = "::")
    if (!setequal(test_key(z), test_key(expected))) stop("Incomplete all-group spatial tests")
    z$Spatial_FDR_significant <- z$Status == "OK" & is.finite(z$FDR_spatial) & z$FDR_spatial < 0.05
    z
  }
  spatial <- check_tests(read_tab("06_Primary_Spatial_Adjusted_Correlations.csv", required), 13L)
  sensitivity <- check_tests(read_tab("07_Spatial_Sensitivity.csv", required), c(8L, 20L))
  spatial <- merge(spatial, consistency, by = "Bacterial_group", all.x = TRUE, sort = FALSE)
  # Post hoc descriptive subset; original all-group FDR values remain unchanged.
  selected <- spatial[spatial$Same_direction_both_patients, , drop = FALSE]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory")
  results <- list(
    "01_Section_Bacterial_Summary" = summaries,
    "02_Patient_Paired_Changes" = paired,
    "03_Between_Patient_Direction_Consistency" = consistency,
    "04_All19_Spatial_Host_Associations" = spatial,
    "05_Spatial_Sensitivity" = sensitivity,
    "06_Consistent_Groups_Exploratory_Spatial" = selected)
  for (nm in names(results)) utils::write.csv(results[[nm]], file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE, na = "")
})
