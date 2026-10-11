# Predictive Modeling
# Pseudobulk elastic net, random forest and XGBoost; patient-level LOPO validation.
# GSE144735 external patients are independent of training but from the same publication.
# Development-only fitting and tuning; six external pairs do not establish clinical validity.

local({

  # Paths and external input download

  sc_root <- "E:/CRC/Results/GSE132465_scRNA/Core_From_Raw"

  output_root <- "E:/CRC/Results/Results/06_Machine Learning"

  if (!dir.exists(sc_root)) {
    stop("Missing scRNA directory: ", sc_root)
  }

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_root)) {
    stop("Cannot create output directory: ", output_root)
  }

  out <- function(name) {
    file.path(output_root, name)
  }

  external_input_root <- file.path(output_root, "GSE144735_Input")
  dir.create(external_input_root, recursive=TRUE, showWarnings=FALSE)
  geo_base <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144735/suppl/"
  annotation_name <- "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"
  counts_name <- "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"
  fetch_geo <- function(name) {
    dest <- file.path(external_input_root, name)
    if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
    tmp <- paste0(dest, ".partial")
    old_timeout <- getOption("timeout"); options(timeout=max(1800,old_timeout))
    on.exit(options(timeout=old_timeout),add=TRUE)
    cat("Downloading GEO input: ",name,"\n",sep="")
    status <- tryCatch(download.file(paste0(geo_base,name),tmp,method="libcurl",mode="wb"),
                       error=function(e) { if(file.exists(tmp)) unlink(tmp); stop(conditionMessage(e)) })
    if(!identical(status,0L) || !file.exists(tmp) || file.info(tmp)$size==0) {
      if(file.exists(tmp)) unlink(tmp)
      stop("Download failed. Place the named GEO file in: ",external_input_root)
    }
    if(!file.rename(tmp,dest)) stop("Cannot finalize downloaded GEO file: ",dest)
    dest
  }
  annotation_file <- fetch_geo(annotation_name)
  count_file <- fetch_geo(counts_name)
  ext_meta <- read.delim(gzfile(annotation_file),check.names=FALSE,stringsAsFactors=FALSE)
  required <- c("Index","Patient","Class","Sample")
  if(!all(required %in% names(ext_meta)) || anyNA(ext_meta[,required])) stop("Invalid KUL3 annotation columns.")
  if(anyDuplicated(ext_meta$Index) || any(!nzchar(as.matrix(ext_meta[,required])))) stop("Invalid KUL3 cell metadata.")
  expected_patients <- c("KUL01","KUL19","KUL21","KUL28","KUL30","KUL31")
  expected_samples <- as.vector(outer(expected_patients,c("-N","-T","-B"),paste0))
  if(!setequal(ext_meta$Patient,expected_patients) || !setequal(ext_meta$Sample,expected_samples)) {
    stop("Unexpected GEO patient/sample IDs; verify input release before continuing.")
  }
  if(any(sub("-[NTB]$","",ext_meta$Sample)!=ext_meta$Patient)) stop("GEO sample/patient mismatch.")

  if(any(ext_meta$Class[grepl("-N$",ext_meta$Sample)]!="Normal") ||
     any(ext_meta$Class[grepl("-T$",ext_meta$Sample)]!="Tumor")) stop("GEO tissue labels inconsistent.")
  cat("GSE144735 deposited cells: ",nrow(ext_meta),"; patients: ",length(unique(ext_meta$Patient)),"\n",sep="")

  # Packages

  packages <- c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "ggplot2",
    "glmnet",
    "ranger",
    "xgboost",
    "pROC"
  )

  for (p in packages) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Install package first: ", p)
    }
  }

  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(pROC)

  set.seed(42)

  # Input file

  seurat_file <- file.path(sc_root, "05_annotated.rds")

  if (!file.exists(seurat_file)) {
    stop("Missing: ", seurat_file)
  }

  # Load gse132465

  obj <- readRDS(seurat_file)
  meta <- obj@meta.data

  required_meta <- c("sample_id", "patient_id", "group")

  if (!all(required_meta %in% colnames(meta))) {
    stop(
      "Missing metadata columns: ",
      paste(setdiff(required_meta, colnames(meta)), collapse = ", ")
    )
  }

  meta$sample_id <- as.character(meta$sample_id)
  meta$patient_id <- as.character(meta$patient_id)
  meta$group <- as.character(meta$group)

  patients <- sprintf("SMC%02d", 1:10)

  samples <- as.vector(
    t(
      outer(
        patients,
        c("-N", "-T"),
        paste0
      )
    )
  )

  if (!setequal(unique(meta$patient_id), patients)) {
    stop("Unexpected patient IDs in GSE132465.")
  }

  if (!setequal(unique(meta$sample_id), samples)) {
    stop("Expected 20 paired GSE132465 samples.")
  }

  if (!all(meta$group %in% c("Normal", "Tumor"))) {
    stop("group must contain Normal and Tumor only.")
  }

  expected_sample <- paste0(meta$patient_id, ifelse(meta$group == "Normal", "-N", "-T"))
  if (anyNA(meta[, required_meta]) || any(meta$sample_id != expected_sample)) {
    stop("Invalid sample/patient/tissue pairing.")
  }

  # Define gene panel

  candidate_genes <- unique(
    c(

      "EPCAM", "KRT8", "KRT18", "CEACAM5",

      "CD3D", "CD3E", "CD4", "CD8A",

      "MS4A1", "CD79A",

      "CD68", "LYZ",

      "NKG7", "GNLY",

      "PECAM1", "VWF"
    )
  )

  present_genes <- intersect(candidate_genes, rownames(obj))
  missing_internal_genes <- setdiff(candidate_genes, rownames(obj))

  cat(
    "\nInternal genes available: ",
    length(present_genes),
    "/",
    length(candidate_genes),
    "\n",
    sep = ""
  )

  if (length(missing_internal_genes) > 0) {
    cat(
      "Genes missing from GSE132465: ",
      paste(missing_internal_genes, collapse = ", "),
      "\n",
      sep = ""
    )
  }

  if (length(present_genes) < 2) {
    stop("Too few genes available in GSE132465.")
  }

  # Define gene signatures

  gene_signatures <- list(
    Epithelial = c("EPCAM", "KRT8", "KRT18", "CEACAM5"),
    T_cell = c("CD3D", "CD3E", "CD4", "CD8A"),
    B_cell = c("MS4A1", "CD79A"),
    Myeloid = c("CD68", "LYZ"),
    Cytotoxic = c("NKG7", "GNLY"),
    Endothelial = c("PECAM1", "VWF")
  )

  # Get raw counts from seurat

  get_counts <- function(o) {

    if (inherits(o[["RNA"]], "Assay5")) {

      layers <- SeuratObject::Layers(
        o[["RNA"]],
        search = "^counts"
      )

      if (!length(layers)) {
        stop("RNA counts layer is missing.")
      }

      if (length(layers) > 1) {
        o <- SeuratObject::JoinLayers(
          o,
          assay = "RNA",
          layers = layers,
          new = "counts"
        )
      }

      counts <- SeuratObject::GetAssayData(
        o,
        assay = "RNA",
        layer = "counts"
      )

    } else {

      counts <- SeuratObject::GetAssayData(
        o,
        assay = "RNA",
        slot = "counts"
      )
    }

    as(counts, "dgCMatrix")
  }

  counts <- get_counts(obj)

  if (!identical(colnames(counts), rownames(meta))) {
    stop("Counts and metadata cell order do not match.")
  }

  # Gse132465 pseudobulk counts

  gene_counts <- matrix(
    0,
    nrow = length(present_genes),
    ncol = length(samples),
    dimnames = list(present_genes, samples)
  )

  library_size <- setNames(
    numeric(length(samples)),
    samples
  )

  for (s in samples) {

    cell_ids <- rownames(meta)[meta$sample_id == s]

    if (!length(cell_ids)) {
      stop("No cells found for sample: ", s)
    }

    gene_counts[, s] <- Matrix::rowSums(
      counts[
        present_genes,
        cell_ids,
        drop = FALSE
      ]
    )

    library_size[s] <- sum(
      counts[
        ,
        cell_ids,
        drop = FALSE
      ]
    )
  }

  if (any(library_size <= 0)) {
    stop("Invalid pseudobulk library size.")
  }

  # Internal log2 cpm

  internal_cpm <- sweep(
    gene_counts,
    2,
    library_size,
    "/"
  ) * 1e6

  internal_log_cpm <- log2(internal_cpm + 1)

  # Internal gene features

  internal_gene_features <- as.data.frame(
    t(internal_log_cpm),
    check.names = FALSE
  )

  internal_gene_features$sample <- rownames(internal_gene_features)
  rownames(internal_gene_features) <- NULL

  raw_gene_columns <- setdiff(
    colnames(internal_gene_features),
    "sample"
  )

  colnames(internal_gene_features)[
    match(raw_gene_columns, colnames(internal_gene_features))
  ] <- paste0("Gene__", raw_gene_columns)

  # Internal signature features

  signature_table <- data.frame(
    sample = samples,
    stringsAsFactors = FALSE
  )

  for (sig in names(gene_signatures)) {

    genes_i <- intersect(
      gene_signatures[[sig]],
      present_genes
    )

    if (!length(genes_i)) {
      next
    }

    feature_name <- paste0("Signature__", sig)

    signature_table[[feature_name]] <- colMeans(
      internal_log_cpm[
        genes_i,
        samples,
        drop = FALSE
      ]
    )
  }

  # Internal ml dataset

  sample_info <- unique(
    meta[
      ,
      c("sample_id", "patient_id", "group"),
      drop = FALSE
    ]
  )

  colnames(sample_info) <- c(
    "sample",
    "patient",
    "group"
  )

  sample_info <- sample_info[
    match(samples, sample_info$sample),
    ,
    drop = FALSE
  ]

  if (anyNA(sample_info$sample)) {
    stop("Failed to match one or more internal samples.")
  }

  ml_data <- merge(
    sample_info,
    internal_gene_features,
    by = "sample",
    all.x = TRUE,
    sort = FALSE
  )

  ml_data <- merge(
    ml_data,
    signature_table,
    by = "sample",
    all.x = TRUE,
    sort = FALSE
  )

  ml_data <- ml_data[
    match(samples, ml_data$sample),
    ,
    drop = FALSE
  ]

  ml_data$target <- ifelse(
    ml_data$group == "Tumor",
    1L,
    0L
  )

  if (nrow(ml_data) != 20) {
    stop("Internal ML dataset must contain exactly 20 samples.")
  }

  if (anyNA(ml_data)) {
    stop("Missing values detected in internal ML dataset.")
  }

  write.csv(
    ml_data,
    out("05_Internal_ML_Dataset.csv"),
    row.names = FALSE
  )

  # Internal feature matrix

  feature_names <- grep(
    "^Gene__|^Signature__",
    colnames(ml_data),
    value = TRUE
  )

  X_internal <- as.matrix(
    ml_data[
      ,
      feature_names,
      drop = FALSE
    ]
  )

  storage.mode(X_internal) <- "double"

  y_internal <- ml_data$target

  if (anyNA(X_internal) || any(!is.finite(X_internal))) {
    stop("Invalid values in internal feature matrix.")
  }

  # Preprocessing function

  prepare_fold <- function(
    X_train,
    X_test
  ) {

    train_sd <- apply(
      X_train,
      2,
      sd
    )

    keep <- is.finite(train_sd) & train_sd > 0

    if (!any(keep)) {
      stop("No non-zero-variance features remain in this fold.")
    }

    X_train <- X_train[, keep, drop = FALSE]
    X_test <- X_test[, keep, drop = FALSE]

    mu <- colMeans(X_train)
    sigma <- apply(X_train, 2, sd)

    X_train_scaled <- sweep(
      X_train,
      2,
      mu,
      "-"
    )

    X_train_scaled <- sweep(
      X_train_scaled,
      2,
      sigma,
      "/"
    )

    X_test_scaled <- sweep(
      X_test,
      2,
      mu,
      "-"
    )

    X_test_scaled <- sweep(
      X_test_scaled,
      2,
      sigma,
      "/"
    )

    list(
      train = X_train_scaled,
      test = X_test_scaled,
      features = colnames(X_train_scaled),
      mean = mu,
      sd = sigma
    )
  }

  # Performance function

  evaluate_model <- function(
    truth,
    probability,
    model_name,
    dataset_name
  ) {

    roc_object <- pROC::roc(
      response = truth,
      predictor = probability,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    auc_value <- as.numeric(
      pROC::auc(roc_object)
    )

    predicted <- ifelse(
      probability >= 0.5,
      1L,
      0L
    )

    TP <- sum(predicted == 1 & truth == 1)
    TN <- sum(predicted == 0 & truth == 0)
    FP <- sum(predicted == 1 & truth == 0)
    FN <- sum(predicted == 0 & truth == 1)

    sensitivity <- ifelse(
      TP + FN > 0,
      TP / (TP + FN),
      NA_real_
    )

    specificity <- ifelse(
      TN + FP > 0,
      TN / (TN + FP),
      NA_real_
    )

    accuracy <- (TP + TN) / length(truth)

    balanced_accuracy <- mean(
      c(sensitivity, specificity),
      na.rm = TRUE
    )

    data.frame(
      Dataset = dataset_name,
      Model = model_name,
      ROC_AUC = auc_value,
      Accuracy = accuracy,
      Balanced_Accuracy = balanced_accuracy,
      Sensitivity = sensitivity,
      Specificity = specificity,
      TP = TP,
      TN = TN,
      FP = FP,
      FN = FN,
      stringsAsFactors = FALSE
    )
  }

  # Internal lopo-cv

  internal_predictions <- list()

  for (i in seq_along(patients)) {

    held_patient <- patients[i]

    cat(
      "\nInternal LOPO fold ",
      i,
      "/",
      length(patients),
      " — held-out: ",
      held_patient,
      "\n",
      sep = ""
    )

    test_idx <- which(
      ml_data$patient == held_patient
    )

    train_idx <- setdiff(
      seq_len(nrow(ml_data)),
      test_idx
    )

    X_train <- X_internal[
      train_idx,
      ,
      drop = FALSE
    ]

    X_test <- X_internal[
      test_idx,
      ,
      drop = FALSE
    ]

    y_train <- y_internal[train_idx]
    y_test <- y_internal[test_idx]

    prep <- prepare_fold(
      X_train,
      X_test
    )

    X_train_s <- prep$train
    X_test_s <- prep$test

    training_patients <- ml_data$patient[train_idx]

    lambda_grid <- 10^seq(2, -5, length.out = 100)
    inner_patients <- sort(unique(training_patients))
    inner_deviance <- matrix(
      NA_real_, nrow = length(inner_patients),
      ncol = length(lambda_grid)
    )

    for (j in seq_along(inner_patients)) {
      inner_valid <- which(training_patients == inner_patients[j])
      inner_train <- setdiff(seq_len(nrow(X_train)), inner_valid)
      inner_prep <- prepare_fold(
        X_train[inner_train, , drop = FALSE],
        X_train[inner_valid, , drop = FALSE]
      )
      if (ncol(inner_prep$train) < 2) {
        stop("Elastic Net requires at least two variable features.")
      }
      inner_fit <- glmnet::glmnet(
        x = inner_prep$train,
        y = y_train[inner_train],
        family = "binomial", alpha = 0.5,
        lambda = lambda_grid, standardize = FALSE
      )
      inner_probability <- as.matrix(predict(
        inner_fit, newx = inner_prep$test,
        s = lambda_grid, type = "response"
      ))
      if (any(!is.finite(inner_probability))) {
        stop("Non-finite inner Elastic Net predictions.")
      }
      inner_probability <- pmin(pmax(inner_probability, 1e-15), 1 - 1e-15)
      inner_truth <- matrix(
        y_train[inner_valid], nrow = length(inner_valid),
        ncol = length(lambda_grid)
      )
      inner_deviance[j, ] <- colMeans(-2 * (
        inner_truth * log(inner_probability) +
          (1 - inner_truth) * log1p(-inner_probability)
      ))
    }

    mean_inner_deviance <- colMeans(inner_deviance)
    selected_lambda <- lambda_grid[which.min(mean_inner_deviance)]
    cat("Elastic Net selected lambda: ", selected_lambda, "\n", sep = "")

    if (ncol(X_train_s) < 2) {
      stop("Elastic Net requires at least two variable features.")
    }
    elastic_fit <- glmnet::glmnet(
      x = X_train_s, y = y_train,
      family = "binomial", alpha = 0.5,
      lambda = lambda_grid, standardize = FALSE
    )
    pred_elastic <- as.numeric(predict(
      elastic_fit, newx = X_test_s,
      s = selected_lambda, type = "response"
    ))

    rf_train <- data.frame(
      target = factor(
        y_train,
        levels = c(0, 1),
        labels = c("Normal", "Tumor")
      ),
      as.data.frame(
        X_train_s,
        check.names = FALSE
      )
    )

    rf_test <- as.data.frame(
      X_test_s,
      check.names = FALSE
    )

    rf_model <- ranger::ranger(
      target ~ .,
      data = rf_train,
      probability = TRUE,
      num.trees = 1000,
      mtry = max(
        1,
        floor(
          sqrt(
            ncol(X_train_s)
          )
        )
      ),
      min.node.size = 2,
      seed = 42
    )

    pred_rf <- predict(
      rf_model,
      data = rf_test
    )$predictions[, "Tumor"]

    dtrain <- xgb.DMatrix(
      data = X_train_s,
      label = y_train
    )

    dtest <- xgb.DMatrix(
      data = X_test_s
    )

    xgb_model <- xgboost::xgb.train(
      params = list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        eta = 0.05,
        max_depth = 2,
        min_child_weight = 2,
        subsample = 0.8,
        colsample_bytree = 0.8,
        lambda = 1,
        alpha = 0.1,
        nthread = 1
      ),
      data = dtrain,
      nrounds = 75,
      verbose = 0
    )

    pred_xgb <- as.numeric(
      predict(
        xgb_model,
        dtest
      )
    )

    internal_predictions[[i]] <- data.frame(
      patient = ml_data$patient[test_idx],
      sample = ml_data$sample[test_idx],
      observed = y_test,
      observed_label = ifelse(
        y_test == 1,
        "Tumor",
        "Normal"
      ),
      ElasticNet = pred_elastic,
      RandomForest = pred_rf,
      XGBoost = pred_xgb,
      stringsAsFactors = FALSE
    )
  }

  internal_pred <- do.call(
    rbind,
    internal_predictions
  )

  rownames(internal_pred) <- NULL

  write.csv(
    internal_pred,
    out("05_Internal_LOPO_Predictions.csv"),
    row.names = FALSE
  )

  internal_performance <- rbind(
    evaluate_model(
      internal_pred$observed,
      internal_pred$ElasticNet,
      "Elastic Net",
      "GSE132465 LOPO"
    ),
    evaluate_model(
      internal_pred$observed,
      internal_pred$RandomForest,
      "Random Forest",
      "GSE132465 LOPO"
    ),
    evaluate_model(
      internal_pred$observed,
      internal_pred$XGBoost,
      "XGBoost",
      "GSE132465 LOPO"
    )
  )

  write.csv(
    internal_performance,
    out("05_Internal_LOPO_Performance.csv"),
    row.names = FALSE
  )

  # Final models: development cohort only

  cat("\nFitting final models on GSE132465 only...\n")
  tune_final <- function(X_train, y_train, training_patients) {
    lambda_grid <- 10^seq(2, -5, length.out = 100)
    ids <- sort(unique(training_patients))
    losses <- matrix(NA_real_, length(ids), length(lambda_grid))
    for (j in seq_along(ids)) {
      va <- which(training_patients == ids[j]); tr <- setdiff(seq_len(nrow(X_train)), va)
      pp <- prepare_fold(X_train[tr,,drop=FALSE], X_train[va,,drop=FALSE])
      if (ncol(pp$train) < 2) stop("Too few variable features for final tuning.")
      fit <- glmnet::glmnet(pp$train, y_train[tr], family="binomial", alpha=.5,
                            lambda=lambda_grid, standardize=FALSE)
      pr <- as.matrix(predict(fit, newx=pp$test, s=lambda_grid, type="response"))
      if (any(!is.finite(pr))) stop("Invalid final tuning predictions.")
      pr <- pmin(pmax(pr, 1e-15), 1-1e-15)
      yy <- matrix(y_train[va], length(va), length(lambda_grid))
      losses[j,] <- colMeans(-2*(yy*log(pr)+(1-yy)*log1p(-pr)))
    }
    list(lambda=lambda_grid[which.min(colMeans(losses))], grid=lambda_grid,
         tuning=data.frame(lambda=lambda_grid, mean_deviance=colMeans(losses)))
  }
  final_tuning <- tune_final(X_internal, y_internal, ml_data$patient)
  final_prep <- prepare_fold(X_internal, X_internal)
  xf <- final_prep$train
  final_elastic <- glmnet::glmnet(xf, y_internal, family="binomial", alpha=.5,
                                  lambda=final_tuning$grid, standardize=FALSE)
  final_rf <- ranger::ranger(
    x=as.data.frame(xf, check.names=FALSE),
    y=factor(y_internal, levels=0:1, labels=c("Normal","Tumor")),
    probability=TRUE, num.trees=1000, mtry=max(1,floor(sqrt(ncol(xf)))),
    min.node.size=2, seed=42)
  final_xgb <- xgboost::xgb.train(
    params=list(objective="binary:logistic", eval_metric="logloss", eta=.05,
                max_depth=2, min_child_weight=2, subsample=.8, colsample_bytree=.8,
                lambda=1, alpha=.1, nthread=1, seed=42),
    data=xgboost::xgb.DMatrix(xf, label=y_internal), nrounds=75, verbose=0)
  saveRDS(list(ElasticNet=final_elastic, RandomForest=final_rf,
               lambda=final_tuning$lambda, preprocessing=final_prep[c("features","mean","sd")],
               candidate_genes=candidate_genes, gene_signatures=gene_signatures,
               training_samples=ml_data[,c("sample","patient","group")],
               threshold=.5), out("06_Final_Models_and_Preprocessing.rds"))
  xgboost::xgb.save(final_xgb, out("06_Final_XGBoost.ubj"))
  write.csv(final_tuning$tuning, out("06_Final_ElasticNet_Tuning.csv"), row.names=FALSE)
  writeLines(capture.output(sessionInfo()), out("06_ML_sessionInfo.txt"))

  read_matrix_header <- function(con) {
    line <- readLines(con,n=1,warn=FALSE)
    if(length(line)!=1) stop("Empty raw UMI file.")
    head <- strsplit(line,"\t",fixed=TRUE)[[1]]
    head <- sub('^"(.*)"$','\\1',head)
    if(head[1]!="Index" || length(head)<2) stop("Unexpected raw UMI matrix header.")
    cells <- head[-1]
    if(anyDuplicated(cells) || !setequal(cells,ext_meta$Index)) stop("UMI/annotation cells do not match.")
    cells
  }
  stream_umi <- function(callback) {
    con <- gzfile(count_file,"rt");on.exit(close(con),add=TRUE)
    cells <- read_matrix_header(con)
    row_id <- 0L
    repeat {
      chunk <- readLines(con,n=100,warn=FALSE)
      if(!length(chunk)) break
      for(line in chunk) {
        row_id <- row_id+1L
        parts <- strsplit(line,"\t",fixed=TRUE)[[1]]
        if(length(parts)!=length(cells)+1L) stop("Bad UMI row length at row ",row_id)
        gene <- sub('^"(.*)"$','\\1',parts[1])
        v <- suppressWarnings(as.numeric(parts[-1]))
        if(!nzchar(gene) || any(!is.finite(v) | v<0 | v!=floor(v))) stop("Invalid raw UMI row: ",row_id)
        callback(gene,v,row_id,cells)
      }
      if(row_id %% 1000L==0L) cat("  UMI genes scanned: ",row_id,"\n",sep="")
    }
    row_id
  }
  con <- gzfile(count_file,"rt");matrix_cells <- read_matrix_header(con);close(con)
  ext_meta <- ext_meta[match(matrix_cells,ext_meta$Index),,drop=FALSE]
  if(any(ext_meta$Patient %in% ml_data$patient) || any(ext_meta$Index %in% colnames(obj))) {
    stop("Observed overlap with development identifiers.")
  }
  ncell <- length(matrix_cells)
  total_umi <- numeric(ncell); detected <- numeric(ncell); mito_umi <- numeric(ncell)
  keep_gene <- logical();gene_ids <- character()
  cat("\nExternal UMI pass 1: gene filtering and cell QC...\n")
  nr <- stream_umi(function(g,v,i,cells) {
    gene_ids[i] <<- g
    keep_gene[i] <<- sum(v>0)>=3
    if(keep_gene[i]) {
      total_umi <<- total_umi+v
      detected <<- detected+(v>0)
      if(grepl("^MT-",g)) mito_umi <<- mito_umi+v
    }
  })
  if(anyDuplicated(gene_ids)) stop("Duplicate gene symbols in GEO raw matrix; resolve before analysis.")
  if(!any(grepl("^MT-",gene_ids[keep_gene]))) stop("No mitochondrial gene symbols found for QC.")
  pct_mt <- ifelse(total_umi>0,100*mito_umi/total_umi,Inf)
  cell_keep <- detected>200 & detected<6000 & pct_mt<20 & total_umi>0
  qc <- data.frame(cell=matrix_cells,patient=ext_meta$Patient,sample=ext_meta$Sample,
                   deposited_class=ext_meta$Class,nFeature=detected,nCount=total_umi,
                   percent_mt=pct_mt,retained=cell_keep)
  write.csv(qc,out("06_GSE144735_Cell_QC.csv"),row.names=FALSE)
  if(!any(cell_keep)) stop("No external cells passed fixed QC.")
  if(!all(present_genes %in% gene_ids[keep_gene])) stop("Development panel genes missing after external gene QC: ",
                                                        paste(setdiff(present_genes,gene_ids[keep_gene]),collapse=", "))
  specimen_ids <- sort(unique(ext_meta$Sample))
  sample_indices <- setNames(lapply(specimen_ids,function(s) which(cell_keep & ext_meta$Sample==s)),specimen_ids)
  nc <- lengths(sample_indices)
  if(any(nc==0)) stop("An external specimen has no cells after QC: ",paste(names(nc)[nc==0],collapse=", "))
  ext_counts <- matrix(0,length(present_genes),length(specimen_ids),dimnames=list(present_genes,specimen_ids))
  ext_lib <- setNames(numeric(length(specimen_ids)),specimen_ids)
  cat("\nExternal UMI pass 2: specimen pseudobulk aggregation...\n")
  nr2 <- stream_umi(function(g,v,i,cells) {
    if(!identical(g,gene_ids[i])) stop("Raw input changed between passes.")
    if(keep_gene[i]) {
      sums <- vapply(sample_indices,function(idx) sum(v[idx]),numeric(1))
      ext_lib <<- ext_lib+sums
      if(g %in% present_genes) ext_counts[g,] <<- sums
    }
  })
  if(nr2!=nr || any(ext_lib<=0)) stop("Invalid external pseudobulk libraries.")
  expected_lib <- vapply(sample_indices,function(idx) sum(total_umi[idx]),numeric(1))
  if(any(ext_lib!=expected_lib)) stop("Two-pass library totals disagree.")
  sm <- data.frame(sample=specimen_ids,patient=sub("-[NTB]$","",specimen_ids),
                   region=sub("^.*-","",specimen_ids),cells=as.integer(nc),library_size=as.numeric(ext_lib))
  sm$group <- ifelse(sm$region=="N","Normal","Tumor")
  sm$comparison <- ifelse(sm$region=="B","Border sensitivity","Primary core vs normal")
  write.csv(sm,out("06_GSE144735_Specimen_Summary.csv"),row.names=FALSE)
  write.csv(data.frame(gene=gene_ids,retained=keep_gene),out("06_GSE144735_Gene_QC.csv"),row.names=FALSE)
  cat("External retained cells: ",sum(cell_keep),"/",ncell,"; genes: ",sum(keep_gene),"/",nr,"\n",sep="")
  print(sm,row.names=FALSE)

  # Identical features and development-only scaling

  logcpm <- log2(sweep(ext_counts,2,ext_lib,"/")*1e6+1)
  ex <- as.data.frame(t(logcpm),check.names=FALSE)
  names(ex) <- paste0("Gene__",names(ex))
  for(sig in names(gene_signatures)) {
    gg <- intersect(gene_signatures[[sig]],present_genes)
    if(length(gg)) ex[[paste0("Signature__",sig)]] <- colMeans(logcpm[gg,,drop=FALSE])
  }
  if(!all(final_prep$features %in% names(ex))) stop("External model features missing.")
  xx <- as.matrix(ex[,final_prep$features,drop=FALSE])
  xx <- sweep(sweep(xx,2,final_prep$mean,"-"),2,final_prep$sd,"/")
  if(any(!is.finite(xx))) stop("Non-finite external scaled features.")
  write.csv(cbind(sm,ex),out("06_GSE144735_External_Features.csv"),row.names=FALSE)
  shift <- do.call(rbind,lapply(final_prep$features,function(f) {
    a <- X_internal[,f];b <- as.matrix(ex)[,f]
    data.frame(feature=f,development_min=min(a),development_median=median(a),development_max=max(a),
               external_min=min(b),external_median=median(b),external_max=max(b),
               external_outside_development_range=sum(b<min(a)|b>max(a)))
  }))
  write.csv(shift,out("06_GSE144735_Feature_Range_Audit.csv"),row.names=FALSE)

  # Frozen-model external predictions and evaluation

  ep <- data.frame(sample=sm$sample,patient=sm$patient,group=sm$group,region=sm$region,
                   observed=as.integer(sm$group=="Tumor"),
                   ElasticNet=as.numeric(predict(final_elastic,newx=xx,s=final_tuning$lambda,type="response")),
                   RandomForest=as.numeric(predict(final_rf,data=as.data.frame(xx,check.names=FALSE))$predictions[,"Tumor"]),
                   XGBoost=as.numeric(predict(final_xgb,xgboost::xgb.DMatrix(xx))))
  models <- c("ElasticNet","RandomForest","XGBoost")
  if(any(!is.finite(as.matrix(ep[,models]))) || any(as.matrix(ep[,models])<0 | as.matrix(ep[,models])>1)) {
    stop("Invalid external probabilities.")
  }
  write.csv(ep,out("06_GSE144735_All_Region_Predictions.csv"),row.names=FALSE)
  primary <- ep[ep$region %in% c("N","T"),,drop=FALSE]
  border <- ep[ep$region %in% c("N","B"),,drop=FALSE]
  if(nrow(primary)!=12 || nrow(border)!=12 || any(table(primary$patient,primary$group)!=1)) {
    stop("Expected six independent matched core/normal pairs after QC.")
  }
  external_metrics <- function(d,label) do.call(rbind,lapply(models,function(m) evaluate_model(d$observed,d[[m]],m,label)))
  primary_perf <- external_metrics(primary,"GSE144735 core vs normal")
  border_perf <- external_metrics(border,"GSE144735 border vs normal (secondary)")
  write.csv(primary,out("06_GSE144735_External_Predictions.csv"),row.names=FALSE)
  write.csv(primary_perf,out("06_GSE144735_External_Performance.csv"),row.names=FALSE)
  write.csv(border_perf,out("06_GSE144735_Border_Performance.csv"),row.names=FALSE)
  combined <- rbind(internal_performance,primary_perf,border_perf)
  write.csv(combined,out("06_GSE144735_Internal_External_Performance.csv"),row.names=FALSE)
  pair_rows <- do.call(rbind,lapply(models,function(m) do.call(rbind,lapply(expected_patients,function(p) {
    z <- primary[primary$patient==p,,drop=FALSE];n <- z[[m]][z$observed==0];t <- z[[m]][z$observed==1]
    data.frame(model=m,patient=p,normal_probability=n,tumor_probability=t,
               difference=t-n,both_correct=n<.5 & t>=.5,tumor_above_normal=t>n)
  }))))
  write.csv(pair_rows,out("06_GSE144735_Paired_Predictions.csv"),row.names=FALSE)
  cat("\nINTERNAL AND EXTERNAL PERFORMANCE (threshold=0.5)\n");print(combined,row.names=FALSE)
  cat("\nEXTERNAL PATIENT PAIRS\n");print(pair_rows,row.names=FALSE)

  plot_roc <- function(d,title,filename) {
    rows <- lapply(models,function(m) {
      ro <- pROC::roc(d$observed,d[[m]],levels=c(0,1),direction="<",quiet=TRUE)
      data.frame(FPR=1-ro$specificities,TPR=ro$sensitivities,
                 Model=paste0(m," (AUC=",sprintf("%.3f",as.numeric(pROC::auc(ro))),")"))
    })
    z <- do.call(rbind,rows)
    g <- ggplot(z,aes(FPR,TPR,color=Model))+geom_path(linewidth=.8)+
      geom_abline(slope=1,intercept=0,linetype=2,color="grey60")+
      coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_classic()+
      labs(title=title,x="False positive rate",y="True positive rate")
    ggsave(out(filename),g,width=7,height=5,dpi=600,bg="white")
  }
  plot_roc(internal_pred,"GSE132465: patient-held-out predictions","06_Internal_ROC.png")
  plot_roc(primary,"GSE144735: independent core tumor vs normal","06_GSE144735_External_ROC.png")
  plot_roc(border,"GSE144735: border vs normal (secondary)","06_GSE144735_Border_ROC.png")
  provenance <- data.frame(file=c(annotation_file,count_file),md5=unname(tools::md5sum(c(annotation_file,count_file))),
                           downloaded_source=paste0(geo_base,c(annotation_name,counts_name)))
  write.csv(provenance,out("06_GSE144735_Input_Provenance.csv"),row.names=FALSE)
  writeLines(c(
    "Development: GSE132465, ten Korean matched patients; fixed marker panel; nested patient-level tuning.",
    "External: GSE144735 (KUL3), six Belgian patients; distinct patients but from the same publication.",
    "Sources: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735",
    "Primary external comparison: six core tumor and six matched normal specimens; all QC-passing cells.",
    "Secondary comparison: six border and the same six normal specimens; not an additional independent cohort.",
    "Publisher-processed raw UMI matrix used, not natural-log TPM; deposited cell selection precedes our QC.",
    "Genes detected in >=3 deposited cells; nFeature >200 and <6000; mitochondrial UMI percentage <20.",
    "Pseudobulk sums all retained cells per specimen; denominator is all genes retained by the fixed gene filter.",
    "Feature construction: log2(CPM+1) for fixed genes; signature scores are means of those log values.",
    "All means/SDs, variable feature masks, tuning and fitted models estimated using development data only.",
    "No external fitting, threshold selection, label-based gene selection, joint integration or batch correction.",
    "Probability threshold 0.5; AUC direction fixed; no flipping predictions when AUC <0.5.",
    "Six external patient pairs are a small cohort; point estimates are uncertain and not clinical validation.",
    "No confidence intervals or comparative significance tests computed; uncertainty must respect patient pairing.",
    "A similar assay does not guarantee successful model transfer; dataset and sampling effects can remain."
  ),out("06_GSE144735_Analysis_Notes.txt"))
  writeLines(capture.output(sessionInfo()),out("06_GSE144735_sessionInfo.txt"))
  cat("\nGSE144735 EXTERNAL WORKFLOW COMPLETED\nOutput folder: ",output_root,"\n",sep="")
})
