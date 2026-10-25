# 01. Configure input and output folders
input_root <- "E:/CRC/Results/Results/02_Clustering_Markers_Annotation"
output_root <- "E:/CRC/Results/Results/03_Cell-Type and Subtype Composition Analysis"

local({
  # Package versions, QC and reclustering can change numerical results/cluster IDs.
  # Uses the existing processed project files; raw GEO import is not rerun.
  # 02. Paths and output directory
  root <- input_root
  work <- output_root
  dir.create(work,recursive=TRUE,showWarnings=FALSE)
  if (!dir.exists(work)) stop("Cannot create output directory: ",work)
  # 03. Required packages
  for(p in c("Seurat","SeuratObject","ggplot2","patchwork"))
    if(!requireNamespace(p,quietly=TRUE)) stop("Install package: ",p)
  library(ggplot2)
  # 04. Input-loading helper
  read_obj <- function(name) {
    f <- file.path(root,name); if(!file.exists(f)) stop("Missing: ",f)
    o <- readRDS(f)
    if(!identical(rownames(o@meta.data),colnames(o))) stop("Metadata mismatch")
    o
  }
  # 05. Metadata and pairing validation
  check_meta <- function(o) {
    m <- o@meta.data
    need <- c("sample_id","patient_id","group","analysis_cell_type")
    if(!all(need %in% names(m)) || anyNA(m[,need])) stop("Metadata invalid")
    if(anyDuplicated(rownames(m)) || any(vapply(m[,need,drop=FALSE],function(z) any(!nzchar(trimws(as.character(z)))),logical(1)))) stop("Blank or duplicated metadata")
    patients <- sprintf("SMC%02d",1:10)
    samples <- as.vector(t(outer(patients,c("-N","-T"),paste0)))
    if(!setequal(m$sample_id,samples) || !setequal(m$patient_id,patients) ||
       !all(m$group %in% c("Normal","Tumor")) ||
       any(m$sample_id!=paste0(m$patient_id,ifelse(m$group=="Normal","-N","-T")))) stop("Pairing invalid")
    m
  }
  # 06. CSV and PNG saving helpers
  write_tab <- function(z,n) write.csv(z,file.path(work,n),row.names=FALSE,na="")
  save_plot <- function(p,n,w=7.09,h=5.5) {
    ggsave(file.path(work,paste0(n,".png")),plot=p,width=w,height=h,
           units="in",device="png",dpi=600,bg="white")
  }
  # 07. Broad composition statistics and plotting
  paired_composition <- function(d, levels, prefix) {
    d$cell_type <- factor(d$cell_type,levels=levels)
    tests <- do.call(rbind,lapply(levels,function(tp) {
      z <- d[as.character(d$cell_type)==tp,]
      nn <- z$percentage[z$group=="Normal"]; names(nn) <- z$patient[z$group=="Normal"]
      tt <- z$percentage[z$group=="Tumor"]; names(tt) <- z$patient[z$group=="Tumor"]
      pp <- intersect(names(nn),names(tt)); delta <- tt[pp]-nn[pp]
      pv <- if(all(delta==0)) 1 else wilcox.test(tt[pp],nn[pp],paired=TRUE,exact=FALSE,correct=TRUE)$p.value
      data.frame(cell_type=tp,n_pairs=length(pp),median_paired_change=median(delta),p_value=pv)
    }))
    tests$p_adj_BH <- p.adjust(tests$p_value,"BH")
    write_tab(d,paste0(prefix,"_Sample_percentages.csv")); write_tab(tests,paste0(prefix,"_Paired_tests.csv"))
    d$group <- factor(d$group,levels=c("Normal","Tumor"))
    p <- ggplot(d,aes(group,percentage)) + geom_line(aes(group=patient),color="#BDBDBD",linewidth=.3) +
      geom_boxplot(aes(fill=group),outlier.shape=NA,width=.45,alpha=.7) + geom_point(size=1.5) +
      facet_wrap(~cell_type,scales="free_y",ncol=3) +
      scale_fill_manual(values=c(Normal="#BDBDBD",Tumor="#208C91")) +
      labs(x=NULL,y="Cells in sample (%)") + theme_classic(base_family="Arial") + theme(legend.position="none")
    save_plot(p,paste0(prefix,"_Composition"),9,6)
    print(tests,row.names=FALSE)
  }
  # 08. Run broad cell-type composition
  o <- read_obj("05_annotated.rds"); m <- check_meta(o)
  levels <- c("B cells","Epithelial cells","Mast cells","Myeloids","Stromal cells","T cells")
  tab <- table(factor(m$sample_id,levels=as.vector(t(outer(sprintf("SMC%02d",1:10),c("-N","-T"),paste0)))),factor(m$analysis_cell_type,levels=levels))
  d <- as.data.frame(tab); names(d) <- c("sample","cell_type","n_cells")
  if(any(rowSums(tab)==0)) stop("Empty sample")
  d$percentage <- 100*d$n_cells/as.numeric(rowSums(tab)[as.character(d$sample)])
  d$patient <- sub("-[NT]$","",d$sample); d$group <- ifelse(grepl("-N$",d$sample),"Normal","Tumor")
  paired_composition(d,levels,"03_Broad")
  
})

local({
  
  # ============================================================
  # 09. Paths
  # ============================================================
  
  root <- output_root
  
  if (!dir.exists(root)) {
    stop("Directory does not exist: ", root)
  }
  
  out <- function(name) {
    file.path(root, name)
  }
  
  # ============================================================
  # 10. Packages
  # ============================================================
  
  packages <- c(
    "Seurat",
    "ggplot2",
    "patchwork"
  )
  
  for (p in packages) {
    
    if (!requireNamespace(p, quietly = TRUE)) {
      
      stop(
        "Install package first: ",
        p
      )
    }
  }
  
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  
  set.seed(42)
  
  # ============================================================
  # 11. Load annotated Seurat object
  # ============================================================
  
  input_file <- file.path(input_root, "05_annotated.rds")
  
  if (!file.exists(input_file)) {
    stop("Missing input file: ", input_file)
  }
  
  obj <- readRDS(input_file)
  
  meta <- obj@meta.data
  
  if (!identical(
    rownames(meta),
    colnames(obj)
  )) {
    
    stop(
      "Metadata row names do not match Seurat cell names."
    )
  }
  
  # ============================================================
  # 12. Required metadata
  # ============================================================
  
  required_columns <- c(
    "sample_id",
    "patient_id",
    "group",
    "analysis_cell_type",
    "GEO_CellSubtype"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(meta)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      "Missing metadata columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  # Convert to character
  
  meta$sample_id <- as.character(
    meta$sample_id
  )
  
  meta$patient_id <- as.character(
    meta$patient_id
  )
  
  meta$group <- as.character(
    meta$group
  )
  
  meta$analysis_cell_type <- as.character(
    meta$analysis_cell_type
  )
  
  meta$GEO_CellSubtype <- as.character(
    meta$GEO_CellSubtype
  )
  
  # ============================================================
  # 13. Validate patients and pairing
  # ============================================================
  
  patients <- sprintf(
    "SMC%02d",
    1:10
  )
  
  samples <- as.vector(
    t(
      outer(
        patients,
        c("-N", "-T"),
        paste0
      )
    )
  )
  
  if (!setequal(
    unique(meta$patient_id),
    patients
  )) {
    
    stop(
      "Unexpected patient IDs in metadata."
    )
  }
  
  if (!setequal(
    unique(meta$sample_id),
    samples
  )) {
    
    stop(
      "Expected 20 matched samples SMC01-N/T ... SMC10-N/T."
    )
  }
  
  if (!all(
    meta$group %in%
    c(
      "Normal",
      "Tumor"
    )
  )) {
    
    stop(
      "group must contain only Normal and Tumor."
    )
  }
  
  expected_sample <- paste0(
    meta$patient_id,
    ifelse(
      meta$group == "Normal",
      "-N",
      "-T"
    )
  )
  
  if (any(
    meta$sample_id != expected_sample
  )) {
    
    stop(
      "Patient / sample / group pairing is inconsistent."
    )
  }
  
  # ============================================================
  # 14. Exact GEO broad cell-type labels
  # ============================================================
  
  geo_broad <- c(
    "B cells",
    "Epithelial cells",
    "Mast cells",
    "Myeloids",
    "Stromal cells",
    "T cells"
  )
  
  display_broad <- c(
    "B",
    "Epithelial",
    "Mast",
    "Myeloid",
    "Stromal",
    "T"
  )
  
  names(display_broad) <- geo_broad
  
  # Check actual values present in metadata
  
  observed_broad <- sort(
    unique(
      meta$analysis_cell_type
    )
  )
  
  unexpected <- setdiff(
    observed_broad,
    geo_broad
  )
  
  if (length(unexpected) > 0) {
    
    stop(
      paste0(
        "Unexpected analysis_cell_type labels: ",
        paste(
          unexpected,
          collapse = ", "
        )
      )
    )
  }
  
  missing_broad <- setdiff(
    geo_broad,
    observed_broad
  )
  
  if (length(missing_broad) > 0) {
    
    stop(
      paste0(
        "These expected broad cell types are missing: ",
        paste(
          missing_broad,
          collapse = ", "
        )
      )
    )
  }
  
  # Add short plotting name
  
  meta$broad_type <- unname(
    display_broad[
      meta$analysis_cell_type
    ]
  )
  
  broad_order <- unname(
    display_broad
  )
  
  # ============================================================
  # 15. Clean subtype
  # ============================================================
  
  meta$subtype <- trimws(
    meta$GEO_CellSubtype
  )
  
  # Convert invalid labels to NA
  
  invalid_subtype <- (
    is.na(meta$subtype) |
      meta$subtype == "" |
      tolower(meta$subtype) %in%
      c(
        "unknown",
        "unassigned",
        "na",
        "n/a",
        "none"
      )
  )
  
  meta$subtype[
    invalid_subtype
  ] <- NA_character_
  
  # ============================================================
  # 16. Total cells in each sample
  # ============================================================
  
  sample_total <- table(
    factor(
      meta$sample_id,
      levels = samples
    )
  )
  
  if (any(
    sample_total == 0
  )) {
    
    stop(
      "One or more samples have zero cells."
    )
  }
  
  # ============================================================
  # 17. Broad cell-type composition
  # ============================================================
  
  broad_table <- table(
    
    factor(
      meta$sample_id,
      levels = samples
    ),
    
    factor(
      meta$broad_type,
      levels = broad_order
    )
    
  )
  
  broad_df <- as.data.frame(
    broad_table,
    stringsAsFactors = FALSE
  )
  
  colnames(
    broad_df
  ) <- c(
    "sample",
    "feature",
    "n_cells"
  )
  
  broad_df$sample <- as.character(
    broad_df$sample
  )
  
  broad_df$feature <- as.character(
    broad_df$feature
  )
  
  broad_df$percentage <-
    
    100 *
    broad_df$n_cells /
    as.numeric(
      sample_total[
        broad_df$sample
      ]
    )
  
  broad_df$patient <- sub(
    "-[NT]$",
    "",
    broad_df$sample
  )
  
  broad_df$group <- ifelse(
    grepl(
      "-N$",
      broad_df$sample
    ),
    "Normal",
    "Tumor"
  )
  
  broad_df$parent <- "Broad"
  
  broad_df$level <- "Broad"
  
  # ============================================================
  # 18. Determine subtype hierarchy
  # ============================================================
  
  sub_meta <- meta[
    !is.na(meta$subtype),
    ,
    drop = FALSE
  ]
  
  subtype_hierarchy <- unique(
    sub_meta[
      ,
      c(
        "broad_type",
        "subtype"
      )
    ]
  )
  
  subtype_hierarchy <- subtype_hierarchy[
    order(
      match(
        subtype_hierarchy$broad_type,
        broad_order
      ),
      subtype_hierarchy$subtype
    ),
    ,
    drop = FALSE
  ]
  
  rownames(
    subtype_hierarchy
  ) <- NULL
  
  if (nrow(subtype_hierarchy) == 0) {
    
    stop(
      "No valid GEO cell subtypes were found."
    )
  }
  
  # ============================================================
  # 19. Calculate subtype percentages
  # ============================================================
  
  subtype_list <- list()
  
  for (
    i in seq_len(
      nrow(
        subtype_hierarchy
      )
    )
  ) {
    
    parent_i <-
      subtype_hierarchy$broad_type[i]
    
    subtype_i <-
      subtype_hierarchy$subtype[i]
    
    cells_i <- meta[
      meta$broad_type == parent_i &
        !is.na(meta$subtype) &
        meta$subtype == subtype_i,
      ,
      drop = FALSE
    ]
    
    counts_i <- table(
      factor(
        cells_i$sample_id,
        levels = samples
      )
    )
    
    subtype_list[[i]] <- data.frame(
      
      sample = samples,
      
      feature = subtype_i,
      
      n_cells =
        as.integer(
          counts_i
        ),
      
      percentage =
        100 *
        as.integer(
          counts_i
        ) /
        as.numeric(
          sample_total[
            samples
          ]
        ),
      
      patient =
        sub(
          "-[NT]$",
          "",
          samples
        ),
      
      group =
        ifelse(
          grepl(
            "-N$",
            samples
          ),
          "Normal",
          "Tumor"
        ),
      
      parent =
        parent_i,
      
      level =
        "Subtype",
      
      stringsAsFactors =
        FALSE
      
    )
    
  }
  
  subtype_df <- do.call(
    rbind,
    subtype_list
  )
  
  rownames(
    subtype_df
  ) <- NULL
  
  # ============================================================
  # 20. Combine broad + subtype
  # ============================================================
  
  composition <- rbind(
    
    broad_df[
      ,
      c(
        "sample",
        "feature",
        "n_cells",
        "percentage",
        "patient",
        "group",
        "parent",
        "level"
      )
    ],
    
    subtype_df[
      ,
      c(
        "sample",
        "feature",
        "n_cells",
        "percentage",
        "patient",
        "group",
        "parent",
        "level"
      )
    ]
    
  )
  
  composition$group <- factor(
    composition$group,
    levels = c(
      "Normal",
      "Tumor"
    )
  )
  
  composition$patient <- factor(
    composition$patient,
    levels = patients
  )
  
  # ============================================================
  # 21. Statistical function
  # ============================================================
  
  paired_test <- function(z) {
    
    normal <- z$percentage[
      z$group == "Normal"
    ]
    
    names(normal) <- as.character(
      z$patient[
        z$group == "Normal"
      ]
    )
    
    tumor <- z$percentage[
      z$group == "Tumor"
    ]
    
    names(tumor) <- as.character(
      z$patient[
        z$group == "Tumor"
      ]
    )
    
    paired_patients <- intersect(
      names(normal),
      names(tumor)
    )
    
    normal <- normal[
      paired_patients
    ]
    
    tumor <- tumor[
      paired_patients
    ]
    
    delta <- tumor - normal
    
    if (
      length(paired_patients) < 2
    ) {
      
      p_value <- NA_real_
      
    } else if (
      all(delta == 0)
    ) {
      
      p_value <- 1
      
    } else {
      
      p_value <- suppressWarnings(
        
        wilcox.test(
          
          tumor,
          
          normal,
          
          paired = TRUE,
          
          exact = FALSE,
          
          correct = TRUE
          
        )$p.value
        
      )
      
    }
    
    data.frame(
      
      n_pairs =
        length(
          paired_patients
        ),
      
      median_normal =
        median(
          normal,
          na.rm = TRUE
        ),
      
      median_tumor =
        median(
          tumor,
          na.rm = TRUE
        ),
      
      median_paired_change =
        median(
          delta,
          na.rm = TRUE
        ),
      
      p_value =
        p_value
      
    )
    
  }
  
  # ============================================================
  # 22. Run test for all populations
  # ============================================================
  
  population_keys <- unique(
    
    composition[
      ,
      c(
        "level",
        "parent",
        "feature"
      )
    ]
    
  )
  
  test_results <- vector(
    "list",
    nrow(
      population_keys
    )
  )
  
  for (
    i in seq_len(
      nrow(
        population_keys
      )
    )
  ) {
    
    current_level <-
      population_keys$level[i]
    
    current_parent <-
      population_keys$parent[i]
    
    current_feature <-
      population_keys$feature[i]
    
    z <- composition[
      
      composition$level ==
        current_level &
        
        composition$parent ==
        current_parent &
        
        composition$feature ==
        current_feature,
      
      ,
      drop = FALSE
      
    ]
    
    current_test <- paired_test(
      z
    )
    
    test_results[[i]] <- data.frame(
      
      level =
        current_level,
      
      parent =
        current_parent,
      
      feature =
        current_feature,
      
      current_test,
      
      stringsAsFactors =
        FALSE
      
    )
    
  }
  
  tests <- do.call(
    rbind,
    test_results
  )
  
  rownames(
    tests
  ) <- NULL
  
  # ============================================================
  # 23. Benjamini-Hochberg FDR
  # ============================================================
  
  tests$FDR <- p.adjust(
    tests$p_value,
    method = "BH"
  )
  
  
  # Exact conditional signed-rank sensitivity: enumerate all signs of
  # NONZERO paired differences, with average ranks for tied magnitudes.
  # This still assumes independent patients and sign symmetry under the null.
  exact_signed_rank <- function(delta) {
    stopifnot(all(is.finite(delta)))
    delta <- delta[delta != 0]
    n <- length(delta)
    if (!n) return(1)
    if (n > 20) stop("Enumeration limited to 20 nonzero pairs")
    ranks <- rank(abs(delta), ties.method = "average")
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n),
                                   KEEP.OUT.ATTRS = FALSE))
    distribution <- drop(signs %*% ranks)
    observed <- sum(sign(delta) * ranks)
    mean(abs(distribution) >= abs(observed) - 1e-10)
  }
  sensitivity <- tests[, c("level", "parent", "feature", "n_pairs",
                           "median_normal", "median_tumor",
                           "median_paired_change", "p_value", "FDR")]
  sensitivity$n_nonzero <- integer(nrow(tests))
  sensitivity$n_zero <- integer(nrow(tests))
  sensitivity$tied_nonzero_magnitudes <- logical(nrow(tests))
  sensitivity$p_exact <- numeric(nrow(tests))
  for (i in seq_len(nrow(tests))) {
    z <- composition[composition$level == tests$level[i] &
                       composition$parent == tests$parent[i] &
                       composition$feature == tests$feature[i], ]
    nn <- z[z$group == "Normal", ]
    tt <- z[z$group == "Tumor", ]
    stopifnot(nrow(nn) == 10, nrow(tt) == 10,
              !anyDuplicated(nn$patient), !anyDuplicated(tt$patient),
              setequal(as.character(nn$patient), patients),
              setequal(as.character(tt$patient), patients))
    nn <- nn[match(patients, as.character(nn$patient)), ]
    tt <- tt[match(patients, as.character(tt$patient)), ]
    delta <- tt$percentage - nn$percentage
    nonzero <- delta[delta != 0]
    sensitivity$n_nonzero[i] <- length(nonzero)
    sensitivity$n_zero[i] <- sum(delta == 0)
    sensitivity$tied_nonzero_magnitudes[i] <- anyDuplicated(abs(nonzero)) > 0
    sensitivity$p_exact[i] <- exact_signed_rank(delta)
  }
  stopifnot(nrow(sensitivity) == 41, all(sensitivity$n_pairs == 10))
  sensitivity$FDR_exact <- p.adjust(sensitivity$p_exact, method = "BH")
  sensitivity$significant_primary <- sensitivity$FDR < 0.05
  sensitivity$significant_exact <- sensitivity$FDR_exact < 0.05
  sensitivity$status_changed <- sensitivity$significant_primary !=
    sensitivity$significant_exact
  write.csv(sensitivity,
            out("05_Composition_Exact_Sensitivity.csv"), row.names = FALSE)
  cat("\nEXACT SIGNED-RANK SENSITIVITY (joint BH over 41 tests)\n")
  print(sensitivity, row.names = FALSE)
  cat("Primary significant:", sum(sensitivity$significant_primary),
      "| Exact significant:", sum(sensitivity$significant_exact),
      "| Status changes:", sum(sensitivity$status_changed), "\n")
  if (any(sensitivity$status_changed))
    print(sensitivity[sensitivity$status_changed, ], row.names = FALSE)
  
  # ============================================================
  # 24. p-value formatting
  # ============================================================
  
  format_stat <- function(x) {
    
    if (is.na(x)) {
      
      return(
        "NA"
      )
      
    }
    
    if (x < 0.001) {
      
      return(
        formatC(
          x,
          format = "e",
          digits = 2
        )
      )
      
    }
    
    formatC(
      x,
      format = "f",
      digits = 3
    )
    
  }
  
  tests$label <- paste0(
    
    "p = ",
    vapply(
      tests$p_value,
      format_stat,
      character(1)
    ),
    
    "\nFDR = ",
    
    vapply(
      tests$FDR,
      format_stat,
      character(1)
    )
    
  )
  
  # ============================================================
  # 25. Save tables
  # ============================================================
  
  write.csv(
    
    composition,
    
    out(
      "04_CellType_Subcluster_Sample_Percentages.csv"
    ),
    
    row.names = FALSE,
    
    na = ""
    
  )
  
  write.csv(
    
    tests,
    
    out(
      "04_CellType_Subcluster_Paired_Wilcoxon_FDR.csv"
    ),
    
    row.names = FALSE,
    
    na = ""
    
  )
  
  write.csv(
    
    subtype_hierarchy,
    
    out(
      "04_CellType_Subcluster_Hierarchy.csv"
    ),
    
    row.names = FALSE,
    
    na = ""
    
  )
  
  # ============================================================
  # 26. Colors
  # ============================================================
  
  normal_color <- "#BDBDBD"
  
  tumor_color <- "#208C91"
  
  # ============================================================
  # 27. Save plot function
  # ============================================================
  
  save_plot <- function(plot_object, filename, width, height) {
    ggsave(filename=out(paste0(filename,".png")),plot=plot_object,
           width=width,height=height,units="in",device="png",dpi=600,bg="white")
  }
  
  # ============================================================
  # 28. Plot function
  # ============================================================
  
  make_composition_plot <- function(
    plot_data,
    plot_tests,
    feature_order,
    title_text,
    columns = 6
  ) {
    
    if (length(feature_order) == 0) {
      
      return(NULL)
      
    }
    
    plot_data <- plot_data[
      plot_data$feature %in%
        feature_order,
      ,
      drop = FALSE
    ]
    
    plot_tests <- plot_tests[
      plot_tests$feature %in%
        feature_order,
      ,
      drop = FALSE
    ]
    
    plot_data$feature <- factor(
      plot_data$feature,
      levels = feature_order
    )
    
    plot_tests$feature <- factor(
      plot_tests$feature,
      levels = feature_order
    )
    
    # Maximum Y for every facet
    
    max_y <- aggregate(
      
      percentage ~ feature,
      
      data =
        plot_data,
      
      FUN =
        function(x) {
          max(
            x,
            na.rm = TRUE
          )
        }
      
    )
    
    colnames(
      max_y
    )[2] <- "max_y"
    
    plot_tests <- merge(
      
      plot_tests,
      
      max_y,
      
      by =
        "feature",
      
      all.x =
        TRUE,
      
      sort =
        FALSE
      
    )
    
    # Position statistical label above boxplots
    
    plot_tests$label_y <- ifelse(
      
      plot_tests$max_y > 0,
      
      plot_tests$max_y *
        1.18,
      
      0.25
      
    )
    
    p <- ggplot(
      
      plot_data,
      
      aes(
        x = group,
        y = percentage
      )
      
    ) +
      
      # Connect matched patient samples
      
      geom_line(
        
        aes(
          group = patient
        ),
        
        color =
          "#BDBDBD",
        
        linewidth =
          0.35,
        
        alpha =
          0.8
        
      ) +
      
      # Boxplot
      
      geom_boxplot(
        
        aes(
          fill = group
        ),
        
        width =
          0.48,
        
        linewidth =
          0.45,
        
        outlier.shape =
          NA,
        
        alpha =
          0.75
        
      ) +
      
      # Patient points
      
      geom_point(
        
        aes(
          fill = group
        ),
        
        shape =
          21,
        
        size =
          1.8,
        
        stroke =
          0.3,
        
        color =
          "black"
        
      ) +
      
      # Statistical annotation
      
      geom_text(
        
        data =
          plot_tests,
        
        aes(
          x = 1.5,
          y = label_y,
          label = label
        ),
        
        inherit.aes =
          FALSE,
        
        size =
          2.5,
        
        lineheight =
          0.9,
        
        vjust =
          0
        
      ) +
      
      facet_wrap(
        
        ~feature,
        
        scales =
          "free_y",
        
        ncol =
          columns
        
      ) +
      
      scale_fill_manual(
        
        values = c(
          
          Normal =
            normal_color,
          
          Tumor =
            tumor_color
          
        )
        
      ) +
      
      scale_y_continuous(
        
        expand =
          expansion(
            mult = c(
              0.03,
              0.35
            )
          )
        
      ) +
      
      labs(
        
        x =
          NULL,
        
        y =
          "Cells in sample (%)",
        
        title =
          title_text
        
      ) +
      
      theme_classic(
        base_family = "Arial",
        base_size = 9
      ) +
      
      theme(
        
        legend.position =
          "none",
        
        axis.text.x =
          element_text(
            color = "black",
            size = 8
          ),
        
        axis.text.y =
          element_text(
            color = "black",
            size = 8
          ),
        
        axis.title.y =
          element_text(
            size = 9
          ),
        
        strip.background =
          element_blank(),
        
        strip.text =
          element_text(
            face = "bold",
            size = 9
          ),
        
        plot.title =
          element_text(
            face = "bold",
            size = 10,
            hjust = 0
          ),
        
        panel.spacing.x =
          grid::unit(
            0.6,
            "lines"
          ),
        
        panel.spacing.y =
          grid::unit(
            0.8,
            "lines"
          )
        
      )
    
    return(p)
    
  }
  
  # ============================================================
  # 29. Broad cell-type plot
  # ============================================================
  
  broad_data <- composition[
    composition$level == "Broad",
    ,
    drop = FALSE
  ]
  
  broad_tests <- tests[
    tests$level == "Broad",
    ,
    drop = FALSE
  ]
  
  p_broad <- make_composition_plot(
    
    plot_data =
      broad_data,
    
    plot_tests =
      broad_tests,
    
    feature_order =
      broad_order,
    
    title_text =
      "Broad cell types",
    
    columns =
      6
    
  )
  
  # ============================================================
  # 30. Subcluster plots
  # ============================================================
  
  subtype_plots <- list()
  
  for (
    parent_i in broad_order
  ) {
    
    subtype_order <- subtype_hierarchy$subtype[
      subtype_hierarchy$broad_type ==
        parent_i
    ]
    
    subtype_order <- unique(
      subtype_order
    )
    
    if (
      length(
        subtype_order
      ) == 0
    ) {
      
      next
      
    }
    
    current_data <- composition[
      
      composition$level ==
        "Subtype" &
        
        composition$parent ==
        parent_i,
      
      ,
      drop = FALSE
      
    ]
    
    current_tests <- tests[
      
      tests$level ==
        "Subtype" &
        
        tests$parent ==
        parent_i,
      
      ,
      drop = FALSE
      
    ]
    
    subtype_plots[[parent_i]] <-
      
      make_composition_plot(
        
        plot_data =
          current_data,
        
        plot_tests =
          current_tests,
        
        feature_order =
          subtype_order,
        
        title_text =
          paste0(
            parent_i,
            " subclusters"
          ),
        
        columns =
          min(
            6,
            length(
              subtype_order
            )
          )
        
      )
    
  }
  
  # ============================================================
  # 31. Save broad plot separately
  # ============================================================
  
  save_plot(
    
    plot_object =
      p_broad,
    
    filename =
      "04_Broad_CellTypes_Composition",
    
    width =
      14,
    
    height =
      3.4
    
  )
  
  # ============================================================
  # 32. Save every subtype group separately
  # ============================================================
  
  for (
    parent_i in names(
      subtype_plots
    )
  ) {
    
    n_sub <- length(
      
      unique(
        
        subtype_hierarchy$subtype[
          subtype_hierarchy$broad_type ==
            parent_i
        ]
        
      )
      
    )
    
    n_rows <- ceiling(
      n_sub / 6
    )
    
    safe_name <- gsub(
      "[^A-Za-z0-9]",
      "_",
      parent_i
    )
    
    save_plot(
      
      plot_object =
        subtype_plots[[parent_i]],
      
      filename =
        paste0(
          "04_",
          safe_name,
          "_Subclusters_Composition"
        ),
      
      width =
        14,
      
      height =
        max(
          3.4,
          3.0 * n_rows
        )
      
    )
    
  }
  
  # ============================================================
  # 33. Construct complete multi-row panel
  # ============================================================
  
  all_plots <- c(
    
    list(
      p_broad
    ),
    
    subtype_plots
    
  )
  
  # Calculate appropriate height
  
  number_subtype_rows <- 0
  
  for (
    parent_i in names(
      subtype_plots
    )
  ) {
    
    n_sub <- length(
      
      unique(
        
        subtype_hierarchy$subtype[
          subtype_hierarchy$broad_type ==
            parent_i
        ]
        
      )
      
    )
    
    number_subtype_rows <-
      number_subtype_rows +
      ceiling(
        n_sub / 6
      )
    
  }
  
  final_height <-
    3.5 +
    (
      3.0 *
        number_subtype_rows
    )
  
  final_plot <- wrap_plots(
    
    all_plots,
    
    ncol =
      1
    
  ) +
    
    plot_annotation(
      
      title =
        paste0(
          "Cell-type and subcluster composition ",
          "in matched CRC Normal and Tumor samples"
        ),
      
      subtitle =
        paste0(
          "Paired Wilcoxon signed-rank test; ",
          "Benjamini-Hochberg FDR correction"
        ),
      
      theme =
        theme(
          
          plot.title =
            element_text(
              family = "Arial",
              face = "bold",
              size = 13
            ),
          
          plot.subtitle =
            element_text(
              family = "Arial",
              size = 9
            )
          
        )
      
    )
  
  # ============================================================
  # 34. Save final figure
  # ============================================================
  
  save_plot(
    
    plot_object =
      final_plot,
    
    filename =
      "04_CellType_Subcluster_Composition",
    
    width =
      14,
    
    height =
      final_height
    
  )
  
  # ============================================================
  # 35. Additional QC tables
  # ============================================================
  
  # Number of cells in every broad cell type
  
  broad_cell_counts <- as.data.frame(
    
    table(
      meta$broad_type
    )
    
  )
  
  colnames(
    broad_cell_counts
  ) <- c(
    "Broad_Cell_Type",
    "Cells"
  )
  
  write.csv(
    
    broad_cell_counts,
    
    out(
      "04_Broad_CellType_Total_Cell_Counts.csv"
    ),
    
    row.names = FALSE
    
  )
  
  # Number of cells in each subtype
  
  subtype_cell_counts <- as.data.frame(
    
    table(
      meta$broad_type,
      meta$subtype,
      useNA = "no"
    )
    
  )
  
  colnames(
    subtype_cell_counts
  ) <- c(
    "Broad_Cell_Type",
    "Subtype",
    "Cells"
  )
  
  subtype_cell_counts <- subtype_cell_counts[
    subtype_cell_counts$Cells > 0,
    ,
    drop = FALSE
  ]
  
  write.csv(
    
    subtype_cell_counts,
    
    out(
      "04_Subcluster_Total_Cell_Counts.csv"
    ),
    
    row.names = FALSE
    
  )
  
  # ============================================================
  # 36. Console summary
  # ============================================================
  
  cat(
    "\n=============================================\n"
  )
  
  cat(
    "CELL-TYPE AND SUBTYPE COMPOSITION ANALYSIS COMPLETED\n"
  )
  
  cat(
    "=============================================\n\n"
  )
  
  cat(
    "Input object:",
    input_file,
    "\n"
  )
  
  cat(
    "Total cells:",
    nrow(meta),
    "\n"
  )
  
  cat(
    "Patients:",
    length(
      unique(
        meta$patient_id
      )
    ),
    "\n"
  )
  
  cat(
    "Samples:",
    length(
      unique(
        meta$sample_id
      )
    ),
    "\n"
  )
  
  cat(
    "Broad cell types:",
    length(
      broad_order
    ),
    "\n"
  )
  
  cat(
    "Subclusters:",
    nrow(
      subtype_hierarchy
    ),
    "\n\n"
  )
  
  cat(
    "Main figure:\n",
    out(
      "04_CellType_Subcluster_Composition.png"
    ),
    "\n\n"
  )
  
  cat(
    "Statistics:\n",
    out(
      "04_CellType_Subcluster_Paired_Wilcoxon_FDR.csv"
    ),
    "\n\n"
  )
  
  cat(
    "Percentages:\n",
    out(
      "04_CellType_Subcluster_Sample_Percentages.csv"
    ),
    "\n\n"
  )
  
  cat(
    "Hierarchy:\n",
    out(
      "04_CellType_Subcluster_Hierarchy.csv"
    ),
    "\n"
  )
  
})

invisible(
  gc()
)
