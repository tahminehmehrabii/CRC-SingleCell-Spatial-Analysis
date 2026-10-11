# Paired Bacterial Analysis
# Paired section summaries retain all 19 groups and existing spatial FDR families.
# Two patients support descriptive direction comparisons, not tumor/normal significance.

local({
  # Paths
  input_dir <- 'E:/CRC/Results/Results/08-Bacteria_Host_Spatial_Association'
  output_dir <- 'E:/CRC/Results/Results/09-Paired_Bacterial_CRC_Association'
  sections <- c('CRC_03_Normal','CRC_03_Tumor','CRC_08_Normal','CRC_08_Tumor')
  patients <- c('CRC_03','CRC_08')
  targets <- c('EPCAM','Signature__Epithelial','Myeloids')
  read_tab <- function(name,required) {
    path <- file.path(input_dir,name)
    if(!file.exists(path)) stop('Run Stage 08 first. Missing: ',path)
    d <- read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
    if(!all(required %in% names(d))) stop('Unexpected columns: ',name)
    d
  }
  key <- function(d) paste(d$section,d$barcode,sep='::')
  old <- options(width=230,digits=7,max.print=100000)
  on.exit(options(old),add=TRUE)
  # Check audited probe counts and complete zero-inclusive spot tables
  audit <- read_tab('02_Bacterial_Extraction_Check.csv',c('Section','Exact_total_agreement'))
  if(nrow(audit)!=4L || anyDuplicated(audit$Section) ||
     !setequal(audit$Section,sections) ||
     !all(tolower(as.character(audit$Exact_total_agreement)) %in% c('true','1')))
    stop('Bacterial extraction audit failed.')
  taxa <- read_tab('03_Bacterial_Groups_All_QC_Spots.csv',
                  c('section','barcode','bacterial_group','bacterial_count'))
  ids <- c('section','barcode','bacterial_group')
  if(anyNA(taxa[,ids]) || any(vapply(taxa[,ids],function(x) any(!nzchar(trimws(x))),logical(1))))
    stop('Missing/blank spot or group IDs.')
  groups <- sort(unique(taxa$bacterial_group))
  if(length(groups)!=19L || !setequal(taxa$section,sections)) stop('Expected 19 groups and four sections.')
  if(anyDuplicated(paste(key(taxa),taxa$bacterial_group,sep='::'))) stop('Duplicate group/spot rows.')
  v <- taxa$bacterial_count
  if(!is.numeric(v) || any(!is.finite(v) | v<0 | v!=floor(v))) stop('Invalid bacterial counts.')
  totals <- tapply(v,key(taxa),sum)
  for(ss in sections) {
    z <- taxa[taxa$section==ss,]
    reference <- unique(z$barcode)
    for(g in groups) if(!setequal(z$barcode[z$bacterial_group==g],reference))
      stop('Incomplete zero-inclusive table: ',ss,' / ',g)
  }
  # Section summaries. mean counts, rather than totals, compare unequal spot counts
  summaries <- do.call(rbind,lapply(sections,function(ss) do.call(rbind,lapply(groups,function(g) {
    z <- taxa[taxa$section==ss & taxa$bacterial_group==g,]
    x <- z$bacterial_count; total <- as.numeric(totals[key(z)])
    data.frame(Patient=sub('_(Normal|Tumor)$','',ss),Tissue=sub('^CRC_[0-9]+_','',ss),
      Section=ss,Bacterial_group=g,QC_spots=length(x),Total_probe_counts=sum(x),
      Mean_counts_per_QC_spot=mean(x),Median_counts_per_QC_spot=median(x),
      Positive_spots=sum(x>0),Detection_percent=100*mean(x>0),
      Fraction_of_all_bacterial_counts=if(sum(total)>0) sum(x)/sum(total) else NA_real_)
  }))))
  # Matched changes; no spot-level tumor-versus-normal p values
  paired <- do.call(rbind,lapply(patients,function(pt) do.call(rbind,lapply(groups,function(g) {
    z <- summaries[summaries$Patient==pt & summaries$Bacterial_group==g,]
    n <- z[z$Tissue=='Normal',]; t <- z[z$Tissue=='Tumor',]
    delta <- t$Mean_counts_per_QC_spot-n$Mean_counts_per_QC_spot
    data.frame(Patient=pt,Bacterial_group=g,
      Normal_mean_counts=n$Mean_counts_per_QC_spot,Tumor_mean_counts=t$Mean_counts_per_QC_spot,
      Mean_count_change=delta,
      Mean_count_ratio=if(n$Mean_counts_per_QC_spot>0) t$Mean_counts_per_QC_spot/n$Mean_counts_per_QC_spot else NA_real_,
      Direction=if(delta>0) 'Higher_in_tumor' else if(delta<0) 'Lower_in_tumor' else 'Unchanged',
      Normal_detection_percent=n$Detection_percent,Tumor_detection_percent=t$Detection_percent,
      Detection_change_percentage_points=t$Detection_percent-n$Detection_percent,
      Bacterial_fraction_change=t$Fraction_of_all_bacterial_counts-n$Fraction_of_all_bacterial_counts)
  }))))
  consistency <- do.call(rbind,lapply(groups,function(g) {
    z <- paired[paired$Bacterial_group==g,]
    data.frame(Bacterial_group=g,Patients=nrow(z),Higher_in_tumor=sum(z$Mean_count_change>0),
      Lower_in_tumor=sum(z$Mean_count_change<0),
      Same_direction_both_patients=all(z$Mean_count_change>0)||all(z$Mean_count_change<0),
      Median_paired_mean_count_change=median(z$Mean_count_change),
      Detection_increased_both=all(z$Detection_change_percentage_points>0),
      Detection_decreased_both=all(z$Detection_change_percentage_points<0))
  }))
  # Reuse all-group spatial tests, retaining their original adjusted p values
  req <- c('Family','Section','Bacterial_group','Host_target','Nclass','Rho','Status','P_spatial','FDR_spatial')
  primary <- read_tab('06_Primary_Spatial_Adjusted_Correlations.csv',req)
  sensitivity <- read_tab('07_Spatial_Sensitivity.csv',req)
  check_tests <- function(d,classes) {
    id <- paste(d$Family,d$Section,d$Bacterial_group,d$Host_target,d$Nclass,sep='::')
    if(anyDuplicated(id) || anyNA(d$Nclass) || !all(d$Nclass %in% classes)) stop('Invalid spatial test identifiers.')
    for(nm in c('P_spatial','FDR_spatial')) if(!is.numeric(d[[nm]]) ||
      any(d[[nm]]<0 | d[[nm]]>1,na.rm=TRUE)) stop('Invalid P/FDR values.')
    if(!is.numeric(d$Rho) || any(abs(d$Rho)>1+1e-8,na.rm=TRUE)) stop('Invalid rho.')
    z <- d[d$Family=='Taxon_SameSpot' & d$Host_target %in% targets,]
    expected <- expand.grid(Section=sections,Bacterial_group=groups,Host_target=targets,Nclass=classes,
                            stringsAsFactors=FALSE)
    k <- function(x) paste(x$Section,x$Bacterial_group,x$Host_target,x$Nclass,sep='::')
    if(!setequal(k(z),k(expected))) stop('Incomplete all-group spatial tests.')
    z$Spatial_FDR_significant <- z$Status=='OK' & is.finite(z$FDR_spatial) & z$FDR_spatial<.05
    z
  }
  spatial <- check_tests(primary,13L)
  sensitivity <- check_tests(sensitivity,c(8L,20L))
  spatial <- merge(spatial,consistency,by='Bacterial_group',all.x=TRUE,sort=FALSE)
  selected <- spatial[spatial$Same_direction_both_patients,]

  notes <- data.frame(Note=c(
    'Two patients, each with one normal and one tumor section. No patient-level significance is claimed.',
    'Mean probe counts and detection percentages include every retained QC spot, including zeros.',
    'Raw probe signal depends on capture efficiency, depth, probe design and tissue sampling; it is not absolute bacterial abundance.',
    'Bacterial fraction uses the total of all bacterial probe groups and is compositional, not host-library normalization.',
    'Direction consistency is descriptive with two pairs; selected groups are exploratory and not validated biomarkers.',
    'Spatial tests are imported from Stage 08: rank-based SpatialPack modified.ttest; primary 13 distance classes, sensitivity 8 and 20.',
    'Original all-group BH FDR values are retained; no correction restricted to selected groups is performed.',
    'Invalid/low-information spatial tests are not nonsignificant evidence: see Status and missing P/FDR.',
    'Spatial associations within sections do not establish CRC causality, viable bacteria, species identity or patient-level replication.',
    'No models or spatial tests are refitted. This script requires completed corrected Stage 08 outputs.'))
  # Export and print all manuscript tables
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(output_dir)) stop('Cannot create output directory.')
  results <- list('01_Section_Bacterial_Summary'=summaries,'02_Patient_Paired_Changes'=paired,
    '03_Between_Patient_Direction_Consistency'=consistency,
    '04_All19_Spatial_Host_Associations'=spatial,'05_Spatial_Sensitivity'=sensitivity,
    '06_Consistent_Groups_Exploratory_Spatial'=selected,'07_Methods_and_Limits'=notes)
  for(nm in names(results)) {
    write.csv(results[[nm]],file.path(output_dir,paste0(nm,'.csv')),row.names=FALSE,na='')
    cat('\n',strrep('=',80),'\n',nm,'\n',sep=''); print(results[[nm]],row.names=FALSE)
  }
  cat('\nGroups with the same raw-mean change direction in both patients: ',
      sum(consistency$Same_direction_both_patients),' / 19\n',sep='')
  cat('Results saved in: ',output_dir,'\n',sep='')
  cat('No tumor-normal significance claim is supported by two patients alone.\n')
})
