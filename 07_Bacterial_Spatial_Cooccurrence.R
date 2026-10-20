# STAGE 07: exploratory all-group bacterial spatial co-occurrence.
# Run AFTER completed stages 04 and 05. No upstream outputs are changed.
# Outputs and resumable per-test caches are confined to a new Stage 07 folder.
# Raw spatial P values are approximate rank adaptations, NOT exact Spearman tests.
# Host-depth-adjusted coefficients and tumor-normal contrasts are descriptive only.
local({
  S04 <- 'E:/CRC/Results/Results/04-Spatial_Analysis_of_Selected_Features'
  S05 <- 'E:/CRC/Results/Results/05-Bacteria_Host_Spatial_Association'
  RAW <- 'E:/CRC/Results/GSE334323_ST'
  REF <- 'E:/CRC/Results/Results/02_Clustering_Markers_Annotation/05_annotated.rds'
  OUT <- 'E:/CRC/Results/Results/07-Bacterial_Spatial_Cooccurrence'
  SECTIONS <- c('CRC_03_Normal','CRC_03_Tumor','CRC_08_Normal','CRC_08_Tumor')
  GSM <- setNames(c('GSM9785428','GSM9785429','GSM9785432','GSM9785433'),SECTIONS)
  GROUPS <- sort(c('Alistipes','Bacteroides_Phocaeicola','Butyricimonas',
    'Campylobacter','Coprobacter','Dialister','Enterobacteriaceae',
    'Fusobacterium','Helicobacter','Intestinimonas','Lachnospiraceae',
    'Leptotrichia','Neisseria','Parabacteroides','Parvimonas',
    'Peptostreptococcus','Porphyromonas','Prevotellaceae','Streptococcaceae'))
  CLASSES <- c(13L,8L,20L)
  MIN_N <- 30L; MIN_POS <- 5L; MIN_COUNTS <- 10L
  REV <- 'stage07-v1'
  pkgs <- c('Matrix','SeuratObject','SpatialPack','digest','ggplot2')
  for(p in pkgs) if(!requireNamespace(p,quietly=TRUE)) stop('Missing package: ',p)
  old <- options(width=190);on.exit(options(old),add=TRUE)
  dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(OUT)) stop('Cannot create Stage 07 folder')
  CACHE <- file.path(OUT,'Cache');dir.create(CACHE,showWarnings=FALSE)
  COMPLETE <- file.path(OUT,'Stage07_Complete.rds')
  if(file.exists(COMPLETE)) unlink(COMPLETE)
  used <- character(); input_hashes <- character()
  record <- function(path) {
    if(!file.exists(path)) stop('Missing input: ',path)
    canonical <- normalizePath(path,winslash='/',mustWork=TRUE)
    if(!canonical %in% names(input_hashes)) {
      input_hashes[canonical] <<- unname(tools::md5sum(canonical))
    } else if(!identical(unname(input_hashes[canonical]),unname(tools::md5sum(canonical)))) {
      stop('Source changed during analysis: ',canonical)
    }
    used <<- unique(c(used,canonical))
    path
  }
  read_csv <- function(path,required) {
    z <- read.csv(record(path),stringsAsFactors=FALSE,check.names=FALSE)
    if(!all(required %in% names(z))) stop('Missing columns: ',path)
    z
  }
  verify <- function(root,name) {
    z <- readRDS(record(file.path(root,name)))
    if(!identical(z$revision,'core-v1')) stop('Wrong upstream revision: ',name)
    paths <- file.path(root,z$files)
    for(p in paths) record(p)
    if(!identical(z$md5,unname(tools::md5sum(paths)))) stop('Upstream checksum failure: ',name)
  }
  verify(S04,'Core_Spatial_Complete.rds')
  verify(S05,'Core_Association_Complete.rds')
  b <- read_csv(file.path(S04,'06_Bacterial_signal_All_spots_MATCHED.csv'),
                c('section','barcode','X','Y','bacterial_total_raw'))
  taxa <- read_csv(file.path(S05,'03_Bacterial_Groups_All_QC_Spots.csv'),
                c('section','barcode','X','Y','bacterial_group','bacterial_count'))
  qc <- read_csv(file.path(S04,'04_QC_summary_MATCHED.csv'),
                c('section','host_features','median_host_UMI','median_host_genes'))
  id <- function(z) paste(z$section,z$barcode,sep='::')
  if(anyNA(b[,c('section','barcode','X','Y','bacterial_total_raw')]) ||
     anyDuplicated(id(b)) || !setequal(b$section,SECTIONS) ||
     !setequal(taxa$section,SECTIONS) || !setequal(taxa$bacterial_group,GROUPS) ||
     anyDuplicated(paste(id(taxa),taxa$bacterial_group,sep='::')) ||
     anyNA(taxa[,c('section','barcode','X','Y','bacterial_group','bacterial_count')]))
    stop('Invalid bacterial spot/group identifiers')
  if(any(!is.finite(taxa$bacterial_count)|taxa$bacterial_count<0|
         taxa$bacterial_count!=floor(taxa$bacterial_count))) stop('Invalid bacterial counts')
  if(any(!is.finite(as.matrix(b[,c('X','Y','bacterial_total_raw')])))) stop('Invalid spatial inputs')
  # Reference is used only to reproduce the Stage 04 host feature mask.
  cat('Loading reference gene names (no single-cell model is rerun)...\n')
  reference <- readRDS(record(REF))
  host_symbols <- rownames(reference[['RNA']]);rm(reference);invisible(gc())
  if(length(host_symbols)<1000L) stop('Reference host gene names unavailable')
  raw_candidates <- unique(c(
    list.files(RAW,recursive=TRUE,full.names=TRUE),
    list.files(file.path(S04,'Matched_RAW_Extracted'),recursive=TRUE,full.names=TRUE),
    list.files(file.path(S05,'Input_Cache'),recursive=TRUE,full.names=TRUE),
    list.files(file.path(OUT,'Raw_Cache'),recursive=TRUE,full.names=TRUE)))
  locate <- function(ss,kind,extract=TRUE) {
    suffix <- switch(kind,features='features\\.tsv(\\.gz)?$',
                     barcodes='barcodes\\.tsv(\\.gz)?$',matrix='matrix\\.mtx(\\.gz)?$')
    hits <- raw_candidates[(grepl(GSM[[ss]],basename(raw_candidates),fixed=TRUE)|
                             grepl(ss,basename(raw_candidates),fixed=TRUE)) &
                            grepl(suffix,basename(raw_candidates),ignore.case=TRUE)]
    if(!length(hits) && extract) {
      archive <- record(file.path(RAW,'GSE334323_RAW.tar'))
      entries <- utils::untar(archive,list=TRUE)
      entries <- entries[grepl(GSM[[ss]],entries,fixed=TRUE)]
      if(!length(entries)||any(grepl('(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))',entries)))
        stop('Missing or unsafe archive entries: ',ss)
      where <- file.path(OUT,'Raw_Cache');dir.create(where,showWarnings=FALSE)
      utils::untar(archive,files=entries,exdir=where)
      raw_candidates <<- unique(c(raw_candidates,list.files(where,recursive=TRUE,full.names=TRUE)))
      return(locate(ss,kind,FALSE))
    }
    if(!length(hits)) stop('Cannot locate ',ss,' ',kind)
    if(length(hits)>1L && length(unique(unname(tools::md5sum(hits))))>1L)
      stop('Conflicting raw copies: ',ss,' ',kind)
    record(hits[1])
  }
  tsv <- function(path) read.delim(path,header=FALSE,quote='',comment.char='',
                                 stringsAsFactors=FALSE,check.names=FALSE)
  datasets <- list();host_audit <- list()
  for(ss in SECTIONS) {
    cat('Aligning all QC spots and reconstructing host depth: ',ss,'\n',sep='')
    bs <- b[b$section==ss,,drop=FALSE]
    if(anyDuplicated(bs[,c('X','Y')])) stop('Duplicate coordinates: ',ss)
    if(nrow(bs)<MIN_N) stop('Too few QC spots')
    z <- taxa[taxa$section==ss,,drop=FALSE]
    g <- vapply(GROUPS,function(nm) {
      a <- z[z$bacterial_group==nm,,drop=FALSE]
      if(nrow(a)!=nrow(bs)||!setequal(a$barcode,bs$barcode)) stop('Incomplete QC grid: ',ss,' ',nm)
      a <- a[match(bs$barcode,a$barcode),,drop=FALSE]
      if(any(abs(a$X-bs$X)>1e-6|abs(a$Y-bs$Y)>1e-6)) stop('Coordinate mismatch')
      a$bacterial_count
    },numeric(nrow(bs)))
    dimnames(g) <- list(bs$barcode,GROUPS)
    if(any(rowSums(g)!=bs$bacterial_total_raw)) stop('Bacterial total disagreement')
    ft <- tsv(locate(ss,'features'));bars <- as.character(tsv(locate(ss,'barcodes'))[[1]])
    mm <- locate(ss,'matrix');con <- if(grepl('\\.gz$',mm)) gzfile(mm,'rt') else file(mm,'rt')
    mat <- tryCatch(Matrix::readMM(con),finally=close(con))
    if(nrow(mat)==length(bars)&&ncol(mat)==nrow(ft)) mat <- t(mat)
    mat <- as(mat,'CsparseMatrix')
    if(nrow(mat)!=nrow(ft)||ncol(mat)!=length(bars)||anyDuplicated(bars)||ncol(ft)<2L)
      stop('Raw matrix dimension mismatch')
    if(any(!is.finite(mat@x)|mat@x<0|mat@x!=floor(mat@x))) stop('Invalid raw matrix')
    fid <- as.character(ft[[1]]);sym <- as.character(ft[[2]])
    missing <- is.na(sym)|!nzchar(sym);sym[missing] <- fid[missing]
    host <- (grepl('^ENSG[0-9]+(\\.[0-9]+)?$',fid)|sym %in% host_symbols) &
      !grepl('_probe[0-9]+$',sym,ignore.case=TRUE) & !is.na(sym) & nzchar(sym)
    if(sum(host)<1000L) stop('Too few human host features')
    ix <- match(bs$barcode,bars);if(anyNA(ix)) stop('QC barcode missing from raw matrix')
    symbols <- unique(sym[host])
    collapse <- Matrix::sparseMatrix(i=match(sym[host],symbols),j=seq_len(sum(host)),
                                    x=1,dims=c(length(symbols),sum(host)))
    hc <- collapse %*% mat[host,ix,drop=FALSE]
    depth <- as.numeric(Matrix::colSums(hc));genes <- as.numeric(Matrix::colSums(hc>0))
    q <- qc[qc$section==ss,,drop=FALSE]
    if(nrow(q)!=1L||length(symbols)!=q$host_features||
       abs(median(depth)-q$median_host_UMI)>1e-6||abs(median(genes)-q$median_host_genes)>1e-6||
       any(depth<100|genes<50)) stop('Reconstructed host QC differs from Stage 04: ',ss)
    bs$host_counts <- depth;bs$host_genes <- genes
    # Rank residuals control for a 3-df smooth function of log host depth.
    # Only descriptive: no spatial P is assigned to these residual coefficients.
    design <- cbind(Intercept=1,splines::ns(log1p(depth),df=3))
    ranks <- apply(g,2,rank,ties.method='average')
    residuals <- stats::lm.fit(design,ranks)$residuals
    dimnames(residuals) <- dimnames(g)
    datasets[[ss]] <- list(spots=bs,counts=g,depth_rank_residuals=residuals)
    host_audit[[ss]] <- data.frame(Section=ss,QC_spots=nrow(bs),Host_features=length(symbols),
      Median_host_counts=median(depth),Median_host_genes=median(genes),
      Bacterial_total=sum(g),Stage04_agreement=TRUE)
    rm(mat,hc,collapse);invisible(gc())
  }
  safe_cor <- function(x,y,method='spearman') {
    if(length(unique(x))<2L||length(unique(y))<2L) return(NA_real_)
    suppressWarnings(stats::cor(x,y,method=method))
  }
  technical <- do.call(rbind,lapply(SECTIONS,function(ss) {
    a <- datasets[[ss]];g <- a$counts;h <- a$spots$host_counts
    signals <- cbind(Total_bacterial_signal=rowSums(g),g)
    do.call(rbind,lapply(colnames(signals),function(nm) {
      x <- signals[,nm]
      data.frame(Section=ss,Patient=sub('_(Normal|Tumor)$','',ss),
        Tissue=sub('^CRC_[0-9]+_','',ss),Group=nm,N_spots=length(x),
        Raw_total=sum(x),Raw_mean=mean(x),Raw_median=median(x),
        Positive_spots=sum(x>0),Detection_percent=100*mean(x>0),
        Host_depth_rho=safe_cor(x,h),
        Mean_spot_signal_per_10000_host=mean(10000*x/h),
        Pooled_signal_per_10000_host=10000*sum(x)/sum(h))
    }))
  }))
  paired <- do.call(rbind,lapply(c('CRC_03','CRC_08'),function(pt) {
    do.call(rbind,lapply(unique(technical$Group),function(g) {
      n <- technical[technical$Patient==pt&technical$Group==g&technical$Tissue=='Normal',]
      t <- technical[technical$Patient==pt&technical$Group==g&technical$Tissue=='Tumor',]
      data.frame(Patient=pt,Group=g,Raw_mean_change=t$Raw_mean-n$Raw_mean,
        Host_scaled_mean_change=t$Mean_spot_signal_per_10000_host-n$Mean_spot_signal_per_10000_host,
        Pooled_host_scaled_change=t$Pooled_signal_per_10000_host-n$Pooled_signal_per_10000_host,
        Detection_change_pp=t$Detection_percent-n$Detection_percent,
        Raw_and_host_scaled_same_sign=sign(t$Raw_mean-n$Raw_mean)==
          sign(t$Mean_spot_signal_per_10000_host-n$Mean_spot_signal_per_10000_host))
    }))
  }))
  pairs <- t(utils::combn(GROUPS,2));colnames(pairs) <- c('Group1','Group2')
  coefficients <- do.call(rbind,lapply(SECTIONS,function(ss) {
    a <- datasets[[ss]]
    do.call(rbind,lapply(seq_len(nrow(pairs)),function(i) {
      x <- a$counts[,pairs[i,1]];y <- a$counts[,pairs[i,2]]
      p1 <- x>0;p2 <- y>0;both <- sum(p1&p2);expected <- sum(p1)*sum(p2)/length(x)
      union <- sum(p1|p2)
      data.frame(Section=ss,Group1=pairs[i,1],Group2=pairs[i,2],N_spots=length(x),
        Counts1=sum(x),Counts2=sum(y),Positive1=sum(p1),Positive2=sum(p2),
        Eligible=length(x)>=MIN_N&&sum(x)>=MIN_COUNTS&&sum(y)>=MIN_COUNTS&&
          sum(p1)>=MIN_POS&&sum(p2)>=MIN_POS&&length(unique(x))>1&&length(unique(y))>1,
        Rho_raw=safe_cor(x,y),
        Rho_host_scaled=safe_cor(x/a$spots$host_counts,y/a$spots$host_counts),
        Partial_rank_depth_coefficient=if(length(unique(x))<2L||length(unique(y))<2L) NA_real_ else
          safe_cor(a$depth_rank_residuals[,pairs[i,1]],a$depth_rank_residuals[,pairs[i,2]],method='pearson'),
        Both_detected_spots=both,Expected_overlap_independence=expected,
        Overlap_ratio=if(expected>0) both/expected else NA_real_,
        Jaccard=if(union>0) both/union else NA_real_)
    }))
  }))
  if(nrow(coefficients)!=684L) stop('Expected 684 all-group comparisons')
  spatial_test <- function(x,y,coords,nclass,eligible) {
    ans <- list(P_spatial=NA_real_,ESS=NA_real_,DF=NA_real_,Status='Low_information',Message='')
    if(!eligible) return(ans)
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(SpatialPack::modified.ttest(
      rank(x,ties.method='average'),rank(y,ties.method='average'),coords,nclass=nclass),
      warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
    if(inherits(fit,'error')) {ans$Status <- 'Error';ans$Message <- conditionMessage(fit);return(ans)}
    if(length(fit$ESS)==1) ans$ESS <- as.numeric(fit$ESS)
    if(length(fit$dof)==1) ans$DF <- as.numeric(fit$dof)
    if(length(warnings)) {ans$Status <- 'Warning';ans$Message <- paste(unique(warnings),collapse=' | ');return(ans)}
    if(!is.finite(ans$ESS)||ans$ESS<=2||!is.finite(ans$DF)||ans$DF<=0||
       length(fit$p.value)!=1L||!is.finite(fit$p.value)||fit$p.value<0||fit$p.value>1) {
      ans$Status <- 'Invalid_spatial_estimate';return(ans)
    }
    ans$Status <- 'OK';ans$P_spatial <- as.numeric(fit$p.value);ans
  }
  cat('Spatial inference: 684 comparisons per setting; 2052 planned total.\n')
  cat('Each completed test is cached in Stage 07; rerun this script to resume.\n')
  spatial <- list();job <- 0L
  for(nclass in CLASSES) {
    block <- vector('list',nrow(coefficients))
    for(i in seq_len(nrow(coefficients))) {
      z <- coefficients[i,,drop=FALSE];a <- datasets[[z$Section]]
      x <- a$counts[,z$Group1];y <- a$counts[,z$Group2];xy <- as.matrix(a$spots[,c('X','Y')])
      signature <- digest::digest(list(REV,x,y,xy,nclass,MIN_N,MIN_POS,MIN_COUNTS,
        as.character(utils::packageVersion('SpatialPack'))),algo='sha256')
      cache_file <- file.path(CACHE,paste0(signature,'.rds'))
      saved <- if(file.exists(cache_file)) tryCatch(readRDS(cache_file),error=function(e)NULL) else NULL
      if(is.null(saved)||!identical(saved$signature,signature)||
         !identical(saved$result_hash,digest::digest(saved$result,algo='sha256'))) {
        result <- spatial_test(x,y,xy,nclass,z$Eligible)
        saved <- list(signature=signature,result=result,result_hash=digest::digest(result,algo='sha256'))
        tmp <- paste0(cache_file,'.writing');saveRDS(saved,tmp)
        if(!file.copy(tmp,cache_file,overwrite=TRUE)) stop('Cannot write Stage 07 cache')
        unlink(tmp)
      }
      block[[i]] <- cbind(z[,c('Section','Group1','Group2','N_spots','Rho_raw')],
                          Nclass=nclass,as.data.frame(saved$result,stringsAsFactors=FALSE))
      job <- job+1L
      if(job%%10L==0L) cat('Completed/cached tests: ',job,' / 2052\n',sep='')
    }
    ans <- do.call(rbind,block)
    ok <- ans$Status=='OK' & is.finite(ans$P_spatial)
    ans$FDR_spatial <- NA_real_
    ans$FDR_spatial[ok] <- stats::p.adjust(ans$P_spatial[ok],method='BH',n=684L)
    ans$Positive_spatial <- ok & ans$FDR_spatial<.05 & ans$Rho_raw>0
    ans$Negative_spatial <- ok & ans$FDR_spatial<.05 & ans$Rho_raw<0
    spatial[[as.character(nclass)]] <- ans
  }
  all_tests <- do.call(rbind,spatial);rownames(all_tests) <- NULL
  primary <- spatial[['13']]
  stability <- coefficients
  for(nclass in CLASSES) {
    z <- spatial[[as.character(nclass)]]
    stopifnot(identical(z[,c('Section','Group1','Group2')],coefficients[,c('Section','Group1','Group2')]))
    stability[[paste0('Positive_FDR_',nclass)]] <- z$Positive_spatial
    stability[[paste0('Negative_FDR_',nclass)]] <- z$Negative_spatial
  }
  stability$Positive_all_settings <- Reduce(`&`,stability[paste0('Positive_FDR_',CLASSES)])
  stability$Negative_all_settings <- Reduce(`&`,stability[paste0('Negative_FDR_',CLASSES)])
  # Descriptive recurrence counts only; no tumor-normal differential network test.
  recurrence <- do.call(rbind,lapply(seq_len(nrow(pairs)),function(i) {
    z <- stability[stability$Group1==pairs[i,1]&stability$Group2==pairs[i,2],]
    t <- z[grepl('_Tumor$',z$Section),];n <- z[grepl('_Normal$',z$Section),]
    data.frame(Group1=pairs[i,1],Group2=pairs[i,2],
      Tumor_positive_primary=sum(t$Positive_FDR_13),Normal_positive_primary=sum(n$Positive_FDR_13),
      Tumor_positive_all_settings=sum(t$Positive_all_settings),
      Normal_positive_all_settings=sum(n$Positive_all_settings),
      Tumor_depth_partial_positive=sum(is.finite(t$Partial_rank_depth_coefficient)&t$Partial_rank_depth_coefficient>0),
      Tumor_low_information=sum(!t$Eligible),Normal_low_information=sum(!n$Eligible))
  }))
  notes <- data.frame(Item=c('Design','Family','Primary method','Distance sensitivity','Eligibility',
      'Host mask','Host scaling','Depth adjustment','Overlap','Recurrence','Taxonomy','Limits','Reference'),
    Description=c(
      'Post hoc exploratory reanalysis motivated by earlier findings; four sections from two paired patients; all 19 deposited probe groups, including Bacteroides_Phocaeicola. This is not independent confirmation.',
      '171 unordered pairs x 4 sections = 684 comparisons per distance setting; BH across full 684 including unavailable tests. Separate adjustment for each setting.',
      'Raw count Spearman coefficient; average ranks passed to SpatialPack modified.ttest with original pixel coordinates, 13 distance classes. Approximate rank adaptation, not exact Spearman inference. Warnings/errors/invalid ESS or df yield missing P/FDR.',
      'Eight and twenty classes; settings are correlated sensitivity checks, not independent replication.',
      'At least 30 spots; each group >=10 total counts and >=5 positive spots, and variation in both groups. All QC spots including zeros are used; low information is not evidence of absence.',
      'Reproduce Stage 04: Ensembl human IDs or reference gene symbols, excluding custom probe suffixes; collapse host symbols and verify QC summaries.',
      'Per-spot counts / host counts x 10000 and pooled count/host count rates are sensitivity measures, not calibrated abundance or full capture correction. Shared denominators can induce correlations; host RNA also reflects biology.',
      'Regress each group average-rank vector on natural spline of log1p(host counts), df=3, within section; Pearson correlation between residuals is descriptive partial-rank association. No adjusted P or FDR claimed.',
      'Detection co-occurrence uses counts>0; observed/independence expected overlap and Jaccard are descriptive. Zeros, depth and spatial structure affect them; no Fisher test.',
      'Counts of positive associations in normal/tumor sections are descriptive, not differential association tests or validated communities.',
      'Deposited labels span taxonomic ranks. Probe cross-reactivity and specificity are not independently validated here. Shared spots do not establish species identity, physical contact, biofilms, cooperation or viability.',
      'No causal inference, validated biomarker, proof of novelty or patient-level significance from two patients. Depth sensitivity does not remove all technical effects.',
      'SpatialPack modified.ttest documentation: https://search.r-project.org/CRAN/refmans/SpatialPack/html/modified.ttest.html'
    ))
  outputs <- list('01_Host_Depth_Extraction_Audit'=do.call(rbind,host_audit),
    '02_Section_Depth_Sensitivity'=technical,'03_Paired_Depth_Sensitivity'=paired,
    '04_All19_Pair_Coefficients'=coefficients,'05_Primary_Spatial_Tests'=primary,
    '06_All_Spatial_Settings'=all_tests,'07_Pair_Stability'=stability,
    '08_Descriptive_Recurrence'=recurrence,'09_Methods_and_Limits'=notes)
  for(nm in names(outputs)) write.csv(outputs[[nm]],file.path(OUT,paste0(nm,'.csv')),row.names=FALSE,na='')
  # Bundle permits read-only checker to recompute coefficients and depth diagnostics.
  manifest <- data.frame(Path=used,MD5=unname(tools::md5sum(used)),stringsAsFactors=FALSE)
  if(anyNA(manifest$MD5)||!identical(manifest$MD5,unname(input_hashes[used])))
    stop('Source missing or changed during analysis; no completion record written')
  write.csv(manifest,file.path(OUT,'10_Input_Manifest.csv'),row.names=FALSE)
  saveRDS(list(revision=REV,groups=GROUPS,sections=SECTIONS,classes=CLASSES,
    min_n=MIN_N,min_pos=MIN_POS,min_counts=MIN_COUNTS,datasets=datasets,
    pairs=pairs,inputs=manifest),file.path(OUT,'11_Checked_Input_Bundle.rds'))
  # Complete, unselected pair heatmaps, with X for unavailable spatial inference.
  pdat <- coefficients
  pdat$Primary_status <- primary$Status
  for(kind in c('Rho_raw','Partial_rank_depth_coefficient')) {
    plot <- ggplot2::ggplot(pdat,ggplot2::aes(x=Group1,y=Group2,fill=.data[[kind]]))+
      ggplot2::geom_tile()+ggplot2::facet_wrap(~Section,ncol=2)+
      ggplot2::scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',
                                   limits=c(-1,1),na.value='grey80')+
      ggplot2::theme_bw(base_size=10)+
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=.5))+
      ggplot2::labs(x=NULL,y=NULL,title=paste('All 19 probe groups:',kind),
        subtitle=if(kind=='Rho_raw') 'X: no reportable primary spatial test; all coefficients shown' else
          'Descriptive depth-adjusted coefficients; no P values or FDR',fill='Coefficient')
    if(kind=='Rho_raw') plot <- plot+ggplot2::geom_text(
      data=pdat[pdat$Primary_status!='OK',],label='x',size=2,color='grey35')
    ggplot2::ggsave(file.path(OUT,paste0('12_',kind,'_Heatmap.png')),plot,
                   width=15,height=13,dpi=300,limitsize=FALSE)
  }
  writeLines(capture.output(sessionInfo()),file.path(OUT,'13_sessionInfo.txt'))
  files <- c(paste0(names(outputs),'.csv'),'10_Input_Manifest.csv','11_Checked_Input_Bundle.rds',
             '12_Rho_raw_Heatmap.png','12_Partial_rank_depth_coefficient_Heatmap.png','13_sessionInfo.txt')
  saveRDS(list(revision=REV,files=files,md5=unname(tools::md5sum(file.path(OUT,files))),
               completed=Sys.time()),COMPLETE)
  cat('\nStage 07 COMPLETE. New output folder: ',OUT,'\n',sep='')
  cat('Run 07_Checker.R and send its console output.\n')
})
