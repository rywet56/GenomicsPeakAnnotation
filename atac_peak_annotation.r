library(GenomicRanges)
library(annotate)
library(org.At.tair.db)

############################################################
# Process GTF file and Add introns to GTF file
############################################################

# # Get GTF file
# ##############
# # Read in the raw gtf file and process it so that the columns make sense
# path <- "/Users/manuel/Downloads/Araport11_GTF_genes_transposons.Mar202021.gtf"
# gtf_df <- read.table(file = path)
# gtf_df$V9 <- NULL; gtf_df$V11 <- NULL; gtf_df$V12 <- NULL; gtf_df$V14 <- NULL
# col_names <- c("chr", "src", "feature", "start", "end", "score", "strand", "readingframe", "transcript_id", "gene_id")
# colnames(gtf_df) <- col_names
# ##############
#
# # Add introns to GTF DataFrame
# ##############
# gtf_introns_df <- rbind(gtf_df, intron_df)
# path <- '/Volumes/work_storg/xiaocai_project/domain_atac_seq/gtf_introns.csv'
# write.csv(x = gtf_introns_df, file = path)
# ##############


############################################################
# All functions for annotation
############################################################

find_tair <- function(tair_id, mtx){
  # used to find whether any annotated column contains a gene and if so, which
  # column(s) and row(s) are annotated with that gene
  tair_cols <- grep(pattern = "tair", x = colnames(mtx))
  tair_cols <- colnames(mtx)[tair_cols]
  cols_idx_dic <- NULL
  for(colname in tair_cols){
    found <- mtx[, colname] %in% tair_id
    if(sum(found) > 0){
      found_idxs <- which(found == TRUE)
      cols_idx_dic[[colname]] <- found_idxs
    }
  }
  return(cols_idx_dic)
}

get_introns <- function(gtf_df){
  transcripts <- unique(gtf_df$transcript_id)
  intron_df <- NULL
  for(transcript in transcripts){
    d <- gtf_df[gtf_df$transcript_id == transcript,]
    d <- d[d$feature == 'exon',]
    no_introns <- (nrow(d)-1)
    # test whether there are any exons and if there are at least one intron
    if(dim(d)[1] > 0 & no_introns > 0){
      row <- as.character(d[1,])
      intron_mtx_temp <- NULL
      for(i in 1:no_introns){
        # create cp of dummy row
        row_cp <- row
        # fill row with intron information
        row_cp[3] <- 'intron'
        row_cp[4] <- d[i, "end"]
        row_cp[5] <- d[i+1, "start"]
        # add row to matrix of introns for transcript
        intron_mtx_temp <- rbind(intron_mtx_temp, row_cp)
      }
      intron_df <- rbind(intron_df, intron_mtx_temp)
    }
  }
  rownames(intron_df) <- (nrow(gtf_df) + 1):(nrow(gtf_df) + nrow(intron_df))
  intron_df <- as.data.frame(intron_df)
  colnames(intron_df) <- colnames(gtf_df)

  return(intron_df)
}

get_genes <- function(gtf_df){
  # select genomic features from GTF
  genes <- gtf_df$feature == 'gene' | gtf_df$feature == 'transposable_element_gene'
  feature_df <- gtf_df[genes, c('chr', 'start', 'end', 'strand', 'gene_id')]
  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_promoters <- function(gtf_df, prom_region=c(-1000, 200)){

  # define TSSs
  gtf_genes <- gtf_df[gtf_df$feature == 'gene', ]
  pos_strand <- gtf_genes$strand == "+"
  tss <- rep(0, length(pos_strand))
  tss[pos_strand] <- (gtf_genes$start)[pos_strand]
  tss[!pos_strand] <- (gtf_genes$end)[!pos_strand]

  # define promoter regions
  ps <- rep(0, length(tss)); pe <- rep(0, length(tss))
  ps[pos_strand] <- tss[pos_strand] + prom_region[1]
  pe[pos_strand] <- tss[pos_strand] + prom_region[2]
  ps[!pos_strand] <- tss[!pos_strand] - prom_region[2]
  pe[!pos_strand] <- tss[!pos_strand] - prom_region[1]
  prom_df <- data.frame(tss=tss, start=ps, end=pe, strand=gtf_genes$strand, chr=gtf_genes$chr, tair=gtf_genes$gene_id)
  
  # define summit
  SUMMIT <- (prom_df$start + prom_df$end)/2
  prom_df[,'SUMMIT'] <- SUMMIT

  return(prom_df)
}

get_promoter_downstream <- function(gtf_df, prom_region=c(0, 1000)){
  
  # define TES
  gtf_genes <- gtf_df[gtf_df$feature == 'gene', ] # maybe replace with get_genes() and include transposable genes
  pos_strand <- gtf_genes$strand == "+"
  tes <- rep(0, length(pos_strand))
  tes[pos_strand] <- (gtf_genes$end)[pos_strand]
  tes[!pos_strand] <- (gtf_genes$start)[!pos_strand]
  
  # define promoter regions
  ps <- rep(0, length(tes)); pe <- rep(0, length(tes))
  ps[pos_strand] <- tes[pos_strand] + prom_region[1]
  pe[pos_strand] <- tes[pos_strand] + prom_region[2]
  ps[!pos_strand] <- tes[!pos_strand] - prom_region[2]
  pe[!pos_strand] <- tes[!pos_strand] - prom_region[1]
  prom_df <- data.frame(tes=tes, start=ps, end=pe, strand=gtf_genes$strand, chr=gtf_genes$chr, tair=gtf_genes$gene_id)
  
  # define summit
  SUMMIT <- (prom_df$start + prom_df$end)/2
  prom_df[,'SUMMIT'] <- SUMMIT
  
  return(prom_df)
}

get_5UTR <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "five_prime_UTR", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_3UTR <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "three_prime_UTR", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_CDS <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "CDS", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_Exon <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "exon", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_Intron <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "intron", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_lncRNA <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "lnc_RNA", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_ncRNA <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "nc_RNA", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_snoRNA <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "snoRNA", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}

get_miRNA <- function(gtf_df){

  # select genomic features from GTF
  feature_df <- gtf_df[gtf_df$feature == "miRNA", c('chr', 'start', 'end', 'strand', 'gene_id')]

  # define feature summit
  SUMMIT <- (feature_df$start + feature_df$end)/2
  feature_df[,'SUMMIT'] <- SUMMIT

  return(feature_df)
}


get_overlaps <- function(peaks_df, regions_df, genomic_feature, summit=TRUE,
                         unique_assignment=TRUE, next_gene_range = c(-2000, 2000),
                         add_tair=TRUE, add_target_index=FALSE){
  
  # peaks_df <- peaks_df
  # regions_df <- regions_df
  # genomic_feature <- "prom"
  
  # calculate peak summit
  SUMMIT <- (peaks_df$start + peaks_df$end)/2
  peaks_df[,'SUMMIT'] <- SUMMIT
  SUMMIT <- (regions_df$start + regions_df$end)/2
  regions_df[,'SUMMIT'] <- SUMMIT

  # use peak summit to calculate overlap
  if(summit){
    peaks_df$start_save <- peaks_df$start
    peaks_df$end_save <- peaks_df$end

    peaks_df$start <- peaks_df$SUMMIT
    peaks_df$end <- peaks_df$SUMMIT
  }
  if(genomic_feature == 'next_gene'){
    # span region around peak summit
    peaks_df$start <- peaks_df$SUMMIT + next_gene_range[1]
    peaks_df$end <- peaks_df$SUMMIT + next_gene_range[2]
  }

  # transform to GR object
  peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(peaks_df)
  regions_gr <- GenomicRanges::makeGRangesFromDataFrame(regions_df)
  
  peaks_df$start <- round(peaks_df$start, 0)
  peaks_df$end <- round(peaks_df$end, 0)
  # regions_gr[1:5,]
  
  # find overlaps of peaks with TSS regions
  res <- IRanges::findOverlaps(query = peaks_gr, subject = regions_gr)
  res <- data.frame(from=from(res), to=to(res))

  # find unique assignment based on distance of peak summit and region summit
  # res is a df which maps regions in peak_df to regions in regions_df
  # one region in peak_df can be mapped to more than one region in regions_df
  # this may make sense if one region is a promoter for two genes (regions in regions_df)
  # in this case the annotated peaks_df will contain one more row
  if(unique_assignment){
    # find distance to middle of peak
    ps <- peaks_df[res$from, 'SUMMIT']
    rs <- regions_df[res$to, 'SUMMIT']
    res[,'distance'] <- abs(ps-rs)
    # keep peak summit that is closest to region summit if there are more than one query-target pairs
    choosen_region <- NULL
    for(i in unique(res$from)){
      reg <- res[res$from == i, ]
      min_row_number <- which.min(reg$distance)
      min_row_name <- rownames(reg)[min_row_number]
      choosen_region <- c(choosen_region, min_row_name)
    }
    # subset overlaps
    res <- res[choosen_region, ]
  }

  # get regions for which no hits were found
  rows_hit <- unique(res$from); rows_hit <- 1:nrow(peaks_df) %in% rows_hit
  # get part of peaks_df for which no hits were found
  pd_nohits <- peaks_df[!rows_hit,]
  pd_nohits[, genomic_feature] <- rep(FALSE, nrow(pd_nohits))
  # get part of peaks_df for which hits were found
  pd_hits <- peaks_df[res$from,] # res$from could be replaced by rows_hit
  pd_hits[, genomic_feature] <- rep(TRUE, nrow(pd_hits))

  if(add_target_index){
    print("Adding row index of target...")
    # add row index of regions_df that was connected to peak in peaks_df
    pd_hits[, "hit_row_idx"] <- res$to
    pd_nohits[, "hit_row_idx"] <- -1
  }
  
  # if feature is promoter, add the genes that regions are promoter for
  if(genomic_feature %in% c("promoter_1kb", "promoter_2kb", "promoter_3kb")){
    gene_type <- paste0(genomic_feature, "_tair")
    pd_nohits[, gene_type] <- rep("", nrow(pd_nohits))
    pd_hits[, gene_type] <- regions_df[res$to, "tair"]
    
    # add distance to TSS of gene
    dist_name <- paste0('dist_tss_', genomic_feature)
    rd <- regions_df[res$to, ]
    strand <- rd[, "strand"]
    strand <- replace(strand, strand=="+", 1)
    strand <- replace(strand, strand=="-", -1)
    strand <- as.numeric(strand)
    pd_nohits[, dist_name] <- rep(NA, nrow(pd_nohits))
    pd_hits[strand==1, dist_name] <- pd_hits[strand==1, 'SUMMIT'] - rd[strand==1, "tss"]
    pd_hits[strand==-1, dist_name] <- rd[strand==-1, "tss"] - pd_hits[strand==-1, 'SUMMIT']
  }
  else if(genomic_feature == "promoter_1kb_downstream"){
    gene_type <- paste0(genomic_feature, "_tair")
    pd_nohits[, gene_type] <- rep("", nrow(pd_nohits))
    pd_hits[, gene_type] <- regions_df[res$to, "tair"]
    # add distance to TES of gene
    dist_name <- paste0('dist_tes_', genomic_feature)
    rd <- regions_df[res$to, ]
    strand <- rd[, "strand"]
    strand <- replace(strand, strand=="+", 1)
    strand <- replace(strand, strand=="-", -1)
    strand <- as.numeric(strand)
    pd_nohits[, dist_name] <- rep(NA, nrow(pd_nohits))
    pd_hits[strand==1, dist_name] <- pd_hits[strand==1, 'SUMMIT'] - rd[strand==1, "tes"]
    pd_hits[strand==-1, dist_name] <- rd[strand==-1, "tes"] - pd_hits[strand==-1, 'SUMMIT']
  }
  # add gene name if next_gene is requested
  else if(genomic_feature == 'next_gene'){
    # add next closest gene
    pd_nohits[, 'next_gene_tair'] <- rep("", nrow(pd_nohits))
    pd_hits[, 'next_gene_tair'] <- regions_df[res$to, 'gene_id']
    # add distance to TSS of next closest gene
    dist_name <- paste0('dist_summit_', genomic_feature)
    rd <- regions_df[res$to, ]
    strand <- rd[, "strand"]
    strand <- replace(strand, strand=="+", 1)
    strand <- replace(strand, strand=="-", -1)
    strand <- as.numeric(strand)
    pd_nohits[, dist_name] <- rep(NA, nrow(pd_nohits))
    pd_hits[strand==1, dist_name] <- pd_hits[strand==1, 'SUMMIT'] - rd[strand==1, 'SUMMIT']
    pd_hits[strand==-1, dist_name] <- rd[strand==-1, 'SUMMIT'] - pd_hits[strand==-1, 'SUMMIT']
  }
  else {
    if(add_tair & "gene_id" %in% colnames(regions_df)){
      gene_type <- paste0(genomic_feature, "_tair")
      pd_nohits[, gene_type] <- rep("", nrow(pd_nohits))
      pd_hits[, gene_type] <- regions_df[res$to, "gene_id"]
    }
  }
  # combine separated DataFrames
  peaks_df <- as.data.frame(rbind(pd_nohits, pd_hits))
  
  # put back original start and end of peak if summit was used
  if(summit){
    peaks_df$start <- peaks_df$start_save
    peaks_df$end <- peaks_df$end_save
    peaks_df$start_save <- NULL
    peaks_df$end_save <- NULL
  }

  return(peaks_df)
}

annotate_peaks_internal <- function(peaks_df, gtf_df, genomic_feature, next_gene_range=c(-2000, 2000),
                                    prom_1kb_region=c(-1000, 0), prom_2kb_region=c(-2000, -1000), prom_3kb_region=c(-3000, -2000),
                                    prom_dowstream_region=c(0, 1000), unique_promoter_assignment=TRUE){

  # Find whether peaks overlap with genomic features
  if(genomic_feature == 'next_gene'){
    regions_df <- get_genes(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature,
                            next_gene_range = next_gene_range, summit = TRUE, unique_assignment = TRUE)
                            # force to use summit instead of range to find overlaps
                            # necessary, otherwise the way I defined finding next genes does not work
  }
  else if(genomic_feature == 'genic'){
    regions_df <- get_genes(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'promoter_1kb'){
    regions_df <- get_promoters(gtf_df = gtf_df, prom_region = prom_1kb_region)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature,
                            unique_assignment=unique_promoter_assignment, summit=TRUE)
  }
  else if(genomic_feature == 'promoter_2kb'){
    regions_df <- get_promoters(gtf_df = gtf_df, prom_region = prom_2kb_region)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature,
                            unique_assignment=unique_promoter_assignment, summit=TRUE)
  }
  else if(genomic_feature == 'promoter_3kb'){
    regions_df <- get_promoters(gtf_df = gtf_df, prom_region = prom_3kb_region)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature,
                            unique_assignment=unique_promoter_assignment, summit=TRUE)
  }
  else if(genomic_feature == 'promoter_1kb_downstream'){
    regions_df <- get_promoter_downstream(gtf_df = gtf_df, prom_region=prom_dowstream_region)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature,
                            unique_assignment=unique_promoter_assignment, summit=TRUE)
  }
  else if(genomic_feature == 'UTR5'){
    regions_df <- get_5UTR(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'UTR3'){
    regions_df <- get_3UTR(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'CDS'){
    regions_df <- get_CDS(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'Exon'){
    regions_df <- get_Exon(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'Intron'){
    regions_df <- get_Intron(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'lncRNA'){
    regions_df <- get_lncRNA(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'ncRNA'){
    regions_df <- get_ncRNA(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'snoRNA'){
    regions_df <- get_snoRNA(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }
  else if(genomic_feature == 'miRNA'){
    regions_df <- get_miRNA(gtf_df = gtf_df)
    anno_df <- get_overlaps(peaks_df = peaks_df, regions_df = regions_df, genomic_feature = genomic_feature)
  }

  return(anno_df)
}

annotate_peaks <- function(peaks_df, gtf_df, summit=TRUE, unique_assignment=TRUE, annotation_priority=NULL,
                           next_gene=TRUE, next_gene_range = c(2000, 2000),
                           prom_1kb_region=c(-1000, 0), prom_2kb_region=c(-2000, -1000), prom_3kb_region=c(-3000, -2000),
                           prom_dowstream_region=c(0, 1000), translate_tair=TRUE,
                           unique_promoter_assignment=FALSE){

  # set annotation_priority if user does not provide input
  if(is.null(annotation_priority)){
    annotation_priority <- c('promoter_1kb', 'promoter_2kb', 'promoter_3kb', 'promoter_1kb_downstream', 
                                   'UTR5', 'Intron', 'Exon', 'UTR3', 'CDS',
                                   'intergenic', 'lncRNA', 'snoRNA', 'miRNA', 'genic')
  }
  if('intergenic' %in% annotation_priority){
    # remove the feature intergenic since there is no function finding it and keeping it leads to error
    # it will be added at the end
    index <- which(annotation_priority == 'intergenic')
    an_pr <- annotation_priority[-index]
  }else{
    an_pr <- annotation_priority
  }

  # go through all features and test whether peaks belong to any of those features
  message("Annotating peaks with selected genomic features ...")
  for(feature in an_pr){
    message(paste0("|-- ", feature))
    peaks_df <- annotate_peaks_internal(peaks_df = peaks_df, gtf_df = gtf_df,
                                        genomic_feature = feature, 
                                        prom_1kb_region=prom_1kb_region, prom_2kb_region=prom_2kb_region, prom_3kb_region=prom_3kb_region,
                                        prom_dowstream_region=prom_dowstream_region,
                                        unique_promoter_assignment = unique_promoter_assignment)
    if(translate_tair){
      # translate tair IDs
      gene_type_tair <- paste0(feature, "_tair")
      gene_type_symbol <- paste0(feature, "_symb")
      peaks_df[, gene_type_symbol] <- tair_to_symbol(peaks_df[, gene_type_tair])
    }
  }

  if(next_gene){
    message("|-- next_gene")
    peaks_df <- annotate_peaks_internal(peaks_df = peaks_df, gtf_df = gtf_df, genomic_feature = 'next_gene',
                                        next_gene_range = next_gene_range)
    if(translate_tair){
      # translate tair IDs of closest genes if user requested finding of closest genes to peaks
      peaks_df$next_gene_symbol <- tair_to_symbol(peaks_df$next_gene_tair)
    }
  }

  # classify all regions that have not been annotated as 'intergenic' (if the user wants to classify intergenic regions)
  if('intergenic' %in% annotation_priority){
    message("|-- intergenic")
    intergenic_regions <- rowSums(peaks_df[, an_pr]) == 0
    peaks_df[, 'intergenic'] <- intergenic_regions
    # adding ID of closest gene
    peaks_df[, 'intergenic_tair'] <- rep("", nrow(peaks_df))
    peaks_df[, 'intergenic_symbol'] <- rep("", nrow(peaks_df))
    if(sum(peaks_df$intergenic) > 0){ # only attempt to add tair and symbol of there are intergeneic regions
      peaks_df[intergenic_regions, 'intergenic_tair'] <- peaks_df[intergenic_regions, 'next_gene_tair']
      if(translate_tair){
        peaks_df[intergenic_regions, 'intergenic_symbol'] <- peaks_df[intergenic_regions, 'next_gene_symbol']
      }
    }
  }

  # find final annotation for each feature
  message("|-- Final Annotation")
  annotation_priority_rev <- base::rev(annotation_priority)
  final_annotation <- rep("", nrow(peaks_df))
  for(feature in annotation_priority_rev){
    final_annotation[peaks_df[, feature]] <- feature
  }
  peaks_df[,'annotation'] <- final_annotation
  
  message("... Genomic Feature Annotation Complete")
  
  return(peaks_df)
}

tair_to_symbol <- function(genes){
  symbols <- NULL
  for(gene in genes){
    k <- annotate::getSYMBOL(gene, data='org.At.tair.db')
    symb <- as.character(k[length(k)])
    if(is.na(symb)){symbols <- c(symbols, '')}
    else{symbols <- c(symbols, symb)}
  }
  symbols <- toupper(symbols)
  return(symbols)
}
