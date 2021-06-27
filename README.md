# GenomicsPeakAnnotation

This repositroy contains functions that are used to annotate eg. ChIP-seq or ATAC-seq data.


```R
ap <- c('promoter_1kb', 'promoter_2kb', 'promoter_3kb', 'promoter_1kb_downstream', 
        'UTR5', 'Intron', 'Exon', 'UTR3', 'CDS',
        'intergenic', 'lncRNA', 'snoRNA', 'miRNA', 'genic')
```
```R
anno <- annotate_peaks(peaks_df = peaks_df, gtf_df = gtf_df, next_gene = TRUE,
                       summit = TRUE, unique_assignment = TRUE, next_gene_range = c(-10000, 10000),
                       translate_tair = TRUE, annotation_priority = ap,
                       prom_1kb_region=c(-1000, 0), prom_2kb_region=c(-2000, -1000), prom_3kb_region=c(-3000, -2000),
                       prom_dowstream_region=c(0, 1000),
                       unique_promoter_assignment=FALSE)
```
peaks_df: [DataFrame]  
containing regions (peaks) with columns specifying chromosome, start and end of the peak. 

gtf_df: [DataFrame]  
containing regions (genomic features) with columns specifying chromosome, start and end of the genomic feature  

next_gene: [boolean]  
specifying whether the next closest gene to peak should be found  

next_gene_range: [numeric vector]  
specifying the range around the peak that should be searched for genomic features   

summit: [boolean]  
specifying whether the peak summit (=TRUE) or the whole peak range (=FALSE) should be used to find genomic features   

unique_assignment: [boolean] 
specifying whether peaks should be assigned to only one genomic feature (=TRUE) or to more than one genomic feature (=FALSE)   

unique_promoter_assignment: [boolean] 
specifying whether peaks should be assigned to only one promoter (=TRUE) or to more than one promoter (=FALSE)   

translate_tair: [boolean]   
specifying whether genes of genomic features that have been associated with a peak, should be translated from the gene symbol to the TAIR ID  

annotation_priority: [character vector]   
specifying the priority of genomic features, which is used to obtain final annotation of peaks  

prom_1kb_region: [numeric vector]   
specifying the range of base pairs upstream/around the gene transcription start site (TSS) that should be used to define proximal promoters  

prom_2kb_region: [numeric vector]   
specifying the range of base pairs upstream/around the gene transcription start site (TSS) that should be used to define medium-proximal promoters  

prom_3kb_region: [numeric vector]   
specifying the range of base pairs upstream/around the gene transcription start site (TSS) that should be used to define distal promoters  

prom_dowstream_region: [numeric vector]   
specifying the range of base pairs upstream/around the gene transcription termination site (TTS) that should be used to define proximal downstream promoters  


