#!/usr/bin/env Rscript
# Explicit TxDb SQLite and organism OrgDb package; no implicit human reference.
suppressPackageStartupMessages({library(ChIPseeker);library(GenomicFeatures);library(AnnotationDbi)})
a<-commandArgs(TRUE);if(length(a)!=4)stop('master.bed TxDb.sqlite OrgDb_package NEW_output.tsv')
if(file.exists(a[4]))stop('Output exists')
b<-read.delim(a[1],header=FALSE);stopifnot(ncol(b)>=4,!anyDuplicated(b[[4]]),all(b[[2]]>=0),all(b[[3]]>b[[2]]))
g<-GenomicRanges::GRanges(b[[1]],IRanges::IRanges(b[[2]]+1,b[[3]]));names(g)<-b[[4]];S4Vectors::mcols(g)$peak_id<-as.character(b[[4]])
txdb<-loadDb(a[2]);z<-as.data.frame(annotatePeak(g,TxDb=txdb,tssRegion=c(-2000,2000),annoDb=a[3]))
stopifnot('peak_id' %in% names(z),!anyDuplicated(z$peak_id),setequal(z$peak_id,b[[4]]))
z<-z[match(b[[4]],z$peak_id),,drop=FALSE];z$start<-z$start-1
stopifnot(identical(as.character(z$seqnames),as.character(b[[1]])),all(z$start==b[[2]]),all(z$end==b[[3]]))
write.table(z,a[4],sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=paste0(a[4],'.sessionInfo.txt'))
