#!/usr/bin/env Rscript
# Offline TERM2GENE table permits GO/KEGG/Reactome/Hallmark snapshots.
suppressPackageStartupMessages(library(clusterProfiler))
a<-commandArgs(TRUE);if(length(a)!=4)stop('target_gene_ids.txt tested_universe.txt term2gene.tsv NEW_output.tsv')
if(file.exists(a[4]))stop('Output exists')
g<-unique(readLines(a[1]));u<-unique(readLines(a[2]));t<-read.delim(a[3],colClasses='character')
stopifnot(all(g %in% u),ncol(t)==2)
z<-enricher(gene=g,universe=u,TERM2GENE=t,pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500,pAdjustMethod='BH')
write.table(as.data.frame(z),a[4],sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=paste0(a[4],'.sessionInfo.txt'))
