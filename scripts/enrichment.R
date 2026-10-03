#!/usr/bin/env Rscript
# Offline TERM2GENE with an explicit tested-gene universe and auditable mapping losses.
suppressPackageStartupMessages(library(clusterProfiler))
a<-commandArgs(TRUE);if(length(a)!=4)stop('target_gene_ids.txt tested_universe.txt term2gene.tsv NEW_output.tsv')
if(file.exists(a[4]))stop('Output exists')
ids<-function(path) { x<-trimws(readLines(path));unique(x[nzchar(x) & !x %in% c('NA','NaN')]) }
g<-ids(a[1]);u<-ids(a[2]);t<-read.delim(a[3],colClasses='character',check.names=FALSE)
if(!length(u) || !all(g %in% u) || ncol(t)!=2 || anyNA(t) || any(!nzchar(t[[1]])) || any(!nzchar(t[[2]])))stop('Invalid target/universe/TERM2GENE contract')
t<-unique(t);mapped_u<-intersect(u,t[[2]]);mapped_g<-intersect(g,mapped_u)
if(!length(mapped_u))stop('No universe genes map to TERM2GENE; check identifier namespace')
dir.create(dirname(a[4]),recursive=TRUE,showWarnings=FALSE)
audit<-data.frame(gene_id=u,in_target=u %in% g,in_term2gene=u %in% mapped_u)
write.table(audit,paste0(a[4],'.mapping.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# enricher uses the intersection of supplied universe with database annotations.
if(length(mapped_g)) {
 z<-enricher(gene=mapped_g,universe=mapped_u,TERM2GENE=t,pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500,pAdjustMethod='BH')
 result<-as.data.frame(z)
} else result<-data.frame()
if(!ncol(result))result<-data.frame(ID=character(),Description=character(),GeneRatio=character(),BgRatio=character(),pvalue=numeric(),p.adjust=numeric(),qvalue=numeric(),geneID=character(),Count=integer())
write.table(result,a[4],sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(paste('State:',if(!length(mapped_g)) 'SKIPPED_NO_MAPPED_TARGET' else if(!nrow(result)) 'NO_ELIGIBLE_TERMS' else 'COMPUTATIONAL_PASS'),paste('Input target:',length(g)),paste('Input universe:',length(u)),paste('Database-mapped target:',length(mapped_g)),paste('Effective universe:',length(mapped_u)),'Gene set sizes: 10..500; BH adjustment; identical identifier namespace must be reviewed.'),paste0(a[4],'.audit.txt'))
capture.output(sessionInfo(),file=paste0(a[4],'.sessionInfo.txt'))
