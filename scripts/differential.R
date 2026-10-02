#!/usr/bin/env Rscript
# Usage: Rscript differential.R counts.tsv metadata.tsv numerator denominator out normalization [spikein.tsv]
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(TRUE)
alpha <- as.numeric(Sys.getenv('CUTTAG_ALPHA','0.05'))
lfc <- as.numeric(Sys.getenv('CUTTAG_LFC','1'))
min_count <- as.numeric(Sys.getenv('CUTTAG_MIN_COUNT','10'))
stopifnot(is.finite(alpha),alpha>0,alpha<1,is.finite(lfc),lfc>=0,is.finite(min_count),min_count>=0)
if(length(a)<6) stop('counts metadata numerator denominator NEW_output conventional|spikein [spikein_table]')
x <- read.delim(a[1],check.names=FALSE); m <- read.delim(a[2],check.names=FALSE)
stopifnot(all(c('peak_id','chrom','start','end') %in% names(x)),all(c('sample_id','biological_sample_id','condition','replicates_confirmed') %in% names(m)))
m <- m[m$condition %in% a[3:4],,drop=FALSE]
stopifnot(!anyDuplicated(m$sample_id),!anyDuplicated(m$biological_sample_id),all(m$replicates_confirmed==TRUE),all(table(factor(m$condition,levels=a[3:4]))>=2),!anyDuplicated(x$peak_id),all(m$sample_id %in% names(x)))
k <- as.matrix(x[,m$sample_id,drop=FALSE]);stopifnot(all(is.finite(k)),all(k>=0),all(k==floor(k)),max(k)<=.Machine$integer.max)
storage.mode(k)<-'integer';rownames(k)<-x$peak_id;rownames(m)<-m$sample_id
m$condition<-factor(m$condition,levels=c(a[4],a[3])); keep<-rowSums(k)>=min_count
if(!any(keep)) stop('No peaks pass count filter')
dds<-DESeqDataSetFromMatrix(k[keep,,drop=FALSE],m,~condition)
if(a[6]=='spikein') {
 if(length(a)!=7) stop('Spike-in table required; calibration must be experimentally accepted')
 s<-read.delim(a[7]);s<-s[match(m$sample_id,s$sample_id),];stopifnot(identical(s$sample_id,m$sample_id),all(is.finite(s$spikein_fragments)),all(s$spikein_fragments>0),all(s$calibration_accepted==TRUE))
 sizeFactors(dds)<-s$spikein_fragments/exp(mean(log(s$spikein_fragments)))
} else if(a[6]!='conventional') stop('Unknown normalization')
if(dir.exists(a[5])) stop('Output exists; choose new model directory')
dir.create(a[5],recursive=TRUE)
write_tsv<-function(z,n) write.table(z,file.path(a[5],n),sep='\t',quote=FALSE,row.names=FALSE)
write_tsv(data.frame(peak_id=x$peak_id,retained=keep),'filter_decisions.tsv')
write_tsv(m,'metadata.tsv')
dds<-DESeq(dds,quiet=TRUE);res<-as.data.frame(results(dds,contrast=c('condition',a[3],a[4]),alpha=alpha))
z<-cbind(x[match(rownames(res),x$peak_id),1:4],res)
z$direction<-'not_significant';sig<-!is.na(z$padj)&z$padj<alpha&abs(z$log2FoldChange)>=lfc
z$direction[sig & z$log2FoldChange>0]<-'gain';z$direction[sig & z$log2FoldChange<0]<-'loss'
write_tsv(z,'complete.tsv');write_tsv(z[z$direction=='gain',],'gain.tsv');write_tsv(z[z$direction=='loss',],'loss.tsv')
write_tsv(data.frame(sample_id=colnames(dds),size_factor=sizeFactors(dds)),'size_factors.tsv')
write_tsv(data.frame(peak_id=rownames(dds),counts(dds,normalized=TRUE),check.names=FALSE),'normalized_counts.tsv')
saveRDS(dds,file.path(a[5],'model.rds'));capture.output(sessionInfo(),file=file.path(a[5],'sessionInfo.txt'))
pdf(file.path(a[5],'MA.pdf'));plotMA(results(dds,contrast=c('condition',a[3],a[4])),main=paste(a[3],'/',a[4]));dev.off()
writeLines(c(paste('Contrast:',a[3],'/',a[4]),paste('Normalization:',a[6]),paste('Design: ~ condition; min_total_count=',min_count,'; FDR<',alpha,'; abs(log2FC)>=',lfc),'Peak signal differences do not constitute RNA expression differences.'),file.path(a[5],'METHODS.txt'))
