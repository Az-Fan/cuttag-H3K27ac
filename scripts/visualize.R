#!/usr/bin/env Rscript
# Descriptive diagnostics; uninformative data produce explicit panels instead of fabricated axes.
suppressPackageStartupMessages(library(DESeq2))
a<-commandArgs(TRUE);if(length(a)!=2)stop('model_directory NEW_output_directory')
if(dir.exists(a[2]))stop('Output exists')
dds<-readRDS(file.path(a[1],'model.rds'));z<-read.delim(file.path(a[1],'complete.tsv'))
v<-log2(counts(dds,normalized=TRUE)+1);dir.create(a[2],recursive=TRUE)
if(!all(is.finite(v)))stop('Nonfinite normalized counts')
correlations<-suppressWarnings(cor(v,method='spearman'))
write.table(correlations,file.path(a[2],'spearman.tsv'),sep='\t',quote=FALSE)
variances<-apply(v,1,var);top<-head(order(variances,decreasing=TRUE),500);top<-top[is.finite(variances[top]) & variances[top]>0]
p<-if(length(top)>=2 && ncol(v)>=2)prcomp(t(v[top,,drop=FALSE]),center=TRUE,scale.=FALSE) else NULL
if(!is.null(p))write.table(data.frame(sample_id=rownames(p$x),p$x),file.path(a[2],'PCA.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
notes<-character()
pdf(file.path(a[2],'diagnostics.pdf'),width=8,height=7)
if(!is.null(p) && ncol(p$x)>=2) {
 plot(p$x[,1],p$x[,2],pch=19,xlab='PC1',ylab='PC2',main='Descriptive PCA: log2(normalized counts + 1)');text(p$x[,1],p$x[,2],rownames(p$x),pos=3,cex=.7)
} else {plot.new();text(.5,.5,'PCA unavailable: fewer than two variable features');notes<-c(notes,'PCA unavailable: fewer than two variable features')}
valid<-is.finite(z$log2FoldChange)&is.finite(z$padj)&z$padj>=0&z$padj<=1
if(any(valid)) {
 x<-z[valid,,drop=FALSE];plot(x$log2FoldChange,pmin(100,-log10(pmax(x$padj,1e-100))),pch=16,cex=.4,col=ifelse(x$direction=='gain','#D55E00',ifelse(x$direction=='loss','#0072B2','grey70')),xlab='log2 fold change',ylab='-log10 FDR (capped 100)')
} else {plot.new();text(.5,.5,'Volcano unavailable: no finite adjusted tests');notes<-c(notes,'Volcano unavailable: no finite adjusted tests')}
sel<-head(order(z$padj,na.last=NA),100);ids<-intersect(z$peak_id[sel],rownames(v));mat<-v[ids,,drop=FALSE]
if(nrow(mat))mat<-mat[apply(mat,1,sd)>0,,drop=FALSE]
if(nrow(mat)>=2 && ncol(mat)>=2)heatmap(t(scale(t(mat))),scale='none',labRow=NA,main='Top tested peaks: row Z scores') else {plot.new();text(.5,.5,'Heatmap unavailable: too few variable tested peaks');notes<-c(notes,'Heatmap unavailable: too few variable tested peaks')}
if(anyNA(correlations))notes<-c(notes,'Some sample correlations are undefined for constant profiles')
dev.off();writeLines(notes,file.path(a[2],'unavailable_diagnostics.txt'));capture.output(sessionInfo(),file=file.path(a[2],'sessionInfo.txt'))
