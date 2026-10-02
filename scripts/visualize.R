#!/usr/bin/env Rscript
# Uses DESeq2 saved model and complete peak results, preserving statistical normalization.
suppressPackageStartupMessages(library(DESeq2))
a<-commandArgs(TRUE);if(length(a)!=2)stop('model_directory NEW_output_directory')
if(dir.exists(a[2]))stop('Output exists')
dds<-readRDS(file.path(a[1],'model.rds'));z<-read.delim(file.path(a[1],'complete.tsv'))
v<-log2(counts(dds,normalized=TRUE)+1);dir.create(a[2],recursive=TRUE)
write.table(cor(v,method='spearman'),file.path(a[2],'spearman.tsv'),sep='\t',quote=FALSE)
top<-head(order(apply(v,1,var),decreasing=TRUE),500);p<-prcomp(t(v[top,,drop=FALSE]),center=TRUE,scale.=FALSE)
write.table(data.frame(sample_id=rownames(p$x),p$x),file.path(a[2],'PCA.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
pdf(file.path(a[2],'diagnostics.pdf'),width=8,height=7)
plot(p$x[,1],p$x[,2],pch=19,xlab='PC1',ylab='PC2',main='Descriptive PCA: log2(normalized counts + 1)');text(p$x[,1],p$x[,2],rownames(p$x),pos=3,cex=.7)
plot(z$log2FoldChange,pmin(100,-log10(pmax(z$padj,1e-100))),pch=16,cex=.4,col=ifelse(z$direction=='gain','#D55E00',ifelse(z$direction=='loss','#0072B2','grey70')),xlab='log2 fold change',ylab='-log10 FDR (capped 100)')
sel<-head(order(z$padj,na.last=NA),100);ids<-intersect(z$peak_id[sel],rownames(v));mat<-v[ids,,drop=FALSE];mat<-mat[apply(mat,1,sd)>0,,drop=FALSE]
if(nrow(mat)>=2)heatmap(t(scale(t(mat))),scale='none',labRow=NA,main='Top tested peaks: row Z scores')
dev.off();capture.output(sessionInfo(),file=file.path(a[2],'sessionInfo.txt'))
