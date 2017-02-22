#title: Olfactory Receptors and Chromatin Regulators/TFs in the Neuronal Lineage
#author: Diya Das
#date: January 30th 2017

library(clusterExperiment)
library(scales)
library(Biobase)

expt_str <- "oeHBCdiff"
clust_dir <- file.path("../output/clust", expt_str)
viz_dir <- file.path("../output/viz", expt_str)
load(file.path(clust_dir, paste0(expt_str, "_slingshot.Rda")))
load(file.path(clust_dir, paste0(expt_str, "_lineageData.Rda")))
load("../data/oeHBCdiff_RSEM_eSet.Rda")

load(file.path(clust_dir,paste0(expt_str,"_cmmerged.Rda")))
ors <- rownames(RSEM_eSet)[grep("^Olfr[^ps]*$",rownames(RSEM_eSet),ignore.case=TRUE)]
or_names <- sapply(ors, function(x) strsplit(x,"\\.")[[1]][1])
or_counts <- data.frame(assayData(RSEM_eSet)$tpm_table[ ors,colnames(nlm)])
or_counts[is.na(or_counts)] <- 0
or_sums <- sapply(or_counts, function(x) tapply(x, or_names,sum))

## Plots
### For an analysis of ORs using Kallisto TPMs, it is sufficient to replace "RSEM" with "kallisto" throughout this script.
png(file.path(viz_dir, paste0(expt_str, "ORchrom_RSEM_",Sys.Date(),".png")), height=11, width=20, units="cm", res=600)
par(mar=c(0,2,1,1),family = 'Helvetica', cex.lab=0.6, cex.axis=0.6, oma=c(3,0,0,0), mgp=c(1,0.2,0),xpd=NA, lwd=0.5, cex=0.5)
layout(cbind(c(1,1,1,2,2,2,3,3,3),c(4:12),c(13:19,0,0),c(20:25,0,0,0)),
       heights=rep(0.5,9))
plot(1:ncol(or_sums),or_sums[1,],pch=16,xlab='',ylim=c(1,max(or_sums)),ylab='TPM', col=alpha(colpalN[nclus.labels],0.5), xaxt='n', yaxt='n', cex=0.5)
axis(2, lwd = 0.5, tck=-0.01)
text(par()$usr[1],par()$usr[4],"OR", adj = c(-0.1,1.2), cex=0.6)
invisible(lapply(2:nrow(or_sums), function(x){
  points(1:ncol(or_sums),or_sums[x,],pch=16,xlab='',ylab='', col=alpha(colpalN[nclus.labels],0.5), cex=0.5)
}
))
nvec = c(100, 5000 )
for (i in 1:length(nvec)){
  n_thresh = nvec[i]
  ORnum <- sapply(1:ncol(or_sums), function(x) sum(or_sums[,x] > n_thresh, na.rm=TRUE))
  plot(1:length(ORnum), ORnum, pch=16, ylab="#ORs", col=alpha(colpalN[nclus.labels],0.5),xaxt='n',xlab='', yaxt='n', cex=0.5)
  text(par()$usr[1],par()$usr[4],paste('>',n_thresh,'TPM'), adj = c(-0.1,1.2), cex=0.6)
  axis(2, lwd = 0.5, tck=-0.01*2)
}
axis(1, lwd = 0.5, tck=-0.01*2)
title(xlab="Developmental Order")
ngenes = c("Ehmt1", "Setdb1","Setdb2","Kdm1a", "Lbr","Lhx2","Ebf1","Ebf3", "Ebf4", "Gnai2","Gnas","Gnal","Gnb1", "Gng8","Gng13","Gng3","Adcy3","Ano2","Cnga2","Cnga4","Atf5","Rtp1")
for (i in 1:length(ngenes)){
  gene <- ngenes[i]
  gene_exp <- 2^(nlm[gene,])-1
  plot(1:ncol(nlm),gene_exp, pch=16,xlab='',ylim=c(1,ceiling(max(max(gene_exp)))),ylab='Counts', col=alpha(colpalN[nclus.labels],0.5),xaxt='n',yaxt='n', cex=0.5)
  text(par()$usr[1],par()$usr[4],gene, adj = c(-0.1,1.2), cex=0.6) #0.6
  axis(2, lwd = 0.5, tck=-0.01*2)
  if (gene %in% c("Ebf4","Gng3","Rtp1")){
    axis(1, lwd = 0.5, tck=-0.01*2)
    title(xlab="Developmental Order")
  }
}
dev.off()



