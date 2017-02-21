
#title: Differentially Expressed Transcription Factor Heatmaps
#author: Russell Fletcher and Diya Das
#date: January 27th 2017

library(NMF)
library(clusterExperiment)
NMF::nmf.options(grid.patch=TRUE)
expt_str <- "oeHBCdiff"

clust_dir <- file.path("../output/clust", expt_str)
viz_dir <- file.path("../output/viz", expt_str)
DE_dir <- file.path("../output/DE", expt_str)

load(file.path(DE_dir,"NL_SL_oneVall500DE_genes.Rda"))
load(file.path(DE_dir, paste0(expt_str, "_NL_geneCl_final.Rda")))
load(file.path(DE_dir, paste0(expt_str, "_SL_geneCl_final.Rda")))
load(file.path(clust_dir, paste0(expt_str,"_lineagedata.Rda")))
TFgenes <- intersect(unlist(read.table("../ref/ATFDB_mm_TF.txt")), rownames(nlm))

NLgClList <- c("m1","m2","m5","m10","m16","m7","m13","m6","m9","m3","m11","m12","m17","m4","m8","m14","m15")
SLgClList <- c("m1","m3","m7","m6","m11","m2","m4","m5","m8","m9","m10","m12")

NL_gclus <- data.frame(genes=c(colnames(cegNL)),clus=factor(c(primaryClusterNamed(cegNL)), levels=NLgClList))
NL_gclus <- NL_gclus[complete.cases(NL_gclus),] 

SL_gclus <- data.frame(genes=c(colnames(cegSL)),clus=factor(c(primaryClusterNamed(cegSL)), levels=SLgClList))
SL_gclus <- SL_gclus[complete.cases(SL_gclus),]

nl_detfs <- intersect(NL_DEgenes, TFgenes)
sl_detfs <- intersect(SL_DEgenes, TFgenes)

nl_var_tf <- nlm[nl_detfs,]
sl_var_tf <- slm[sl_detfs,]
breakv <- c(min(nl_var_tf), seq(0, quantile(nl_var_tf[nl_var_tf > 0], .98, na.rm = TRUE), length = 50), max(nl_var_tf))

NL_tfs <- unlist(sapply(levels(NL_gclus$clus), function(x){
    intersect(rownames(nl_var_tf) , NL_gclus$genes[NL_gclus$clus==x])}))
pdf(file=file.path(DE_dir, paste0(expt_str, "_DE_OneVall_HM_Neur.pdf")), width=11, height=8.5)
NMF::aheatmap(nl_var_tf[NL_tfs,], color=seqPal5, Colv=NA, Rowv=NA, annCol = data.frame(Clusters=droplevels(nclus.labels)), annColors = list(Clusters=colpalN), breaks = breakv)
dev.off()

SL_tfs <- unlist(sapply(levels(SL_gclus$clus), function(x){
  intersect(rownames(sl_var_tf) , SL_gclus$genes[SL_gclus$clus==x])}))
pdf(file=file.path(DE_dir, paste0(expt_str, "_DE_OneVall_HM_Sus.pdf")), width=11, height=8.5)
NMF::aheatmap(sl_var_tf[SL_tfs,], color=seqPal5, Colv=NA, Rowv=NA, annCol = data.frame(Clusters=droplevels(sclus.labels)), annColors = list(Clusters=colpalS), breaks=breakv)
dev.off()

