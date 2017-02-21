#title: Differentially Expressed Transcription Factor Heatmaps
#author: Diya Das
#date: January 30th 2017

library(clusterExperiment)
expt_str <- "oeHBCdiff"
esh <- gsub("Expt","E",expt_str)
DE_dir <- file.path("../output/DE", expt_str)
clust_dir <- file.path("../output/clust", expt_str)

load(file.path(clust_dir, paste0(esh,"_cmmerged.Rda")))
load(file.path(clust_dir, paste0(esh, "_slingshot.Rda")))

pairsDE <- getBestFeatures(transform(cmobj), primaryCluster(cmobj), contrastType="Pairs",number=Inf,p.value=0.05)
write.table(pairsDE, file.path(DE_dir,"all_DE_pairs.txt"),quote=FALSE,sep="\t")

get_tuple <- function(lineages, n) rbind(do.call(rbind,lapply(lineages[n], function(x) data.frame(cbind(x[-1],x[-length(x)]), stringsAsFactors = F))))

contrast_tuple <- unique(get_tuple(lineages, 1:3))

contrast_name <- data.frame(contrast = paste0("X",contrast_tuple$X1,"-X",contrast_tuple$X2), names= c("H1-rH", "H2-H1", "G-H2", "I1-G", "I2-I1", "I3-I2", "iO-I3", "mO-iO", "MV1-GBC","MV2-MV1", "iS-H2", "mS-iS"))

plotContrasts <- function(DE, contrast_tuple, contrast_name, fname, pval_thresh){
  par(family='URWHelvetica', cex.lab=0.5, cex.axis=0.5, oma=c(0,0,0,0), mgp=c(1.2,0.15,0),xpd=F, lwd=0.5, pty="s", mar=c(0,2,1,1))
  write(c("Contrast name", "Positive Change", "Negative change", paste("max_pval <", pval_thresh)), file=fname, sep = "\t", append=FALSE, ncolumns = 4)
  tmplist <- apply(contrast_tuple, 1, function(x){
    desired_contrast <- paste0("X",x[1],"-X",x[2])
    opposite_contrast <- paste0("X",x[2],"-X",x[1])
    if (desired_contrast %in% levels(DE$Contrast)){
      cont_sub <<- subset(DE, abs(logFC) > 1 & Contrast==desired_contrast)
    } else if (opposite_contrast %in% levels(DE$Contrast)) {
      cont_sub <<- subset(DE, abs(logFC) > 1 & Contrast==opposite_contrast)
      cont_sub$logFC <<- -cont_sub$logFC
    } else {
      stop('stop something wrong')
    }
  
    # plotting
    pdf(file.path(DE_dir,paste0(expt_str, "_", contrast_name$names[contrast_name$contrast==desired_contrast],"_volcano.pdf")), width=2, height = 2)
    par(cex.lab=0.5, cex.axis=0.5, oma=c(0,0,0,0), mgp=c(1.2,0.15,0),xpd=F, lwd=0.5, pty="s", mar=c(1,1,1,1))
    with(cont_sub, plot(logFC, -log10(adj.P.Val),xlim=c(-12,12), pch=20, col="grey", ylim=c(0,100),cex=0.5, xaxt='n',yaxt='n',xlab='',ylab=''))
    axis(1, lwd = 0.5, tck=-0.01*2, labels=FALSE);axis(2, lwd = 0.5, tck=-0.01*2, labels=FALSE)
    with(subset(cont_sub, adj.P.Val<pval_thresh & abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="black", cex=0.5))
    
    # label number of DE genes
    neg_num <- with(subset(cont_sub, adj.P.Val<pval_thresh & logFC < -1), length(Feature))
    pos_num <- with(subset(cont_sub, adj.P.Val<pval_thresh & logFC > 1), length(Feature))
    usr <- par( "usr" )
    text( usr[ 1 ], usr[ 4 ], neg_num, adj = c( -0.5, 2 ), col = "dodgerblue4", font=4)
    text( usr[ 2 ], usr[ 4 ], pos_num, adj = c( 1.5, 2 ), col = "firebrick4", font=4)
    
    dev.off()

    with(subset(cont_sub, adj.P.Val<pval_thresh & abs(logFC)>1), write(c(contrast_name$names[contrast_name$contrast==desired_contrast], sum(logFC>0), sum(logFC<0), max(adj.P.Val)), ncolumns=4, file=fname, sep = "\t", append=TRUE))
  })
}

pval_thresh=1e-2

plotContrasts(pairsDE, contrast_tuple, contrast_name, fname=file.path(DE_dir,paste0(expt_str,"_pairsDE_tuple_p",pval_thresh,"_FC1_2.txt")), pval_thresh)
