rm(list=ls()); options(getClass.msg=FALSE)

library(biomaRt)
library(parallel)
library(limma)
library(optparse)
library(clusterExperiment)
printf <- function(...) cat(sprintf(...))

option_list <- list(
  make_option("--nrot", default="1", type="double")
)

opt <- parse_args(OptionParser(option_list=option_list))
nrot <- opt$nrot

load("../output/clust/oeHBCdiff/MSigDBselectedGeneSets.Rdata")
load("../output/clust/oeHBCdiff/oeHBCdiff_finalClusterObject.Rda")

lognorm <- transform(cmobj2)
lognorm <- lognorm[bm[,1],]
cont <- clusterContrasts(cmobj2, contrastType="OneAgainstAll")$contrastMatrix
X <- as.factor(primaryCluster(cmobj2))
design <- model.matrix(~X - 1)
head(design)

coresToUse <- detectCores()-1
lapply(list(msigdbSets), 
       function(geneSets) {
         sprintf("Testing romer for gene set family %s (%d gene sets)\n", geneSets$description, length(geneSets$idx))
         
         message(system.time(res <- mclapply(seq_len(NCOL(cont)), function(i) romer(lognorm, geneSets$idx, design, cont[,i],nrot=nrot), mc.cores = coresToUse)))
         names(res) <- colnames(cont)
         
         
         message(length(res))
         lapply(seq_along(res), function(i) {message(typeof(res[[i]]));message(names(res[[i]])) })           
         res <- lapply(seq_along(res), 
                       function(i) { ord = order(res[[i]][, "Mixed"]); res[[i]][ord, ] } )
         
         system("mkdir -p ../output/oeHBCdiff/romer")
         save(res, design, cont, file=sprintf("../output/oeHBCdiff/romer/romer_res_%s.rda", geneSets$description))
         
         lapply(seq_along(res), 
                function(i) {
                  contrast_name = gsub("\\)|\\(|\\+", '', gsub('/', '_', colnames(cont)[i]), perl=TRUE)
                  
                  write.table(res[[i]], file=sprintf("../output/oeHBCdiff/romer/%s_contrast%d_%s.txt", geneSets$description, i, contrast_name), quote=FALSE, sep='\t')
                  
                  
                  #more user-friendly output
                  stopifnot(all(colnames(res[[i]]) == c("NGenes", "Up", "Down",  "Mixed")))
                  res[[i]] = as.data.frame(cbind(res[[i]], NA))
                  colnames(res[[i]])[length(colnames(res[[i]]))] =  "genes_in_gene_set"
                  for(j in seq_len(nrow(res[[i]]))) {
                    curGenes = geneSets$genes[[rownames(res[[i]])[j]]]
                    res[[i]][j, "genes_in_gene_set"] = paste(curGenes, collapse = '; ')
                  }
                  
                  colnames(res[[i]]) = c("#genes_in_set", "pval_upregulated", "pval_downregulated", "pval_DE", "genes_in_gene_set")
                  write.csv(res[[i]], file=sprintf("../output/oeHBCdiff/romer/%s_contrast%d_%s.csv", geneSets$description, i, contrast_name), row.names = TRUE)
                })
         
       })
