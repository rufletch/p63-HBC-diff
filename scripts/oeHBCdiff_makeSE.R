rm(list=ls()); options(getClass.msg=FALSE)
library(scone)
library(BiocParallel)
library(MASS)
library(RUVSeq)
library(clusterExperiment)
library(optparse)
library(SummarizedExperiment)

option_list <- list(
  make_option("--expt", default="", type="character", help="full form, e.g. Expt1"),
  make_option("--ncores", default="1", type="double"),
  make_option("--nrm", type="character")
)

opt <- parse_args(OptionParser(option_list=option_list))
expt_str <- opt$expt

register(MulticoreParam(workers = opt$ncores))

clust_dir <- paste0("../output/clust/",expt_str)

scone_out <- get(load(file.path(clust_dir, paste0(expt_str,"_scone.Rda")))[1])

load(file.path(clust_dir,paste0(expt_str,"_filtdata.Rda")))

nrm <- opt$nrm
nrmstr <- gsub(",","_",gsub("_k=|_","",nrm))
scone_out <- scone_out$normalized_data[[nrm]]
se <- SummarizedExperiment(list(counts = expm1(scone_out)),
                           colData=data.frame(batch=batch, expt=expt,  nreads=qc$NREADS, ralign=qc$RALIGN, pctmrna=qc$PCT_MRNA_BASES))
save(se, file = file.path(clust_dir,paste0(expt_str,"_",nrmstr,"_se.Rda")))
