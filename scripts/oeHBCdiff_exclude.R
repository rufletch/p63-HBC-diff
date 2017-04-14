#title: Excluding Contaminant Cells
#author: Diya Das and Russell Fletcher
#date: January 30th 2017

rm(list=ls()); options(getClass.msg=FALSE)

library(clusterExperiment)
library(optparse)

option_list <- list(
  make_option("--expt", default="", type="character", help="full form, e.g. Expt1"),
  make_option("--ncores", default="1", type="double"),
  make_option("--nrm", type="character")
)

opt <- parse_args(OptionParser(option_list=option_list))
expt_str <- opt$expt
normstr <- opt$nrm

BiocParallel::register(MulticoreParam(workers = opt$ncores))

clust_dir <- paste0("../output/clust/",expt_str)

load(file.path(clust_dir,paste0(expt_str,"_",normstr,"_se.Rda")))

counts <- assay(se)

omp_cyp <- counts["Omp",] > 200 & (counts["Cyp2g1",] > 200 | counts["Cyp1a2",] > 200)
table(colData(se)$expt[omp_cyp]) 

omp_cyp_list <- colnames(counts)[omp_cyp]

reg3g_list <- colnames(counts)[counts["Reg3g",]>=100]

save(reg3g_list,omp_cyp_list,file=file.path("../ref",paste0(expt_str,"_",normstr,"_exclude.Rda")))
