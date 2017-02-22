rm(list=ls()); options(getClass.msg=FALSE)
library(scone)
library(BiocParallel)
library(optparse)

option_list <- list(
  make_option("--expt", default="", type="character", help="full form, e.g. Expt1"),
  make_option("--ncores", default="1", type="double")
)

opt <- parse_args(OptionParser(option_list=option_list))
expt_str <- opt$expt

register(MulticoreParam(workers = opt$ncores))

clust_dir <- paste0("../output/clust/",expt_str)

set.seed(1999)
load(file.path(clust_dir, paste0(expt_str,"_filtdata.Rda")))
expt <- droplevels(expt)
batch <- droplevels(batch)

hk <- read.table(file.path("../ref", "hklist.txt"))
hk <- intersect(rownames(counts), hk[,1])
del <- read.table(file.path("../ref", "oeHBCdiff_de.txt"))
del <- intersect(rownames(counts), del[,1])
cc <- read.table(file.path("../ref","cell_cycle.txt"))
cc <- intersect(rownames(counts), cc[,1])

# Generate Scores and Ranking
print(system.time({
scone_out <- scone(counts, imputation=list(none=impute_null), impute_args=list(0), return_norm="in_memory", scaling=list(none=identity, fq=FQT_FN, deseq=DESEQ_FN, tmm=TMM_FN), k_ruv=3, k_qc=3, ruv_negcon=hk, qc=as.matrix(qc), adjust_bio="yes", bio=expt, adjust_batch="yes", batch=batch,run=TRUE, evaluate=TRUE, eval_negcon=cc, eval_poscon=del, eval_kclust = 10:12)
}))

save(scone_out, file = file.path(clust_dir,paste0(expt_str,"_scone.Rda")))