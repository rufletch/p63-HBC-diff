rm(list=ls()); options(getClass.msg=FALSE)
library(BiocParallel)
library(clusterExperiment)
library(optparse)

option_list <- list(
  make_option("--expt", default="", type="character", help="full form, e.g. Expt1"),
  make_option("--ncores", default="1", type="double"),
  make_option("--nrm", type="character")
)

opt <- parse_args(OptionParser(option_list=option_list))
expt_str <- opt$expt

register(MulticoreParam(workers = opt$ncores))

clust_dir <- paste0("../output/clust/",expt_str)

seed=927501

nrmstr <- opt$nrm
load(file.path(out_dir,paste0(expt_str,"_",nrmstr,"_se.Rda")))

beta <- 0.8
minSize <- 1

cmobj <- clusterMany(se, dimReduce="PCA", nPCADims=50, clusterFunction="hierarchical01", 
alphas=0.1, betas=beta, minSizes=minSize, subsample=TRUE, sequential=TRUE, ncores=detectCores(), 
subsampleArgs=list(resamp.num=100,ncores=1,clusterFunction="kmeans",clusterArgs=list(nstart=10)),
seqArgs=list(k.min=3, top.can=5), random.seed=seed, run=TRUE, ks = 4:15, verbose=TRUE, isCount=TRUE)

save(cmobj,beta,minSize, seed, file=file.path(out_dir, paste0(expt_str,"_", nrmstr,"_cm.Rda")))

cm2 <- combineMany(cmobj, proportion=0.7, whichClusters = clusterLabels(cmobj))
cm2 <- makeDendrogram(cm2,dimReduce="var",ndims=5000, ignoreUnassignedVar=TRUE)
  
save(cm2, file=file.path(clust_dir,paste0(normstr,"_",clname,"_cons.Rda")))