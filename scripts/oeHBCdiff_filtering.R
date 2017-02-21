rm(list=ls()); options(getClass.msg=FALSE)
library(scone)
library(Biobase)
library(BiocParallel)
library(optparse)
library(RColorBrewer)

colpal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(9, "Set3"), brewer.pal(6, "Dark2"))

option_list <- list(
  make_option("--expt", default="", type="character", help="full form, e.g. Expt1"),
  make_option("--norm", default=NULL, type="character", help="normalization for excluded samples, if given")
)

opt <- parse_args(OptionParser(option_list=option_list))
expt_str <- opt$expt

register(SerialParam())

out_dir <- paste0("../output/clust/",expt_str)

# Load data and extract single cells
load("../data/oeHBCdiff_Cufflinks_eSet.Rda")
eSet <- Cufflinks_eSet

#####-------- Filter out contaminants 
### 1. first run: use NULL to choose get normalized counts from all samples so thresholds can be set
### 2. second run, run with populated excluded_samples_list

excluded_samples_list <- NULL
if (!is.null(opt$norm)){
  excluded_samples_list <<- load(paste0("../ref/",expt_str,"_",opt$norm,"_exclude.Rda"))
  message("using sample to exclude list")
}

message(paste(dim(eSet)[1], "genes,", dim(eSet)[2], "samples"))

# Exclude cells that are known to be contaminants (nonsensory epithelium and neuron/sus doublets)
if(length(excluded_samples_list) > 0){
  seqIDsToExclude=vector()
  for (i in seq_along(excluded_samples_list)){
    print(excluded_samples_list[i])
    seqIDsToExclude <- append(seqIDsToExclude,get(excluded_samples_list[i]))
    print(length(seqIDsToExclude))
  }
  desiredSamples <- !(toupper(phenoData(eSet)$sample_sequencing_id) %in% toupper(seqIDsToExclude))
  eSet <<- eSet[, desiredSamples]
}
message(paste("Dimensions after dropping contaminants:", dim(eSet)[1], "genes,", dim(eSet)[2], "samples"))

# Drop redundant and uninformative QC metrics
qc <- protocolData(eSet)@data[,-(6:9)] #these QC metrics apply only to paired end sequencing so we ignore them
qc <- subset(qc, select = -MEDIAN_5PRIME_TO_3PRIME_BIAS) # redundant, as we include 5PRIME_BIAS and 3PRIME_BIAS

#####-----Gene and sample filtering

# Sample filtering: all genes 0 or NA, all QC metrics NA
is.failed <-
  apply(is.na(exprs(eSet)) | (exprs(eSet) == 0), 2, all) | #all genes 0 or NA
  apply(is.na(pData(protocolData(eSet))), 1, all) #all QC metrics NA
print(paste(sum(is.failed),"samples failed pipeline:"),quote = F)
write.table(colnames(eSet)[is.failed],file = paste0("../output/EDA/",expt_str,"/failed_preproc_list.txt"),row.names=F, col.names=F,quote = F)

eSet <- eSet[,!is.failed]
print(sprintf("Removed %d failed samples.", sum(is.failed)))

# Filtering of Transcripts: ERCCs removed
geneidx <- !grepl("^ERCC-", featureData(eSet)$Gene_Symbol) # drop ERCCs
eSet <- eSet[geneidx,]

message(paste("After removing ERCC:",dim(eSet)[1], "genes,", dim(eSet)[2], "samples"))

# Select only transcripts above 0 TPM in at least one sample and no NaN in any sample
is.expressed.sc <- rowMeans(exprs(eSet)) > 0
eSet <- eSet[which(is.expressed.sc),]

message(paste("After removing 0TPM:",dim(eSet)[1], "genes,", dim(eSet)[2], "samples"))

# Remove transcripts with NA counts
is.counts.na <- apply(assayData(eSet)$counts_table, 1, function(x) any(is.na(x)))
eSet <- eSet[which(!is.counts.na),]
print("Removed undetected transcripts.")

message(paste("After removing NA counts:",dim(eSet)[1], "genes,", dim(eSet)[2], "samples"))

ct <- assayData(eSet)$counts_table
qc <- qc[colnames(ct),]

message(paste("Dimensions after dropping undetected transcripts:", dim(eSet)[1], "genes,", dim(eSet)[2], "samples"))

#####----- Sample/Gene Filtering

# Params for gene filtering
thresh_fail <-  40  # Tophat counts
num_fail <- 5  # cells

hk=intersect(read.table("../ref/hklist.txt")$V1,rownames(ct))

# Initial Gene Filtering
init.gf.vec <- rowSums(ct > thresh_fail) > num_fail

message(paste("Number of genes in initial gene filter:", sum(init.gf.vec)))

# Metric based filtering
pdf(file=file.path("../output/EDA",expt_str,paste0(expt_str,"_mfilt.pdf")), width=11, height=8)
mfilt <- metric_sample_filter(ct,mixture = FALSE, plot = TRUE,hist_breaks = 20,zcut = 4, gene_filter = NULL,nreads = qc$NREADS,ralign = qc$RALIGN, hard_ralign=85, suff_ralign = NULL, suff_nreads = 100000, pos_controls = NULL)
dev.off()

mfilt <- !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)
mfilt.gf.vec <- rowSums(ct[,mfilt] > thresh_fail) > num_fail

qc <- qc[mfilt,]
eSet <- eSet[mfilt.gf.vec,mfilt] 

eSet.gf.vec <- rowSums(assayData(eSet)$counts_table > thresh_fail) > num_fail
eSet <- eSet[eSet.gf.vec,]

message(paste("Dimensions of eSet post counts-in-cells:",dim(eSet)[1], "genes,", dim(eSet)[2], "samples"))

expt <- as.character(eSet$MD_expt_condition)
expt[grep("K5ERRY_UI_72HPT",expt)] <- "K5ERRY_UI"
expt[grep("K5ERRY_UI_96HPT",expt)] <- "K5ERRY_UI"
expt <- factor(expt, levels=c("K5ERRY_UI", "K5ERP63CKO_UI_24HPT", "K5ERP63CKO_UI_48HPT", "K5ERP63CKO_UI_96HPT", "K5ERP63CKO_UI_7DPT", "K5ERP63CKO_UI_14DPT", "SOX2EGFP+ICAM-F3-SRB1-")) #to order levels
batch <- eSet$MD_c1_run_id

counts <- assayData(eSet)$counts_table

save(counts, qc, batch, expt, file = file.path(out_dir, paste0(expt_str,"_filtdata.Rda")))