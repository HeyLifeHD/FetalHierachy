#libraries
library(methrix)
library(BSgenome.Mmusculus.UCSC.mm10)

#directories
input.dir <- "/omics/groups/OE0219/internal/FetalHierachy/WGBS/methylDackel/methylationCalls/"
output.dir <- "/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/"
dir.create(output.dir,recursive=TRUE)
#functions
# collapse_samples <- function(m, col=NULL){

#   if (!(col %in% colnames(m@colData))){
#     stop("The provided column name is not in the dataset. Please provide a name that is in the colData. ")
#   }
#   m <- m[!is.na(m@colData[,col])]
#   groupby <- factor(m@colData[,col])
#   sp <- split(seq(along = groupby), groupby)

#   if (methrix:::is_h5(m)){
#     stop("Not implemented yet.")
#     # new_cov <- sapply(sp, function(i) rowSums(m@assays[[2]][,i, drop = FALSE], na.rm=T))
#     # new_beta <- m@assays[[1]]*m@assays[[2]]
#     # new_beta <- sapply(sp, function(i) rowSums(new_beta[,i, drop = FALSE], na.rm=T))
#     # new_beta <- new_beta/new_cov
#   } else {
#     new_cov <- sapply(sp, function(i) matrixStats::rowSums2(m@assays@data@listData$cov[,i, drop = FALSE], na.rm=TRUE))
#     new_beta <- m@assays@data@listData$beta*m@assays@data@listData$cov
#     new_beta <- sapply(sp, function(i) matrixStats::rowSums2(new_beta[,i, drop = FALSE], na.rm=TRUE))
#     new_beta <- new_beta/new_cov
#   }

#   new_beta[new_cov==0] <- NA
#   new_cov[new_cov==0] <- NA
#   combined <- m[,unlist(lapply(sp, '[', 1))]
#   dimnames(combined)=dimnames(new_cov)

#   combined@assays@data@listData$beta <- new_beta
#   combined@assays@data@listData$cov <- new_cov

#   combined
# }
#load reference CpGs
mm10_cpgs = methrix::extract_CPGs(ref_genome = "BSgenome.Mmusculus.UCSC.mm10")

#load begraph paths
bdg_PBAT_files = dir(path = input.dir , recursive=TRUE,pattern = "chr.bedgraph$", full.names = TRUE)
length(bdg_PBAT_files)
#prepare sample sheet
#twgbs
SampleID_Lib<-sapply(strsplit(bdg_PBAT_files, "/", fixed=TRUE), "[",12)
SampleID_Lib <- gsub("fetal-liver", "fetalLiver",SampleID_Lib)
SampleID_Lib <- gsub("rep-", "",SampleID_Lib)
Tissue <- sapply(strsplit(SampleID_Lib, "_", fixed=TRUE), "[",1)
Celltype <-  sapply(strsplit(SampleID_Lib, "_", fixed=TRUE), "[",2)
Origin <-  sapply(strsplit(SampleID_Lib, "_", fixed=TRUE), "[",3)
Genotype <-sapply(strsplit(SampleID_Lib, "_", fixed=TRUE), "[",4)
Age <-sapply(strsplit(SampleID_Lib, "_", fixed=TRUE), "[",6)
Age<- gsub("-", "_", Age)
Replicate <-sapply(strsplit(SampleID_Lib, "_", fixed=TRUE), "[",7)
Replicate <- gsub("rep-", "",Replicate)
sample_anno <- data.frame(SampleID=paste0(Tissue, "_", Age, "_",Replicate), group=paste0(Tissue, "_",Age) ,Tissue=Tissue,
    Celltype=Celltype, Origin=Origin, Genotype=Genotype, Age=Age, Replicate=Replicate)
rownames(sample_anno)<- sample_anno$SampleID
saveRDS(sample_anno, file.path(output.dir, "sample_anno.rds"))
sample_anno <- readRDS(file.path(output.dir, "sample_anno.rds"))

#read begraphs
meth = methrix::read_bedgraphs(files = bdg_PBAT_files, ref_cpgs = mm10_cpgs, 
    chr_idx = 1, start_idx = 2, beta_idx=4 ,M_idx = 5, U_idx = 6,
    stranded = TRUE, collapse_strands = TRUE,
    coldata=sample_anno,
    n_threads=3 #, h5=TRUE, h5_dir=file.path(output.dir, "190923_methrix")
    )
dir.create(file.path(output.dir ,"methrix"))
saveRDS(meth, file.path(output.dir ,"methrix", "methrix_original.rds"))
#collapse data
#meth_collapsed <- collapse_samples(meth, col="SampleID")
#saveRDS(meth_collapsed, file.path(output.dir ,"methrix",  "methrix.rds"))

#QC report
#meth <- readRDS( file.path(output.dir ,"methrix",  "methrix.rds"))
dir.create(file.path(output.dir, "methrix",  "QC_report_noFilter"))
methrix::methrix_report(meth = meth, recal_stats=TRUE,output_dir = file.path(output.dir, "methrix",  "QC_report_noFilter"))

# 1.2: mask low coverage CpGs
meth_fil <- meth
meth_fil <- methrix::mask_methrix(m = meth_fil, low_count=1, high_quantile = 0.99)
saveRDS(meth_fil, file.path(output.dir ,"methrix",  "methrix_masked.rds"))
# 1.3: filter out uncovered CpGs
meth_fil <- methrix::remove_uncovered(m = meth_fil)
saveRDS(meth_fil, file.path(output.dir , "methrix_filter.rds"))
#QC report
dir.create(file.path(output.dir, "methrix", "QC_report_Filter"))
methrix::methrix_report(meth = meth_fil, output_dir = file.path(output.dir, "methrix", "QC_report_Filter"),recal_stats=TxRUE)

#export as bsseq
bsseq_all <-methrix::methrix2bsseq(meth)
bsseq_all_sub <-methrix::methrix2bsseq(meth_fil)
dir.create(file.path(output.dir, "bsseq"))
saveRDS(bsseq_all_sub, file.path(output.dir , "bsseq", "bsseq_all_fil.rds"))
saveRDS(bsseq_all, file.path(output.dir , "bsseq", "bsseq_all.rds"))
