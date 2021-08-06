##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##Prepare Homer input and run in comman line

#Libraries
library(DSS)
library(bsseq)
library(doParallel)
library(ChIPseeker)
library(foreach)
library(SummarizedExperiment)
library(doMC)
library(rtracklayer)
library(HDF5Array)
library(org.Mm.eg.db)
library(ggpubr)
library(randomcoloR)
library(RColorBrewer)
library(VennDiagram)
library(ChIPpeakAnno)
library(dendextend)

#Directories
odcf.dir <-"/omics/groups/OE0219/internal/FetalHierachy/WGBS/methylDackel/methylationCalls/"
input.dir <-"/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/"
analysis.dir <- "/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/DMRAnal"
datasets.dir <- "c010-datasets/Internal/COPD/enrichment_databases/"

#load data
#bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "sig_dmrs_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "sig_dmrs_anno_reduced.rds"))

#assign direction
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
#export lists
dmrs_final_df<-list()
for(i in names(dmrs_final)){
    dmrs_final_df[[i]] <- as.data.frame(dmrs_final[[i]])
    dmrs_final_df[[i]]$peak_id <- 1:nrow(dmrs_final_df[[i]])

    dir.create(file.path(analysis.dir, "homer",i, "all"), recursive=TRUE)
    dir.create(file.path(analysis.dir, "homer",i, "hypo"), recursive=TRUE)
    dir.create(file.path(analysis.dir, "homer",i, "hyper"), recursive=TRUE)
    #all
    write.table(data.frame(chr=dmrs_final_df[[i]]$seqnames, start=dmrs_final_df[[i]]$start, end=dmrs_final_df[[i]]$end, id=dmrs_final_df[[i]]$peak_id, 
        notUsed=NA, strand="+"),
        file.path(analysis.dir,  "homer", i, "all", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    #hypo
    write.table(data.frame(chr=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hypo",]$seqnames, 
        start=dmrs_final_df[[i]][dmrs_final[[i]]$direction=="hypo",]$start, end=dmrs_final_df[[i]][dmrs_final[[i]]$direction=="hypo",]$end, 
        id=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hypo",]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer", i,"hypo", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    #hyper
    write.table(data.frame(chr=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$seqnames, 
        start=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$start, end=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$end, 
        id=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer", i,"hyper", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
}
dmrs_red_df <- as.data.frame(dmrs_red)
dmrs_red_df$peak_id <- 1:nrow(dmrs_red_df)
write.table(data.frame(chr=dmrs_red_df$seqnames, 
        start=dmrs_red_df$start, end=dmrs_red_df$end, 
        id=dmrs_red_df$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer", paste0("DMRs_bg.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

#run in command line
cd /omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/DMRAnal/homer
chmod 777 -R ./
conda activate homer2

for file in `ls */*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} mm10 ${path} -size given -preparsedDir ${path}/ -p 10
    echo ${path}
    echo $file
done


#with custom bg
#cd /omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/DMRAnal/homer
chmod 777 -R ./
#conda activate homer2

for file in `ls */*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} mm10 ${path}_customBG -size given -preparsedDir ${path}/ -p 10 -bg /omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/DMRAnal/homer/DMRs_bg.bed
    echo ${path}
done
