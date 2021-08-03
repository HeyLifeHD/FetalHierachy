##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##Run Lola enrichment

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
library(LOLA)

#Directories
odcf.dir <-"/omics/groups/OE0219/internal/FetalHierachy/WGBS/methylDackel/methylationCalls/"
input.dir <-"/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/"
analysis.dir <- "/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/DMRAnal"
datasets.dir <- "c010-datasets/Internal/COPD/enrichment_databases/"

#load data
#bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "sig_dmrs_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "sig_dmrs_anno_reduced.rds"))

#load LOLA
datasets.dir <- "c010-datasets/Internal/COPD/enrichment_databases/"

#load LOLA
datasets.dir <- "c010-datasets/Internal/COPD/enrichment_databases/"
regionDB_genomicRegions <- loadRegionDB(file.path(datasets.dir, "mm10_genomic_regions"))
regionDB_homer <- loadRegionDB(file.path(datasets.dir, "homer/homer_lola"))
regionDB_chipSeq <- loadRegionDB(file.path(datasets.dir, "lola_chipseq"))
regionDB_msigdb <- loadRegionDB(file.path(datasets.dir, "msigdb"))

#Create Lists
results_Core <- list()
results_genomicRegions <- list()
results_homer <- list()
results_chipSeq <- list()
results_msigdb <- list() 

#Run Enrichment
dir.create(file.path(analysis.dir, "LOLA"))
for (i in names(dmrs_final)){
#stratify open and closed
hypo<- dmrs_final[[i]][ dmrs_final[[i]]$direction == "hypo",]
print(length(hypo))
hyper <- dmrs_final[[i]][ dmrs_final[[i]]$direction =="hyper",]
print(length(hyper))
userSets<- list(hypo, hyper)
names(userSets)<- c("hypo","hyper" )
#set  Universe
userUnisverse <-dmrs_red
#Run analysis
#results_Core[[i]]= runLOLA(userSets, userUniverse, regionDB_Core, cores=4)
results_genomicRegions[[i]]= runLOLA(userSets, userUnisverse, regionDB_genomicRegions, cores=3)
results_homer[[i]]= runLOLA(userSets, userUnisverse, regionDB_homer, cores=3)
results_chipSeq[[i]]= runLOLA(userSets, userUnisverse, regionDB_chipSeq, cores=3)
results_msigdb[[i]]= runLOLA(userSets, userUnisverse, regionDB_msigdb, cores=3)
print(i)
}
dir.create(file.path(analysis.dir, "LOLA"))
saveRDS(results_genomicRegions, file.path(analysis.dir, "LOLA", "results_genomicRegions.rds"))
saveRDS(results_homer, file.path(analysis.dir, "LOLA", "results_homer.rds"))
saveRDS(results_chipSeq, file.path(analysis.dir, "LOLA", "results_chipSeq.rds"))
saveRDS(results_msigdb, file.path(analysis.dir, "LOLA", "results_msigdb.rds"))

results_genomicRegions<- readRDS(file.path(analysis.dir, "LOLA", "results_genomicRegions.rds"))
results_homer<- readRDS(file.path(analysis.dir, "LOLA", "results_homer.rds"))
results_chipSeq <- readRDS(file.path(analysis.dir, "LOLA", "results_chipSeq.rds"))
results_msigdb <- readRDS(file.path(analysis.dir, "LOLA", "results_msigdb.rds"))
#Plot genomic regions in bubble plot
#function 
label_func <- function(x){
    breaks <- x
    breaks[breaks>=200] <- ">=200"
    breaks
}
#subset multi-cell
results <- lapply(results_genomicRegions, function(x) x[grep("multi-cell", x$filename),])
#plotting
dir.create(file.path(analysis.dir, "LOLA","genomicRegions"))
g<- list(NULL)

library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
col <- col[c(1,3)]

for(i in names(results)){
  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,"_"),`[`, 2)
combined_data$filename<-sapply(strsplit(combined_data$filename,"fe-"),`[`, 2)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir, "LOLA","genomicRegions", paste0("EnrichLOLA_",names(results[i]), "_genomicRegions.pdf")))
print(g[i])
dev.off()
}



#Plot homer results in bubble plot
#subset multi-cell
results <- results_homer
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="hypo",]$filename, 20), head(x[x$userSet=="hyper",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
dir.create(file.path(analysis.dir, "LOLA","homer"))
g<- list(NULL)

for(i in names(results)){
  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir, "LOLA","homer", paste0("EnrichLOLA_",names(results[i]), "_homer.pdf")))
print(g[i])
dev.off()
}



#Plot regionDB_chipSeq results in bubble plot
#subset multi-cell
results <- results_chipSeq
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="hypo",]$filename, 20), head(x[x$userSet=="hyper",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
dir.create(file.path(analysis.dir, "LOLA","ChipSeq"))
g<- list(NULL)

for(i in names(results)){
  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir, "LOLA","ChipSeq", paste0("EnrichLOLA_",names(results[i]), "_ChipSeq.pdf")), width=20)
print(g[i])
dev.off()
}




#For MsigDB
results <- results_msigdb
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="hypo",]$filename, 20), head(x[x$userSet=="hyper",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
dir.create(file.path(analysis.dir, "LOLA","MSigDB"))
g<- list(NULL)

for(i in names(results)){
  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir, "LOLA","MSigDB", paste0("EnrichLOLA_",names(results[i]), "_MSigDB.pdf")))
print(g[i])
dev.off()
}


#For MsigDB only hallmarks
results <- results_msigdb
results <- lapply(results, function(x){
    x <- x[x$collection =="hallmarks",]
    temp <- unique(c(head(x[x$userSet=="hypo",]$filename, 20), head(x[x$userSet=="hyper",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
dir.create(file.path(analysis.dir, "LOLA","MSigDBHallmarks"))
g<- list(NULL)

for(i in names(results)){
  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)
#adjust filename 
combined_data$filename <- gsub("_", " ", combined_data$filename)
combined_data$filename <- gsub("hallmark ", "", combined_data$filename)
combined_data$filename <- stringr::str_to_sentence(combined_data$filename)
temp <- split(combined_data, combined_data$filename)
temp <- unlist(lapply(temp, function(x){
    x <- mean(x$pValueLog)
    x
}))
temp <- temp[order(temp, decreasing=TRUE)]
combined_data$filename <- as.factor(combined_data$filename)
combined_data$filename <- ordered(combined_data$filename  , names(temp))

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir, "LOLA","MSigDBHallmarks", paste0("EnrichLOLAHallmarks_",names(results[i]), "_MSigDB_orderSign.pdf")))
print(g[i])
dev.off()
temp <- g[[i]]
ggsave(plot=temp, height=7, width=7,filename=file.path(analysis.dir, "LOLA","MSigDBHallmarks", paste0("EnrichLOLAHallmarks_",names(results[i]), "_MSigDB_orderSign_new.pdf")),useDingbats=FALSE)

}
