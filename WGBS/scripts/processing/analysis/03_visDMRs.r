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
library(pheatmap)

#Directories
odcf.dir <-"/omics/groups/OE0219/internal/FetalHierachy/WGBS/methylDackel/methylationCalls/"
input.dir <-"/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/"
analysis.dir <- "/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/DMRAnal"

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "sig_dmrs_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "sig_dmrs_anno_reduced.rds"))

#colors
col = list(
    group= c("bonemarrow_p2"="#0058b4", "bonemarrow_p5"="#2188c9", 
    "fetalLiver_e12_5"="#fbbb25", "fetalLiver_e14_5"= "#a41220","fetalLiver_5"="#fca349", 
    "fetalLiver_p2"="#ff6b36", "fetalLiver_p5"="#e34e2e"),
    Tissue=c(bonemarrow="#ababab", fetalLiver="#99a637"),
    Age=c("p2"="#252525", "p5"="#737373", "e12_5"="#9babcf", "e14_5"="#99a637"))

#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
dmrs_final$all <- dmrs_red
mcols(dmrs_final$all)$direction <- "hypo"

#subset dmrs with a abs(difference) >0.1 
#dmrs_final_sub <- lapply(dmrs_final, function(x){
#    x <- x[abs(x$diff.Methy)>0.3,]
#    x
#})
#lapply(dmrs_final_sub, function(x)table(x$direction))

#export text files
for(i in names(dmrs_final)){
        dir.create(file.path(analysis.dir, i, "visualization"), recursive=TRUE)
    write.table(as.data.frame(dmrs_final[[i]]),file.path(analysis.dir,i, paste0(i,"_","dmrs_final.txt")),row.names = FALSE, quote=FALSE, sep="\t")
}

#width of dmrs
for(i in names(dmrs_final)){
    temp <- width(dmrs_final[[i]])
    pdf(file.path(analysis.dir, i, "visualization", paste0("DMR_width_histo",i, ".pdf")), height=4, width=4)
    print(gghistogram(as.data.frame(temp), x="temp", fill="grey", add_density=TRUE,add="mean",bins=50, rug=TRUE, xlab="Width of DMRs [bp]", ylab="# of DMRs")+xscale("log2"))
    dev.off()
    print(mean(temp))
    dmrs_final[[i]]$comparison <- i
}
width <- lapply(dmrs_final, function(x){
    width <- as.data.frame(width(x))
    width$comparison <- x$comparison
    width
})
width$all <- NULL
width <- do.call("rbind", width)
colnames(width)<-c("width", "comparison")
pdf(file.path(analysis.dir,"all","visualization", "DMr_width_allComparison.pdf"), height=4, width=7)
gghistogram(width, x="width", color="comparison", add="mean",add_density=TRUE,bins=50, rug=TRUE, xlab="Width of DMRs [bp]", ylab="# of DMRs", palette="jco",legend="right") +xscale("log2")
dev.off()

#pca of all original dmrs
#PCA
for(i in names(dmrs_final)){
    
    colnames( mcols(dmrs_final[[i]])) <- gsub(".", "-",colnames( mcols(dmrs_final[[i]])), fixed=TRUE )
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    
    meth_dmr <- meth_dmr[complete.cases(meth_dmr),]
    ir.pca <- prcomp(t(as.matrix(meth_dmr)),
                 center = T,
                 scale. = F) 
    
    print(i)
    print(summary(ir.pca))
    pc <- ir.pca$x
    pc<- as.data.frame(cbind(pc , pData(bsseq_all)))


    for(j in names(col)){
        pdf(file.path(analysis.dir, i, "visualization", paste0("PC12_",j,".pdf")),height = 5, width = 5)
        print(ggscatter(pc, x="PC1", y="PC2",
                color = j, shape = "Replicate",
                ellipse = F , mean.point = FALSE,palette= col[[j]],
                star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
                "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
                theme(legend.position="right",legend.title = element_text(, size=10, 
                                            face="bold")))
        dev.off()
        pdf(file.path(analysis.dir, i, "visualization",paste0("PC23_",j,".pdf")),height = 5, width = 5)
        print(ggscatter(pc, x="PC2", y="PC3",
                color = j, shape = "Replicate",
                ellipse = F , mean.point = FALSE,palette= col[[j]],
                star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
                "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
                theme(legend.position="right",legend.title = element_text(, size=10, 
                                            face="bold")))
        dev.off()
        pdf(file.path(analysis.dir, i, "visualization", paste0("PC34_",j,".pdf")),height = 5, width = 5)
        print(ggscatter(pc, x="PC3", y="PC4",
                color = j, shape = "Replicate",
                ellipse = F , mean.point = FALSE,palette= col[[j]],
                star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,3]*100,2), 
                "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
                theme(legend.position="right",legend.title = element_text(, size=10, 
                                            face="bold")))
        dev.off()
    print(i)
}
}


#Sample Clustering
for( i in names(dmrs_final)){
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    drld <- dist(t(as.matrix(meth_dmr)))
    hrld<- hclust(drld)
    dend<- hrld%>% as.dendrogram 
    anno <- pData(bsseq_all)
    #col1 <- hrld$labels
    #names(col1)<- anno[col1,]$group
    #col1 <- col1[order.dendrogram(dend)]
    #col1 <- ifelse(names(col1)=="TAM_tumor", "#EFC00099", ifelse(names(col1)=="BMDM_tumor","#0073C299",  ifelse(names(col1)=="TAM_healthy","#CD534CFF","#86868699")))
    dend <- dend %>% 
    #set("branches_k_color", k=4, c("#EFC00099", "#0073C299", "#CD534CFF","#86868699" )) %>%
    set("branches_lwd", 2) %>%
    #set("labels_colors",col1) %>% 
    set("labels_cex", .6 )%>%
    set("leaves_pch", 19)%>% 
    set("leaves_cex", 1.5)#%>% 
    #set("leaves_col", col1)
    
    pdf(file.path(analysis.dir, i, "visualization", "Clustering_DMR_meth.pdf"), height = 6, width = 5)
    print(dend %>% plot)
    dev.off()

}


#Heatmap of DMRs
#get contrasts
set.seed(42)
annovst <- as.data.frame(colData(bsseq_all))[, names(col)] 

for(i in names(dmrs_final)){
    #all dmrs in each comparison with all samples
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    rownames(meth_dmr)<- paste0("Dmr.",1:nrow(meth_dmr),"_", dmrs_final[[i]]$SYMBOL)
    meth_dmr <- meth_dmr[complete.cases(meth_dmr),]
    # pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
    #     clustering_distance_rows= "manhattan", clustering_method ="ward.D2",
    #     scale="row",fontsize_row=5,  annotation_color=col,
    #     filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_rowScale_WarD2_Manhatten.pdf"))
    pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        clustering_distance_rows= "manhattan", clustering_method ="ward.D2",
        scale="none",fontsize_row=5,  annotation_color=col,
        filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_noScale_WarD2_Manhatten.pdf"))  
    # pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
    #     clustering_distance_rows= "manhattan", clustering_method ="ward.D2", clustering_distance_cols= "manhattan",
    #     scale="row",fontsize_row=5,  annotation_color=col,
    #     filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_rowScale_WarD2_Manhatten_both.pdf"))
    pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        clustering_distance_rows= "manhattan", clustering_method ="ward.D2",clustering_distance_cols= "manhattan", 
        scale="none",fontsize_row=5,  annotation_color=col,
        filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_noScale_WarD2_Manhatten_both.pdf"))
    print(i)
}

#location of dmrs
library(circlize)
for(i in names(dmrs_final))[!names(dmrs_final) %in% c("all")]){
    dmrs_final_df <- as.data.frame(dmrs_final[[i]])
    pdf(file.path(analysis.dir,i,"visualization", paste0("DMRlocation_CirclePlot",".pdf")))
    circos.initializeWithIdeogram(species = "mm10", chromosome.index = paste0("chr", c(1:19, "X", "Y")))
    #bed_list_en = list(dmrs_final_df[[i]][dmrs_final_df[[i]]$direction =="hypo",], dmrs_final_df[[i]][dmrs_final_df[[i]]$direction =="hyper",])
    #circos.genomicRainfall(bed_list_en, pch = 16, cex = 0.8, col = c("#0000FF80", "#FF000080"))
    circos.genomicDensity(dmrs_final_df[dmrs_final_df$direction =="hypo",], col = c("#FF000080"), track.height = 0.1)
    circos.genomicDensity(dmrs_final_df[dmrs_final_df$direction =="hyper",], col = c("#0000FF80"), track.height = 0.1)
    dev.off()
    circos.clear()
    print(i)
}

# #pie plot of dmr dsistribution
# #define colors
# #col <- brewer.pal(3,"RdBu")
# col <- c("red","blue")
# #For loop
# for (comp in names(dmrs_final)){
# # calculate Distribution
# table <- as.data.frame(table(dmrs_final[[comp]]$direction))
# labs <- paste0(table$Var1,"\n(", round((table$Freq/sum(table$Freq)*100)),"%)\n", table$Freq)
# pdf(file.path(analysis.dir,comp,"visualization",paste0("Direction_Pie_DMRs.pdf")), height=3.5, width=3.5)
# print(ggpie(table, "Freq", fill="Var1", palette=col, label=labs, lab.pos = "in", main="DMR Analysis",submain=comp,  lab.font = c(5, "bold", "white")) + rremove("legend"))
# dev.off()
# #stratified barplot
# dis <- as.data.frame(dmrs_final[[comp]]$diff.Methy)
# colnames(dis)<- "dis"
# df <- data.frame( ">-100"=sum(dis< (-100)),  "-10 to -100"=sum(dis<(-10) & dis>(-100)), "-1 to -10"=sum(dis<(-1) & dis>(-10)), "0 to -1"=sum(dis<0 & dis>(-1)),
# "0 to 1"=sum(dis>0 & dis<1), "1 to 10"=sum(dis>1 & dis<10), "10 to 100"=sum(dis>10 & dis<100), ">100"=sum(dis<100))
# colnames(df)<- c( ">-100","-10 to -100",  "-1 to -10","0 to -1", "0 to 1" , "1 to 10", "10 to 100", ">100")
# df <- data.frame(distance= colnames(df), quant=t(df)[,1])
# df$quant <- df$quant/sum(df$quant)
# pdf(file.path(analysis.dir,comp,"visualization", paste0("Distance_to_TSS_barplot.pdf")), height=3.5)
# print(ggbarplot(df, x="distance", y="quant",  fill= "gray", ylab="Number of DMRs", main="CpG analysis")+rotate_x_text(angle = 90))
# dev.off()
# #Difference in Methylation
# # df <- cut(dis$dis, breaks=c((-0.8), (-0.6) ,(-0.4), (-0.2),0, 0.2,0.4,0.6,0.8), labels=c("-0.8 to -0.6","-0.6 to -0.4", "-0.4 to -0.2","-0.2 to -0","0 to 0.2","0.2 to 0.4","0.4 to 0.6","0.6 to 0.9" ))
# # df<-as.data.frame(t(table(df)))
# # df$Var1 <- c(rep("hypo",4), rep("hyper",4))
# # pdf(file.path(analysis.dir,comp,"visualization",paste0("DifferenceMeth_stratiefied.pdf")), height=3.5, width=3.5)
# # print(ggbarplot(df ,x="df", y="Freq", fill="Var1", palette=col, xlab="Difference in methylation", ylab="Number of DMRs") +rotate_x_text(angle = 45)+ rremove("legend"))
# # dev.off()
# }