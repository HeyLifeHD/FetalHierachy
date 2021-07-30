##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##Exploratory analysis

#libraries
library(bsseq)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggpubr)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library("EnsDb.Mmusculus.v79")
library(RnBeads.mm10)

#directories
odcf.dir <-"/omics/groups/OE0219/internal/FetalHierachy/WGBS/methylDackel/methylationCalls/"
input.dir <-"/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/"
analysis.dir <- "/omics/groups/OE0219/internal/FetalHierachy/WGBS/210729_analysis/ExplorAnal"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all.rds"))

#colors
col = list(
    group= c("bonemarrow_hsc_pr-fetal_wt_p2"="#0058b4", "bonemarrow_hsc_pr-fetal_wt_p5"="#2188c9", 
    "fetalLiver_hsc_pr-fetal_wt_e12-5"="#fbbb25", "fetalLiver_hsc_pr-fetal_wt_e14-5"="#fca349", 
    "fetalLiver_hsc_pr-fetal_wt_p2"="#ff6b36", "fetalLiver_hsc_pr-fetal_wt_p5"="#e34e2e"),
    Tissue=c(bonemarrow="#ababab", fetalLiver="#99a637"),
    Age=c("p2"="#252525", "p5"="#737373", "e12-5"="#9babcf", "e14-5"="#99a637"))

#Plot average Methylation 
## extract the methylation values
pheno <- colData(bsseq_all)
meth_per_cpg <- bsseq::getMeth(bsseq_all, type = "raw")
meth_per_cpg <- as.data.frame(meth_per_cpg)
av.meth_per_cpg<- data.frame(AverageMethylation=colMeans(meth_per_cpg, na.rm = TRUE))
av.meth_per_cpg <- as.data.frame(cbind(av.meth_per_cpg, pheno))
av.meth_per_cpg

#plotting
for(i in names(col)){
    #prepare plotting
    #compare_means(AverageMethylation ~ group,  data = av.meth_per_cpg, method = "t.test")
    med <- aggregate(av.meth_per_cpg[,"AverageMethylation",drop=FALSE], list(av.meth_per_cpg[,colnames(av.meth_per_cpg)==i]), median)
    ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
    #plotting
    pdf(file.path(analysis.dir, paste0("AvMeth_",i,".pdf")), height=4, width=4)
    print(ggpubr::ggboxplot(av.meth_per_cpg, x=i, y="AverageMethylation", color=i,ylim=c(0,1), order=ord,
        palette =col[[i]],
        add = "jitter")  +  
        #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
        #stat_compare_means(label.y = .25) +
        rremove("legend")+
        rremove("xlab")+
        rotate_x_text(angle = 45))
    dev.off()
}

#Plot average Methylation of CpGislans
#cpg locations
cgi <- unlist(rnb.get.annotation("cpgislands",assembly="mm10"))

## extract the methylation values
meth_per_cpg_cpgi <- bsseq::getMeth(bsseq_all, type = "raw", what="perRegion", regions=cgi)
av.meth_per_cpgi <- data.frame(AverageMethylation=colMeans(meth_per_cpg_cpgi, na.rm = TRUE))
av.meth_per_cpgi <- as.data.frame(cbind(av.meth_per_cpgi , pheno))
av.meth_per_cpgi 

#plotting
for(i in names(col)){
    #prepare plotting
    #compare_means(AverageMethylation ~ group,  data = av.meth_per_cpg, method = "t.test")
    med <- aggregate(av.meth_per_cpgi [,"AverageMethylation",drop=FALSE], list(av.meth_per_cpgi [,colnames(av.meth_per_cpgi )==i]), median)
    ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
    #plotting
    pdf(file.path(analysis.dir, paste0("AvMethCpGisland_",i,".pdf")), height=4, width=4)
    print(ggpubr::ggboxplot(av.meth_per_cpgi, x=i, y="AverageMethylation", color=i,ylim=c(0,1), order=ord,ylab="Average CpGisland methylation",
        palette =col[[i]],
        add = "jitter")  +  
        #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
        #stat_compare_means(label.y = .25) +
        rremove("legend")+
        rremove("xlab")+
        rotate_x_text(angle = 45))
    dev.off()
}

#Heatmap of most variable CpGs
#meth_per_cpg <- bsseq::getMeth(bsseq_all, type = "raw")
topVarCpGs<- head(order(rowVars(as.matrix(meth_per_cpg)), decreasing=TRUE),20000)
meth_per_cpg<- as.data.frame(meth_per_cpg)
matrld <- meth_per_cpg[topVarCpGs,]
annovst <- as.data.frame(colData(bsseq_all))[, names(col)] 

pheatmap(matrld,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,
    show_colnames=FALSE, scale="row",fontsize_row=5,  annotation_colors=col,
    filename=file.path(analysis.dir,"Heatmap20000mvCpG.pdf"))
pheatmap(matrld,  annotation_col=as.data.frame(annovst), annotation_colors=col,
    show_rownames=FALSE,show_colnames=FALSE, scale="none",fontsize_row=5,  
    filename=file.path(analysis.dir,"Heatmap20000mvCpG_noScale.pdf"))



#Heatmap of most variable average Promoter Methylation
#Heatmap on Prtomoter
ENS<-EnsDb.Mmusculus.v79
seqlevelsStyle(ENS) <- "UCSC"
GENE = genes(ENS)

#1500 upstream and 500 downstream as standard
#20000mv
Promoter<- promoters(GENE, upstream = 1500, downstream = 500)
PromoterMeth <- bsseq::getMeth(bsseq_all, regions=Promoter, what="perRegion", type="raw")
row.names(PromoterMeth)<-Promoter$symbol
PromoterMeth<- as.matrix(PromoterMeth)
topVarGenesRld<- head(order(rowVars(PromoterMeth), decreasing=TRUE),20000)
anno <- Promoter
anno<-as.data.frame(anno)
annotop <- anno[topVarGenesRld,]
matrld <- PromoterMeth[topVarGenesRld,]
#plot
pheatmap(matrld,  annotation_col=as.data.frame(annovst),
    show_rownames=FALSE,show_colnames=FALSE, scale="row",fontsize_row=5, annotation_colors=col,
    filename=file.path(analysis.dir,"Heatmap_20000mvPromoterMeth.pdf"))
pheatmap(matrld,  annotation_col=as.data.frame(annovst),
    show_rownames=FALSE,show_colnames=FALSE, scale="none",fontsize_row=5, annotation_colors=col,  
    filename=file.path(analysis.dir,"Heatmap_20000mvPromoterMeth_noScale.pdf"))
#1000 mv
topVarGenesRld<- head(order(rowVars(PromoterMeth), decreasing=TRUE),1000)
anno <- Promoter
anno<-as.data.frame(anno)
annotop <- anno[topVarGenesRld,]
matrld <- PromoterMeth[topVarGenesRld,]
pheatmap(matrld,  annotation_col=as.data.frame(annovst), annotation_colors=col,
    show_rownames=FALSE,show_colnames=FALSE, scale="row",fontsize_row=5, 
    filename=file.path(analysis.dir,"Heatmap_1000mvPromoterMeth.pdf"))
pheatmap(matrld,  annotation_col=as.data.frame(annovst), annotation_colors=col,
    show_rownames=FALSE,show_colnames=FALSE, scale="none",fontsize_row=5,  
    filename=file.path(analysis.dir,"Heatmap_1000mvPromoterMeth_noScale.pdf"))


#PCA mv Promoter
#meth_per_cpg <- as.data.frame(bsseq::getMeth(bsseq_all, type = "raw"))
topVarGenesRld<- head(order(rowVars(PromoterMeth), decreasing=TRUE),20000)
anno <- Promoter
anno<-as.data.frame(anno)
annotop <- anno[topVarGenesRld,]
matrld <- PromoterMeth[topVarGenesRld,]
ir.pca <- prcomp(t(matrld),
                 center = TRUE,
                 scale. = TRUE) 
summary(ir.pca)
x <- ir.pca$x
x<- as.data.frame(cbind(x , pheno))

for(i in names(col)){
pdf(file.path(analysis.dir, paste0("PC12_20000Promoter_subsetted_",i,".pdf")),height = 5, width = 5)
print(ggscatter(x, x="PC1", y="PC2",
          color = i, shape = "Replicate",
          ellipse = F , mean.point = FALSE,palette= col[[i]],
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold")))
dev.off()
pdf(file.path(analysis.dir, paste0("PC23_20000Promoter_subsetted_",i,".pdf")),height = 5, width = 5)
print(ggscatter(x, x="PC2", y="PC3",
          color = i, shape = "Replicate",
          ellipse = F , mean.point = FALSE,palette= col[[i]],
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold")))
dev.off()
pdf(file.path(analysis.dir, paste0("PC34_20000Promoter_subsetted_",i,".pdf")),height = 5, width = 5)
print(ggscatter(x, x="PC3", y="PC4",
          color = i, shape = "Replicate",
          ellipse = F , mean.point = FALSE,palette= col[[i]],
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,3]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold")))
dev.off()
print(i)
}


#PCA mv CpG
topVarCpGs<- head(order(rowVars(as.matrix(meth_per_cpg)), decreasing=TRUE),20000)
ir.pca <- prcomp(t(meth_per_cpg[topVarCpGs,]),
                 center = TRUE,
                 scale. = TRUE) 
summary(ir.pca)
x <- ir.pca$x
x<- as.data.frame(cbind(x , pheno))


for(i in names(col)){
pdf(file.path(analysis.dir, paste0("PC12_20000mvCpGs_subsetted_",i,".pdf")),height = 5, width = 5)
print(ggscatter(x, x="PC1", y="PC2",
          color = i, shape = "Replicate",
          ellipse = F , mean.point = FALSE,palette= col[[i]],
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold")))
dev.off()
pdf(file.path(analysis.dir, paste0("PC23_20000mvCpGs_subsetted_",i,".pdf")),height = 5, width = 5)
print(ggscatter(x, x="PC2", y="PC3",
          color = i, shape = "Replicate",
          ellipse = F , mean.point = FALSE,palette= col[[i]],
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold")))
dev.off()
pdf(file.path(analysis.dir, paste0("PC34_20000mvCpGs_subsetted_",i,".pdf")),height = 5, width = 5)
print(ggscatter(x, x="PC3", y="PC4",
          color = i, shape = "Replicate",
          ellipse = F , mean.point = FALSE,palette= col[[i]],
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,3]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold")))
dev.off()
print(i)
}



#Donor and Genotype
pdf(file.path(analysis.dir, "PC12_20000mvCpGs_subsetted_DonorGenotype.pdf"),height = 5, width = 5)
ggscatter(x, x="PC1", y="PC2",
          color = "Donor", shape = "Genotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC23_20000mvCpGs_subsetted_DonorGenotyp.pdf"),height = 5, width = 5)
ggscatter(x, x="PC2", y="PC3",
          color = "Donor", shape = "Genotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC34_20000mvCpGs_subsetted_DonorGenotyp.pdf"),height = 5, width = 5)
ggscatter(x, x="PC3", y="PC4",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
          star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
          "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
#Donor and Epigenotype
pdf(file.path(analysis.dir, "PC12_20000mvCpGs_subsetted_DonorEpigenotype.pdf"),height = 5, width = 5)
ggscatter(x, x="PC1", y="PC2",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC23_20000mvCpGs_subsetted_DonorEpigenotype.pdf"),height = 5, width = 5)
ggscatter(x, x="PC2", y="PC3",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC34_20000mvCpGs_subsetted_DonorEpigenotype.pdf"),height = 5, width = 5)
ggscatter(x, x="PC3", y="PC4",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
          star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
          "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
#Celltype and tissue
pdf(file.path(analysis.dir, "PC12_20000Promoter_subsetted_CelltypeTissue.pdf"),height = 5, width = 5)
ggscatter(x, x="PC1", y="PC2",
          color = "Celltype", shape = "Tissue",label="Patient",repel=TRUE,
          ellipse = F , mean.point = FALSE,palette= c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC23_20000Promoter_subsetted_CelltypeTissue.pdf"),height = 5, width = 5)
ggscatter(x, x="PC2", y="PC3",
          color = "Celltype", shape = "Tissue",label="Patient",repel=TRUE,
          ellipse = F , mean.point = FALSE,palette= c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC34_20000Promoter_subsetted_CelltypeTissue.pdf"),height = 5, width = 5)
ggscatter(x, x="PC3", y="PC4",
          color = "Celltype", shape = "Tissue",label="Patient",repel=TRUE,
          ellipse = F , mean.point = FALSE,palette= c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),
          star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
          "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()

pdf(file.path(analysis.dir, "PC12_20000mvCpG_subsetted_DonorGenotype.pdf"),height = 5, width = 5)
ggscatter(x, x="PC1", y="PC2",
          color = "Donor", shape = "Genotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#f99f58", 
                    D213 = "#ee5529", D360 = "#dc263f", D124 = "#de5332", D123 = "#c33126", 
                    adult ="#c2c2c2", cordblood = "#a9a9a9"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC23_20000mvCpG_subsetted_DonorGenotyp.pdf"),height = 5, width = 5)
ggscatter(x, x="PC2", y="PC3",
          color = "Donor", shape = "Genotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#f99f58", 
                    D213 = "#ee5529", D360 = "#dc263f", D124 = "#de5332", D123 = "#c33126", 
                    adult ="#c2c2c2", cordblood = "#a9a9a9"),
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC34_20000mvCpG_subsetted_DonorGenotyp.pdf"),height = 5, width = 5)
ggscatter(x, x="PC3", y="PC4",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#f99f58", 
                    D213 = "#ee5529", D360 = "#dc263f", D124 = "#de5332", D123 = "#c33126", 
                    adult ="#c2c2c2", cordblood = "#a9a9a9"),
          star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
          "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()

#Sample Clustering
topVarCpGs<- head(order(rowVars(as.matrix(meth_per_cpg)), decreasing=TRUE),200000)
dissimilarity <- 1 - cor(meth_per_cpg[topVarCpGs,], use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
hrld<- hclust(distance)
dend<- hrld%>% as.dendrogram 
#col1 <- hrld$labels
#names(col1)<- anno[col1,]$group
#col1 <- col1[order.dendrogram(dend)]
#col1 <- ifelse(names(col1)=="TAM_tumor", "#EFC00099", ifelse(names(col1)=="BMDM_tumor","#0073C299",  ifelse(names(col1)=="TAM_healthy","#CD534CFF","#86868699")))
#dend <- dend %>% 
#set("branches_k_color", k=4, c("#EFC00099", "#0073C299", "#CD534CFF","#86868699" )) %>%
#set("branches_col",col1) %>%
#set("branches_lwd", 2) %>%
#set("labels_colors",col1) %>% 
#set("labels_cex", .6 )%>%
#set("leaves_pch", 19)%>% 
#set("leaves_cex", 1.5)%>% 
#set("leaves_col", col1)
pdf(file.path(analysis.dir, "Clustering_correlation_200000CpGs.pdf"), height = 7, width = 5)
dend %>% plot
dev.off()

#export bigwig
meth_per_cpg_all <- as.data.frame(bsseq::getMeth(bsseq_all, type = "raw"))
library(BSgenome.Mmusculus.UCSC.mm10)
mm10 <- BSgenome.Mmusculus.UCSC.mm10
anno <- granges(bsseq_all)
dir.create(file.path(analysis.dir,"bw"))
for(i in rownames(pheno)){
    temp<- anno 
    score(temp) <- meth_per_cpg_all[,i]
    temp<- temp[!is.na(score(temp)),]
    seqlengths(temp)<-seqlengths(mm10)[names(seqlengths(temp))]
    export.bw(temp, file.path(analysis.dir,"bw",paste0(i, ".bw")))
    print(i)
}




