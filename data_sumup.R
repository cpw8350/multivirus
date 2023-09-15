
{
  library(BiocManager)
  library(hgu133a2.db)
  library(org.Hs.eg.db)
  library(illuminaHumanv4.db)
  library(GEOquery)
  library(tidyr)
  library(sva)
  library(dplyr)
  library(limma)
  library(openxlsx)
  library(openxlsx)
  library(hgu133a2.db)
  library(hgu95av2.db)
  library(org.Hs.eg.db)
  library(illuminaHumanv4.db)
  library(GEOquery)
  library(tidyr)
  library(sva)
  library(bladderbatch)
  library(dplyr)
  library(limma)
  library(ggplot2)
  library(reshape)
  library(reshape2)
  library(RColorBrewer)
  library(ggsci)
  library(WGCNA)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(futile.logger)
  library(VennDiagram) 
  library(tibble)
}


#data load
rm(list=ls())
gc()
options(stringsAsFactors = FALSE)
setwd("D:/susceptibility/data/GEO/数据预处理/RSV")  
exprs00_RSV<-read.csv('exprs00_RSV.csv',row.names = 1)  #64*11644
exprs00_RSV <- exprs00_RSV %>% rownames_to_column("symbol")
RSV_group <- read.csv("group_RSV.csv",row.names = 1)

setwd("D:/susceptibility/data/GEO/数据预处理/HRV/sym+infect")  
exprs00_HRV<-read.csv('exprs00_HRV.csv',row.names = 1)  #64*11644
HRV_group <- read.csv("group_HRV.csv",row.names = 1)

setwd("D:/susceptibility/data/GEO/数据预处理/H1N1/作者给的信息")
exprs00_h1n1<-read.csv('exprs00_h1n1.csv',row.names = 1)  #64*11644
h1n1_group <- read.csv("group_h1n1.csv",row.names = 1)
#expression data
exprs<-Reduce(function(x,y) merge(x,y,by='symbol'),list(exprs00_RSV,exprs00_h1n1,exprs00_HRV),accumulate = FALSE)
row.names(exprs)<-exprs$symbol
exprs<-exprs[,-1]

#group
batch<-data.frame(subject=colnames(exprs),batch=factor(c(rep("RSV",20),rep("H1N1",65),rep("HRV",32))))
group1 <- rbind(RSV_group,h1n1_group,HRV_group)
write.csv(group1,file="sym_group.csv",row.names = T)
#scale
library('limma')
exprs_n<-normalizeBetweenArrays(exprs)   ##use the limma to do the chip normalization
#boxplot(exprs00_h1n1_n,outline=FALSE,notch=T,col=group_h1n1,las=2,main='the distribution of H1N1_group after correction')
exprs_1 <- as.data.frame(exprs_n)
exprs_1$gene <- rownames(exprs_1)
library('reshape')
exprs_2 <- melt(exprs_1,id="gene")
colnames(exprs_2)<- c("gene","subject","expression")
exprs00_2 <- merge(exprs_2,group1,by="subject")

setwd("D:/susceptibility/data/GEO/数据预处理")

#batch effect
remove_batch<-function(exprs,group_infor,dataname){
  library(sva)
  exprs<-as.matrix(exprs)
  group_infor<-group_infor[colnames(exprs),]
  mod0<-model.matrix(~as.factor(group),data = group_infor)
  combat_edata = ComBat(dat=exprs, batch=group_infor$batch, mod=mod0, par.prior=TRUE, prior.plots=FALSE)
  #write.csv(combat_edata, paste("D:/susceptibility/data/GEO/数据预处理/",dataname,"_combat_data没有标准化.csv",sep = ""))
}
remove_batch(exprs,group1,"all")  


#pca
mytheme<-theme_bw()+theme(legend.position="right",
                          #panel.border=element_blank(),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          plot.title=element_text(size=20,
                                                  colour="black",
                          ), #family="CA"),
                          axis.title =element_text(size=20,
                                                   colour="black",
                          ),#family="CA"),
                          legend.text=element_text(size=20,colour="black",
                          ),#family="CA"),
                          legend.key=element_blank(),
                          axis.text=element_text(size=20,colour="black",
                          ),#family="CA"),
                          strip.text=element_text(size=20,colour="#085A9C",
                          ),#family="CA"),
                          strip.background=element_blank())

do_pca<-function(exprs,group,picturename){
  pca<-prcomp(t(exprs),scale= T)
  pca.var<-pca$sdev  
  pca.var.per<-round(pca.var/sum(pca.var)*100 ,1) 
  pca.data<-data.frame(subject=row.names(pca$x),x=pca$x[,1],y=pca$x[,2]) 
  pca.data<- merge(pca.data,group,by="subject")
  pp<- ggplot(data=pca.data,aes(x=x,y=y))+#,shape=group
    geom_point(aes(colour=batch),size=5)+mytheme+ scale_color_lancet()+
    xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
    ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
  pdf(paste("D:/susceptibility/data/GEO/数据预处理/",picturename,"batch.pdf",sep = "_"),height=8,width=8)
  print(pp)
  dev.off()
  png(paste("D:/susceptibility/data/GEO/数据预处理/",picturename,"batch.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
  pp<- ggplot(data=pca.data,aes(x=x,y=y))+geom_point(aes(colour=group),size=5)+ scale_color_lancet()+mytheme+
    xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
    ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
  pdf(paste("D:/susceptibility/data/GEO/数据预处理/",picturename,"group.pdf",sep = "_"),height=8,width=8)
  print(pp)
  dev.off()
  png(paste("D:/susceptibility/data/GEO/数据预处理/",picturename,"group.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
}
##apply the do_pca function to check the performance of removing batch effects
do_pca(exprs,group1,"group_before")
do_pca(exprs_combat,group1,"combat")



#################################################################################################################
#wgcna
library(tibble)
setwd("D:/susceptibility/data/GEO/数据预处理/集合数据")
exprs_combat<-read.csv('all_combat_data.csv',row.names = 1)

exprs_combat <- t(exprs_combat)
exprs_combat <- as.data.frame(exprs_combat)

{
library(WGCNA)
library(stringr)
library(ggplot2)
library(grid)
library(futile.logger)
library(VennDiagram) 
}
mytheme<-theme_bw()+theme(legend.position="right",
                          #panel.border=element_blank(),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          plot.title=element_text(size=16,
                                                  colour="black",
                          ), #family="CA"),
                          axis.title =element_text(size=16,
                                                   colour="black",
                          ),#family="CA"),
                          legend.text=element_text(size=16,colour="black",
                          ),#family="CA"),
                          legend.key=element_blank(),
                          axis.text=element_text(size=16,colour="black",
                          ),#family="CA"),
                          strip.text=element_text(size=16,colour="#085A9C",
                          ),#family="CA"),
                          strip.background=element_blank())

###all data  expression data: sym1 group: sym_group
gsg = goodSamplesGenes(exprs_combat, verbose = 3)
gsg$allOK;
#==========================================[ data prepare ]===========================================


exprs<-apply(exprs_combat,1,var)
exprs<-as.data.frame(t(exprs_combat[which(exprs>quantile(exprs, probs = seq(0, 1, 0.2))[5]),]))  ##2329*65
Traits <- subset(group1,select=c("gseID","group"))
Traits<- Traits[rownames(exprs),]
sampletree<-hclust(dist(exprs),method = 'average')

traitColors = numbers2colors(as.numeric(as.factor(Traits$group)), signed = TRUE,centered=TRUE)
plotDendroAndColors(sampletree, 
                    traitColors,
                    groupLabels = Traits$type,
                    cex.dendroLabels = 0.8,
                    marAll = c(1,4,3,1),
                    cex.rowText = 0.02,
                    main = "Sample dendrogram and trait heatmap")

#===================================[ build network ]===============================
#soft_threshold
powers<-c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(exprs, powerVector = powers,verbose = 5,)  
sizeGrWindow(9, 5)
par(mfrow = c(1,2)) 
cex1 = 0.90
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

sft$powerEstimate
####One-step building blocks:14 modules

cor <- WGCNA::cor
net = blockwiseModules(exprs, power =sft$powerEstimate,maxBlockSize = 5000, TOMType = "unsigned", minModuleSize = 15,reassignThreshold = 0, 
                       mergeCutHeight = 0.2,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",verbose = 5,deepSplit = 4) #sft$powerEstimate 
table(net$colors) ##minModuleSize=20:12;15:16
cor<-stats::cor
mergedColors = labels2colors(net$colors)
#write.csv(table(mergedColors),'mergedColors.csv')
number_module <- data.frame(module=names(table(mergedColors)),number=table(mergedColors))
number_module$module_color <- sort(unique(mergedColors[net$blockGenes[[1]]]))

pp<- ggplot(number_module,aes(x=module,y=number.Freq,fill=module))+
  geom_bar(stat = "identity",width = 0.8, alpha = 0.8)+
  scale_fill_manual(values = number_module$module_color)+mytheme+
  geom_text(aes(label = number.Freq),color="#E6E8FA", position = position_dodge(0.1))+
  theme(axis.text.x=element_text(size=10,angle = 90))+
  theme(legend.text = element_text(size = 15))
print(pp)

#dendrogram
gene_color <-mergedColors[net$blockGenes[[1]]]
plotDendroAndColors(net$dendrograms[[1]],gene_color, "Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
###save variable
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = "wgcna_median-networkConstruction-auto.RData")

#=======================================[ Screening modules related to phenotype ]===================================

nGenes = ncol(exprs)
nSamples = nrow(exprs)
design=model.matrix(~0+ as.factor(Traits$group))
colnames(design)=c('asym','sym')
##The first principal component of each module
MEs0 = moduleEigengenes(exprs, moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor = cor(MEs, design,use = "p",method = "spearman")
moduleTraitCor1<-moduleTraitCor[order(moduleTraitCor[,"sym"],decreasing = TRUE),]
moduleTraitPvalue = corPvalueStudent(moduleTraitCor1, nSamples)
#signif保留小数
textMatrix = paste(signif(moduleTraitCor1, 2), "\n(",signif(moduleTraitPvalue, 2), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor1)
sizeGrWindow(9,7)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor1,
               xLabels = colnames(design),
               yLabels = rownames(moduleTraitCor1),
               ySymbols = rownames(moduleTraitCor1),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab = 0.5,
               cex.text = 0.5,
               zlim = c(-1,1),
               width = 0.2,
               main = paste("Module-trait relationships"))

MET = orderMEs(cbind(MEs, design))
#################
library(reshape)
library(reshape2)
cor<-as.data.frame(moduleTraitCor1)
cor$Module <- rownames(cor)
cor_melt <-melt(cor,id="Module")
colnames(cor_melt)<-c("Module","group","R2")

pp<-as.data.frame(moduleTraitPvalue)
pp$Module <- rownames(pp)
pp_melt <-melt(pp,id="Module")
colnames(pp_melt)<-c("Module","group","Pvalue")
cor_pp <-cbind(cor_melt,pp_melt)
cor_pp<-cor_pp[,-c(4,5)]
cor_pp<-cor_pp[order(cor_pp$Module),]
cor_pp$mm_color <- rep(number_module$module,each=2) 
ggplot(cor_pp,aes(R2,-log10(Pvalue),color=mm_color)) + 
  geom_point(aes(size=-log10(Pvalue),shape=group)) +
  labs(x="R2",y="- Log10Pvalue",fill="") +
  scale_color_manual(values =unique(cor_pp$mm_color) )+
  theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
        panel.background = element_blank(),
        panel.border=element_rect(fill='transparent', color='black'), 
        axis.ticks.x = element_blank(),
        axis.title.x  = element_text(size = 16),
        axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        legend.background = element_blank())+
  geom_hline(aes(yintercept = -log10(0.05)),colour="grey",linetype=2,size=1)
library(pheatmap)
pp<-as.data.frame(moduleTraitCor1)
pp<-pp[order(pp$asym),]
anno_row <- data.frame("Module"=rownames(moduleTraitCor1))
rownames(anno_row)<-rownames(pp)

ann_colors <- list(Module=c(MEbrown="brown",MEmagenta="magenta",MEyellow="yellow",
                            MEblue="blue",MEblack="black",MEpurple="purple",MEcyan="cyan",MEgreen="green",
                            MEgreenyellow='greenyellow',MElightyellow="lightyellow",MEroyalblue="royalblue",
                            MEmidnightblue="midnightblue",MEgrey60="grey60",MElightgreen="lightgreen",
                            MEred="red",MEgreen="green" ,MElightcyan="lightcyan",
                            MEgrey="grey",MEpink= "pink", MEturquoise="turquoise",
                            MEtan="tan",MEsalmon="salmon"))

ann_colors <- list(Module=c(MEturquoise="turquoise",MEblue="blue",MEmagenta="#FF80FFFF",MEblack="black",
                            MElightcyan="lightcyan",MEyellow="yellow",MEmidnightblue="midnightblue",MEgreenyellow='greenyellow',
                            MEpurple="purple",MEgrey="grey",MEcyan="cyan",MEbrown="brown",MEtan="tan",MEsalmon="salmon",
                            MEgreen="green",MEred="red",MEpink= "pink"))


pheatmap(pp,cluster_rows=F,
         show_rownames = FALSE,
         cluster_cols=F,
         display_numbers=T,
         number_format="%.2f",
         border="white",
         fontsize_number=8,
         fontsize_col = 8,
         fontsize_row = 8,
         #angle_col = 90,
         legend_breaks=c(-0.5,-0.2,0,0.2,0.5),
         legend_labels=c("-0.5","-0.2","0","0.2","0.5"),
         #color = mycol,
         fontsize=8,
         annotation_row=anno_row,
         annotation_colors=ann_colors)

#
g <- colnames(exprs)
a <- rbind(g,mergedColors)
a <- t(a) %>% as.data.frame()
a1 <- a[a$mergedColors=='red',]
a2 <- a[a$mergedColors=='blue',]
a3 <- a[a$mergedColors=='greenyellow',]
a4 <- a[a$mergedColors=='midnightblue',]
b <- rbind(a1,a2,a3,a4)

write.csv(b,'modules.csv')

b$g
diff1 <- b$g
gene.df1 <- bitr(diff1,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
count(gene.df1,ENTREZID) %>% dplyr::filter(n>1) %>% glimpse()
gene1 <- gene.df1$ENTREZID
ego_BP <- enrichGO(gene = gene1,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
ego_result_BP <- as.data.frame(ego_BP)[1:20, ]
ego_result_BP <- ego_result_BP[order(ego_result_BP$p.adjust,decreasing=F),]
options(repr.plot.width =15, repr.plot.height =100)
ggplot(data=ego_result_BP, aes(x=reorder(Description,-p.adjust),y=-log10(p.adjust))) + #横纵轴取值
  geom_bar(stat="identity", width=0.8,color='#66C3A5',fill='#66C3A5') + #柱状图的宽度，可以自己设置
  #scale_fill_manual(values = 'green') + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("-Log10Padjust") + 
  labs(title = "The Most Enriched GO Terms")+
  theme(axis.text = element_text(size = 15))


#####
exprs <- read.csv('all_combat_data.csv',row.names = 1)
exprs.gnames<-rownames(exprs)
group2<-group1[colnames(exprs),]

coldata <- data.frame(row.name=colnames(exprs),condition=exprs.cl)

exprs.cl <- as.numeric(as.factor(group1$group)) #group信息数字化
exprs.cl[exprs.cl==1]<-0
exprs.cl[exprs.cl==2]<-1
exprs.cl
library(RankProd)
RP.adv.out <- RP.advance(exprs,exprs.cl,group1$batch, logged=TRUE,na.rm=FALSE,  gene.names=exprs.gnames,
                         plot = FALSE,rand = 123,MinNumOfValidPairs =1)     ###class1:asym,class2:sym

#
degs<-topGene(RP.adv.out,cutoff = 0.1,method="pfp",logged=TRUE,logbase=2,gene.names=exprs.gnames) 
degs<-topGene(RP.adv.out,cutoff = 1,method="pfp",logged=TRUE,logbase=2,gene.names=exprs.gnames) 

down<-degs$Table1
up<-degs$Table2
#
commgene <- intersect(rownames(up),rownames(down))

#
up<-up[!rownames(up)%in%commgene,]  
down<-down[!rownames(down)%in%commgene,]






















