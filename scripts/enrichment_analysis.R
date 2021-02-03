
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(ggplot2)

keytypes(org.Mm.eg.db)
genelist_down<- bitr(gene,fromType="ENSEMBL",toType=c("ENTREZID"),OrgDb=org.Mm.eg.db)
x <- enrichPathway(gene=genelist_down$ENTREZID, pvalueCutoff = 0.01, readable=TRUE, organism = "mouse")

t<-read.csv('KO_21vsWT_21.txt', header = T, sep = '\t')
t<-t[,c(1:8)]
t <- na.omit(t)
t<-t[t$padj < 0.01,]

ranking <- -log10(t$pvalue)*sign(t$log2FoldChange)
ranks_rnaseq <- as.data.frame(cbind (t$ensemblID,ranking))

t_up <- na.omit(t[t$log2FoldChange > 1, 'ensemblID'])
t_down <- na.omit(t[t$log2FoldChange < -1, 'ensemblID'])

gene <- t$ensemblID
gene <- na.omit(gene)


ego <- enrichGO(gene          = gene,
                OrgDb         = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego_sim <-simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)



p2<-barplot(x, showCategory=10,font.size = 7, x='GeneRatio')+scale_size(range=c(2, 9))+
  scale_x_discrete(labels=function(y) stringr::str_wrap(y,width=25))
p2

p3<-barplot(ego_sim, showCategory=10,font.size = 7, x='GeneRatio')+scale_size(range=c(2, 9))+
  scale_x_discrete(labels=function(y) stringr::str_wrap(y,width=25))
p3



dotplot(ego, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(x, showCategory=30) + ggtitle("dotplot for ORA")
emapplot(ego, node_scale=1.5,layout="kk") 
emapplot(x, node_scale=1.5,layout="kk")+scale_size(range=c(2, 9))
cowplot::plot_grid(p3, p2, nrow=2, labels=LETTERS[1:2])



########
suppressPackageStartupMessages(library(optparse))

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(crosstalk))
suppressPackageStartupMessages(library(DT))
database_list <- list(mouse='org.Mm.eg.db', human='org.Hg.eg.db')
species <- tolower('Mouse')

suppressPackageStartupMessages(library(database_list[[species]], character.only = TRUE))



t<-read.csv('KO_21vsWT_21.txt', header = T, sep = '\t')
t<-t[,c(1:8)]
t <- na.omit(t)#ignore cutoff sig

t_p0.01<-t[t$padj < 0.01,]

down <- t_p0.01[t_p0.01$log2FoldChange < 0, 'ensemblID']

ego <- enrichGO(gene          = down,
                OrgDb         = database_list[[species]],
                keyType = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = FALSE)

ego <-simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
gene_ont <- ego[,c(1:7, 9)]


ego_MF <- enrichGO(gene          = down,
                OrgDb         = database_list[[species]],
                keyType = "ENSEMBL",
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = FALSE)
ego_MF <-simplify(ego_MF, cutoff=0.7, by="p.adjust", select_fun=min)
geneMF_ont <- ego_MF[,c(1:7, 9)]

ego_CC <- enrichGO(gene          = down,
                   OrgDb         = database_list[[species]],
                   keyType = "ENSEMBL",
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = FALSE)
ego_CC <-simplify(ego_CC, cutoff=0.7, by="p.adjust", select_fun=min)
geneCC_ont <- ego_CC[,c(1:7, 9)]



p2<-barplot(ego_MF, showCategory=20,font.size = 9)+scale_size(range=c(2, 9))+
  scale_x_discrete(labels=function(y) stringr::str_wrap(y,width=25))
p2

###### GSEA

#use ReactomePA
#requires Entrez gene ID as input
#convert ensembl ID to entrez gene ID
genelist<- bitr(t$ensemblID,fromType="ENSEMBL",toType=c("ENTREZID"),OrgDb=database_list[[species]])

#some input gene ID might fail to map to entrezID
#taking the IDs that have been successfully mapped
mapped_genes<- t[which(genelist$ENSEMBL %in% t$ensemblID), ]

#calculate ranking score
ranking <- -log10(mapped_genes$padj)*sign(mapped_genes$log2FoldChange)

#create ranked gene list
ranks_rnaseq <- as.data.frame(cbind (genelist$ENTREZID,ranking))
geneList = as.numeric(ranks_rnaseq[,2])
names(geneList) = as.character(ranks_rnaseq[,1])
geneList = sort(geneList, decreasing = TRUE)





#write.table(ranks_rnaseq, 'KO_21vsWT_21.rnk',col.name=FALSE,sep="\t",row.names=FALSE,quote=FALSE)
#ranks_rnaseq<-ranks_rnaseq[-2379,]


y <- gsePathway(geneList, 
                organism = 'mouse',
                pvalueCutoff = 1,
                minGSSize = 15,
                pAdjustMethod = "BH", 
                verbose = FALSE)
yyy<-y[(y[,'qvalues'] < 0.25),]
yy<-as.data.frame(y)

#number of the total gene set being analyzed
nrow(y)

# number of gene sets have false discovery rate of less than 25%
# table for these gene sets
sum(y[,'qvalues'] < 0.25)

down_cc<-bscols(
  datatable(y[y[,'qvalues'] < 0.25,], 
            escape = F, selection = "none", extensions='Scroller',
            options = list(scrollY = 350, scrollX = TRUE, deferRender=T, scroller=T))
)

#enrichment map
p4 <- emapplot(y, node_scale=1.5,layout="kk", showCategory = 30, color = "qvalues") 

edox <- setReadable(y, 'org.Mm.eg.db', 'ENTREZID')
heatplot(edox, showCategory = 20, foldChange=mapped_genes$log2FoldChange)




genelist2<- bitr(t$ensemblID,fromType="ENSEMBL",toType=c("SYMBOL"),OrgDb=database_list[[species]])

#some input gene ID might fail to map to entrezID
#taking the IDs that have been successfully mapped
mapped_genes2<- t[which(genelist$ENSEMBL %in% t$ensemblID), ]

#calculate ranking score
ranking2 <- -log10(mapped_genes2$padj)*sign(mapped_genes2$log2FoldChange)

#create ranked gene list
ranks_rnaseq2 <- as.data.frame(cbind (genelist2$SYMBOL,ranking2))
geneList2 = as.numeric(ranks_rnaseq2[,2])
names(geneList2) = as.character(ranks_rnaseq2[,1])
geneList2 = sort(geneList2, decreasing = TRUE)


c5 <- read.gmt('Mouse_Human_Reactome_January_13_2021_symbol.gmt')
y2 <- GSEA(geneList2, TERM2GENE=c5,minGSSize = 15, pvalueCutoff = 1)
y22<-y2[y2[,'qvalues'] < 0.25,]
emapplot(y22, node_scale=1.5,layout="kk", showCategory = 30, color = "qvalues") 
newID <- unlist(lapply(y2[,'ID'], function(x){unlist(strsplit(x, '%'))[1]}))
y2$ID <- newID
