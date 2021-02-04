#'---
#'title: "Enrichment analysis"
#'author: Lauren
#'output: 
#'  html_document: 
#'    toc: yes
#'    code_folding: hide
#'    highlight: tango
#'    theme: lumen
#'    df_print: paged
#'---
#' This is an example of one way to write one R script


#'
#' load packages
#'
#+ message=FALSE, warning=FALSE
#suppressPackageStartupMessages(library(optparse))
#suppressPackageStartupMessages(library(clusterProfiler))
#suppressPackageStartupMessages(library(enrichplot))
#suppressPackageStartupMessages(library(ReactomePA))
#suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(crosstalk))
#suppressPackageStartupMessages(library(DT))

#option_list = list(
#  make_option(c('-s','--species'), type="character", help="species in this analysis"),
#  make_option(c('-r','--de_result'), type="character", help="path to the DE result file"),
#  make_option(c("-d", '--directory'), type="character",help="directory to output html result"),
# make_option(c("-g", '--gsea'), action="store_true", default = FALSE, help="add this flag for GSEA")
#) 

#opt_parser <- OptionParser(option_list=option_list)
#opt <- parse_args(opt_parser)

#currentpath <- paste0(getwd(),'/scripts/enrichment_analysis.R')
#gsea_analysis <- opt$gsea
#output_dir <- opt$directory
#dir.create(output_dir,mode="0770")#create the output dir
#output_name <- gsub('.txt', '-enrichment_report', rev(unlist(strsplit(opt$de_result, '/')))[1], )

#database_list <- list(mouse='org.Mm.eg.db', human='org.Hg.eg.db')
#species <- tolower(opt$species)

#suppressPackageStartupMessages(library(database_list[[species]], character.only = TRUE))

#de_output <- read.csv(opt$de_result, header = T, sep = '\t')
#de_output <- na.omit(de_output)

#de_down <- de_output[de_output$regulation == 'Down',]
#de_up <- de_output[de_output$regulation == 'Up',]

#de_down_genelist <- de_down$ensemblID
#de_up_genelist <- de_up$ensemblID

#'
#' # Over Representation Analysis (GO term enrichment)
#' 
#' GO BP - Biological processes
#' 
#' GO CC - Cellular components
#' 
#' GO MF - Molecular functions
#'
#'## Down-regulated DEGs {.tabset}
#'
#' Some descriptions of the table
#'
#'### Enriched BP terms
#+ message=FALSE, warning=FALSE
#downGO_bp <- enrichGO(gene          = de_down_genelist,
#                OrgDb         = database_list[[species]],
#                keyType = "ENSEMBL",
#                ont           = "BP",
#                pAdjustMethod = "BH",
#                pvalueCutoff  = 0.01,
#                qvalueCutoff  = 0.05,
#                readable      = FALSE)

#downGO_bp <-simplify(downGO_bp, cutoff=0.7, by="p.adjust", select_fun=min)

bscols(
  datatable(downGO_bp[,c(1:7, 9)], 
            escape = F, selection = "none", extensions='Scroller',
            options = list(scrollY = 350, scrollX = TRUE, deferRender=T, scroller=T))
)

#'### Enriched CC terms
#+ message=FALSE, warning=FALSE
#downGO_cc <- enrichGO(gene          = de_down_genelist,
#                      OrgDb         = database_list[[species]],
#                      keyType = "ENSEMBL",
#                      ont           = "CC",
#                      pAdjustMethod = "BH",
#                      pvalueCutoff  = 0.01,
#                      qvalueCutoff  = 0.05,
#                      readable      = FALSE)

#downGO_cc <-simplify(downGO_cc, cutoff=0.7, by="p.adjust", select_fun=min)

bscols(
  datatable(downGO_cc[,c(1:7, 9)], 
            escape = F, selection = "none", extensions='Scroller',
            options = list(scrollY = 350, scrollX = TRUE, deferRender=T, scroller=T))
)


#'### Enriched MF terms
#+ message=FALSE, warning=FALSE
#downGO_mf <- enrichGO(gene          = de_down_genelist,OrgDb         = database_list[[species]],keyType = "ENSEMBL",
#                      ont           = "MF",pAdjustMethod = "BH",
#                      pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable      = FALSE)

#downGO_mf <-simplify(downGO_mf, cutoff=0.7, by="p.adjust", select_fun=min)

bscols(
  datatable(downGO_mf[,c(1:7, 9)], 
            escape = F, selection = "none", extensions='Scroller',
            options = list(scrollY = 350, scrollX = TRUE, deferRender=T, scroller=T))
)

#'
#'### Top 20 enriched GO_BP terms for the down-regulated genes
#' 
#' Some descriptions here
#' 
#+ message=FALSE, warning=FALSE, fig.width=8, fig.height=7
barplot(downGO_bp, showCategory=20,font.size = 9)+scale_size(range=c(2, 9))+
  scale_x_discrete(labels=function(y) stringr::str_wrap(y,width=25))
#' 
#'


#'## Up-regulated DEGs {.tabset}
#'
#' Some descriptions here
#'
#'### Enriched BP terms
#+ message=FALSE, warning=FALSE
#upGO_bp <- enrichGO(gene          = de_up_genelist,
#                      OrgDb         = database_list[[species]],
#                      keyType = "ENSEMBL",
#                      ont           = "BP",
#                      pAdjustMethod = "BH",
#                      pvalueCutoff  = 0.01,
#                      qvalueCutoff  = 0.05,
#                      readable      = FALSE)

#upGO_bp <-simplify(upGO_bp, cutoff=0.7, by="p.adjust", select_fun=min)

bscols(
  datatable(upGO_bp[,c(1:7, 9)], 
            escape = F, selection = "none", extensions='Scroller',
            options = list(scrollY = 350, scrollX = TRUE, deferRender=T, scroller=T))
)

#'### Enriched CC terms
#+ message=FALSE, warning=FALSE
#upGO_cc <- enrichGO(gene          = de_up_genelist,
#                      OrgDb         = database_list[[species]],
#                      keyType = "ENSEMBL",
#                      ont           = "CC",
#                     pAdjustMethod = "BH",
#                      pvalueCutoff  = 0.01,
#                      qvalueCutoff  = 0.05,
#                      readable      = FALSE)

#upGO_cc <-simplify(upGO_cc, cutoff=0.7, by="p.adjust", select_fun=min)

bscols(
  datatable(upGO_cc[,c(1:7, 9)], 
            escape = F, selection = "none", extensions='Scroller',
            options = list(scrollY = 350, scrollX = TRUE, deferRender=T, scroller=T))
)

#'### Enriched MF terms
#+ message=FALSE, warning=FALSE
#upGO_mf <- enrichGO(gene          = de_up_genelist,
#                    OrgDb         = database_list[[species]],
#                    keyType = "ENSEMBL",
#                    ont           = "MF",
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = 0.01,
#                    qvalueCutoff  = 0.05,
#                    readable      = FALSE)

#upGO_mf <-simplify(upGO_mf, cutoff=0.7, by="p.adjust", select_fun=min)

bscols(
  datatable(upGO_mf[,c(1:7, 9)], 
            escape = F, selection = "none", extensions='Scroller',
            options = list(scrollY = 350, scrollX = TRUE, deferRender=T, scroller=T))
)
#'
#'### Top 20 enriched GO_BP terms for the up-regulated genes
#' 
#' Some descriptions here
#' 
#+ message=FALSE, warning=FALSE, fig.width=8, fig.height=7
barplot(upGO_bp, showCategory=20,font.size = 9)+scale_size(range=c(2, 9))+
  scale_x_discrete(labels=function(y) stringr::str_wrap(y,width=25))
#' 
#'


#+ results='asis', echo=FALSE
# if true, hide
if (!gsea_analysis) {cat("<!---")}

#' # Gene set enrichment analysis 
#' 
#' Gene sets information from Reactome pathway database

#+ message=FALSE, warning=FALSE, eval=gsea_analysis
#use ReactomePA
#requires Entrez gene ID as input
#convert ensembl ID to entrez gene ID
#genelist<- bitr(de_output$ensemblID,fromType="ENSEMBL",toType=c("ENTREZID"),OrgDb=database_list[[species]])

#some input gene ID might fail to map to entrezID
#taking the IDs that have been successfully mapped
#mapped_genes<- de_output[which(genelist$ENSEMBL %in% de_output$ensemblID), ]

#calculate ranking score
#ranking <- -log10(mapped_genes$padj)*sign(mapped_genes$log2FoldChange)

#create ranked gene list
#ranks_rnaseq <- as.data.frame(cbind (genelist$ENTREZID,ranking))
#geneList = as.numeric(ranks_rnaseq[,2])
#names(geneList) = as.character(ranks_rnaseq[,1])
#geneList = sort(geneList, decreasing = TRUE)

y <- gsePathway(geneList, 
                organism = species,
                pvalueCutoff = 1,
                minGSSize = 15,
                pAdjustMethod = "BH", 
                verbose = FALSE)



#' 
#' ## table for all analyzed gene sets 
#' 
#' The table has been reordered decreasing by on the FDR qvalue 
#' 
#' More descriptions here on the columns
#' 
#+ message=FALSE, warning=FALSE, fig.width=8, fig.height=7, eval=gsea_analysis
yy<-as.data.frame(y)
yy <-yy[order(yy$qvalues),]
bscols(
  datatable(yy, 
            escape = F, selection = "none", extensions='Scroller',
            options = list(scrollY = 350, scrollX = TRUE, deferRender=T, scroller=T))
)

#' ## enrichment map for top 30 enriched gene sets
#' 
#' Some descriptions of the map here
#' 
#+ message=FALSE, warning=FALSE, fig.width=8, fig.height=7,eval=gsea_analysis
#emapplot(y, node_scale=1.5,layout="kk", showCategory = 30, color = "qvalues") 

#+ results='asis', echo=FALSE
#if true, hide <- default
if (!gsea_analysis) {cat("-->")}

# /* spin and knit this file to html
knitr::spin(hair = currentpath, knit = FALSE)
rmarkdown::render(currentpath, output_file = output_name, output_format = "html_document", output_dir = output_dir)
# */
